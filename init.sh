#!/usr/bin/env bash
set -euo pipefail

ROOT=$(pwd)

# Making needed directories
mkdir -p data/models

if [[ "$#" -gt 0 ]]
then
    
    if [[ "$#" -gt 1 ]]
    then
        echo "
This script either runs without an argument, or with the --data flag.
You supplied $# arguments.
        "
        exit 1;
    fi
    
    if [[ "$1" != "--data" ]]
    then
        echo "Unknown argument $1"
        exit 1;
    fi
    
    echo "

####################
# Downloading data #
####################

    "
    mkdir -p "$ROOT/temp_data" && cd "$ROOT/temp_data"
    
    echo "Downloading data files"
    # Whole CHM13 genome
    wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.3_CHM13_T2T_v1.1/GCA_009914755.3_CHM13_T2T_v1.1_genomic.fna.gz"
    
    # Drosophila genome
    wget "https://obj.umiacs.umd.edu/marbl_publications/hicanu/dmel_hifi_40x.fasta.gz"#reads
    wget "http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.35_FB2020_04/fasta/dmel-all-chromosome-r6.35.fasta.gz"
    
    # Preprocessing data
    echo "Pre-Processing Whole human genome"
    gunzip -c "./GCA_009914755.3_CHM13_T2T_v1.1_genomic.fna.gz" | \
    "$ROOT/tools/fastatools/fastatools" upper | \
    "$ROOT/tools/fastatools/fastatools" rename --regex '\.[0-9] Homo.*' > \
    "$ROOT/data/whole_human_genome.fa"
    
    echo "Pre-processing drosophila genome"
    gunzip -c "dmel-all-chromosome-r6.35.fasta.gz" | \
    "$ROOT/tools/fastatools/fastatools" subset -s 1 -e 7 > \
    "$ROOT/data/whole_drosophila_genome.fa"
    
    mv "dmel_hifi_40x.fasta.gz" "$ROOT/data/real_drosophila_reads.fa"
    
    # Move tandemtools genome
    cp "$ROOT/tools/TandemTools/test_data/simulated_del.fasta" "$ROOT/data/tandemtools.ref.fa"
    
    cd "$ROOT" && rm -rf "$ROOT/temp_data"
fi

exit 0

echo "

######################
# Setting up NanoSim #
######################

"
ln -s "$ROOT/tools/NanoSim/src/simulator.py" "$ROOT/bin/nanosim.py"
echo "Extracting nanosim pre-trained model"
tar -xzf "$ROOT/tools/NanoSim/pre-trained_models/human_NA12878_DNA_FAB49712_guppy_flipflop.tar.gz" --directory "$ROOT/data/models"



echo "

##################
# Building PBSIM #
##################

"
cd "$ROOT/tools/pbsim2"
./configure && make
cd "$ROOT"
ln -s "$ROOT/tools/pbsim2/src/pbsim" "$ROOT/bin/pbsim"
ln -s "$ROOT/tools/pbsim2/data/P6C4.model" "$ROOT/data/models"




echo "

#####################
# Building Minimap2 #
#####################

"
if [[ "$OSTYPE" == "darwin"* ]] #MacOs
then
    cd "$ROOT/tools/minimap2"
    make
    cd "$ROOT"
    ln -s "$ROOT/tools/minimap2/minimap2" "$ROOT/bin"
else #Linux
    echo "Downloading minimap2 release"
    cd "$ROOT/tools"
    wget "https://github.com/lh3/minimap2/releases/download/v2.22/minimap2-2.22_x64-linux.tar.bz2"
    tar -xjvf "minimap2-2.22_x64-linux.tar.bz2"
    cd "$ROOT"
    ln -s "$ROOT/tools/minimap2-2.22_x64-linux/minimap2" "$ROOT/bin"
fi



echo "

######################
# Building Winnowmap #
######################

"
cd "$ROOT/tools/Winnowmap"
if [[ "$OSTYPE" == "darwin"* ]] #MacOS
then
    if command -v gcc-10 &> /dev/null
    then
        GCC="gcc-10"
        GPP="g++-10"
    elif command -v gcc-11 &> /dev/null
    then
        GCC="gcc-11"
        GPP="g++-11"
    else
        echo "You need at least gcc 10 to compile Winnowmap"
        echo "Using Homebrew: brew install gcc@10"
        exit 1
    fi
    make CC=$GCC CXX=$GPP
else #Linux
    make
fi
cd "$ROOT"
ln -s "$ROOT/tools/Winnowmap/bin/meryl" "$ROOT/bin/meryl"
ln -s "$ROOT/tools/Winnowmap/bin/winnowmap" "$ROOT/bin/winnowmap"


echo "

######################
# Building Bedtools #
######################

"
if [[ "$OSTYPE" == "darwin"* ]]
then
    cd "$ROOT/tools/bedtools2"
    make
    cd "$ROOT"
    ln -s "$ROOT/tools/bedtools2/bin/bedtools" "$ROOT/bin"
else
    echo "Downloading bedtools release"
    cd "$ROOT/tools"
    wget "https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary"
    chmod +x "bedtools.static.binary"
    cd "$ROOT"
    ln -s "$ROOT/tools/bedtools.static.binary" "$ROOT/bin/bedtools"
fi


echo "

#############################
# Setting up Go executables #
#############################

"
# Setup rename_sequences and reduce_sequences
if [[ "$OSTYPE" == "darwin"* ]] #MacOS
then
    echo "Extracting MacOS reduce and rename binaries"
    xz -cdk "$ROOT/tools/reduce_sequences/bin/macos_amd64.xz" > "$ROOT/bin/reduce_sequences"
    xz -cdk "$ROOT/tools/rename_sequences/bin/macos_amd64.xz" > "$ROOT/bin/rename_sequences"
else #Linux
    echo "Extracting Linux reduce and rename binaries"
    xz -cdk "$ROOT/tools/reduce_sequences/bin/linux_amd64.xz" > "$ROOT/bin/reduce_sequences"
    xz -cdk "$ROOT/tools/rename_sequences/bin/linux_amd64.xz" > "$ROOT/bin/rename_sequences"
fi
chmod +x "$ROOT/bin/reduce_sequences"
chmod +x "$ROOT/bin/rename_sequences"


echo "Succesfully initialized pipeline environment!"
