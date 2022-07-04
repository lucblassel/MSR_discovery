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

##########################
# Setting up bigBedToBed #
##########################

"
if [[ ! -f "$ROOT/bin/bigBedToBed" ]]; then
    if [[ "$OSTYPE" == "darwin"* ]]; then # MacOS
        wget "https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigBedToBed"
    else # Linux
        wget "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed"
    fi
    chmod +x "bigBedToBed"
    mv "bigBedToBed" "$ROOT/bin/bigBedToBed"
else
    echo "bigBedToBed is already present"
    echo ""
fi


echo "

##########################
# Setting up FastaTools #
##########################

"
if [[ ! -f "$ROOT/bin/fastatools" ]]; then
    ln -s "$ROOT/tools/fastatools/fastatools" "$ROOT/bin/fastatools"
else
    echo "FastaTools is already present"
    echo ""
fi


echo "

######################
# Setting up NanoSim #
######################

"
if [[ ! -f "$ROOT/bin/nanosim.py" ]]; then
    ln -s "$ROOT/tools/NanoSim/src/simulator.py" "$ROOT/bin/nanosim.py"
else
    echo "Nanosim executable already present."
fi
if [[ ! -d "$ROOT/data/models/human_NA12878_DNA_FAB49712_guppy_flipflop" ]]; then
    echo "Extracting nanosim pre-trained model"
    tar -xzf "$ROOT/tools/NanoSim/pre-trained_models/human_NA12878_DNA_FAB49712_guppy_flipflop.tar.gz" --directory "$ROOT/data/models"
else
    echo "Nanosim model already present."
fi



echo "

##################
# Building PBSIM #
##################

"
if [[ ! -f "$ROOT/bin/pbsim" ]]; then
    cd "$ROOT/tools/pbsim2"
    ./configure && make
    cd "$ROOT"
    cp "$ROOT/tools/pbsim2/src/pbsim" "$ROOT/bin/pbsim"
else
    echo "pbsim binary already present."
fi
if [[ ! -f "$ROOT/data/models/P6C4.model" ]]; then
    cp "$ROOT/tools/pbsim2/data/P6C4.model" "$ROOT/data/models"
else
    echo "Pbsim model already present"
fi


echo "

#####################
# Building Minimap2 #
#####################

"
if [[ ! -f "$ROOT/bin/minimap2" ]]; then
    if [[ "$OSTYPE" == "darwin"* ]]; then # MacOs
        cd "$ROOT/tools/minimap2"
        make
        cd "$ROOT"
        cp "$ROOT/tools/minimap2/minimap2" "$ROOT/bin"
    else # Linux
        echo "Downloading minimap2 release"
        cd "$ROOT/tools"
        wget "https://github.com/lh3/minimap2/releases/download/v2.22/minimap2-2.22_x64-linux.tar.bz2"
        tar -xjvf "minimap2-2.22_x64-linux.tar.bz2"
        rm "minimap2-2.22_x64-linux.tar.bz2"
        cd "$ROOT"
        cp "$ROOT/tools/minimap2-2.22_x64-linux/minimap2" "$ROOT/bin"
    fi
else
    echo "minimap2 binary already present."
fi


echo "

#####################
# Building Samtools #
#####################

"
if [[ ! -f "$ROOT/bin/samtools" ]]; then
    cd "$ROOT/tools"
    wget "https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2"
    tar -xjvf "samtools-1.12.tar.bz2"
    rm "samtools-1.12.tar.bz2"
    cd "samtools-1.12"
    autoheader
    autoconf -Wno-syntax
    ./configure --without-curses --disable-lzma --disable-bz2
    make
    cp samtools "$ROOT/bin/samtools"
else
    echo "samtools binary already present."
fi

echo "

######################
# Building Winnowmap #
######################

"
if [[ ! -f "$ROOT/bin/winnowmap" || ! -f "$ROOT/bin/meryl" ]]; then
    cd "$ROOT/tools/Winnowmap"
    if [[ "$OSTYPE" == "darwin"* ]]; then # MacOS
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
    else # Linux
        make
    fi
    cd "$ROOT"
    cp "$ROOT/tools/Winnowmap/bin/meryl" "$ROOT/bin/meryl"
    cp "$ROOT/tools/Winnowmap/bin/winnowmap" "$ROOT/bin/winnowmap"
else
    echo "winnowmap and meryl binaries already present."
fi

echo "

#####################
# Building Bedtools #
#####################

"
if [[ ! -f "$ROOT/bin/bedtools" ]]; then
    if [[ "$OSTYPE" == "darwin"* ]]; then # MacOS
        cd "$ROOT/tools/bedtools2"
        make
        cd "$ROOT"
        cp "$ROOT/tools/bedtools2/bin/bedtools" "$ROOT/bin"
    else # Linux
        echo "Downloading bedtools release"
        cd "$ROOT/tools"
        wget "https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary"
        chmod +x "bedtools.static.binary"
        cd "$ROOT"
        cp "$ROOT/tools/bedtools.static.binary" "$ROOT/bin/bedtools"
    fi
else
    echo "bedtools binary already present."
fi

echo "

###############################
# Setting up reduce_sequences #
###############################

"
if [[ ! -f "$ROOT/bin/reduce_sequences" ]]; then
    if [[ "$OSTYPE" == "darwin"* ]]; then # MacOS
        xz -cdk "$ROOT/tools/reduce_sequences/bin/macos_amd64.xz" > "$ROOT/bin/reduce_sequences"
    else # Linux
        xz -cdk "$ROOT/tools/reduce_sequences/bin/linux_amd64.xz" > "$ROOT/bin/reduce_sequences"
    fi
    chmod +x "$ROOT/bin/reduce_sequences"
else
    echo "reduce_sequences binary already present."
fi

echo "

###############################
# Setting up rename_sequences #
###############################

"
if [[ ! -f "$ROOT/bin/rename_sequences" ]]; then
    if [[ "$OSTYPE" == "darwin"* ]]; then # MacOS
        xz -cdk "$ROOT/tools/rename_sequences/bin/macos_amd64.xz" > "$ROOT/bin/rename_sequences"
    else # Linux
        xz -cdk "$ROOT/tools/rename_sequences/bin/linux_amd64.xz" > "$ROOT/bin/rename_sequences"
    fi
    chmod +x "$ROOT/bin/rename_sequences"
else
    echo "rename_sequences binary already present."
    echo ""
fi


echo "Succesfully initialized pipeline environment!"


    echo "

####################
# Downloading data #
####################

    "
    mkdir -p "$ROOT/temp_data" && cd "$ROOT/temp_data"

    if [[ ! -f "$ROOT/data/whole_human_genome.fa" ]]; then
        wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.3_T2T-CHM13v1.1/GCA_009914755.3_T2T-CHM13v1.1_genomic.fna.gz"
        echo "Pre-Processing Whole human genome"
        gunzip -c "./GCA_009914755.3_T2T-CHM13v1.1_genomic.fna.gz" | \
        "$ROOT/tools/fastatools/fastatools" upper | \
        "$ROOT/tools/fastatools/fastatools" rename --file "$ROOT/tools/renamer.txt" > \
        "$ROOT/data/whole_human_genome.fa"
    else
        echo "Reference Human genome already downloaded and processed."
    fi

    if [[ ! -f "$ROOT/data/whole_drosophila_genome.fa" ]]; then
        wget "http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.35_FB2020_04/fasta/dmel-all-chromosome-r6.35.fasta.gz"
        echo "Pre-processing drosophila genome"
        gunzip -c "dmel-all-chromosome-r6.35.fasta.gz" | \
        "$ROOT/tools/fastatools/fastatools" subset -s 1 -e 7 > \
        "$ROOT/data/whole_drosophila_genome.fa"
    else
        echo "Reference Drosophila genome already downloaded and processed."
    fi

    if [[ ! -f "$ROOT/data/whole_ecoli_genome.fa" ]]; then
       wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=nucleotide&id=U00096.2&retmode=text&rettype=fasta" -O "genome.fa"
       "$ROOT/tools/fastatools/fastatools" rename -r "\.2.*genome" -s " " -i "genome.fa" > "$ROOT/data/whole_ecoli_genome.fa"
       rm "genome.fa"
    else
        echo "Reference E. coli genome already downloaded"
    fi

    if [[ ! -f "$ROOT/data/tandemtools.ref.fa" ]]; then
        cp "$ROOT/tools/TandemTools/test_data/simulated_del.fasta" "$ROOT/data/tandemtools.ref.fa"
    else
        echo "Already processed simulated centromeric sequence."
    fi

    if [[ ! -f "$ROOT/data/chm13.repeats.bed" ]]; then
        wget "https://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.1/rmsk/rmsk.bigBed"
        "$ROOT/bin/bigBedToBed" ./rmsk.bigBed "$ROOT/data/chm13.repeats.bed"
    else
        echo "Already processed RepeatMasker track"
    fi

    cd "$ROOT" && rm -rf "$ROOT/temp_data"
fi
