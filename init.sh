#!/usr/bin/env bash
set -euo pipefail

ROOT=$(pwd)

# Making needed directories
mkdir -p data/models


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
if [[ "$OSTYPE" == "darwin"* ]]
then
    cd "$ROOT/tools/minimap2"
    make
    cd "$ROOT"
    ln -s "$ROOT/tools/minimap2/minimap2" "$ROOT/bin"
else
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
