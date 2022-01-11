# Tools used in the HPC project

## Simulation tools

- [x] Pbsim: git checkout _(submodule)_
- [x] Nanosim: git checkout _(submodule)_

## Mappers

- [x] Minimap: version + git checkout _(include release folder)_
- [x] Winnowmap: version (git checkout ?) + meryl _(submodule)_
- [ ] TandemTools: git checkout (for reference sequence)

## Evaluation tools

- [ ] paftoolsCustom.js: include in subdir (forked from minimap...)
- [ ] bedtools intersect

## Pipelining

- [ ] Nextflow _(mention as dependency, don't include in repo)_

## Custom tools

### RF generation

- [x] makeAllUniqueFunction.py _(in bin/)_

### Sequence manipulation

- [ ] Library: reduction-functions github (keep in README)
- [x] Reducer: reduceDelete*linux_amd64 *(submodule)\_
- [x] Renamer: applyOffsets*linux_amd64 *(submodule)\_

## Plotting and results

- [ ] mapeval gathering
- [ ] mapevalGrapher
- [ ] Rmarkdown

# Data used

### Sequence references

- [ ] CHM13 T2T version..
- [ ] TandemTools reference (git checkout)
- [ ] drosophila reference

### Simulation models

- [ ] nanosim model (git checkout)
- [ ] pbsim model (git checkout)

### repeats

- [ ] repeat masker

### real datasets

- [ ] HG002 Nanopore
- [ ] HG002 PacBio

# Pipelines & script

## Init script

- [x] Init script / pipeline that:
  - [x] Builds the different tools in submodules
  - [x] copies / links executables to bin folder
  - [x] Decompress pre-trained models

This script works for MacOS and Linux in theory, with amd64 CPU architectures.
It will build:

- Winnowmap on MacOS/Linux
- PBSim on MacOS/Linux
- Minimap2 on MacOS
  It will download and setup prebuilt binaries for:
- Minimap2 on Linux
- Go reduce_sequences and rename_sequences on MacOS/Linux

Dependencies:

- You will [need zlib development libraries](https://github.com/lh3/minimap2#installation) to build minimap2 on MacOs.
- You will need gcc 10 or 11 to build winnowmap on macos, you can install it with homebrew: `brew install gcc@10` or `gcc@11`.

### description

## Pipelines

- [ ] MSR selection pipeline (Write it, low coverage all MSRs)
- [ ] Function evaluation pipeline (selected MSRs):
  - [ ] Whole genome
  - [ ] Drosophila
  - [ ] Tandemtools simulated centromere

## Figure generation

- [ ] Plotting RMarkdown.
