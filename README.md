# Mapping-friendly Sequence Reduction function discovery

## Introduction

This is the repository linked to the ["Mapping-friendly sequence reductions: going beyond homopolymer compression"]() paper.  
It contains the necessary pipelines, tools and information in order to rerun the analyses performed during the project as well as explore this subject further.

## Quick start

This is made to work on Linux or MacOS systems. You must first clone this repository with the `--recursive` flag in order to also get the tools included as submodules. You can then run the init script that builds the tools that need to be compiled and makes sure the various data files and pre-trained models are in the right place, the `--data` flag ensure the reference genomes are downloaded as well. You should then be able to run the pipelines.

```shell
git clone --recursive git@github.com:lucblassel/MSR_discovery.git
cd MSR_discovery
./init.sh --data
nextflow run <pipeline-name>
```

## Dependencies

You will need the following dependencies available on your system in order to run the pipelines:

- gcc >= 10 to compile [winnowmap](https://github.com/marbl/Winnowmap) on MacOs. You can install it with HomeBrew on MacOs: `brew install gcc@10`
- zlib development file to compile [minimap](https://github.com/lh3/minimap2)
- [k8 javascript shell](https://github.com/attractivechaos/k8) to run our fork of [paftools.js](https://github.com/lh3/minimap2/blob/master/misc/README.md)
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) in order to execute the pipelines

## Pipelines

All the pipelines are implemented in [NextFlow](https://www.nextflow.io). You can run them by using: `nextflow run <pipeline-name>`.

### MSR Selection

This pipeline generate all the unique MSRs _(as described in the paper)_ and evaluated each of them on a low coverage human genome simulated dataset.  
The 10 best MSRs w.r.t error rate and 10 best w.r.t number of reads mapped are selected for further analysis and are used as input for further pipelines.

### MSR In depth evaluation

The previously selected MSRs are evaluated on a wider range of use cases.  
For a given reference, a set of reads is simulated with a coverage of approximately 1.5. The reads and reference are transformed with each MSR, then the transformed reads are mapped to transformed reference. The resulting mapping is then evaluated.  
This process is done for each possible combination of:

- Reference:
  - T2T CHM13 v1.1 whole human genome reference
  - TandemTools simulated centromeric reference
  - Whole Drosophila _melanogaster_ reference
- Simulator:
  - NanoSim with R94 model
  - PBSim with P6C4 model
- Mapper:
  - minimap2
  - winnowmap

Each possible mapping is evaluated using the `mapeval` command in our fork of `paftools.js`.  
Additionaly, a mapping is evaluated on a subset of reads that are mapped to repeated regions of the genome.
A single `.csv` file with all the evaluation is produced and can be used to generate plots and tables.

### Plotting and Tables

This has been implemented as an Rmarkdown notebook.

## Included tools

We have made the choice of including all the tools used in this project as git submodules when possible.

- Read simulation
  - Nanosim
  - PbSim2
- Mapping
  - Minimap2
  - Winnowmap
- Read manipulation:
  - fastatools
  - lucblassel/rename_sequences
  - lucblassel/reduce_sequences
- result file manipulation:
  - our custom fork of paftools
  - bedtools
  - bigBedToBed

The [init.sh](./init.sh) script will set up the tools before you run the pipeline, it currently only supports amd64 architectures on Linux and MacOS. It will build:

- Winnowmap on MacOS/Linux
- PBSim on MacOS/Linux
- Minimap2 on MacOS
- Bedtools on MacOS

It will download and setup prebuilt binaries for:

- Minimap2 on Linux
- Bedtools on Linux
- Go reduce_sequences and rename_sequences on MacOS/Linux

If you specify the `--data` flag it will also download the reference datasets used in the analysis to the correct directories, and pre-process it for the pipelines to work.

# Tools used in the HPC project

## Simulation tools

- [x] Pbsim: git checkout _(submodule)_
- [x] Nanosim: git checkout _(submodule)_

## Mappers

- [x] Minimap: version + git checkout _(include release folder)_
- [x] Winnowmap: version (git checkout ?) + meryl _(submodule)_
- [x] TandemTools: git checkout (for reference sequence)

## Evaluation tools

- [x] paftoolsCustom.js: include in subdir (forked from minimap 2.18...) depends on [k8 javascript shell](https://github.com/attractivechaos/k8)
- [ ] bedtools intersect

## Pipelining

- [x] Nextflow _(mention as dependency, don't include in repo)_

## Custom tools

### RF generation

- [x] makeAllUniqueFunction.py _(in bin/)_

### Sequence manipulation

- [ ] Library: reduction-functions github (keep in README)
- [x] Reducer: `reduceDelete_linux_amd64` _(submodule)_
- [x] Renamer: `applyOffsets_linux_amd64` _(submodule)_
- [x] fastatools for renaming references and other fasta manipulation

## Plotting and results

- [ ] mapeval gathering
- [ ] mapevalGrapher
- [ ] Rmarkdown

# Data used

### Sequence references

- [x] [CHM13 T2T v1.1](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13.draft_v1.1.fasta.gz)
- [x] TandemTools reference (git checkout)
- [ ] drosophila reference
- [ ] Add option to DL to init script

### Simulation models

- [x] nanosim model (git checkout)
- [x] pbsim model (git checkout)

### repeats

- [x] [repeat masker](https://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.1/rmsk/rmsk.bigBed) _(CHM13 v1.1)_

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

- [x] MSR selection pipeline (Write it, low coverage all MSRs)
- [x] Function evaluation pipeline (selected MSRs):
  - [x] Whole genome
  - [x] Drosophila
  - [x] Tandemtools simulated centromere

## Figure generation

- [x] Plotting RMarkdown.
