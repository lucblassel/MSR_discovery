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
A sample config file [`nextflow.config`](./nextflow.config) is available for you to fill out and adapt to your HPC environment. 

### SSR evaluation
```shell
nextflow run msr_selection.nf -resume
```
This pipeline executes the following steps: 
  - generate all SSRs as described in the paper and saves them to the `data/SSRs` directory
  - evaluate each SSR using simulated reads and saves the resulting mapping `paf` file and `mapeval` ouput file to the `results/SSR_eval/[SSR]` directory
  - generate a gathered `csv` file of `mapeval` outputs for all evaluated SSRs in the `results/SSR_eval/evaluations.csv` file. 

Using this file you can select the top MSRs using the following command: 
```shell
bin/selectMSRs.py \
  --csv results/SSR_eval/evaluations.csv \
  --top 20 \
  --dir data/SSRs/
```
This will perform the selection method described in the paper and a subdirectory in the `data/SSRs` directory called `MSRs` with symbolic links to the relevant SSR `json` files. It will also create a text file which lists the best MSR in each selection category.  
If the `--dir` flag is not specified this script will create a `json` file in the current directory listing the selected MSRs as well as the best in each category. 
### MSR In depth evaluation [TODO]

The previously selected MSRs in `data/SSRs/MSRs` are evaluated on a wider range of use cases.  
For a given reference, a set of reads is simulated with a coverage of approximately 1.5. The reads and reference are transformed with each MSR, then the transformed reads are mapped to transformed reference. The resulting mapping is then evaluated.  
This process is done for each possible combination of:

- Reference:
  - T2T CHM13 v1.1 whole human genome reference
  - TandemTools simulated centromeric reference
  - Whole Drosophila _melanogaster_ reference
  - Whole Escherischia _coli_ reference
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

All the necessary result data as well as an R-markdown file are available in the `figure_generation` directory. With this you should be able to reproduce all the figures and tables used in the article (both in the main text and supplementary material).

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

If you already have binaries for the aforementioned tools, you can place them in the `bin` directory and the [`init.sh`](./init.sh) script will skip those.
