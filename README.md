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

- [ ] Init script / pipeline that:
  - [ ] Builds the different tools in submodules
  - [ ] copies / links executables to bin folder
  - [ ] Decompress pre-trained models

## Pipelines

- [ ] MSR selection pipeline (Write it, low coverage all MSRs)
- [ ] Function evaluation pipeline (selected MSRs):
  - [ ] Whole genome
  - [ ] Drosophila
  - [ ] Tandemtools simulated centromere

## Figure generation

- [ ] Plotting RMarkdown.
