# CARP-snakemake
![](https://img.shields.io/badge/BASH-282F37) ![](https://img.shields.io/badge/snakemake-038F71) ![](https://img.shields.io/badge/Singularity-40B12E)  ![](https://img.shields.io/badge/docker-0291E7)
## Description
CARP (Comprehensive ab initio Repeat Pipeline) is a pipeline for repeat detection and annotation. This project aims to produce a reliable, scalable, and reproducible workflow as a series of chained rules for CARP.
Reducing the painful process like installing dependencies and managing the computational resources.

More information on CARP refers to:
- [Superior ab initio identification, annotation and characterisation of TEs and segmental duplications from genome assemblies](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0193588)
- [CARP manual](https://github.com/carp-te/carp-documentation/blob/master/CARP.pdf) (Updated in 2018-10-26)

Singularity is employed as a container that suits both HPC and local. 
AB-BLAST replaces the WU-BLAST because it is the updated version of WU-BLAST.

## Contributing


## Instructions
### Server
```bash
# 2023-12-23 Complete Krishna block
$ cd /your/Snakefile/location/
$ snakemake --cores 8 output/1_krishna/scalefold17.gff
```

```bash
# 2023-12-24 Complete igor block
$ cd /your/Snakefile/location/
$ snakemake --cores 8 output/2_igor/scalefold17_krishna.json
```

```bash
# 2024-01-04 gffer code block
$ cd /your/Snakefile/location/
$ snakemake --cores 8 output/3_seqer/scalefold17_krishna.igor.gff
```

```bash
# 2024-01-05 repeat_annotation block
$ snakemake --cores 8 output/4_rep_anno/scalefold17_ConsensusSequences.fa
# Update the existing code block output files directory
```

```bash
# 2024-01-05 repeat_annotation block
$ snakemake --cores 8 output/4_rep_anno/scalefold17_ConsensusSequences.fa
# 2024-01-15 Finished
```


## Developing log
