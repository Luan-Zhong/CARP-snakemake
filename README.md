# CARP-snakemake

<img src="https://github.com/Luan-Zhong/CARP-snakemake/assets/37817639/55e9023b-45a3-487a-ac40-57cd4abb9ffa" width="400">



CARP (Comprehensive *ab initio* Repeat Pipeline) is a pipeline for repeat detection and annotation. This project aims to produce a reliable, scalable, and reproducible workflow as a series of chained rules for CARP.
Reducing the painful process like installing dependencies and managing the computational resources.

## Advantages
- One command for repeat analysis
- No painful steps in dependencies
- 

## Description


More information on CARP refers to:
- [Superior ab initio identification, annotation and characterisation of TEs and segmental duplications from genome assemblies](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0193588)
- [CARP manual](https://github.com/carp-te/carp-documentation/blob/master/CARP.pdf) (Updated in 2018-10-26)
- [CARP github page](https://github.com/carp-te/carp-documentation)

Singularity is employed as a container that suits both HPC and local. 
AB-BLAST replaces the WU-BLAST because it is the updated version of WU-BLAST.


## Pipeline overview
### Minimal input
- `fasta` file 

### Process and tools
- Pairwise alignment
  - [krishna](https://github.com/biogo/examples/tree/master/krishna)
  - [matrix](https://github.com/biogo/examples/tree/master/krishna/matrix)
  - [bundle](https://github.com/biogo/examples/tree/master/bundle): split genome sequence files into smaller files
- Report repeat feature family groupings in JSON format
  - [igor](https://github.com/biogo/examples/tree/master/igor)
- Generate consensus sequences from genome intervals
  - [seqer](https://github.com/biogo/examples/tree/master/igor/seqer)
- Annotate consensus sequences with repeat families
  - [censor V4.2.29](https://www.girinst.org/downloads/software/censor/)
  - [Repbase library 28.12 (vertebrate)](http://www.girinst.org/repbase/update/browse.php?type=All&format=FASTA&autonomous=on&nonautonomous=on&simple=on&division=Vertebrata&letter=A)
  - [AB-BLAST](https://blast.advbiocomp.com/)
- Classify consensus sequences
  - [ClassifyConsensusSequences.java](https://github.com/carp-te/carp-documentation/blob/master/code/ClassifyConsensusSequences.java)
- Filter sequences
  [AB-BLAST](https://blast.advbiocomp.com/)
  - Identify potential protein sequences
    - [uniprot protein dataset](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz)
    
  - Identify GB_TE sequences
  - Identify potential retrovirus sequences

  
  - [EDirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/#chapter6.Getting_Started)
  - [260118_GB_TE.fa](esearch -db protein -query "reverse transcriptase or transposon or repetitive element or RNA-directed DNA polymerase or pol protein or non-LTR retrotransposon or mobile element or retroelement or polyprotein or retrovirus or (group-specific antigen gag) or polymerase (pol)" | efetch -format fasta > 260118_GB_TE.fa)
  - [Retrovirus datasets](https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=11632)


## Quick start
### Set up conda environment
```bash
$ cd /your/Snakefile/location/
```
### Example dataset

### Running the pipeline
```bash
$ snakemake --use-conda --cores 8
```

### Server
```bash
# 2023-12-23 Complete Krishna block
$ cd /your/Snakefile/location/
$ snakemake --cores 8 output/1_krishna/scalefold17.gff
# with singularity
$ snakemake --cores 8 --use-singularity output/scalefold17.gff
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
# 2024-01-15 Finished
```
## Benchmarks

Data from CARP manual.
Analysis runs on a machine with 512GB RAM, running Red Hat Linux.

| Genome             | Krishna Threads | Genome DB Size | Krishna run time (hh:mm) | Igor run time (hh:mm) | Seqer run time (hh:mm) | Consensus sequences size | Censor first run time (hh:mm) | reportJ.sh (hh:mm) | phobos run time (hh:mm) |
|--------------------|-----------------|----------------|--------------------------|-----------------------|------------------------|--------------------------|-------------------------------|---------------------|-------------------------|
| Human              | 8               | 3.0G           | ~200                     | 128:30                | 2:23                   | 38M                      | 5:14                          | 19:30               | <00:10                  |
| Bearded Dragon     | 8               | 1.8G           | ~23                      | 73:11                 | <4                     | 88M                      | 22:21                         | <178                | <00:10                  |
| New Bearded Dragon | -               | -              | -                        | -                     | -                      | 88M                      | 7:13                          | 18:28               | <00:10                  |
| Anolis             | 6               | 1.8G           | 76:52                    | 97:32                 | 2:40                   | 63M                      | 9:42                          | 78                  | <00:10                  |
| Chicken            | 4               | 1017M          | 5                        | <4                    | <1                     | 18M                      | 3:01                          | <24                 | <00:10                  |
| Opossum            | 8               | 3.5G           | ~83                      | 61:48                 | 4:52                   | 60M                      | 13:17                         | 80                  | <01:00                  |
| Platypus           | 8               | 2.0G           | ~99                      | 191:34                | 10:16                  | 162M                     | 17:34                         | 115                 | <01:00                  |
| Echidna            | 8               | 1.9G           | 360                      | 12:43                 | 9:02                   | 59M                      | 18:40                         | 105                 | 01:54                   |



## Citation
If you using CARP snakemake pipeline, please cite:
> **Superior *ab initio* identification, annotation and characterisation of TEs and segmental duplications from genome assemblies**<br/>
> Lu Zeng, R. Daniel Kortschak, Joy M. Raison, Terry Bertozzi, David L. Adelson. *PLOS ONE*. 2018. <br/>
> [doi:https://doi.org/10.1371/journal.pone.0193588](https://doi.org/10.1371/journal.pone.0193588) [[pdf](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0193588&type=printable)]
> [[github page](https://github.com/carp-te/carp-documentation)]
> [[manual (Updated in 2018-10-26)](https://github.com/carp-te/carp-documentation/blob/master/CARP.pdf)]
