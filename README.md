# CARP-snakemake
<img src="https://github.com/Luan-Zhong/CARP-snakemake/assets/37817639/55e9023b-45a3-487a-ac40-57cd4abb9ffa" width="200">
<br/>

## Overview
- What is it and what does it do?

  CARP (Comprehensive *ab initio* Repeat Pipeline) is a pipeline for repeat detection and annotation. This project aims to produce a reliable, scalable, and reproducible workflow as a series of chained rules for CARP.

- Accepted input:
  
  Genome, cotigs in fasta format

- Expected output:

  Identify repeat consensus sequences from an input genome and classify them into repeat families

### What's new?
* One-step dependency installation with a apptainer container
* Support for resuming from checkpoints
* HPC/SLURM support with no root privileges required
* Replace `censor` with `INS`, which does not require `AB-BLAST` or `WU-BLAST`
* Add a summary table showing the coverage of different repeat elements based on the RepBase classification system
* Add optional integrate visualization script to visualized found repeat consensus annotation
* Add optional interactive html reports for easier viewing of results
  <img width="1494" height="787" alt="CARP_snkm_report_screenshot" src="https://github.com/user-attachments/assets/6fa06a56-bffc-4215-bd4a-e78d400a724d" />






## Citation
If you using CARP snakemake pipeline, please cite:
> **Superior *ab initio* identification, annotation and characterisation of TEs and segmental duplications from genome assemblies**<br/>
> Lu Zeng, R. Daniel Kortschak, Joy M. Raison, Terry Bertozzi, David L. Adelson. *PLOS ONE*. 2018. <br/>
> [doi:https://doi.org/10.1371/journal.pone.0193588](https://doi.org/10.1371/journal.pone.0193588) [[pdf](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0193588&type=printable)]
> [[github page](https://github.com/carp-te/carp-documentation)]
> [[manual (Updated in 2018-10-26)](https://github.com/carp-te/carp-documentation/blob/master/CARP.pdf)]

## Workflow overview

### Workflow structure
<img width="316" height="738" alt="carp_snkm" src="https://github.com/user-attachments/assets/d3c1c5ef-053b-46af-b77d-c8695ce96408" />

### Process and tools
| Step | Software | Purpose |
| --- | --- | ---|
| 00_bundle | [bundle](https://github.com/biogo/examples/tree/master/bundle)| split genome sequence files into smaller files |
| 01_krishna | [krishna](https://github.com/biogo/examples/tree/master/krishna)  </br> [matrix](https://github.com/biogo/examples/tree/master/krishna/matrix) | Pairwise alignment |
| 02_igor |  [igor](https://github.com/biogo/examples/tree/master/igor) | Repeat family groupings into json file |
| 03_seqer | [seqer](https://github.com/biogo/examples/tree/master/igor/seqer) | Generate consensus sequences from genome intervals |
| 04_finalize_consensus | `bash cat` | Integrate consensus |
| 05_repeats_in_repbase | `ins-fast` which is motified from [ins](https://github.com/kortschak/ins) | Consensus annotation with RepBase and our known repeat | 
| 06_classify_consensus | ClassifyConsensusSequences_CARP_snkm.java | Generate annotation and grouping consensus into known, partial, and not known |
| 07_filter_sequences | `blastx` `tblastx` | Search not known consensus against the protien databases (Swiss-Prot, Genbank-TE protein database, retrovirus) | 
| 08_get_protein_info | GetProteins_CARP_snkm.java | Generate annotation for consensus based on the protein search result |
| 09_check_ssr | `phobos` v3.3.12 | Identify simple sequence repeat | 
| 10_identify_ssr | IdentifySSRs_CARP_snkm.java | Generate annotation for consensus based on the phobos search result |
| 11_generate_annotated_repeat_lib | GenerateAnnotatedLibrary_CARP_snkm.java | Integrate all the annotation result |
| (optional) 12_genome_annotation | `ins-fast` </br> ins_json_to_map.py </br> ins_map_to_tbl.py| Annotate whole genome and generate summary table |
| (optional) 13_plot_ins_map | plot_ins_map.R | Generate consensus annotation map | 
| (optional) 14_make_report | make_report.py | Generate an interactive html file |


    - [uniprot protein dataset](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz)
  - [EDirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/#chapter6.Getting_Started)
  - [260118_GB_TE.fa](esearch -db protein -query "reverse transcriptase or transposon or repetitive element or RNA-directed DNA polymerase or pol protein or non-LTR retrotransposon or mobile element or retroelement or polyprotein or retrovirus or (group-specific antigen gag) or polymerase (pol)" | efetch -format fasta > 260118_GB_TE.fa)
  - [Retrovirus datasets](https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=11632)

## Installation
### For Linux server
#### Setup conda environment
Prerequesite: [conda](https://www.anaconda.com/docs/getting-started/miniconda/install/linux-install) and [apptainer](https://apptainer.org/docs/admin/main/installation.html) are installed
```bash
$ conda create -n snkm9 -c bioconda snakemake=9.27.0
$ conda activate snkm9
$ snakemake --version
```
NB: `--version` should show 9.27.0, if snakemake version is not shown, please double check your conda installation and snakemake download

#### Clone repo
```bash
$ cd <whereever/you/want/to/put>
$ git clone https://github.com/Luan-Zhong/CARP-snakemake
$ cd CARP-snakemake
```
#### Pull container under the `CARP-snakemake/lib`

```bash
$ apptainer pull <myproject.sif> <docker://ghcr.io/USERNAME/myproject:1.0>
$ mv <myproject.sif> lib
```

#### Copy your input fasta data under `CARP-snakemake/data`

```bash
$ cp <your/data/absolute/path> data
# Make sure data successfully copies and in right size
$ ls -lah data
```

#### Dry run whole pipeline to make sure everything are prepared
```bash
$ snakemake --snakefile Snakefile --use-singularity --printshellcmds --cores 4 -n --config preFix=<your_input_data_prefix>
```

```bash
$ snakemake --snakefile Snakefile --use-singularity --printshellcmds --cores 4 --config preFix=<your_input_data_prefix> -- output/test_input/02_igor/test_input_krishna.gff
```
NB: If your input data is `test_input.fa`, `preFix=<your_input_data_prefix>` will be `preFix=test_input`

#### Actual run
```bash
$ snakemake --snakefile Snakefile --use-singularity --cores 4 --printshellcmds --config preFix=<your_input_data_prefix> --report output/<<your_input_data_prefix>_snkm_statistics_report.html
```

#### Where is the result?
Your result will be under `CARP-snakemake/output/<your_input_data_prefix>`

#### Example: run with example dataset
Dry run
```bash
$ cd CARP-snakemake
$ cp <my/input_fasta/location>/test_input.fa data
$ snakemake --snakefile Snakefile --use-singularity --printshellcmds --cores 4 -n --config preFix=test_input 
```
You will expect
```bash

```
Actual run

```bash
$ snakemake --snakefile Snakefile --use-singularity --cores 4 --printshellcmds --config preFix=test_input --report output/<your_input_data_prefix>/snkm_statistics_report.html
```
After process is completed, you will expect this output structure
```bash
.
└── test_input
    ├── 00_bundle
    │   ├── test_input.fa-0.fa
    │   └── test_input.fa-1.fa
    ├── 01_krishna
    │   ├── krishna-20260922153925-Tue-345216.log
    │   ├── krishna-20260922153925-Tue-345217.log
    │   ├── krishna-20260922153925-Tue-345219.log
    │   ├── test_input.fa-0.gff
    │   ├── test_input.fa-0_test_input.fa-1.gff
    │   └── test_input.fa-1.gff
    ├── 02_igor
    │   ├── test_input_krishna.gff
    │   └── test_input_krishna.json
    ├── 03_seqer
    │   ├── consensus
    │   │   ├── family000000_consensus.fq
    │   │   ├── family000000.mfa
    │   │   ├── family000001_consensus.fq
    │   │   ├── family000001.mfa
    │   │   ├── family......_consensus.fq
    │   │   └── family.......mfa
    │   ├── test_input_krishna.igor.gff
    │   ├── test_input.mfa
    │   └── test_input_temp_ConsensusSequences.fa
    ├── 04_finalize_consensus
    │   ├── test_input_ConsensusSequences.fa
    │   └── test_input_ConsensusSequences.fa.idx
    ├── 05_repeats_in_repbase
    │   ├── censor.350970.log
    │   ├── censor.351400.log
    │   ├── censor.352206.log
    │   ├── censor.355258.log
    │   ├── test_input_ConsensusSequences.fa.aln
    │   ├── test_input_ConsensusSequences.fa.found
    │   ├── test_input_ConsensusSequences.fa.map
    │   └── test_input_ConsensusSequences.fa.masked
    ├── 06_classify_consensus
    │   ├── test_input_check.txt
    │   ├── test_input_ClassifyConsensusSequences.log
    │   ├── test_input_known.txt
    │   ├── test_input_notKnown.fa
    │   ├── test_input_notKnown.fa.gff
    │   └── test_input_partial.txt
    ├── 07_filter_sequences
    │   ├── test_input_notKnown.fa.ervwb
    │   ├── test_input_notKnown.fa.ervwb.gff
    │   ├── test_input_notKnown.fa.spwb
    │   ├── test_input_notKnown.fa.spwb.gff
    │   ├── test_input_notKnown.fa.tewb
    │   └── test_input_notKnown.fa.tewb.gff
    ├── 08_get_protein_info
    │   ├── test_input_GetProteins.log
    │   ├── test_input_notKnownNotProtein.fa
    │   └── test_input_protein.txt
    ├── 09_check_ssr
    │   └── test_input_notKnownNotProtein.phobos
    ├── 10_identify_ssr
    │   ├── test_input_GB_TE.fa
    │   ├── test_input_IdentifySSRs.log
    │   ├── test_input_notKnownNotProteinNotSSR.fa
    │   └── test_input_SSR.txt
    └── 11_generate_annotated_repeat_lib
        ├── test_input_Denovo_TE_Library.fasta
        └── test_input_GenerateAnnotatedLibrary.log
```

### For HPC slurm
#### Setup conda environment in HPC home directory
Prerequesite: conda is installed
```bash
$ conda create -n snkm9 -c bioconda snakemake=9.27.0
$ conda activate snkm9
$ snakemake --version
```
NB: If snakemake version is not shown, please double check your conda installation and snakemake download

#### Clone repo
```bash
$ cd <whereever/you/want/to/put>
$ git clone https://github.com/Luan-Zhong/CARP-snakemake
$ cd CARP-snakemake
```
#### Pull container under the `CARP-snakemake/lib`

```bash
$ apptainer pull <myproject.sif> <docker://ghcr.io/USERNAME/myproject:1.0>
$ mv <myproject.sif> lib
```

#### Copy your input fasta data under `CARP-snakemake/data`

```bash
$ cp <your/data/absolute/path> ./data
# Make sure data successfully copies and in right size
$ ls -lah ./data
```
#### Edit your slurm script
```
  


## Time and Benchmarks

Depends on the genome size and repetitive, large genome and repeatitive genome (e.g.:T2T) will take longer time.

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

