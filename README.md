<img src="https://github.com/Luan-Zhong/CARP-snakemake/assets/37817639/55e9023b-45a3-487a-ac40-57cd4abb9ffa" width="200">

# CARP-snakemake

## Overview
- What is it and what does it do?

  CARP (Comprehensive *ab initio* Repeat Pipeline) is a pipeline for repeat detection and annotation. This project aims to produce a reliable, scalable, and reproducible workflow as a series of chained rules for CARP

- Accepted input:
  
  Genome, cotigs in fasta format

- Expected output:

  Identify repeat consensus sequences from an input genome and classify them into repeat families

### What's new?
| Feature | Previous | CARP-snakemake |
| --- | --- | ---|
| Execution | Manual shell scripts | snakemake managed workflow |
| Cofiguration | Parameters in scripts | Central yaml configuration |
| Dependencies | Manual installation | Appatainer container |
| Recovery | Manual rerun | Automatic reuse of completed outputs; only missing or failed steps are rerun |
| HPC support | Manual adaption for cluster environment | HPC slurm support without requiring root privileges |
| Repeat identification | `censor` requiring `AB-BLAST` or `WU-BLAST` | `INS-fast`, an update version of our in-house `INS` software |
| Coverage summary | No integrated | Summary table of repeat-element coverage based on the RepBase classification system |
| Visualization | No integrated | Optional visualization script for identified repeat consensus annotations |
| Report | Static output files | Optional interactive html reports for easier result exploration |


<p align="center">
  <img
    width="1494"
    height="787"
    alt="CARP_snkm_report_screenshot"
    src="https://github.com/user-attachments/assets/6fa06a56-bffc-4215-bd4a-e78d400a724d"
  >
</p>

<p align="left">
  <strong>Figure 1. Example of CARP-snakemake produced report </strong>
</p>


<p align="center">
  <img
    width="3360"
    height="1392"
    alt="consensus_annotation_map"
    src="https://github.com/user-attachments/assets/58a5c78d-6245-419d-9973-e0d5b132797f"
  >
</p>

<p align="left">
  <strong>Figure 2. Example of consensus annotation map </strong>
</p>

## Workflow overview
<p align="center">
  <img
    width="316"
    height="738"
    alt="Pipeline workflow"
    src="https://github.com/user-attachments/assets/d3c1c5ef-053b-46af-b77d-c8695ce96408"
  >
</p>

<p align="left">
  <strong>Figure 3. Overview of the CARP-snakemake workflow showing the core rules corresponding to the steps in the original CARP pipeline</strong> Newly added optional visualization rules are not shown.
</p>




## Requirement

CARP-snakemake requires the following software on the host system:
- Snakemake 9.27.0
- Apptainer
- SLURM workload manager (optional, required only for execution on an HPC cluster)

All software dependencies required by the pipeline are included in the provided Apptainer image, `src/carp_open_source_v2.sif`. Root privileges are not required to run the prebuilt container.
Pipeline parameters and the input sample prefix must be specified in `config.yaml`.

## Installation
### Setup conda environment
Prerequesite: [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install/linux-install) and [apptainer](https://apptainer.org/docs/admin/main/installation.html) are installed
```bash
conda create -n snkm9 -c bioconda snakemake=9.27.0
conda activate snkm9
snakemake --version
```
NB: `--version` should show 9.27.0, if snakemake version is not shown, please double check your conda installation and snakemake download

### Clone CARP-snakemake
```bash
cd <whereever/you/want/to/put>
git clone https://github.com/Luan-Zhong/CARP-snakemake
cd CARP-snakemake
```
### Pull container under the `CARP-snakemake/src`

```bash
docker pull ghcr.io/luan-zhong/carp_snkm_sif_v2:latest
mv <myproject.sif> src
```

### Download databases
Make sure to move the downloaded databases to the lib directory and update the file names in the "# Reference libraries section" of `config.yaml`
- RepBase FASTA format</br>
  Download RepBase differnt releases from [Dfam Repbase Archive](https://dfam.org/repbase/releases).
- Swissprot library </br>
  Download uniprot protein data set from [UniProt download](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz)
  ```bash
  wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
  ```
- Download GB_TE dataset </br>
  First install EDirect from EDirect download page [Entrez Direct: E-utilities on the Unix Command Line](https://www.ncbi.nlm.nih.gov/books/NBK179288/#chapter6.Getting_Started). Then download the data:
  ```bash
  esearch -db protein -query "reverse transcriptase or transposon or repetitive element or RNA-directed DNA polymerase or pol protein or non-LTR retrotransposon or mobile gag) or polymerase (pol)" | efetch -format fasta > GB_TE.2022-07-14.fa
  ```
- Download retrovirus datasets </br>
  Download retrovirus datasets from your browser: [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Viruses,taxid:11632)

## Quick start
```bash
# Dry run
conda activate snkm9
snakemake --snakefile Snakefile --use-singularity --printshellcmds --configfile config.yaml --cores 4 -n
```

## Input data 

Copy input fasta to data folder

```bash
cp <your/data/absolute/path> data
```

Put your data under data folder, and modify `config.yaml` preFix according to your input data name
```yaml
# =========== Input genome name =================
# Sample prefixes matching FASTA files under data/
# Example
# If your input file is data/test_input.fa
# Your preFix will be

preFix:
  - test_input
```

## Configuration
The pipeline uses a `configure.yaml` file to specify input data, output directories, computational resources, and CARP parameters. Users can also enable or disable optional modules through this file.

The pipeline has two rounds of searches for `ins`. Users can control whether soft masking is enabled in the first round by setting `repeats_in_repbase_first_softmask: false
` or `annotate_genome_first_softmask: true` and  in `configure.yaml`.

All file paths and parameter values should be checked before running the pipeline.

## Usage
* Usage 1: Run the core rules corresponding to the steps in the original CARP pipeline
  ```bash
  snakemake --snakefile Snakefile --use-singularity --printshellcmds --configfile config.yaml --cores 4
  ```

* Usage 2: Run whole piepline with the genome annotation, and visualization steps
  - You can check available consensus name by
    ```bash
    awk '/^>/{sub(/^>/, ""); print $1}' output/test_input/04_finalize_consensus/test_input_ConsensusSequences.fa
    ```
  - Edit the `config.yaml` file. To generate PNG files only for selected consensus sequences, specify their consensus IDs in the configuration file
    ```yaml
    ####################
    # Optional outputs #
    ####################
    # ========= Enable rule annotate_genome =======================
    # Enable rule annotate_genome
    annotate_genome: true
    # Soft-mask matches before searching the Denovo TE library
    annotate_genome_first_softmask: false
    # ========= Enable rule plot_ins_map ==========================
    # plot_ins_map: false
    # Keep an empty list when plot_ins_map is false
    # plot_consensus_ids: []
    # If you want to generate annotation png files for selected consensu
    # First enable plot_ins_map
    plot_ins_map: true
    # Then, list the consensus you select
    # NB: Consensus must exactly match the first field of FASTA headers
    plot_consensus_ids:
      - family000001_consensus
      - family000007_consensus
      - family000025_consensus
    # ========= Enable rule make_report ===========================
    # Generate the final interactive html report
    make_report: true
    ```
  
  - Execute snakemake piepline
    ```bash
    snakemake --snakefile Snakefile --use-singularity --printshellcmds --configfile config.yaml --cores 4
    ```

* Usage 3: Run certain step (Only run krishna as example)
  ```bash
  snakemake --snakefile Snakefile --use-singularity --printshellcmds --configfile config.yaml --cores 4 output/test_input/02_igor/test_input_krishna.gff
  ```
* Usage 4: Rerun from incomplete
  ```bash
  snakemake --snakefile Snakefile --use-singularity --printshellcmds --configfile config.yaml --cores 4 --rerun-incomplete
  ``` 

## Output
Your result will be under `CARP-snakemake/output/<your_input_data_prefix>`

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


## Reproducibility

### Tools and softwares
| Step | Software/Scripts | Purpose |
| --- | --- | ---|
| 00_bundle | `Go v1.22.3` </br> [bundle](https://github.com/biogo/examples/tree/master/bundle)| Split genome sequence files into smaller files |
| 01_krishna | [krishna](https://github.com/biogo/examples/tree/master/krishna)  </br> [matrix](https://github.com/biogo/examples/tree/master/krishna/matrix) | Pairwise alignment |
| 02_igor |  [igor](https://github.com/biogo/examples/tree/master/igor) | Repeat family groupings into json file |
| 03_seqer | [seqer](https://github.com/biogo/examples/tree/master/igor/seqer) | Generate consensus sequences from genome intervals |
| 04_finalize_consensus | `bash cat` | Integrate consensus |
| 05_repeats_in_repbase | `ins-fast` which is motified from [ins](https://github.com/kortschak/ins) | Consensus annotation with RepBase and our known repeat | 
| 06_classify_consensus | `java Temurin JDK 8u504-b01` </br> ClassifyConsensusSequences_CARP_snkm.java | Generate annotation and grouping consensus into known, partial, and not known |
| 07_filter_sequences | ncbi-blast v2.15.0 `blastx` `tblastx` | Search not known consensus against the protien databases (Swiss-Prot, Genbank-TE protein database, retrovirus) | 
| 08_get_protein_info | GetProteins_CARP_snkm.java | Generate annotation for consensus based on the protein search result |
| 09_check_ssr | `phobos` v3.3.12 | Identify simple sequence repeat | 
| 10_identify_ssr | IdentifySSRs_CARP_snkm.java | Generate annotation for consensus based on the phobos search result |
| 11_generate_annotated_repeat_lib | GenerateAnnotatedLibrary_CARP_snkm.java | Integrate all the annotation result |
| (optional) 12_genome_annotation | `ins-fast` </br> ins_json_to_map.py </br> ins_map_to_tbl.py| Annotate whole genome and generate summary table |
| (optional) 13_plot_ins_map | `R 4.5.2` </br> `Gviz 1.54.0`plot_ins_map.R | Generate consensus annotation map | 
| (optional) 14_make_report | `python3`make_report.py | Generate an interactive html file |


### For HPC slurm
#### Setup conda environment in HPC home directory
Prerequesite: conda is installed
```bash
conda create -n snkm9 -c conda-forge -c bioconda --strict-channel-priority snakemake=9.27.0
conda activate snkm9
snakemake --version
```
NB: If snakemake version is not shown, please double check your conda installation and snakemake download

#### Clone repo
```bash
cd <whereever/you/want/to/put>
git clone https://github.com/Luan-Zhong/CARP-snakemake
cd CARP-snakemake
```
#### Pull container under the `CARP-snakemake/lib`

```bash
docker pull ghcr.io/luan-zhong/carp_snkm_sif_v2:latest
mv <myproject.sif> src
```

#### Copy your input fasta data under `CARP-snakemake/data`

```bash
cp <your/data/absolute/path> ./data
```
#### Setup `lib` directory
#### Edit `swimming_carp.sl` to match your hpc account and resource requirements
#### Submit your job
```bash
sbatch swimming_carp.sl
```

## Citation
If you using CARP snakemake pipeline, please cite:
> **Superior *ab initio* identification, annotation and characterisation of TEs and segmental duplications from genome assemblies**<br/>
> Lu Zeng, R. Daniel Kortschak, Joy M. Raison, Terry Bertozzi, David L. Adelson. *PLOS ONE*. 2018. <br/>
> [doi:https://doi.org/10.1371/journal.pone.0193588](https://doi.org/10.1371/journal.pone.0193588) [[pdf](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0193588&type=printable)]
> [[github page](https://github.com/carp-te/carp-documentation)]
> [[manual (Updated in 2018-10-26)](https://github.com/carp-te/carp-documentation/blob/master/CARP.pdf)]

Please also cite Phobos when using this software. This software is free for academic use only:
> **Genome-wide analysis of tandem repeats in Daphnia pulex - a comparative approach**<br/>
> Christoph Mayer, Florian Leese & Ralph Tollrian. *BMC Genomics*. 2010. <br/>
> [https://doi.org/10.1186/1471-2164-11-277](https://doi.org/10.1186/1471-2164-11-277)
> [[webpage](https://www.ruhr-uni-bochum.de/ecoevo/cm/cm_phobos.htm)]

