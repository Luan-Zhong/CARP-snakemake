# SnakeCARP
# This is the snakemake prpeline for CARP
# 2026-09-20
# Version: 3.0
# snakemake version: conda version 9.27.0

# a) Use singularity container 
# b) Now, Outputs all in "output/" directory, see global variables OUTDIR

# Understand snakefile

# a) Use previous block result in next block
# e.g.: use igor json out file as seqer input
# rules.igor.output.json


# Load workflow settings from the project configuration file
configfile: "config.yaml"

####################
# Global variables #
####################

#========= Set input genome (control by config.yaml) =================
# Read sample prefixes from config.yaml
# Example: test_input corresponds to data/test_input.fa
preFix = config.get("preFix")
if preFix is None:
    raise ValueError("[ERROR] preFix is required in config.yaml")
# One sample name or a list of sample names
if isinstance(preFix, str):
    preFix = [preFix]
# Not putting any sample name
if (
    not isinstance(preFix, list)
    or not preFix
    or any(
        not isinstance(prefix, str) or not prefix
        for prefix in preFix
    )
):
    raise ValueError(
        "[ERROR] preFix must contain one or more sample names"
    )


#========= Set input directory (control by config.yaml) =================
INPUT_DIR = config.get("input_dir")

if not isinstance(INPUT_DIR, str) or not INPUT_DIR:
    raise ValueError(
        "[ERROR] input_dir must be a non-empty path in config.yaml"
    )
# Add exactly one trailing slash for later path construction
INPUT_DIR = INPUT_DIR.rstrip("/") + "/"


#========= Set output directory (control by config.yaml) =================
OUTPUT_DIR = config.get("output_dir")
if not isinstance(OUTPUT_DIR, str) or not OUTPUT_DIR:
    raise ValueError(
        "[ERROR] output_dir must be a non-empty path in config.yaml"
    )
if not OUTPUT_DIR.rstrip("/"):
    raise ValueError(
        "[ERROR] output_dir cannot be the filesystem root"
    )
# OUTDIR retains the name used throughout the existing Snakefile
OUTDIR = OUTPUT_DIR.rstrip("/") + "/"


#========= Set threads (control by config.yaml) =================
MAX_THREADS = config.get("max_threads")

if (
    isinstance(MAX_THREADS, bool)
    or not isinstance(MAX_THREADS, int)
    or MAX_THREADS < 1
):
    raise ValueError(
        "[ERROR] max_threads must be an integer greater than zero"
    )

#========= Set maximum genome bundle size (control by config.yaml) =================
BUNDLE_SIZE = config.get("bundle_size")
if (
    isinstance(BUNDLE_SIZE, bool)
    or not isinstance(BUNDLE_SIZE, int)
    or BUNDLE_SIZE < 1
):
    raise ValueError(
        "[ERROR] bundle_size must be an integer greater than zero"
    )

#========= Set Krishna parameters (control by config.yaml) =================
# Minimum alignment identity expressed as a proportion
# Example: 0.94 requires 94% identity
KRISHNA_FILTID = config.get("krishna_filter_identity")

if (
    isinstance(KRISHNA_FILTID, bool)
    or not isinstance(KRISHNA_FILTID, (int, float))
    or not 0 < KRISHNA_FILTID <= 1
):
    raise ValueError(
        "[ERROR] krishna_filter_identity must be greater than 0 and no greater than 1"
    )
# Minimum alignment length expressed in bases
# Example: 400 requires an alignment length of at least 400 bases
KRISHNA_FILTLEN = config.get("krishna_filter_length")

if (
    isinstance(KRISHNA_FILTLEN, bool)
    or not isinstance(KRISHNA_FILTLEN, int)
    or KRISHNA_FILTLEN < 1
):
    raise ValueError(
        "[ERROR] krishna_filter_length must be an integer greater than zero"
    )


#========= Set Seqer parameters (control by config.yaml) =================
# Maximum number of members used for one repeat family
SEQER_MAX_FAMILY_SIZE = config.get("seqer_max_family_size")

if (
    isinstance(SEQER_MAX_FAMILY_SIZE, bool)
    or not isinstance(SEQER_MAX_FAMILY_SIZE, int)
    or SEQER_MAX_FAMILY_SIZE < 1
):
    raise ValueError(
        "[ERROR] seqer_max_family_size must be an integer greater than zero"
    )


# Minimum retained sequence length expressed as a proportion
# Example: 0.95 requires at least 95 percent of the expected length
SEQER_MIN_LENGTH = config.get("seqer_min_length")

if (
    isinstance(SEQER_MIN_LENGTH, bool)
    or not isinstance(SEQER_MIN_LENGTH, (int, float))
    or not 0 < SEQER_MIN_LENGTH <= 1
):
    raise ValueError(
        "[ERROR] seqer_min_length must be greater than 0 and no greater than 1"
    )

#========= Set INS parameters (control by config.yaml) =================

# Control RepBase soft masking before searching the our-known library
REPEATS_IN_REPBASE_FIRST_SOFTMASK = config.get(
    "repeats_in_repbase_first_softmask"
)

# Reject strings such as "true" because the YAML value must be boolean
if not isinstance(REPEATS_IN_REPBASE_FIRST_SOFTMASK, bool):
    raise ValueError(
        "[ERROR] repeats_in_repbase_first_softmask must be true or false in config.yaml"
    )


# Control RepBase soft masking before searching the Denovo TE library
ANNOTATE_GENOME_FIRST_SOFTMASK = config.get(
    "annotate_genome_first_softmask"
)

# Reject strings such as "true" because the YAML value must be boolean
if not isinstance(ANNOTATE_GENOME_FIRST_SOFTMASK, bool):
    raise ValueError(
        "[ERROR] annotate_genome_first_softmask must be true or false in config.yaml"
    )

#============ Reference libraries ==================================

# Return one required library path from config.yaml
# config_key:
#   Name of the YAML setting that contains the file path
#   Example: "repbase_library" reads config["repbase_library"]
def read_library_path(config_key):
    library_path = config.get(config_key)

    if not isinstance(library_path, str) or not library_path:
        raise ValueError(
            f"[ERROR] {config_key} must be a non-empty path in config.yaml"
        )

    return library_path
# Read reference-library paths from config.yaml
REPBASE_LIBRARY = read_library_path("repbase_library")
OUR_KNOWN_LIBRARY = read_library_path("our_known_library")
SWISSPROT_LIBRARY = read_library_path("swissprot_library")
GB_TE_LIBRARY = read_library_path("gb_te_library")
RETROVIRUS_LIBRARY = read_library_path("retrovirus_library")

#========= Enable/disable consensus plot (control by config.yaml) =================
# Control whether stage 12 INS genome annotation is requested
ANNOTATE_GENOME = config.get("annotate_genome")
if not isinstance(ANNOTATE_GENOME, bool):
    raise ValueError(
        "[ERROR] annotate_genome must be true or false in config.yaml"
    )

# Control whether consensus annotation plots are requested
PLOT_INS_MAP = config.get("plot_ins_map")
if not isinstance(PLOT_INS_MAP, bool):
    raise ValueError(
        "[ERROR] plot_ins_map must be true or false in config.yaml"
    )


# Read the consensus IDs selected for plotting
PLOT_CONSENSUS_IDS = config.get("plot_consensus_ids")

# Require a YAML list even when plotting is disabled
if not isinstance(PLOT_CONSENSUS_IDS, list):
    raise ValueError(
        "[ERROR] plot_consensus_ids must be a YAML list"
    )

# Require at least one selected ID when plotting is enabled
if PLOT_INS_MAP and not PLOT_CONSENSUS_IDS:
    raise ValueError(
        "[ERROR] plot_consensus_ids must contain at least one ID when plot_ins_map is true"
    )

# Reject empty values, non-string values, and commas used as separators
if any(
    not isinstance(consensus_id, str)
    or not consensus_id
    or "," in consensus_id
    for consensus_id in PLOT_CONSENSUS_IDS
):
    raise ValueError(
        "[ERROR] plot_consensus_ids must contain non-empty IDs without commas"
    )

# Reject duplicated IDs because each ID should produce one PNG
if len(set(PLOT_CONSENSUS_IDS)) != len(PLOT_CONSENSUS_IDS):
    raise ValueError(
        "[ERROR] plot_consensus_ids must not contain duplicated IDs"
    )

# Convert the YAML list into one comma-separated command-line argument
# Example: ["family1", "family2"] becomes "family1,family2"
PLOT_CONSENSUS_IDS_ARGUMENT = ",".join(PLOT_CONSENSUS_IDS)

# Control whether the final HTML report is requested
MAKE_REPORT = config.get("make_report")
if not isinstance(MAKE_REPORT, bool):
    raise ValueError(
        "[ERROR] make_report must be true or false in config.yaml"
    )


def optional_genome_annotation(wildcards):
    # Return no report dependency when genome annotation is disabled
    if not ANNOTATE_GENOME:
        return []

    # Make the report wait for the stage 12 summary table
    return (
        OUTDIR
        + f"{wildcards.prefix}/12_genome_annotation/{wildcards.prefix}.mfa.tbl"
    )


def optional_consensus_plots(wildcards):
    # Return no report dependency when consensus plotting is disabled
    if not PLOT_INS_MAP:
        return []

    # Make the report wait for the consensus plot directory
    return (OUTDIR + f"{wildcards.prefix}/13_consensus_plots")


# Always request the core pipeline output from stage 11
FINAL_OUTPUTS = expand(
    OUTDIR
    + "{prefix}/11_generate_annotated_repeat_lib/{prefix}_Denovo_TE_Library.fasta",
    prefix=preFix,
)

# Add the stage 12 summary table only when genome annotation is enabled
if ANNOTATE_GENOME:
    FINAL_OUTPUTS.extend(
        expand(
            OUTDIR
            + "{prefix}/12_genome_annotation/{prefix}.mfa.tbl",
            prefix=preFix,
        )
    )

# Add the stage 13 consensus plot directory only when plotting is enabled
if PLOT_INS_MAP:
    FINAL_OUTPUTS.extend(
        expand(
            OUTDIR
            + "{prefix}/13_consensus_plots",
            prefix=preFix,
        )
    )

# Add the HTML report only when report generation is enabled
if MAKE_REPORT:
    FINAL_OUTPUTS.extend(
        expand(
            OUTDIR + "{prefix}/{prefix}_CARP_snkm_report.html",
            prefix=preFix,
        )
    )

################
# Main content #
################

onstart:
    print("#################################\n### Starting CARP... ><(((('> ###\n#################################\n")
    shell("mkdir -p {OUTDIR}")

# Use singularity local contianer
# singularity exec src/carp_open_source_v2.sif
singularity: "src/carp_open_source_v2.sif"


# Set output file inside expand function to run seleted step
rule all:
    input:
        FINAL_OUTPUTS



################
# Rules Proper #
################
# ================== Genome split by bundle =============================
checkpoint bundle:
    input:
        genome = INPUT_DIR + "{prefix}.fa"
    output:
        directory(OUTDIR + "{prefix}/00_bundle")
    params:
        out_dir = OUTDIR + "{prefix}/00_bundle",
        bundle_size = BUNDLE_SIZE
    shell:
        """
        echo "======= Running bundle... ======="

        mkdir -p {params.out_dir}

        genome=$(realpath {input.genome})
        cd {params.out_dir}

        bundle -bundle {params.bundle_size} -in "$genome"

        echo "======= Bundle completed! Output is in {params.out_dir} ======="
        """
        
def bundled_genomes(wildcards):
    bundle_dir = checkpoints.bundle.get(prefix=wildcards.prefix).output[0]
    files = sorted(glob.glob(str(bundle_dir) + "/*.fa"))

    if not files:
        raise ValueError(
            f"[ERROR] No bundled FASTA files found in {bundle_dir}"
        )

    return files


# ===================== Pairwise alignment with krishna ===========================
checkpoint krishna:
    input:
        genomes=bundled_genomes
    output:
        directory(OUTDIR + "{prefix}/01_krishna")
    params:
        out_dir = OUTDIR + "{prefix}/01_krishna/",
        thrd = MAX_THREADS,
        filtid = KRISHNA_FILTID,
        filtlen = KRISHNA_FILTLEN
    shell:
        """
        echo "====================== Krishna: Performing pairwise alignment... ======================"
        mkdir -p {params.out_dir}
        echo "-> Created directory {params.out_dir}"
        
        genomes=$(realpath {input.genomes})
        echo "-> $genomes will be processed with CARP"
        
        cd {params.out_dir}
        echo "-> Towards {params.out_dir}"
        
        echo "-> Running krishna..."
        matrix -threads={params.thrd} -krishnaflags="-tmp=. -threads={params.thrd} -log -filtid={params.filtid} -filtlen={params.filtlen}" $genomes
        echo "====================== Krishna completed! Output in {params.out_dir} ======================"
        """

import glob

def krishna_gffs_for_igor(wildcards):
    krishna_dir = checkpoints.krishna.get(prefix=wildcards.prefix).output[0]

    gffs = sorted(glob.glob(krishna_dir + "/*.gff"))

    if not gffs:
        raise ValueError(
            f"[ERROR] No Krishna GFF files found in {krishna_dir}"
        )

    return gffs

# Repeat family groupings into json file by igor
rule igor:
    input: 
        krishna_gffs = krishna_gffs_for_igor
    output: 
        gff = OUTDIR + "{prefix}/02_igor/{prefix}_krishna.gff",
        json = OUTDIR + "{prefix}/02_igor/{prefix}_krishna.json"
    params:
        out_dir = OUTDIR + "{prefix}/02_igor/",
    shell: 
        """
        echo "====================== igor: Grouping repeat families into JSON file... ======================"
        mkdir -p {params.out_dir}
        echo "-> Created directory {params.out_dir}"
        
        cat {input.krishna_gffs} > {output.gff}
        echo "-> Running igor"
        igor -in {output.gff} -out {output.json}
        echo "====================== igor completed! Output in {params.out_dir} ======================"
        """

# Use seqer to generate consensus sequences from genome intervals
rule seqer:
    input: 
        igor_out = rules.igor.output.json,
        genome = INPUT_DIR + "{prefix}.fa"
    output: 
        gff = OUTDIR + "{prefix}/03_seqer/{prefix}_krishna.igor.gff",
        cat = OUTDIR + "{prefix}/03_seqer/{prefix}.mfa",
        consensus = OUTDIR + "{prefix}/03_seqer/{prefix}_temp_ConsensusSequences.fa"
    params:
        out_dir = OUTDIR + "{prefix}/03_seqer/",
        consensus_dir = OUTDIR + "{prefix}/03_seqer/consensus/",
        maxFam = SEQER_MAX_FAMILY_SIZE,
        minLen = SEQER_MIN_LENGTH,
        thrd = MAX_THREADS
    shell:
        """
        echo "================ seqer: Generating consensus sequences from genome intervals... ================"
        echo "-> Created directories: {params.out_dir} and {params.consensus_dir}"
        mkdir -p {params.out_dir}
        mkdir -p {params.consensus_dir}
        echo "-> Running gffer..."
        gffer < {input.igor_out} > {output.gff}
        cat {input.genome} > {output.cat}
        
        echo "-> Running seqer..."
        seqer -aligner=muscle -dir={params.consensus_dir} -fasta=true -maxFam={params.maxFam} -subsample=true -minLen={params.minLen} -threads={params.thrd} -ref={output.cat} {output.gff}
        echo "-> seqer completed! Output in {params.out_dir}"
        echo "-> Combining consensus..."
        find {params.consensus_dir} -maxdepth 1 -name '[!.]*.fq' -print0 | xargs -r0 cat > {output.consensus}
        echo "====================== Consensus combined! Output in {params.out_dir}consensus ======================"
        """

rule finalize_consensus:
    input:
        OUTDIR + "{prefix}/03_seqer/{prefix}_temp_ConsensusSequences.fa"
    output:
        OUTDIR + "{prefix}/04_finalize_consensus/{prefix}_ConsensusSequences.fa"
    params:
        out_dir = OUTDIR + "{prefix}/04_finalize_consensus"
    shell:
        """
        echo "================ finalize_consensus ================"
        mkdir -p {params.out_dir}
        cat {input} > {output}
        echo "================ finalize_consensus completed! check files ================"
        """

rule repeats_in_repbase:
    input:
        consensus = OUTDIR + "{prefix}/04_finalize_consensus/{prefix}_ConsensusSequences.fa",
        repbase = REPBASE_LIBRARY,
        our_known = OUR_KNOWN_LIBRARY,
        converter = "src/ins_json_to_map.py"
    output:
        json = OUTDIR + "{prefix}/05_repeats_in_repbase/{prefix}_ConsensusSequences.fa.ins.json",
        map = OUTDIR + "{prefix}/05_repeats_in_repbase/{prefix}_ConsensusSequences.fa.map"
    params:
        out_dir = OUTDIR + "{prefix}/05_repeats_in_repbase",
        thrd = MAX_THREADS,
        # Convert the Python boolean to the lowercase value required by ins-fast
        first_softmask = str(
            REPEATS_IN_REPBASE_FIRST_SOFTMASK
        ).lower()
    log:
        OUTDIR + "{prefix}/05_repeats_in_repbase/{prefix}_INS.log"
    shell:
        """
        echo "====================== Identifying repeats in RepBase... ======================"
        mkdir -p {params.out_dir}
        exec > >(tee {log}) 2>&1
        ins-fast \
            -query {input.consensus} \
            -lib {input.repbase} \
            -lib {input.our_known} \
            -mode normal \
            -first-softmask={params.first_softmask} \
            -fast-realign=true \
            -cull=true \
            -pool=true \
            -cores {params.thrd} \
            -json \
            > {output.json}
        #  Convert json file to map which map was generated by censor
        python3 {input.converter} {output.json} {output.map}
        echo "====================== repeats_in_repbase completed! Output in {params.out_dir}======================"
        """

rule classify_consensus:
    input:
        script = "src/ClassifyConsensusSequences_CARP_snkm.java",
        consensus = OUTDIR + "{prefix}/04_finalize_consensus/{prefix}_ConsensusSequences.fa",
        map = OUTDIR + "{prefix}/05_repeats_in_repbase/{prefix}_ConsensusSequences.fa.map",
        repbase = REPBASE_LIBRARY,
        our_known = OUR_KNOWN_LIBRARY
    output:
        known = OUTDIR + "{prefix}/06_classify_consensus/{prefix}_known.txt",
        partial = OUTDIR + "{prefix}/06_classify_consensus/{prefix}_partial.txt",
        check = OUTDIR + "{prefix}/06_classify_consensus/{prefix}_check.txt",
        not_known = OUTDIR + "{prefix}/06_classify_consensus/{prefix}_notKnown.fa",
        not_known_gff = OUTDIR + "{prefix}/06_classify_consensus/{prefix}_notKnown.fa.gff",
        log = OUTDIR + "{prefix}/06_classify_consensus/{prefix}_ClassifyConsensusSequences.log",
    params:
        out_dir = OUTDIR + "{prefix}/06_classify_consensus"
    shell:
        """
        echo "====================== Processing and classifying consensus sequences... ======================"
        mkdir -p {params.out_dir}
        javac -d src {input.script}
        echo "-> Classify consensus sequence"
        java -cp src ClassifyConsensusSequences_CARP_snkm \
            {input.consensus} \
            {input.map} \
            {input.repbase} \
            {input.our_known} \
            {output.known} \
            {output.partial} \
            {output.check} \
            {output.not_known} \
            {output.not_known_gff} \
            > "{output.log}" 2>&1
        echo "====================== Classification completed! Output in {params.out_dir} ======================"
        """

rule filter_sequences:
    input:
        not_known = OUTDIR + "{prefix}/06_classify_consensus/{prefix}_notKnown.fa",
        protein_lib = SWISSPROT_LIBRARY,
        gb_te_lib = GB_TE_LIBRARY,
        retrovirus_lib = RETROVIRUS_LIBRARY
    output:
        protein = OUTDIR + "{prefix}/07_filter_sequences/{prefix}_notKnown.fa.spwb.ncbi",
        protein_gff = OUTDIR + "{prefix}/07_filter_sequences/{prefix}_notKnown.fa.spwb.gff",
        gb_te = OUTDIR + "{prefix}/07_filter_sequences/{prefix}_notKnown.fa.tewb.ncbi",
        gb_te_gff = OUTDIR + "{prefix}/07_filter_sequences/{prefix}_notKnown.fa.tewb.gff",
        retrovirus = OUTDIR + "{prefix}/07_filter_sequences/{prefix}_notKnown.fa.ervwb.ncbi",
        retrovirus_gff = OUTDIR + "{prefix}/07_filter_sequences/{prefix}_notKnown.fa.ervwb.gff"
    params:
        out_dir = OUTDIR + "{prefix}/07_filter_sequences",
        thrd = MAX_THREADS
    log:
        OUTDIR + "{prefix}/07_filter_sequences/logs/{prefix}_filter_sequences.log"
    shell:
        """
        mkdir -p {params.out_dir}/db {params.out_dir}/logs

        # Run log
        exec > >(tee {log}) 2>&1

        echo "====================== Identifying notknowns from protein database... ======================"
        makeblastdb -in "{input.protein_lib}" -dbtype prot -out "{params.out_dir}/db/sprot"

        blastx \
            -db {params.out_dir}/db/sprot \
            -query {input.not_known} \
            -max_hsps 1 \
            -seg no \
            -evalue 0.00001 \
            -num_threads {params.thrd} \
            -max_target_seqs 1 \
            -word_size 2 \
            -outfmt 6 \
            -out {output.protein}

        # Preserve the existing Swiss-Prot ID without adding another sp| prefix
        awk \
            '{{print $1"\t""blast""\t""hit""\t"$7"\t"$8"\t"$11"\t"".""\t"".""\t""Target "$2" "$9" "$10}}' \
            {output.protein} \
            > {params.out_dir}/tmp

        # Swap descending query coordinates
        awk \
            '{{if($4>$5) print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9" "$10" "$11" "$12; else print $0}}' \
            {params.out_dir}/tmp \
            > {output.protein_gff}

        echo "-> Identifying notknowns from GB_TE database... "
        makeblastdb -in "{input.gb_te_lib}" -dbtype prot -out "{params.out_dir}/db/GB_TE.new"

        blastx \
            -db {params.out_dir}/db/GB_TE.new \
            -query {input.not_known} \
            -max_hsps 1 \
            -seg no \
            -evalue 0.00001 \
            -num_threads {params.thrd} \
            -max_target_seqs 1 \
            -word_size 2 \
            -outfmt 6 \
            -out {output.gb_te}

        awk \
            '{{print $1"\t""blast""\t""hit""\t"$7"\t"$8"\t"$11"\t"".""\t"".""\t""Target sp|"$2" "$9" "$10}}' \
            {output.gb_te} \
            > {params.out_dir}/tmp

        awk \
            '{{if($4>$5) print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9" "$10" "$11" "$12; else print $0}}' \
            {params.out_dir}/tmp \
            > {output.gb_te_gff}

        echo "-> Identifying notknowns from retrovirus database... "
        makeblastdb -in "{input.retrovirus_lib}" -dbtype nucl -out "{params.out_dir}/db/retrovirus"

        tblastx \
            -db {params.out_dir}/db/retrovirus \
            -query {input.not_known} \
            -max_hsps 1 \
            -seg no \
            -evalue 0.00001 \
            -num_threads {params.thrd} \
            -max_target_seqs 1 \
            -word_size 2 \
            -outfmt 6 \
            -out {output.retrovirus}

        # Match step 11 annotation keys by removing gi|number| and the final |
        # Example: gi|114146046|gb|DQ838494.1| becomes gb|DQ838494.1
        awk \
            '{{sub(/^gi[|][^|]+[|]/, "", $2); sub(/[|]$/, "", $2); print $1"\t""blast""\t""hit""\t"$7"\t"$8"\t"$11"\t"".""\t"".""\t""Target "$2" "$9" "$10}}' \
            {output.retrovirus} \
            > {params.out_dir}/tmp

        awk \
            '{{if($4>$5) print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9" "$10" "$11" "$12; else print $0}}' \
            {params.out_dir}/tmp \
            > {output.retrovirus_gff}

        echo "====================== Sequence filtering completed! Output in {params.out_dir} ========================"
        """

rule get_protein_info_from_consensus:
    input:
        script = "src/GetProteins_CARP_snkm.java",
        not_known = OUTDIR + "{prefix}/06_classify_consensus/{prefix}_notKnown.fa",
        protein_hits = OUTDIR + "{prefix}/07_filter_sequences/{prefix}_notKnown.fa.spwb.gff"
    output:
        protein = OUTDIR + "{prefix}/08_get_protein_info/{prefix}_protein.txt",
        not_protein = OUTDIR + "{prefix}/08_get_protein_info/{prefix}_notKnownNotProtein.fa",
        log = OUTDIR + "{prefix}/08_get_protein_info/{prefix}_GetProteins.log"
    params:
        out_dir = OUTDIR + "{prefix}/08_get_protein_info"
    shell:
        """
        echo "==================== Retrieving protein information from consensus sequences... ===================="
        mkdir -p {params.out_dir}

        javac -d src {input.script}

        java -cp src GetProteins_CARP_snkm \
            {input.not_known} \
            {input.protein_hits} \
            {output.protein} \
            {output.not_protein} \
            > "{output.log}" 2>&1

        echo "====================== Protein information retrieved! Output in {params.out_dir} ======================"
        """

rule check_ssr:
    input:
        not_protein = OUTDIR + "{prefix}/08_get_protein_info/{prefix}_notKnownNotProtein.fa",
    output:
        out_dir = OUTDIR + "{prefix}/09_check_ssr/{prefix}_notKnownNotProtein.phobos"
    shell:
        """
        echo "====================== Check for simple sequence repeats... ======================"
        phobos-linux-gcc4.1.2 -r 7 --outputFormat 0 --printRepeatSeqMode 0 {input.not_protein} > {output}
        echo "====================== SSR check completed! Output in {output.out_dir} ======================"
        """

rule identify_ssr_from_phobos_output:
    input:
        script = "src/IdentifySSRs_CARP_snkm.java",
        gbte_fa = GB_TE_LIBRARY,
        phobos = OUTDIR + "{prefix}/09_check_ssr/{prefix}_notKnownNotProtein.phobos",
        not_protein = OUTDIR + "{prefix}/08_get_protein_info/{prefix}_notKnownNotProtein.fa"
    output:
        ssr = OUTDIR + "{prefix}/10_identify_ssr/{prefix}_SSR.txt",
        not_ssr = OUTDIR + "{prefix}/10_identify_ssr/{prefix}_notKnownNotProteinNotSSR.fa",
        ssr_log = OUTDIR + "{prefix}/10_identify_ssr/{prefix}_IdentifySSRs.log",
        gbte_out = OUTDIR + "{prefix}/10_identify_ssr/{prefix}_GB_TE.fa"
    params:
        out_dir = OUTDIR + "{prefix}/10_identify_ssr"
    shell:
        """
        echo "================== Identifying SSRs sequences from Phobos output... =================="
        mkdir -p {params.out_dir}

        javac -d src {input.script}

        phobos=$(realpath {input.phobos})
        not_protein=$(realpath {input.not_protein})

        java -cp src IdentifySSRs_CARP_snkm \
            "$phobos" \
            "$not_protein" \
            "{output.ssr}" \
            "{output.not_ssr}" \
            > "{output.ssr_log}" 2>&1

        perl -pe "s/^>/>gi|GBTE|sp|/g" \
            {input.gbte_fa} > {output.gbte_out}

        sed -i 's/ /| /' {output.gbte_out}

        echo "================== Data format conversion completed! Output in {params.out_dir} =================="
        """

rule generate_annotated_repeat_lib:
    input:
        script = "src/GenerateAnnotatedLibrary_CARP_snkm.java",
        consensus = OUTDIR + "{prefix}/04_finalize_consensus/{prefix}_ConsensusSequences.fa",
        consensus_map = OUTDIR + "{prefix}/05_repeats_in_repbase/{prefix}_ConsensusSequences.fa.map",
        te_gff = OUTDIR + "{prefix}/07_filter_sequences/{prefix}_notKnown.fa.tewb.gff",
        gb_te = OUTDIR + "{prefix}/10_identify_ssr/{prefix}_GB_TE.fa",
        erv_gff = OUTDIR + "{prefix}/07_filter_sequences/{prefix}_notKnown.fa.ervwb.gff",
        retrovirus = RETROVIRUS_LIBRARY,
        ssr = OUTDIR + "{prefix}/10_identify_ssr/{prefix}_SSR.txt",
        protein = OUTDIR + "{prefix}/08_get_protein_info/{prefix}_protein.txt",
        known = OUTDIR + "{prefix}/06_classify_consensus/{prefix}_known.txt"
    output:
        library = OUTDIR + "{prefix}/11_generate_annotated_repeat_lib/{prefix}_Denovo_TE_Library.fasta",
        log = OUTDIR + "{prefix}/11_generate_annotated_repeat_lib/{prefix}_GenerateAnnotatedLibrary.log"
    params:
        out_dir = OUTDIR + "{prefix}/11_generate_annotated_repeat_lib"
    shell:
        """
        echo "================== Finally, generating annotated repeat library... =================="
        mkdir -p {params.out_dir}
        exec > >(tee {output.log}) 2>&1

        echo Start time is `date`
        echo Running on `hostname`

        javac -d src {input.script}

        CONSENSUS=$(realpath "{input.consensus}")
        CONSENSUS_MAP=$(realpath "{input.consensus_map}")
        TE_GFF=$(realpath "{input.te_gff}")
        GB_TE=$(realpath "{input.gb_te}")
        ERV_GFF=$(realpath "{input.erv_gff}")
        RETROVIRUS=$(realpath "{input.retrovirus}")
        SSR=$(realpath "{input.ssr}")
        PROTEIN=$(realpath "{input.protein}")
        KNOWN=$(realpath "{input.known}")
        LIBRARY=$(realpath -m "{output.library}")

        cd lib

        java -cp ../src GenerateAnnotatedLibrary_CARP_snkm \
            "${{CONSENSUS}}" \
            "${{CONSENSUS_MAP}}" \
            "${{TE_GFF}}" \
            "${{GB_TE}}" \
            "${{ERV_GFF}}" \
            "${{RETROVIRUS}}" \
            "${{SSR}}" \
            "${{PROTEIN}}" \
            "${{KNOWN}}" \
            "${{LIBRARY}}" 
        echo "================== Check TE library results in {params.out_dir} =================="
        echo End time is `date`
        """
# This is the block to use RepeatMasker to 
# rule annotate_genome:
#     input:
#         repbase = "lib/Vertebrates_RepBase31.09.fa",
#         denovo = OUTDIR + "{prefix}/11_generate_annotated_repeat_lib/{prefix}_Denovo_TE_Library.fasta",
#         mfa_file = OUTDIR + "{prefix}/03_seqer/{prefix}.mfa"
#     output:
#         directory(OUTDIR + "{prefix}/12_genome_annotation")
#     params:
#         thrd = MAX_THREADS,
#         out_dir = OUTDIR + "{prefix}/12_genome_annotation",
#     shell:
#         """
#         echo "================== Annotate genome by RepeatMasker =================="
#         mkdir -p {params.out_dir}
#         combined="{params.out_dir}/{wildcards.prefix}_combined_library2.fa"
#         cat {input.repbase} {input.denovo} > "$combined"
#         env FAMDB_DIR= RepeatMasker -engine rmblast -pa 16 -a -nolow -norna -dir {params.out_dir} -lib "$combined" {input.mfa_file}
#         echo "================== Check genome annotation result in {params.out_dir} ================== "
#         """ 

rule annotate_genome:
    input:
        repbase = REPBASE_LIBRARY,
        denovo = OUTDIR + "{prefix}/11_generate_annotated_repeat_lib/{prefix}_Denovo_TE_Library.fasta",
        mfa_file = OUTDIR + "{prefix}/03_seqer/{prefix}.mfa",
        json_converter = "src/ins_json_to_map.py",
        tbl_converter = "src/ins_map_to_tbl.py",
        classes = "src/repbase_repeat_classes.tsv"
    output:
        query = OUTDIR + "{prefix}/12_genome_annotation/{prefix}.mfa",
        json = OUTDIR + "{prefix}/12_genome_annotation/{prefix}.ins.json",
        map = OUTDIR + "{prefix}/12_genome_annotation/{prefix}.ins.map",
        masked = OUTDIR + "{prefix}/12_genome_annotation/{prefix}.mfa-masked.fasta",
        tbl = OUTDIR + "{prefix}/12_genome_annotation/{prefix}.mfa.tbl"
    params:
        thrd = MAX_THREADS,
        # Convert the Python boolean to the lowercase value required by ins-fast
        first_softmask = str(
            ANNOTATE_GENOME_FIRST_SOFTMASK
        ).lower(),
        out_dir = OUTDIR + "{prefix}/12_genome_annotation"
    log:
        OUTDIR + "{prefix}/12_genome_annotation/logs/{prefix}_INS.log"
    shell:
        """
        mkdir -p "{params.out_dir}/logs"
        exec > >(tee "{log}") 2>&1

        echo Running on `hostname`
        echo Start time is `date`
        trap 'echo End time is `date`' EXIT

        # Stage the query here because INS writes the masked FASTA beside its query
        cp "{input.mfa_file}" "{output.query}"

        echo "================== Annotate genome with INS =================="

        ins-fast \
            -query {output.query} \
            -lib {input.repbase} \
            -lib {input.denovo} \
            -mode normal \
            -first-softmask={params.first_softmask} \
            -fast-realign=true \
            -cull=true \
            -pool=true \
            -cores {params.thrd} \
            -json \
            > "{output.json}"

        # Convert the INS JSON to the 10-column compatibility MAP file
        python3 "{input.json_converter}" "{output.json}" "{output.map}"

        # Summarize genome coverage using the RepBase classification hierarchy
        python3 "{input.tbl_converter}" \
            --map "{output.map}" \
            --genome "{output.query}" \
            --repbase "{input.repbase}" \
            --te_library "{input.denovo}" \
            --classes "{input.classes}" \
            --output "{output.tbl}"

        echo "================== Genome annotation completed in {params.out_dir} =================="
        """

rule plot_ins_map:
    input:
        script = "src/plot_ins_map.R",
        consensus = OUTDIR + "{prefix}/04_finalize_consensus/{prefix}_ConsensusSequences.fa",
        map = OUTDIR + "{prefix}/05_repeats_in_repbase/{prefix}_ConsensusSequences.fa.map",
        protein = OUTDIR + "{prefix}/07_filter_sequences/{prefix}_notKnown.fa.spwb.ncbi",
        gb_te = OUTDIR + "{prefix}/07_filter_sequences/{prefix}_notKnown.fa.tewb.ncbi",
        retrovirus = OUTDIR + "{prefix}/07_filter_sequences/{prefix}_notKnown.fa.ervwb.ncbi"
    output:
        directory(OUTDIR + "{prefix}/13_consensus_plots"),
    params:
        # Pass selected IDs as one comma-separated argument
        consensus_ids = PLOT_CONSENSUS_IDS_ARGUMENT
    log:
        OUTDIR + "{prefix}/13_consensus_plots/{prefix}_plot_ins_map.log"
    shell:
        """
        mkdir -p {output}
        exec > >(tee "{log}") 2>&1

        Rscript {input.script} \
            {input.consensus} \
            {input.map} \
            {input.protein} \
            {input.gb_te} \
            {input.retrovirus} \
            {output} \
            {params.consensus_ids}
        """

rule make_report:
    input:
        script = "src/make_report.py",
        # Always wait until the core pipeline has completed stage 11
        core_library = rules.generate_annotated_repeat_lib.output.library,
        # Wait for stage 12 only when annotate_genome is enabled
        annotation = optional_genome_annotation,
        # Wait for plots only when plot_ins_map is enabled
        consensus_plots = optional_consensus_plots
    output:
        file = OUTDIR + "{prefix}/{prefix}_CARP_snkm_report.html"
    params:
        sample_dir = OUTDIR + "{prefix}"
    shell:
        """
        echo "================== Generating report ... =================="
        python3 {input.script} {params.sample_dir} -o {output.file}
        echo "================== Completed! CARP managed to swim to shore! =================="
        """