# SnakeCARP
# This is the snakemake prpeline for CARP
# 2023-05-30
# Version: 2.0
# a) Use singularity container lib/CARP_02_w_py.sif
# b) Now, Outputs all in "output/" directory, see global variables OUTDIR

# Understand snakefile

# a) Use previous block result in next block
# e.g.: use igor json out file as seqer input
# rules.igor.output.json

# Understand my shortcut
# thrd: threads
# 

####################
# Global variables #
####################

# Set input genome here
# e.g.: genome file "X.fa"
# Set preFix = ["X"]
preFix = ["X"]

# Set some directory
# Don't touch, Unless you want to change output directory name
OUTDIR = "output/"

# Set some resources
# threads
MAX_THREADS = 4


onstart:
    print("#################################\n#### Starting CARP... ><(((('>
 ####\n#################################\n")
    shell("mkdir -p output")

# Use singularity local contianer
singularity: "lib/CARP_02_w_py.sif"


# Set output file inside expand function to run seleted step
# krishna: "output/{prefix}.gff"
# igor: "output/{prefix}_krishna.json"
# seqer: "output/{prefix}_krishna.igor.gff"
# repeats_in_repbase:
# "output/ConsensusSequences.fa.aln",
# "output/ConsensusSequences.fa.found",
# "output/ConsensusSequences.fa.idx",
# "output/ConsensusSequences.fa.map",
# "output/ConsensusSequences.fa.masked",

# classify_consensus: 
# "output/known.txt",
# "output/partial.txt",
# "output/check.txt",
# "output/notKnown.fa",
# "output/notKnown.fa.gff",
# "lib/ClassifyConsensusSequences.class",
# "output/ClassifyConsensusSequences.log"

# filter_sequences:
# "output/notKnown.fa.spwb",
# "output/notKnown.fa.tewb",
# "output/notKnown.fa.ervwb",

rule all:
    input:
        expand("output/{prefix}_krishna.igor.gff", prefix=preFix),
        # filter_sequences:
        "output/notKnown.fa.spwb",
        "output/notKnown.fa.spwb.gff",
        "output/notKnown.fa.tewb",
        "output/notKnown.fa.tewb.gff",
        "output/notKnown.fa.ervwb",
        "output/notKnown.fa.ervwb.gff",
        # get_protein_info_from_consensus
        "output/protein.txt",
        "output/notKnownNotProtein.fa",
        "output/GetProteins.log",
        # check_ssr
        "output/notKnownNotProtein.phobos",
        # identify_ssr_from_phobos_output
        "output/IdentifySSRs.log",
        "output/GB_TE.fa",
        # generate_annotated_repeat_lib
        "output/library/Denovo_TE_Library.fasta",
        "output/GenerateAnnotatedLibrary.log",
        "lib/GenerateAnnotatedLibrary.class"




################
# Rules Proper #
################
# Pairwise alignment with krishna
rule krishna:
    input: 
        "data/{prefix}.fa"
    output:
        OUTDIR + "{prefix}.gff"
    params:
        out_dir = OUTDIR,
        thrd = MAX_THREADS,
        filtid = 0.99,
        filtlen = 400
    shell:
        """
        echo "====================== Krishna: Performing pairwise alignment... ======================"
        matrix -threads={params.thrd} -krishnaflags="-tmp={params.out_dir} -threads={params.thrd} -log -filtid={params.filtid} -filtlen={params.filtlen} -target={input}" {input}
        mv {wildcards.prefix}.gff {params.out_dir}
        mv krishna-*.log {params.out_dir}
        echo "====================== Krishna completed! ======================"
        """

# Repeat family groupings into json file by igor
rule igor:
    input: rules.krishna.output
    output: 
        gff = OUTDIR + "{prefix}_krishna.gff",
        json = OUTDIR + "{prefix}_krishna.json"
    params:
        out_dir = OUTDIR
    shell: 
        """
        echo "====================== igor: Grouping repeat families into JSON file... ======================"
        find {params.out_dir} \
            -maxdepth 1 \
            -name '[!.]*.gff' \
            -print0 \
            |xargs \
            -r0 cat \
            > {output.gff}
        igor -in {output.gff} -out {output.json}
        echo "====================== igor completed! ======================"
        """

# Use seqer to generate consensus sequences from genome intervals
rule seqer:
    input: 
        igor_out = rules.igor.output.json,
        genome = rules.krishna.input
    output: 
        gff = OUTDIR + "{prefix}_krishna.igor.gff",
        cat = OUTDIR + "{prefix}.mfa",
        consensus = OUTDIR + "{prefix}_temp_ConsensusSequences.fa"
    params:
        out_dir = OUTDIR,
        maxFam = 100,
        minLen = 0.95,
        thrd = MAX_THREADS,
    shell:
        """
        echo "================ seqer: Generating consensus sequences from genome intervals... ================"
        gffer < {input.igor_out} > {output.gff}
        cat {input.genome} > {output.cat}
        seqer -aligner=muscle -dir={params.out_dir}/consensus -fasta=true -maxFam={params.maxFam} -subsample=true -minLen={params.minLen} -threads={params.thrd} -ref={output.cat} {output.gff}
        echo "====================== seqer completed! ======================"
        echo "====================== Combining consensus... ======================"
        find {params.out_dir}/consensus -maxdepth 1 -name '[!.]*.fq' -print0 | xargs -r0 cat > {output.consensus}
        echo "====================== Consensus combined! ======================"
        """



# This block is for snakemake prefix rule
rule finalize_consensus:
    input:
        expand(OUTDIR + "{prefix}_temp_ConsensusSequences.fa", prefix=preFix)
    output:
        OUTDIR + "ConsensusSequences.fa"
    shell:
        """
        cat {input} > {output}
        """

rule repeats_in_repbase:
    input:
        OUTDIR + "ConsensusSequences.fa"
    output:
        OUTDIR + "ConsensusSequences.fa.aln",
        OUTDIR + "ConsensusSequences.fa.found",
        OUTDIR + "ConsensusSequences.fa.idx",
        OUTDIR + "ConsensusSequences.fa.map",
        OUTDIR + "ConsensusSequences.fa.masked",
        # "lib/ClassifyConsensusSequences.class"
        # "output/known.txt",
        # "output/partial.txt",
        # "output/check.txt",
        # "output/notKnown.fa",
        # "output/notKnown.fa.gff"
    params:
        out_dir = OUTDIR,
        cpus = 4,
        repbase = "lib/29.04_RepBase_All.fa",
        our_known = "lib/our_known_reps_20130520.fasta",
    shell:
        """
        echo "====================== Identifying repeats in RepBase... ======================"
        censor -bprm cpus={params.cpus} -lib {params.repbase} -lib {params.our_known} {input}
        echo "====================== Identification completed, moving results to outdir...======================"
        mv ConsensusSequences.fa.* censor.*.log {params.out_dir}
        echo "====================== Files moved successfully! ======================"
        """

rule classify_consensus:
    input:
        script = "lib/ClassifyConsensusSequences.java",
        consensus = OUTDIR + "ConsensusSequences.fa",
        map = OUTDIR + "ConsensusSequences.fa.map"
    output:
        OUTDIR + "known.txt",
        OUTDIR + "partial.txt",
        OUTDIR + "check.txt",
        OUTDIR + "notKnown.fa",
        OUTDIR + "notKnown.fa.gff",
        "lib/ClassifyConsensusSequences.class",
        OUTDIR + "ClassifyConsensusSequences.log",
    params:
        out_dir = OUTDIR
    shell:
        """
        echo "====================== Processing and classifying consensus sequences...======================"
        cd lib
        javac ClassifyConsensusSequences.java
        java ClassifyConsensusSequences > ClassifyConsensusSequences.log 2>&1
        cd ..
        mv lib/ClassifyConsensusSequences.log {params.out_dir}
        echo "====================== Classification completed! ======================"
        """

rule filter_sequences:
    input:
        OUTDIR + "notKnown.fa"
    output:
        protein = OUTDIR + "notKnown.fa.spwb",
        protein_gff = OUTDIR + "notKnown.fa.spwb.gff",
        gb_te = OUTDIR + "notKnown.fa.tewb",
        gb_te_gff = OUTDIR + "notKnown.fa.tewb.gff",
        retrovirus = OUTDIR + "notKnown.fa.ervwb",
        retrovirus_gff = OUTDIR + "notKnown.fa.ervwb.gff",
    params:
        out_dir = OUTDIR,
        protein_lib = "lib/sprot/uniprot_sprot.fasta",
        gb_te_lib = "lib/BlastDB/GB_TE.new",
        retrovirus_lib = "lib/BlastDB/all_retrovirus",
        cpus = 4,
        script = "lib/wublastx2gff.py"
    shell:
        """
        echo "====================== Identifying notknowns from protein database... ======================"
        blastx {params.protein_lib} {input} -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus={params.cpus} > {output.protein}
        python {params.script} {output.protein} > {output.protein_gff}
        echo "====================== Identifying notknowns from GB_TE database... ======================"
        blastx {params.gb_te_lib} {input} -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus={params.cpus} > {output.gb_te}
        python {params.script} {output.gb_te} > {output.gb_te_gff}
        echo "====================== Identifying notknowns from retrovirus database... ======================"
        tblastx {params.retrovirus_lib} {input} -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus={params.cpus} > {output.retrovirus}
        python {params.script} {output.retrovirus} > {output.retrovirus_gff}
        echo "====================== Sequence filtering completed! ========================"
        """

rule get_protein_info_from_consensus:
    input:
        "lib/GetProteins.java"
    output:
        OUTDIR + "protein.txt",
        OUTDIR + "notKnownNotProtein.fa",
        "lib/GetProteins.class",
        OUTDIR + "GetProteins.log"
    params:
        out_dir = OUTDIR
    shell:
        """
        echo "==================== Retrieving protein information from consensus sequences... ===================="
        cd lib
        javac GetProteins.java
        java GetProteins > GetProteins.log 2>&1
        cd ..
        mv lib/GetProteins.log {params.out_dir}
        echo "====================== Protein information retrieved! ======================"
        """

rule check_ssr:
    input:
        OUTDIR + "notKnownNotProtein.fa"
    output:
        OUTDIR + "notKnownNotProtein.phobos"
    shell:
        """
        echo "====================== Check for simple sequence repeats... ======================"
        phobos-linux-gcc4.1.2 -r 7 --outputFormat 0 --printRepeatSeqMode 0 {input} > {output}
        echo "====================== Simple sequence repeats check completed! ======================"
        """

rule identify_ssr_from_phobos_output:
    input:
        script = "lib/IdentifySSRs.java",
        gbte_fa = "lib/GB_TE.2022-07-14.fa"
    output:
        ssr_log = OUTDIR + "IdentifySSRs.log",
        gbte_out = OUTDIR + "GB_TE.fa"
    params:
        out_dir = OUTDIR
    shell:
        """
        echo "================== Identifying SSRs sequences from Phobos output... =================="
        cd lib
        javac IdentifySSRs.java
        java IdentifySSRs > IdentifySSRs.log 2>&1
        cd ..
        mv lib/IdentifySSRs.log {params.out_dir}
        echo "============== SSRs sequences identification from Phobos output completed! =============="
        echo "=== GB_TE and retrovirus data are downloaded by using EDirect. Convert the data format before generating the annotated repeat library... ==="
        perl -pe "s/^>/>gi|GBTE|sp|/g" {input.gbte_fa} > {output.gbte_out}
        sed -i 's/ /| /' {output.gbte_out}
        echo "================== Data format conversion completed! =================="
        """

rule generate_annotated_repeat_lib:
    input:
        "lib/GenerateAnnotatedLibrary.java"
    output:
        OUTDIR + "library/Denovo_TE_Library.fasta",
        OUTDIR + "GenerateAnnotatedLibrary.log",
        "lib/GenerateAnnotatedLibrary.class"
    params:
        out_dir = OUTDIR
    shell:
        """
        "================== Finally, generating annotated repeat library... =================="
        mkdir -p output/library
        cd lib
        javac GenerateAnnotatedLibrary.java
        java GenerateAnnotatedLibrary > GenerateAnnotatedLibrary.log 2>&1
        cd ..
        mv lib/GenerateAnnotatedLibrary.log {params.out_dir}
        "================== CARP has completed the annotation. Check your results at output/library/Denovo_TE_Library.fasta! =================="
        """
