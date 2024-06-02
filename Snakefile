# SnakeCARP
# This is the snakemake prpeline for CARP
# 2023-05-30
# Version: 2.0
# a) Use singularity container lib/CARP_01.sif
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
    print("#################################\n#### Starting CARP... |Xx> ####\n#################################\n")
    shell("mkdir -p output")

# Use singularity local contianer
singularity: "lib/CARP_01.sif"


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

rule all:
    input:
        expand("output/{prefix}_krishna.igor.gff", prefix=preFix),
        "output/notKnown.fa.spwb",
        "output/notKnown.fa.tewb",
        "output/notKnown.fa.ervwb",




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
        echo "====================== Krishna Start... ======================"
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
        echo "====================== igor Start... ======================"
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
        thrd = 4,
    shell:
        """
        echo "====================== Seqer Start... ======================"
        gffer < {input.igor_out} > {output.gff}
        cat {input.genome} > {output.cat}
        seqer -aligner=muscle -dir={params.out_dir}/consensus -fasta=true -maxFam={params.maxFam} -subsample=true -minLen={params.minLen} -threads={params.thrd} -ref={output.cat} {output.gff}
        echo "====================== Seqer complete! ======================"
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
        gb_te = OUTDIR + "notKnown.fa.tewb",
        retrovirus = OUTDIR + "notKnown.fa.ervwb"
    params:
        out_dir = OUTDIR,
        protein_lib = "lib/sprot/uniprot_sprot.fasta",
        gb_te_lib = "lib/BlastDB/GB_TE.new",
        retrovirus_lib = "lib/BlastDB/all_retrovirus",
        cpus = 4
    shell:
        """
        echo "====================== Identifying notknown from protein database... ======================"
        blastx {params.protein_lib} {input} -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus={params.cpus} > {output.protein}
        echo "====================== Identifying notknown from GB_TE database... ======================"
        blastx {params.gb_te_lib} {input} -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus={params.cpus} > {output.gb_te}
        echo "====================== Identifying notknown from Retrovirus database... ======================"
        tblastx {params.retrovirus_lib} {input} -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus={params.cpus} > {output.retrovirus}
        """
