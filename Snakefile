# SnakeCARP
# This is the snakemake prpeline for CARP
# 2023-12-23
# Version: 1.0

# Pairwise alignment with krishna
rule krishna:
    input: 
        "fa/{prefix}.fa"
    output:
        "output/1_krishna/{prefix}.gff"
    params:
        krishna_output_dir = "output/1_krishna/",
        threads = 8,
        filtid = 0.94,
        filtlen = 400
    shell:
        """
        mkdir -p {params.krishna_output_dir}
        krishna -threads={params.threads} \
            -tmp={params.krishna_output_dir} \
            -log \
            -filtid={params.filtid} \
            -filtlen={params.filtlen} \
            -query {input} \
            -target {input} \
            -out {output}
        # Move log file to krishna folder
        mv krishna-*.log {params.krishna_output_dir}
        """

# Repeat family groupings into json file by igor
rule igor:
    input: "output/1_krishna/{prefix}.gff"
    output: "output/2_igor/{prefix}_krishna.json"
    params:
        igor_output_dir = "output/2_igor/"
    shell: 
        """
        mkdir -p {params.igor_output_dir}
        find  output/1_krishna \
            -maxdepth 1 \
            -name '[!.]*.gff' \
            -print0 \
            |xargs \
            -r0 cat \
            > {params.igor_output_dir}{wildcards.prefix}_krishna.gff
        igor -in {params.igor_output_dir}{wildcards.prefix}_krishna.gff -out {output}
        """

# Use seqer to generate consensus sequences from genome intervals        
rule seqer:
    input: "output/2_igor/{prefix}_krishna.json"
    output: 
            gff = "output/3_seqer/{prefix}_krishna.igor.gff",
            consensus_dir = directory("output/3_seqer/consensus/{prefix}/")
    params:
        seqer_output_dir = "output/3_seqer/"
    shell:
        """
        mkdir -p {output.consensus_dir}
        # Use gffer convert JSON file to gff file
        gffer < {input} > {output.gff}
        cat fa/{wildcards.prefix}.fa >> {params.seqer_output_dir}{wildcards.prefix}.mfa
        seqer -aligner=muscle \
        -dir={output.consensus_dir} \
        -fasta=true \
        -maxFam=100 \
        -subsample=true \
        -minLen=0.95 \
        -threads=12 \
        -ref={params.seqer_output_dir}{wildcards.prefix}.mfa {output.gff}
        """

rule repeat_annotation:
    input:
        # Find all the .fq files by snakemake glob_wildcards function
        consensus_dir = rules.seqer.output.consensus_dir
    output: 
        'output/4_rep_anno/{prefix}_ConsensusSequences.fa'
    params:
        rep_anno_output_dir = 'output/4_rep_anno/',
        rep_lib = 'libraries/replib/'
    shell:
        """
        mkdir -p {params.rep_anno_output_dir}
        # Find all the .fq files combine as ConsensusSequences.fa
        cat {input.consensus_dir}/[!.]*.fq > {output}
        censor -bprm cpus=8 \
        -lib {params.rep_lib}Vertebrates.fa \
        -lib {params.rep_lib}our_known_reps_20130520.fasta \
        {output}
        mv censor*.log {params.rep_anno_output_dir}
        mv {wildcards.prefix}_ConsensusSequences.fa.* {params.rep_anno_output_dir}
        """
