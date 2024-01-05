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
    output: "output/3_seqer/{prefix}_krishna.igor.gff"
    params:
        seqer_output_dir = "output/3_seqer/"
    shell:
        """
        mkdir -p {params.seqer_output_dir}
        # Use gffer convert JSON file to gff file
        gffer < {input} > {params.seqer_output_dir}{wildcards.prefix}_krishna.igor.gff
        cat fa/{wildcards.prefix}.fa >> {params.seqer_output_dir}{wildcards.prefix}.mfa
        seqer -aligner=muscle \
        -dir=output/3_seqer/consensus/ \
        -fasta=true \
        -maxFam=100 \
        -subsample=true \
        -minLen=0.95 \
        -threads=12 \
        -ref={params.seqer_output_dir}{wildcards.prefix}.mfa {output}
        """

rule repeat_annotation:
    input: 
        'output/3_seqer/consensus/.fa'
    output: 
        'output/4_rep_anno/ConsensusSequences.fa'
    params:
        rep_anno_output_dir = 'output/4_rep_anno/',
        rep_lib = 'libraries/replib/'
    shell:
        """
        mkdir -p {params.rep_anno_output_dir}
        # Check if .fq file exist
        if [ -f 'output/3_seqer/consensus/your_file.fq' ]; then
            echo '### fq file exists ###'
        fi
        cat 'done' > {output}
        # Combine .fq files
        cat output/3_seqer/consensus/*.fq > {output}
        # Use censor to annotate consensus sequences with the Repbase library
        censor -bprm cpus=8 \
        -lib {params.rep_lib}Vertebrates.fa \
        -lib {params.rep_lib}our_known_reps_20130520.fasta \
        {output}
        # Move censor log file and related files to repeat annotation folder
        mv censor.*.log {wildcards.prefix}_ConsensusSequences.fa.* {params.rep_anno_output_dir}
        """
