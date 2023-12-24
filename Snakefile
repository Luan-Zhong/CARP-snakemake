# SnakeCARP
# This is the snakemake prpeline for CARP
# 2023-12-23
# Version: 1.0

# Pairwise alignment with krishna
rule krishna:
    input: 'fa/{prefix}.fa'
    output: 'output/1_krishna/{prefix}.gff'
    params:
        krishna_output_dir = 'output/1_krishna/',
        threads = 8,
        filtid = 0.94,
        filtlen = 400
    shell:
        '''
        mkdir -p {params.krishna_output_dir}
        krishna -threads={params.threads} \
            -tmp=./ \
            -log \
            -filtid={params.filtid} \
            -filtlen={params.filtlen} \
            -query {input} \
            -target {input} \
            -out {output}
        '''

# Repeat family groupings into json file by igor
rule igor:
    input: 'output/1_krishna/{prefix}.gff'
    output: 'output/2_igor/{prefix}_krishna.json'
    params:
        igor_output_dir = 'output/2_igor/'
    shell: 
        '''
        mkdir -p {params.igor_output_dir}
        find  output/1_krishna \
            -maxdepth 1 \
            -name '[!.]*.gff' \
            -print0 \
            |xargs \
            -r0 cat \
            > {params.igor_output_dir}{wildcards.prefix}_krishna.gff
        igor -in {params.igor_output_dir}{wildcards.prefix}_krishna.gff \
        -out {output}
        '''
        
