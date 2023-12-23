# SnakeCARP
# This is the snakemake prpeline for CARP
# 2023-12-23
# Version: 1.0

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

rule igor:
    input: 
