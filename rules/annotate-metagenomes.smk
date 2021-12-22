#!/usr/bin/env python

"""
name: annotate
description: This snakemake subpipeline is part of the MAGs building pipeline
dependency: prepare_phase_2
author: Lucas Paoli (paolil@ethz.ch)
rules:
    - antismash [done]
    - centrifuge [TBA]
"""

rule antismash:
    input:
        assembly = lambda wildcards: INPUT.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
    output:
        marker = touch(magpipe_config.metag_dir + '/{dataset}/{sample}/antismash/{sample}-antismash.done')
    params:
        antismash = magpipe_config.dependencies["antismash"],
        outdir = magpipe_config.metag_dir + '/{dataset}/{sample}/antismash/{sample}-antismash'
    #conda:
    #    magpipe_config.envs_dir + "/antismash_env.yaml" # Currently using the module system
    threads:
        12 # use 4 for big samples (>30G RAM)
    resources:
        load = 4,
        mem = 3000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 1200
    benchmark:
        magpipe_config.metag_dir + '/{dataset}/{sample}/antismash/{sample}-antismash.benchmark'
    log:
        log = magpipe_config.metag_dir + '/{dataset}/{sample}/antismash/{sample}-antismash.log',
        command = magpipe_config.metag_dir + '/{dataset}/{sample}/antismash/{sample}-antismash.command',
        qoutfile = magpipe_config.metag_dir + '/{dataset}/{sample}/antismash/{sample}-antismash.qout',
        qerrfile = magpipe_config.metag_dir + '/{dataset}/{sample}/antismash/{sample}-antismash.qerr'
    shell:
        '''
        command="{params.antismash} {input} --output-dir {params.outdir} --genefinding-tool prodigal-m --cpus {threads}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''
