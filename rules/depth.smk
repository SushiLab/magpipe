#!/usr/bin/env python

"""
name: depth_snake
description: This snakemake subpipeline is part of the MAGs building pipeline
dependency: None
author: Lucas Paoli (paolil@ethz.ch)
rules:
    - depth [done]
    - atomize [TBA]
    - fragmentize [TBA]
    - randomize [TBA]
"""

rule depth:
    input:
        lambda wildcards: ancient(magpipe_config.input.loc[wildcards.dataset, "backmapping"] + '/' + wildcards.sample + '/')
    output:
        #FIXME when running on koch, temporary science files can't be removed and jobs fail
        tmp_bams = temp(directory(magpipe_config.fast_dir + '/{dataset}/backmapping/{sample}{depth_sfx}/')),
        marker = touch(magpipe_config.metag_dir + '/{dataset}/{sample}/depth_files/{sample}{depth_sfx}-depth.done')
    params:
        depthfile = magpipe_config.metag_dir + '/{dataset}/{sample}/depth_files/{sample}{depth_sfx}-depth.tsv' 
    conda:
        magpipe_config.envs_dir + "/metabat2_env.yaml"
    threads:
        16
    resources:
        load = 40,
        mem = 8000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 60 # in minutes
    benchmark:
        magpipe_config.metag_dir + '/{dataset}/{sample}/depth_files/{sample}{depth_sfx}-depth.benchmark'
    log:
        log = magpipe_config.metag_dir + '/{dataset}/{sample}/depth_files/{sample}{depth_sfx}-depth.log',
        command = magpipe_config.metag_dir + '/{dataset}/{sample}/depth_files/{sample}{depth_sfx}-depth.command',
        qoutfile = magpipe_config.metag_dir + '/{dataset}/{sample}/depth_files/{sample}{depth_sfx}-depth.qout',
        qerrfile = magpipe_config.metag_dir + '/{dataset}/{sample}/depth_files/{sample}{depth_sfx}-depth.qerr'
    shell:
        '''
        command="
        export OMP_NUM_THREADS={threads}; rsync {input}/*.bam {output.tmp_bams}; jgi_summarize_bam_contig_depths --outputDepth {params.depthfile} {output.tmp_bams}/*.bam &> {log.log}";
        echo "$command" > {log.command};
        eval "$command"
        '''

