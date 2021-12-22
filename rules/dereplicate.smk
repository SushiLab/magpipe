#!/usr/bin/env python

"""
name: dereplicate_snake
description: This snakemake subpipeline is part of the MAGs building pipeline
dependency: prepare_phase_2
author: Lucas Paoli (paolil@ethz.ch)
rules:
    - drep [done]
    - specI [TBA]
"""

rule drep:
    input:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-export_evaluation.done'
    output:
        touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-drep_{ani}.done')
    params:
        genomes = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-genomes/*.fa',
        checkm_summary = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-evaluation_for_drep.csv',
        drep_output = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-drep_{ani}',
        ani = '{ani}',
        overlap = 0.2
    conda:
        magpipe_config.envs_dir + "/drep_env.yaml"
    threads:
        64 #48
    resources:
        load = 100,
        mem = 3000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 21000 # In minutes, 21000 min = 350 h ~=14.5 days
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-drep_{ani}.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-drep_{ani}.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-drep_{ani}.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-drep_{ani}.qout',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-drep_{ani}.qerr'
    shell:
        '''
        command="dRep dereplicate {params.drep_output} -g {params.genomes} --genomeInfo {params.checkm_summary} -comp 0 -con 1000 -sa {params.ani} -nc {params.overlap} -p {threads}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''

rule prepare_specI:
    input:
        magpipe_config.flag_gene_calling
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI/prepare_specI.done')
    params:
        magpipe_scripts = magpipe_config.scripts,
        annotations = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations',
        output = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI'
    conda:
        magpipe_config.envs_dir + "/pyutils_env.yaml"
    threads:
        1 
    resources:
        load = 4,
        mem = 3000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 180
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI/prepare_specI.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI/prepare_specI.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI/prepare_specI.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI/prepare_specI.qoutfile',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI/prepare_specI.qerrfile'
    shell:
        '''
        command="{params.magpipe_scripts}/prepare_specI.py -a {params.annotations} -o {params.output}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''

rule specI:
    input:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI/prepare_specI.done'
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI.done')
    params:
        specI = magpipe_config.dependencies['specI'],
        workdir = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI',
        mode = "snakemake_lsf"
        #mode = "snakemake_sge"
    # only needs snakemake
    #conda: 
    #    magpipe_config.envs_dir + "/pyutils_env.yaml"
    threads:
        1 
    resources:
        load = 4,
        mem = 5000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 7195
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI/specI.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI/specI.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI/specI.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI/specI.qoutfile',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-specI/specI.qerrfile'
    shell:
        '''
        command="{params.specI} -o {params.workdir} -m {params.mode}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''
