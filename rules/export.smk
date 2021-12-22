#!/usr/bin/env python

"""
name: export
description: This snakemake subpipeline is part of the MAGs building pipeline
dependency: binning_snake, evaluate_snake
author: Lucas Paoli (paolil@ethz.ch)
rules:
    - export_evaluation [done]
    - integrate_external [TBA]
"""

rule export_evaluation:
    input:
        magpipe_config.flag_export_mags
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-export_evaluation.done')
    params:
        scripts = magpipe_config.scripts,
        metag_dir = magpipe_config.metag_dir,
        datasets = config['datasets'],
        method = config['binning'], # this should be a single method, failsafe in place
        datasets_table = magpipe_config.dataset_info,
        methods_table = magpipe_config.method_info,
        completeness = magpipe_config.cpl,
        contamination = magpipe_config.ctn
    conda:
        magpipe_config.envs_dir + "/pyutils_env.yaml"
    threads:
        1
    resources:
        load = 100,
        mem = 5000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 60
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-export_evaluation.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-export_evaluation.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-export_evaluation.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-export_evaluation.qout',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-export_evaluation.qerr'
    shell:
        '''
        command="{params.scripts}/export_evaluation.py --results_dir {params.metag_dir} \
        --datasets {params.datasets} --method {params.method} --output_dir {wildcards.path} \
        --filter --completeness {params.completeness} --contamination {params.contamination} \
        --export_sequences --table_for_drep --datasets_table {params.datasets_table} --methods_table {params.methods_table}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''
