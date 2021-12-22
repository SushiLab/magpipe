#!/usr/bin/env python

"""
name: Documentation
description: This is a dummy rule to use as a template in the pipeline
dependency: None
author: Lucas Paoli (paolil@ethz.ch)
"""

rule documentation:
    input:
        # If the input is the assembly
        lambda wildcards: ancient(INPUT.loc[wildcards.dataset, "backmapping"] + '/' + wildcards.sample + '/')
        # Other option is to call a marker from another rule
        marker = DIR + SUBDIR + NAME + '-dependency.done'
        # And other file input...
    output:
        # Use a marker as output
        marker = touch(DIR + SUBDIR + NAME + '-rule.done')
        # Use ther output files or directory with the necessary stamps, i.e. `temp()` or `directory()`
    params: # Non file parameters
        # e.g. the commandline tool
        tool = DEPENDENCIES["tool"],
        # e.g. the magpipe magpipe scripts
        magpipe_scripts = MAGPIPE_SCRIPTS
        #...
    conda: # the conda envionment, if appropriate
        ENVS_DIR + "/rule_env.yaml"
    threads: # Self explanatory
        16
    resources:
        # Define a load to better manage the worklow, used with `--resources load=XX`
        load = 10,
        # defined in the qsub line or bsub usage string `-R \"rusage[mem={params.mem},scratch={params.scratch}]\"`
        mem = 8000, # unit: Mb * threads
        scratch = 1000, # unit: Mb * threads
        # you can use a callable for the pipelinie to increase the resource if previous attempt failed
        time = 60# unit: minutes or hours:minutes
        # You can also go fancy with something like:
        time = lambda wildcards, attempt: attempt * 60
    benchmark:
        DIR + SUBDIR + NAME + '-rule.benchmark'
    log:
        log = DIR + SUBDIR + NAME + '-rule.log'
        command = DIR + SUBDIR + NAME + '-rule.command'
        qoutfile = DIR + SUBDIR + NAME + '-rule.qout'
        qerrfile = DIR + SUBDIR + NAME + '-rule.qerr'
    shell:
        '''
        command="export OMP_NUM_THREADS={threads}; export MKL_NUM_THREADS={threads}; export GOTO_NUM_THREADS={threads}; # Some tools using OMP/MKL require that to not go over the specified number of threads.
        do this with that &> {log.log}";
        echo "$command" > {log.command};
        eval "$command"
        '''
