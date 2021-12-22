#!/usr/bin/env python

"""
name: evaluate_snake
description: This snakemake subpipeline is part of the MAGs building pipeline
dependency: binning_snake
author: Lucas Paoli (paolil@ethz.ch)
rules:
    - checkm [done]
    - checkm_plots [done]
    - anvio_evaluate [TBA]
    - compare [done]
"""

rule checkm:
    input:
        marker = magpipe_config.metag_dir  + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-bins.done'
    output:
        marker = touch(magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-checkm.done')
    params:
        bins = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/bins/',
        checkm_lineage = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/checkm_lineage/',
        summary = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-checkm.summary'
    conda:
        magpipe_config.envs_dir + "/checkm_env.yaml"
    threads:
        8
    resources:
        load = 10,
        mem = 10000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 240 # in minutes
    benchmark:
        magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-checkm.benchmark'
    log:
        log = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-checkm.log',
        command = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-checkm.command',
        qoutfile = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-checkm.qout',
        qerrfile = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-checkm.qerr'
    shell:
        '''
        command="checkm lineage_wf -x fa {params.bins} {params.checkm_lineage} --threads {threads} -f {params.summary} --tab_table &> {log.log}";
        echo "$command" > {log.command};
        eval "$command"
        '''

rule anvi_evaluate:
    input:
        marker = magpipe_config.metag_dir  + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-bins.done',
        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
        #assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.contigs.min1000.fasta.gz'
    output:
        marker = touch(magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-anvi_evaluate.done'),
        pre_fasta = temp(magpipe_config.fast_dir + '/{dataset}/anvio_db/{sample}/{binning_method}/{sample}.scaffolds.min1000.pre.fasta'), # FIXME test temp outputs
        anvio_fasta = temp(magpipe_config.fast_dir + '/{dataset}/anvio_db/{sample}/{binning_method}/{sample}-{binning_method}.fasta') # FIXME test temp outputs
    params:
        scripts = magpipe_config.scripts,
        anvio_prefix = magpipe_config.fast_dir + '/{dataset}/anvio_db/{sample}/{binning_method}/{sample}-{binning_method}',
        bins_directory = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/bins/',
        collection_name = lambda wildcards: 'COL_' + re.sub("[^0-9a-zA-Z]+", "_", wildcards.sample) + '_' + wildcards.binning_method.replace('.', '_'),
        sample = lambda wildcards: 'SAMPLE_' + re.sub("[^0-9a-zA-Z]+", "_", wildcards.sample),
        summary = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-anvi_evaluate.summary',
        prompt = 'Did not find any bin from the binning results. This can happen be you should probably check that it makes sense.'
    conda:
        magpipe_config.envs_dir + "/anvio_env.yaml"
    threads:
        8#2
    resources:
        load = 4,
        mem = 7500,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 360 # in minutes
    benchmark:
        magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-anvi_evaluate.benchmark'
    log:
        log = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-anvi_evaluate.log',
        command = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-anvi_evaluate.command',
        qoutfile = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-anvi_evaluate.qout',
        qerrfile = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}-{binning_method}-anvi_evaluate.qerr',
    shell:
        '''
        command="
        {params.scripts}/prepare_anvio_collection.py \
        -i {params.bins_directory} \
        -o {params.anvio_prefix}.anvio.collection

        gunzip -c {input.assembly} | sed 's/\s.*$//g' > \
        {output.pre_fasta}

        anvi-script-reformat-fasta \
        {output.pre_fasta} \
        -o {output.anvio_fasta} \
        --keep-ids <(cut -f 1 {params.anvio_prefix}.anvio.collection)
        
        if [ -s {params.anvio_prefix}.anvio.collection ];
        then 
            anvi-gen-contigs-database \
            -f {output.anvio_fasta} \
            -o {params.anvio_prefix}.contigs.db \
            --split-length 0 # Don't split

            # Must create a blank profile to import collection... also must specify a sample name
            anvi-profile \
            -c {params.anvio_prefix}.contigs.db \
            -o {params.anvio_prefix}.blank.profile \
            --blank-profile --skip-hierarchical-clustering --sample-name {params.sample}

            anvi-import-collection \
            -c {params.anvio_prefix}.contigs.db \
            -p {params.anvio_prefix}.blank.profile/PROFILE.db \
            -C {params.collection_name} \
            --contigs-mode \
            {params.anvio_prefix}.anvio.collection

            anvi-run-hmms -c {params.anvio_prefix}.contigs.db -T {threads}

            anvi-estimate-genome-completeness \
            -c {params.anvio_prefix}.contigs.db \
            -p {params.anvio_prefix}.blank.profile/PROFILE.db \
            -C {params.collection_name} \
            -o {params.summary}
        else
            echo {params.prompt}
        fi"
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''

rule checkm_ext_genomes:
    output:
        marker = touch('{external_dataset}/{chunk}-checkm_ext.done')
    params:
        genome_chunk = '{external_dataset}/{chunk}/',
        checkm_output = '{external_dataset}/{chunk}-checkm_output',
        summary = '{external_dataset}/{chunk}-checkm_ext.summary'
    conda:
        magpipe_config.envs_dir + "/checkm_env.yaml"
    threads:
        8
    resources:
        load = 10,
        mem = 10000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 360 # in minutes
    benchmark:
        '{external_dataset}/{chunk}-checkm_ext.benchmark'
    log:
        log = '{external_dataset}/{chunk}-checkm_ext.log',
        command = '{external_dataset}/{chunk}-checkm_ext.command',
        qoutfile = '{external_dataset}/{chunk}-checkm_ext.qout',
        qerrfile = '{external_dataset}/{chunk}-checkm_ext.qerr'
    shell:
        '''
        command="checkm lineage_wf -x fa {params.genome_chunk} {params.checkm_output} --threads {threads} -f {params.summary} --tab_table &> {log.log}";
        echo "$command" > {log.command};
        eval "$command"
        '''

rule anvi_evaluate_ext_genomes:
    output:
        marker = touch('{external_dataset}/anvio_db/{genome}-anvi_evaluate_ext.done')
    params:
        genome = '{external_dataset}/processed_genomes/{genome}.fa',
        contigs_db = '{external_dataset}/anvio_db/{genome}.contigs.db',
        summary = '{external_dataset}/anvio_db/{genome}.summary'
    conda:
        magpipe_config.envs_dir + "/anvio_env.yaml"
    threads:
        2
    resources:
        load = 2,
        mem = 7500,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 30 # in minutes
    benchmark:
        '{external_dataset}/anvio_db/{genome}.benchmark'
    log:
        log = '{external_dataset}/anvio_db/{genome}.log',
        command = '{external_dataset}/anvio_db/{genome}.command',
        qoutfile = '{external_dataset}/anvio_db/{genome}.qout',
        qerrfile = '{external_dataset}/anvio_db/{genome}.qerr'
    shell:
        '''
        command="
        anvi-gen-contigs-database \
        -f {params.genome} \
        -o {params.contigs_db} \
        --split-length 0 # Don't split

        anvi-run-hmms -c {params.contigs_db} -T {threads}

        anvi-estimate-genome-completeness \
        -c {params.contigs_db} \
        -o {params.summary}"
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''
