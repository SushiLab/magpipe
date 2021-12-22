#!/usr/bin/env python

"""
name: identify-elements_genomes 
description: rules to identify the different genetic elements in genomes, e.g. plasmids, phages, prophages.. 
dependency: Primary analyses
author: Lucas Paoli (paolil@ethz.ch)
rules:
    - plasflow
    - plasmidfinder
    - cbar
    - vibrant
    - virsorter
    - deepvirfinder
    - checkv
    - phigaro
    - eukrep
    - ccontig
"""

rule checkv_genome:
    input:
        ancient(magpipe_config.flag_gene_calling)
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/checkv/{MAG}-checkv.done')
    params:
        folder = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/checkv/{MAG}-checkv',
        genome = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-genomes/{MAG}.fa',
        checkv_db = '/nfs/cds/scratch/paolil/databases/checkv-db-v0.6/'
    conda:
        magpipe_config.envs_dir + "/checkv_env.yaml"
    threads:
        8
    resources:
        load = 20,
        mem = 1000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 60
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/checkv/{MAG}-checkv.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/checkv/{MAG}-checkv.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/checkv/{MAG}-checkv.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/checkv/{MAG}-checkv.qoutfile',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/checkv/{MAG}-checkv.qerrfile'
    shell:
        '''
        command="checkv end_to_end {params.genome} {params.folder} -t {threads} -d {params.checkv_db}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''

rule vibrant_genome:
    input:
        ancient(magpipe_config.flag_gene_calling)
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/vibrant/{MAG}-vibrant.done')
    params:
        folder = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/vibrant/{MAG}-vibrant',
        genome = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-genomes/{MAG}.fa',
        vibrant_db = '/nfs/cds/scratch/paolil/databases/vibrant-db/databases/',
        vibrant_files = '/nfs/cds/scratch/paolil/databases/vibrant-db/files/'
    conda:
        magpipe_config.envs_dir + "/vibrant_env.yaml"
    threads:
        8
    resources:
        load = 20,
        mem = 4000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 60
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/vibrant/{MAG}-vibrant.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/vibrant/{MAG}-vibrant.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/vibrant/{MAG}-vibrant.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/vibrant/{MAG}-vibrant.qoutfile',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/vibrant/{MAG}-vibrant.qerrfile'
    shell:
        '''
        command="VIBRANT_run.py -i {params.genome} -folder {params.folder} -d {params.vibrant_db} -m {params.vibrant_files} -t {threads}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''

rule phigaro_genome:
    input:
        ancient(magpipe_config.flag_gene_calling)
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/phigaro/{MAG}-phigaro.done')
    params:
        phigaro = magpipe_config.dependencies["phigaro"],
        folder = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/phigaro/{MAG}-phigaro',
        genome = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-genomes/{MAG}.fa'
    threads:
        8
    resources:
        load = 20,
        mem = 3000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 60
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/phigaro/{MAG}-phigaro.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/phigaro/{MAG}-phigaro.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/phigaro/{MAG}-phigaro.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/phigaro/{MAG}-phigaro.qoutfile',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/phigaro/{MAG}-phigaro.qerrfile'
    shell:
        '''
        command="{params.phigaro} -f {params.genome} -o {params.folder} -e tsv html stdout -p -d -t {threads}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''

rule integron_finder_genome:
    input:
        ancient(magpipe_config.flag_gene_calling)
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/integron_finder/{MAG}-integron_finder.done')
    params:
        integron_finder = magpipe_config.dependencies["integron_finder"], 
        folder = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/integron_finder/{MAG}-integron_finder',
        genome = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-genomes/{MAG}.fa'
    conda:
        magpipe_config.envs_dir + "/integron_finder_env.yaml"
    threads:
        8
    resources:
        load = 20,
        mem = 3000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 60
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/integron_finder/{MAG}-integron_finder.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/integron_finder/{MAG}-integron_finder.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/integron_finder/{MAG}-integron_finder.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/integron_finder/{MAG}-integron_finder.qoutfile',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/integron_finder/{MAG}-integron_finder.qerrfile'
    shell:
        '''
        command="{params.integron_finder} {params.genome} --outdir {params.folder} --keep-tmp --local-max --promoter-attI --func-annot --cpu {threads}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''


rule isescan_genome:
    input:
        ancient(magpipe_config.flag_gene_calling)
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/isescan/{MAG}-isescan.done')
    params:
        outdir = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/isescan/',
        proteome = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/isescan/{MAG}-isescan-proteome',
        hmmres = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/isescan/{MAG}-isescan-hmm',
        genome = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-genomes/{MAG}.fa'
    conda:
        magpipe_config.envs_dir + "/isescan_env.yaml"
    threads:
        16
    resources:
        load = 33,
        mem = 3000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 60
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/isescan/{MAG}-isescan.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/isescan/{MAG}-isescan.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/isescan/{MAG}-isescan.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/isescan/{MAG}-isescan.qoutfile',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/isescan/{MAG}-isescan.qerrfile'
    shell:
        '''
        command="cd {params.outdir}; isescan.py {params.genome} {params.proteome} {params.hmmres} --nthread {threads}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''


# FIXME
#rule plasflow_genome:
#    input:
#        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
#    output:
#        marker = touch(magpipe_config.metag_dir + '/{dataset}/{sample}/plasflow/{sample}-plasflow.done'),
#        unzipped_assembly = temp(magpipe_config.metag_dir + "/{dataset}/{sample}/plasflow/{sample}.scaffolds.min1000.fasta")
#    params:
#        output_table = magpipe_config.metag_dir + "/{dataset}/{sample}/plasflow/{sample}-plasflow.tsv",
#        batch_size = 10000,
#        threshold = 0.7
#    conda:
#        magpipe_config.envs_dir + "/plasflow_env.yaml"
#    threads:
#        #1
#        2
#    resources:
#        load = 10,
#        mem = lambda wildcards, attempt: attempt * 96000, # the memory req is high and quite variable
#        scratch = 500,
#        time = lambda wildcards, attempt: attempt * 60
#    benchmark:
#        magpipe_config.metag_dir + "/{dataset}/{sample}/plasflow/{sample}-plasflow.benchmark"
#    log:
#        log = magpipe_config.metag_dir + "/{dataset}/{sample}/plasflow/{sample}-plasflow.log",
#        command = magpipe_config.metag_dir + "/{dataset}/{sample}/plasflow/{sample}-plasflow.command",
#        qoutfile = magpipe_config.metag_dir + "/{dataset}/{sample}/plasflow/{sample}-plasflow.qout",
#        qerrfile = magpipe_config.metag_dir + "/{dataset}/{sample}/plasflow/{sample}-plasflow.qerr"
#    shell:
#        """
#        command="export OMP_NUM_THREADS={threads}; export MKL_NUM_THREADS={threads}; export GOTO_NUM_THREADS={threads}; gunzip -c {input.assembly} > {output.unzipped_assembly}; PlasFlow.py --input {output.unzipped_assembly} --output {params.output_table} --threshold {params.threshold} --batch_size {params.batch_size} &> {log.log}";
#        echo "$command" > {log.command};
#        eval "$command"
#        """
#
#
#rule virsorter_genome:
#    input:
#        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
#    output:
#        marker = touch(magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}-virsorter.done"),
#        unzipped_assembly = temp(magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}.scaffolds.min1000.fasta")
#    params:
#        virsorter = magpipe_config.dependencies["virsorter"],
#        database = magpipe_config.dependencies["virsorter_database"],
#        output_directory = magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}-virsorter"
#    conda:
#        magpipe_config.envs_dir + "/virsorter_env.yaml"
#    threads:
#        8
#    resources:
#        load = 7,
#        mem = 3000,
#        scratch = 500,
#        time = lambda wildcards, attempt: attempt * 1400 # in min., just below 24h
#    benchmark:
#        magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}-virsorter.benchmark"
#    log:
#        log = magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}-virsorter.log",
#        command = magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}-virsorter.command",
#        qoutfile = magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}-virsorter.qout",
#        qerrfile = magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}-virsorter.qerr"
#    shell: # tested and it works, but if res is empty use --no_c option
#        """
#        command="gunzip -c {input.assembly} > {output.unzipped_assembly}; which perl; {params.virsorter} -f {output.unzipped_assembly} --db 2 --wdir {params.output_directory} --ncpu {threads} --data-dir {params.database} &> {log.log}";
#        echo "$command" > {log.command};
#        eval "$command"
#        """
#
#rule ccontig_genome:
#    input:
#        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
#    output:
#        marker = touch(magpipe_config.metag_dir + "/{dataset}/{sample}/ccontig/{sample}-ccontig.done")
#    params:
#        ccontig = magpipe_config.dependencies["ccontig"],
#        output_table = magpipe_config.metag_dir + "/{dataset}/{sample}/ccontig/{sample}-ccontig.tsv"
#    conda:
#        magpipe_config.envs_dir + "/ccontig_env.yaml"
#    threads:
#        1
#    resources:
#        load = 2,
#        mem = 10000,
#        scratch = 500,
#        time = lambda wildcards, attempt: attempt * 30
#    benchmark:
#        magpipe_config.metag_dir + "/{dataset}/{sample}/ccontig/{sample}-ccontig.benchmark"
#    log:
#        log = magpipe_config.metag_dir + "/{dataset}/{sample}/ccontig/{sample}-ccontig.log",
#        command = magpipe_config.metag_dir + "/{dataset}/{sample}/ccontig/{sample}-ccontig.command",
#        qoutfile = magpipe_config.metag_dir + "/{dataset}/{sample}/ccontig/{sample}-ccontig.qout",
#        qerrfile = magpipe_config.metag_dir + "/{dataset}/{sample}/ccontig/{sample}-ccontig.qerr"
#    shell:
#        """
#        command="{params.ccontig} -i <(gunzip -c {input.assembly}) -o {params.output_table} &> {log.log}";
#        echo "$command" > {log.command};
#        eval "$command"
#        """
#
#
#rule plasmidfinder_genome:
#    input:
#        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
#    output:
#        marker = touch(magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}-plasmidfinder.done"),
#        unzipped_assembly = temp(magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}.scaffolds.min1000.fasta")
#    params:
#        output_directory = magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}-plasmidfinder",
#        database = magpipe_config.dependencies["plasmidfinder_database"],
#        blastn = magpipe_config.dependencies["plasmidfinder_blastdb"],
#        coverage = 0.66,
#        identity = 0.5
#    conda:
#        magpipe_config.envs_dir + "/plasmidfinder_env.yaml"
#    threads:
#        1
#    resources:
#        load = 4,
#        mem = 16000,
#        scratch = 500,
#        time = lambda wildcards, attempt: attempt * 60
#    benchmark:
#        magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}-plasmidfinder.benchmark"
#    log:
#        log = magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}-plasmidfinder.log",
#        command = magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}-plasmidfinder.command",
#        qoutfile = magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}-plasmidfinder.qout",
#        qerrfile = magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}-plasmidfinder.qerr"
#    shell:
#        """
#        command="gunzip -c {input.assembly} > {output.unzipped_assembly}; mkdir -p {params.output_directory}; plasmidfinder.py -i {output.unzipped_assembly} -o {params.output_directory} -p {params.database} -mp {params.blastn} -l {params.coverage} -t {params.identity} --extented_output &> {log.log}";
#        echo "$command" > {log.command};
#        eval "$command"
#        """
#
#
#rule cbar_genome:
#    input:
#        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
#    output:
#        marker = touch(magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}-cbar.done"),
#        unzipped_assembly = temp(magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}.scaffolds.min1000.fasta")
#    params:
#        cbar = magpipe_config.dependencies["cbar"],
#        output_file = magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}-cbar.txt"
#    conda:
#        magpipe_config.envs_dir + "/cbar_env.yaml"
#    threads:
#        #1
#        2
#    resources:
#        load = 10,
#        mem = 64000, #cBar needs enough memory space to create the virtual java enviroment
#        scratch = 500,
#        time = lambda wildcards, attempt: attempt * 30
#    benchmark:
#        magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}-cbar.benchmark"
#    log:
#        log = magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}-cbar.log",
#        command = magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}-cbar.command",
#        qoutfile = magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}-cbar.qout",
#        qerrfile = magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}-cbar.qerr"
#    shell:
#        """
#        command="gunzip -c {input.assembly} > {output.unzipped_assembly}; {params.cbar} {output.unzipped_assembly} {params.output_file} &> {log.log}";
#        echo "$command" > {log.command};
#        eval "$command"
#        """
#
#
#rule deepvirfinder_genome:
#    input:
#        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
#    output:
#        marker = touch(magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}-deepvirfinder.done")
#    params:
#        deepvirfinder = magpipe_config.dependencies["deepvirfinder"],
#        output_directory = magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}-deepvirfinder"
#    conda:
#        magpipe_config.envs_dir + "/deepvirfinder_env.yaml"
#    threads:
#        10
#    resources:
#        load = 10,
#        mem = 8000,
#        scratch = 1000,
#        time = lambda wildcards, attempt: attempt * 1400 # in min., just below 24h
#    benchmark:
#        magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}-deepvirfinder.benchmark"
#    log:
#        log = magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}-deepvirfinder.log",
#        command = magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}-deepvirfinder.command",
#        qoutfile = magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}-deepvirfinder.qout",
#        qerrfile = magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}-deepvirfinder.qerr"
#    shell: # FIXME define THEANO_FLAGS to use gpu, i.e. device=cuda,force_device=True,floatX=float32 the deepvirfinder code has `os.environ['THEANO_FLAGS'] = "floatX=float32,openmp=True"` and `os.environ['THEANO_FLAGS'] = "mode=FAST_RUN,device=gpu0,floatX=float32"`
#        """
#        command="export OMP_NUM_THREADS={threads}; export MKL_NUM_THREADS={threads}; export GOTO_NUM_THREADS={threads}; python {params.deepvirfinder} -i <(gunzip -c {input.assembly}) -o {params.output_directory} -c 1 &> {log.log}";
#        echo "$command" > {log.command};
#        eval "$command"
#        """
#
#rule eukrep_genome:
#    input:
#        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
#    output:
#        marker = touch(magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}-eukrep.done"),
#        unzipped_assembly = temp(magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}.scaffolds.min1000.fasta"),
#        euk = magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}.euk",
#        prok = magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}.prok"
#    params:
#        None
#    conda:
#        magpipe_config.envs_dir + "/eukrep_env.yaml"
#    threads:
#        1
#    resources:
#        load = 2,
#        mem = 8000,
#        scratch = 500,
#        time = lambda wildcards, attempt: attempt * 60 # in min.
#    benchmark:
#        magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}-eukrep.benchmark"
#    log:
#        log = magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}-eukrep.log",
#        command = magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}-eukrep.command",
#        qoutfile = magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}-eukrep.qout",
#        qerrfile = magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}-eukrep.qerr"
#    shell: # using python $(which command) to overwrite the shebang that fails if you install inconda through pip as it's too long...
#        """
#        command="gunzip -c {input.assembly} > {output.unzipped_assembly}; python $(which EukRep) -i {output.unzipped_assembly} -o {output.euk} --prokarya {output.prok} --min 1000 --seq_names -m balanced --tie skip &> {log.log}";
#        echo "$command" > {log.command};
#        eval "$command"
#        """
