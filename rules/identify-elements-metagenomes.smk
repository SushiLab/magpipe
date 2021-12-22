#!/usr/bin/env python

"""
name: plasmids
description: rules to annotate plasmids in metagenomic assemblies
dependency: Primary analyses
author: Daniel Gehrig (gehrigd@biol.ethz.ch) and Lucas Paoli (paolil@ethz.ch)
rules:
    - plasflow
    - virsorter
    - ccontig
    - plasmidfinder
    - cbar
    - deepvirfinder
    - eukrep
"""

rule plasflow:
    input:
        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
    output:
        marker = touch(magpipe_config.metag_dir + '/{dataset}/{sample}/plasflow/{sample}-plasflow.done'),
        unzipped_assembly = temp(magpipe_config.metag_dir + "/{dataset}/{sample}/plasflow/{sample}.scaffolds.min1000.fasta")
    params:
        output_table = magpipe_config.metag_dir + "/{dataset}/{sample}/plasflow/{sample}-plasflow.tsv",
        batch_size = 10000,
        threshold = 0.7
    conda:
        magpipe_config.envs_dir + "/plasflow_env.yaml"
    threads:
        #1
        2
    resources:
        load = 10,
        mem = lambda wildcards, attempt: attempt * 96000, # the memory req is high and quite variable
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 60
    benchmark:
        magpipe_config.metag_dir + "/{dataset}/{sample}/plasflow/{sample}-plasflow.benchmark"
    log:
        log = magpipe_config.metag_dir + "/{dataset}/{sample}/plasflow/{sample}-plasflow.log",
        command = magpipe_config.metag_dir + "/{dataset}/{sample}/plasflow/{sample}-plasflow.command",
        qoutfile = magpipe_config.metag_dir + "/{dataset}/{sample}/plasflow/{sample}-plasflow.qout",
        qerrfile = magpipe_config.metag_dir + "/{dataset}/{sample}/plasflow/{sample}-plasflow.qerr"
    shell:
        """
        command="export OMP_NUM_THREADS={threads}; export MKL_NUM_THREADS={threads}; export GOTO_NUM_THREADS={threads}; gunzip -c {input.assembly} > {output.unzipped_assembly}; PlasFlow.py --input {output.unzipped_assembly} --output {params.output_table} --threshold {params.threshold} --batch_size {params.batch_size} &> {log.log}";
        echo "$command" > {log.command};
        eval "$command"
        """


rule virsorter:
    input:
        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
    output:
        marker = touch(magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}-virsorter.done"),
        unzipped_assembly = temp(magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}.scaffolds.min1000.fasta")
    params:
        virsorter = magpipe_config.dependencies["virsorter"],
        database = magpipe_config.dependencies["virsorter_database"],
        output_directory = magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}-virsorter"
    conda:
        magpipe_config.envs_dir + "/virsorter_env.yaml"
    threads:
        8
    resources:
        load = 7,
        mem = 3000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 1400 # in min., just below 24h
    benchmark:
        magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}-virsorter.benchmark"
    log:
        log = magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}-virsorter.log",
        command = magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}-virsorter.command",
        qoutfile = magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}-virsorter.qout",
        qerrfile = magpipe_config.metag_dir + "/{dataset}/{sample}/virsorter/{sample}-virsorter.qerr"
    shell: # tested and it works, but if res is empty use --no_c option
        """
        command="gunzip -c {input.assembly} > {output.unzipped_assembly}; which perl; {params.virsorter} -f {output.unzipped_assembly} --db 2 --wdir {params.output_directory} --ncpu {threads} --data-dir {params.database} &> {log.log}";
        echo "$command" > {log.command};
        eval "$command"
        """

rule ccontig:
    input:
        #assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.contigs.min1000.fasta.gz'
    output:
        marker = touch(magpipe_config.metag_dir + "/{dataset}/{sample}/ccontig/{sample}-ccontig.done")
    params:
        ccontig = magpipe_config.dependencies["ccontig"],
        output_table = magpipe_config.metag_dir + "/{dataset}/{sample}/ccontig/{sample}-ccontig.tsv"
    conda:
        magpipe_config.envs_dir + "/ccontig_env.yaml"
    threads:
        1
    resources:
        load = 2,
        mem = 10000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 30
    benchmark:
        magpipe_config.metag_dir + "/{dataset}/{sample}/ccontig/{sample}-ccontig.benchmark"
    log:
        log = magpipe_config.metag_dir + "/{dataset}/{sample}/ccontig/{sample}-ccontig.log",
        command = magpipe_config.metag_dir + "/{dataset}/{sample}/ccontig/{sample}-ccontig.command",
        qoutfile = magpipe_config.metag_dir + "/{dataset}/{sample}/ccontig/{sample}-ccontig.qout",
        qerrfile = magpipe_config.metag_dir + "/{dataset}/{sample}/ccontig/{sample}-ccontig.qerr"
    shell:
        """
        command="{params.ccontig} -i <(gunzip -c {input.assembly}) -o {params.output_table} &> {log.log}";
        echo "$command" > {log.command};
        eval "$command"
        """


rule plasmidfinder:
    input:
        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
    output:
        marker = touch(magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}-plasmidfinder.done"),
        unzipped_assembly = temp(magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}.scaffolds.min1000.fasta")
    params:
        output_directory = magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}-plasmidfinder",
        database = magpipe_config.dependencies["plasmidfinder_database"],
        blastn = magpipe_config.dependencies["plasmidfinder_blastdb"],
        coverage = 0.66,
        identity = 0.5
    conda:
        magpipe_config.envs_dir + "/plasmidfinder_env.yaml"
    threads:
        1
    resources:
        load = 4,
        mem = 16000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 60
    benchmark:
        magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}-plasmidfinder.benchmark"
    log:
        log = magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}-plasmidfinder.log",
        command = magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}-plasmidfinder.command",
        qoutfile = magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}-plasmidfinder.qout",
        qerrfile = magpipe_config.metag_dir + "/{dataset}/{sample}/plasmidfinder/{sample}-plasmidfinder.qerr"
    shell:
        """
        command="gunzip -c {input.assembly} > {output.unzipped_assembly}; mkdir -p {params.output_directory}; plasmidfinder.py -i {output.unzipped_assembly} -o {params.output_directory} -p {params.database} -mp {params.blastn} -l {params.coverage} -t {params.identity} --extented_output &> {log.log}";
        echo "$command" > {log.command};
        eval "$command"
        """


rule cbar:
    input:
        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
    output:
        marker = touch(magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}-cbar.done"),
        unzipped_assembly = temp(magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}.scaffolds.min1000.fasta")
    params:
        cbar = magpipe_config.dependencies["cbar"],
        output_file = magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}-cbar.txt"
    conda:
        magpipe_config.envs_dir + "/cbar_env.yaml"
    threads:
        #1
        2
    resources:
        load = 10,
        mem = 64000, #cBar needs enough memory space to create the virtual java enviroment
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 30
    benchmark:
        magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}-cbar.benchmark"
    log:
        log = magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}-cbar.log",
        command = magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}-cbar.command",
        qoutfile = magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}-cbar.qout",
        qerrfile = magpipe_config.metag_dir + "/{dataset}/{sample}/cbar/{sample}-cbar.qerr"
    shell:
        """
        command="gunzip -c {input.assembly} > {output.unzipped_assembly}; {params.cbar} {output.unzipped_assembly} {params.output_file} &> {log.log}";
        echo "$command" > {log.command};
        eval "$command"
        """


rule deepvirfinder:
    input:
        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
    output:
        marker = touch(magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}-deepvirfinder.done"),
        unzipped_assembly = temp(magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}.scaffolds.min1000.fasta"),
        filtered_assembly = temp(magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}.scaffolds.filtered.fasta")
    params:
        deepvirfinder = magpipe_config.dependencies["deepvirfinder"],
        output_directory = magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}-deepvirfinder"
    conda:
        magpipe_config.envs_dir + "/deepvirfinder_env.yaml"
    threads:
        60
    resources:
        load = 10,
        mem = 10000,
        scratch = 1000,
        time = lambda wildcards, attempt: attempt * 1400 # in min., just below 24h
    benchmark:
        magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}-deepvirfinder.benchmark"
    log:
        log = magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}-deepvirfinder.log",
        command = magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}-deepvirfinder.command",
        qoutfile = magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}-deepvirfinder.qout",
        qerrfile = magpipe_config.metag_dir + "/{dataset}/{sample}/deepvirfinder/{sample}-deepvirfinder.qerr"
    shell: # FIXME define THEANO_FLAGS to use gpu, i.e. device=cuda,force_device=True,floatX=float32 the deepvirfinder code has `os.environ['THEANO_FLAGS'] = "floatX=float32,openmp=True"` and `os.environ['THEANO_FLAGS'] = "mode=FAST_RUN,device=gpu0,floatX=float32"`
        '''
        gunzip -c {input.assembly} > {output.unzipped_assembly};
        seqtk subseq {output.unzipped_assembly} <(seqtk comp {output.unzipped_assembly} | awk \'{{if ($2 <= 2100000) {{print $1}} }}\' | cut -f1 -d: ) > {output.filtered_assembly}; # FIXME this is not to have to reinstall the env
        command="export OMP_NUM_THREADS={threads}; export MKL_NUM_THREADS={threads}; export GOTO_NUM_THREADS={threads};
        python {params.deepvirfinder} -i {output.filtered_assembly} -o {params.output_directory} -c 1 &> {log.log}";
        echo "$command" > {log.command};
        eval "$command"
        '''

rule eukrep:
    input:
        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
    output:
        marker = touch(magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}-eukrep.done"),
        unzipped_assembly = temp(magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}.scaffolds.min1000.fasta"),
        euk = magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}.euk",
        prok = magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}.prok"
    params:
        None
    conda:
        magpipe_config.envs_dir + "/eukrep_env.yaml"
    threads:
        1
    resources:
        load = 2,
        mem = 8000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 60 # in min.
    benchmark:
        magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}-eukrep.benchmark"
    log:
        log = magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}-eukrep.log",
        command = magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}-eukrep.command",
        qoutfile = magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}-eukrep.qout",
        qerrfile = magpipe_config.metag_dir + "/{dataset}/{sample}/eukrep/{sample}-eukrep.qerr"
    shell: # using python $(which command) to overwrite the shebang that fails if you install inconda through pip as it's too long...
        """
        command="gunzip -c {input.assembly} > {output.unzipped_assembly}; python $(which EukRep) -i {output.unzipped_assembly} -o {output.euk} --prokarya {output.prok} --min 1000 --seq_names -m balanced --tie skip &> {log.log}";
        echo "$command" > {log.command};
        eval "$command"
        """
