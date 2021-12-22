#!/usr/bin/env python

"""
name: binning_snake
description: This snakemake subpipeline is part of the MAGs building pipeline
dependency: depth_snake
author: Lucas Paoli (paolil@ethz.ch)
rules:
    - metabat2_default [done]
    - metabat2_sensitive [done]
    - vamb [needs testing]
    - maxbin2 [TBA]
    - metaWRAP [TBA]
    - DAStool [TBA?]
"""


def should_we_compute_depth(wildcards, magpipe_config = magpipe_config):
    """
    A simple helper function to know wheter we should skip the depth step or not
    """
    backmapping = magpipe_config.input.loc[wildcards.dataset, "backmapping"]
    if "depth" in backmapping:
        return [] 
    else:
        return ancient(magpipe_config.metag_dir + '/' + '/'.join([wildcards.dataset, wildcards.sample, 'depth_files', wildcards.sample + wildcards.depth_sfx + '-depth.done']))


def where_is_depthfile(wildcards, magpipe_config = magpipe_config):
    """
    A simple helper function to know where is the depthfile 
    """
    backmapping = magpipe_config.input.loc[wildcards.dataset, "backmapping"]
    if "depth" in backmapping:
        return backmapping + '/' + wildcards.sample + '.depth'
    else:
        return magpipe_config.metag_dir + '/' + '/'.join([wildcards.dataset, wildcards.sample, 'depth_files', wildcards.sample + wildcards.depth_sfx + '-depth.tsv'])


rule metabat2:
    input:
        marker = should_we_compute_depth,
        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz'
        #assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.contigs.min1000.fasta.gz'
    output:
        tmp_assembly = temp(magpipe_config.metag_dir + '/{dataset}/{sample}/metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '/metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '.tmp_assembly.fasta'),
        seqtk_comp = temp(magpipe_config.metag_dir + '/{dataset}/{sample}/metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '/metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '.tmp_assembly.seqtkcomp'),
        marker = touch(magpipe_config.metag_dir + '/{dataset}/{sample}/metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '/{sample}-metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '-bins.done')
    params:
        #bins_dir = magpipe_config.metag_dir + '/{dataset}/{sample}/metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '/raw_bins/',
        #prefix = magpipe_config.metag_dir + '/{dataset}/{sample}/metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '/raw_bins/{sample}-metabat2{depth_sfx}' + magpipe_config.methods.suffix[0],
        bins_dir = magpipe_config.metag_dir + '/{dataset}/{sample}/metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '/bins/',
        prefix = magpipe_config.metag_dir + '/{dataset}/{sample}/metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '/bins/{sample}-metabat2{depth_sfx}' + magpipe_config.methods.suffix[0],
        depthfile = where_is_depthfile,
        scripts = magpipe_config.scripts,
        min_bin_length = magpipe_config.methods.min_bin_size[0],
        min_contig_size = magpipe_config.methods.min_contig_size[0],
        max_edges = magpipe_config.methods.max_edges[0],
        min_cv = magpipe_config.methods.min_cv[0]
    conda:
        magpipe_config.envs_dir + "/metabat2_env.yaml"
    threads:
        16
    resources:
        load = 20,
        mem = 2000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 60 # in minutes
    benchmark:
        magpipe_config.metag_dir + '/{dataset}/{sample}/metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '/{sample}-metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '-metabat2.benchmark'
    log:
        log = magpipe_config.metag_dir + '/{dataset}/{sample}/metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '/{sample}-metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '-metabat2.log',
        command = magpipe_config.metag_dir + '/{dataset}/{sample}/metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '/{sample}-metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '-metabat2.command',
        qoutfile = magpipe_config.metag_dir + '/{dataset}/{sample}/metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '/{sample}-metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '-metabat2.qout',
        qerrfile = magpipe_config.metag_dir + '/{dataset}/{sample}/metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '/{sample}-metabat2{depth_sfx}' + magpipe_config.methods.suffix[0] + '-metabat2.qerr'
    shell:
        '''
        command="
        gunzip -c {input.assembly} > {output.tmp_assembly};
        metabat2 -i {output.tmp_assembly} -a {params.depthfile} -o {params.prefix} --minContig {params.min_contig_size} --maxEdges {params.max_edges} -x {params.min_cv} --numThreads {threads} --minClsSize {params.min_bin_length} --saveCls -v &> {log.log};
        seqtk comp {output.tmp_assembly} > {output.seqtk_comp};
        #{params.scripts}/process_metabat2.py -m {params.bins_dir} -c {output.seqtk_comp} -a {output.tmp_assembly} -l {params.min_bin_length}
        ";
        echo "$command" > {log.command};
        eval "$command" &>> {log.log};
        '''

# FIXME: to be tested with the exports, new environmnent etc
rule vamb:
    input:
        assembly = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz',
        bamfiles = lambda wildcards: magpipe_config.input.loc[wildcards.dataset, "backmapping"] + '/' + wildcards.sample + '/'
    output:
        marker = touch(magpipe_config.metag_dir + '/{dataset}/{sample}/vamb{depth_sfx}_c2k/{sample}_vamb{depth_sfx}_c2k-bins.done'),
        fastassembly = temp(magpipe_config.metag_dir + '/{dataset}/{sample}/vamb{depth_sfx}_c2k/{sample}.assembly.fa')
    params:
        vamb = magpipe_config.dependencies["vamb"],
        scripts = magpipe_config.scripts,
        vamb_out = magpipe_config.metag_dir + '/{dataset}/{sample}/vamb{depth_sfx}_c2k/vamb_out/',
        min_contig_size = 2000,
        cluster_file = magpipe_config.metag_dir + '/{dataset}/{sample}/vamb{depth_sfx}_c2k/vamb_out/clusters.tsv',
        bins_out = magpipe_config.metag_dir + '/{dataset}/{sample}/vamb{depth_sfx}_c2k/bins/',
        bins_prefix = '{sample}_vamb{depth_sfx}_c2k'
    conda:
        magpipe_config.envs_dir + "/vamb_env.yaml"
    threads:
        32
    resources:
        load = 20,
        mem = 16000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 500 # in minutes
    benchmark:
        magpipe_config.metag_dir + '/{dataset}/{sample}/vamb{depth_sfx}_c2k/{sample}_vamb{depth_sfx}_c2k-vamb.benchmark'
    log:
        log = magpipe_config.metag_dir + '/{dataset}/{sample}/vamb{depth_sfx}_c2k/{sample}_vamb{depth_sfx}_c2k-vamb.log',
        command = magpipe_config.metag_dir + '/{dataset}/{sample}/vamb{depth_sfx}_c2k/{sample}_vamb{depth_sfx}_c2k-vamb.command',
        qoutfile = magpipe_config.metag_dir + '/{dataset}/{sample}/vamb{depth_sfx}_c2k/{sample}_vamb{depth_sfx}_c2k-vamb.qout',
        qerrfile = magpipe_config.metag_dir + '/{dataset}/{sample}/vamb{depth_sfx}_c2k/{sample}_vamb{depth_sfx}_c2k-vamb.qerr'
    shell:
        '''
        command="
        gunzip -c {input.assembly} > {output.fastassembly};
        export OMP_NUM_THREADS={threads}; export MKL_NUM_THREADS={threads}; export GOTO_NUM_THREADS={threads}; {params.vamb} {params.outdir} {output.fastassembly} {input.bamfiles}/*.bam -m {params.min_contig_size} -p {threads} &> {log.log};
        python {params.scripts}/process_vamb.py -c {params.cluster_file} -a {input.assembly} -o {params.outdir} -p {params.prefix} -m {params.min_bin_size} &>> {log.log}";
        echo "$command" > {log.command};
        eval "$command"
        '''
