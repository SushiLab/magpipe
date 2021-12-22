#!/usr/bin/env python

"""
name: Deprecated
description: Old rules that were written at some point. Might come handy.
dependency: None
author: Lucas Paoli (paolil@ethz.ch)
"""

# Note: If unitem fails, manually check and fix the shebang to #!/usr/bin/env python
rule compare:
    input:
        assembly = lambda wildcards: INPUT.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz',
        bins_1 = magpipe_config.metag_dir + '/{dataset}/{sample}/{method1}/bins/',
        bins_2 = magpipe_config.metag_dir + '/{dataset}/{sample}/{method2}/bins/'
    output:
        marker = touch(magpipe_config.metag_dir + '/{dataset}/{sample}/compare_binning_methods/{sample}_{method1}_vs_{method2}.compare.done')
    params:
        output = magpipe_config.metag_dir + '/{dataset}/{sample}/compare_binning_methods/{sample}_{method1}_vs_{method2}.compare.output.tsv'
    conda:
        magpipe_config.envs_dir + "/unitem_env.yaml"
    threads:
        1
    resources:
        load = 5,
        mem = 4000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 30 # in minutes
    benchmark:
        magpipe_config.metag_dir + '/{dataset}/{sample}/compare_binning_methods/{sample}_{method1}_vs_{method2}.compare.benchmark'
    log:
        log = magpipe_config.metag_dir + '/{dataset}/{sample}/compare_binning_methods/{sample}_{method1}_vs_{method2}.compare.log',
        command = magpipe_config.metag_dir + '/{dataset}/{sample}/compare_binning_methods/{sample}_{method1}_vs_{method2}.compare.command',
        qoutfile = magpipe_config.metag_dir + '/{dataset}/{sample}/compare_binning_methods/{sample}_{method1}_vs_{method2}.compare.qout',
        qerrfile = magpipe_config.metag_dir + '/{dataset}/{sample}/compare_binning_methods/{sample}_{method1}_vs_{method2}.compare.qerr'
    shell:
        '''
        command="unitem compare {input.assembly} {input.bins_1} {input.bins_2} {params.output} -x fa -y fa &> {log.log}";
        echo "$command" > {log.command};
        eval "$command"
        '''

rule plots:
    input:
        marker = magpipe_config.metag_dir  + '/{dataset}/{sample}/{binning_method}/{sample}_{binning_method}.checkm.done',
        assembly = lambda wildcards: INPUT.loc[wildcards.dataset, "assemblies"] + '/' + wildcards.sample + '.scaffolds.min1000.fasta.gz',
        bam = lambda wildcards: INPUT.loc[wildcards.dataset, "bamfiles"] + '/' + wildcards.sample + '/'
    output:
        marker = touch(magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}_{binning_method}.plots.done')
    params:
        checkm_lineage = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/checkm_lineage/',
        bins = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/bins/',
        plots = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/checkm_plots/',
        tetra = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}_{binning_method}.checkm.tetra.tsv',
        coverage = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}_{binning_method}.checkm.coverage.tsv',
        profile = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}_{binning_method}.checkm.profile.tsv'
    conda:
        magpipe_config.envs_dir + "/checkm_env.yaml"
    threads:
        4
    resources:
        load = 5,
        mem = 16000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 180 # in minutes
    benchmark:
        magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}_{binning_method}.plots.benchmark'
    log:
        log = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}_{binning_method}.plots.log',
        command = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}_{binning_method}.plots.command',
        qoutfile = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}_{binning_method}.plots.qout',
        qerrfile = magpipe_config.metag_dir + '/{dataset}/{sample}/{binning_method}/{sample}_{binning_method}.plots.qerr'
    shell:
        '''
        command="checkm tetra <(gunzip -c {input.assembly}) {params.tetra} &> {log.log};
        #checkm coverage -x fa -t {threads} {params.bins} {params.coverage} {input.bam}/*.bam &>> {log.log};
        #checkm profile {params.coverage} --tab_table -f {params.profile} &>> {log.log};
        checkm bin_qa_plot -x fa --image_type pdf {params.checkm_lineage} {params.bins} {params.plots} &>> {log.log};
        checkm dist_plot -x fa --image_type pdf {params.checkm_lineage} {params.bins} {params.plots} {params.tetra} 95 &>> {log.log};
        checkm nx_plot -x fa --image_type pdf {params.bins} {params.plots} &>> {log.log};
        checkm len_plot -x fa --image_type pdf {params.bins} {params.plots} &>> {log.log};
        checkm len_hist -x fa --image_type pdf {params.bins} {params.plots} &>> {log.log};
        checkm marker_plot -x fa --image_type pdf {params.checkm_lineage} {params.bins} {params.plots} &>> {log.log};
        #checkm par_plot -x fa --image_type pdf {params.checkm_lineage} {params.bins} {params.plots} {params.coverage} &>> {log.log}";
        echo "$command" > {log.command};
        eval "$command"
        '''
