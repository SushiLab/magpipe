#!/usr/bin/env python

"""
name: annotate
description: This snakemake subpipeline is part of the MAGs building pipeline
dependency: prepare_phase_2
author: Lucas Paoli (paolil@ethz.ch)
rules:
    - GTDBtk [done] # Issue with qsub: faster to run out of the queue
    - prokka [To test]
    - fetch_MGs [TBA]
    - centrifuge [TBA]
"""

rule gtdbtk:
    input:
        marker = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-export_evaluation.done'
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-gtdbtk.done')
    params:
        gtdbtk = magpipe_config.dependencies["gtdbtk"],
        magpipe_scripts = magpipe_config.scripts,
        extension = 'fa',
        genomes = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-genomes',
        gtdbtk_output = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-gtdbtk',
        eval_table = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-evaluate_summary.tsv',
        processed_prefix = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-dictionaries/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn)
    conda: # FIXME effectively this step is ran manually, just run the rule as a dummy and use the .command file to run it manually in an appropriate conda environment
        magpipe_config.envs_dir + "/gtdbtk_env.yaml"
    threads:
        3
    resources:
        load = 100,
        mem = 200000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 7200 # In min, 7200 min = 120 hours
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-gtdbtk.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-gtdbtk.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-gtdbtk.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-gtdbtk.qout',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-gtdbtk.qerr'
    shell:
        '''
        command="{params.gtdbtk} classify_wf --genome_dir {params.genomes} --out_dir {params.gtdbtk_output} --extension {params.extension} --cpus {threads};
        {params.magpipe_scripts}/process_mags_taxonomy.py -e {params.eval_table} -g {params.gtdbtk_output} -o {params.processed_prefix}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''

rule prokka:
    input:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-gtdbtk.done'
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/prokka/{MAG}-{domain}-prokka.done')
    params:
        prokka = magpipe_config.dependencies["prokka"],
        magpipe_scripts = magpipe_config.scripts,
        mag = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-genomes/{MAG}.fa',
        outdir = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/prokka/{MAG}-{domain}-prokka'
    conda:
        magpipe_config.envs_dir + "/prokka_env.yaml"
    threads:
        4
    resources:
        load = 2,
        mem = 3000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 30
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/prokka/{MAG}-{domain}-prokka.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/prokka/{MAG}-{domain}-prokka.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/prokka/{MAG}-{domain}-prokka.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/prokka/{MAG}-{domain}-prokka.qoutfile',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/prokka/{MAG}-{domain}-prokka.qerrfile'
    shell:
        '''
        command="{params.prokka} {params.mag} --outdir {params.outdir} --prefix {wildcards.MAG}-prokka-original --kingdom {wildcards.domain} --genus {wildcards.domain} --species {wildcards.MAG} --strain unknown --cdsrnaolap --cpus {threads};
        {params.magpipe_scripts}/process_prokka.py -p {params.outdir}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''

rule fetch_mgs:
    input:
        ancient(dynamic('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/prokka/{MAG}-{domain}-prokka.done'))
        #ancient(magpipe_config.flag_gene_calling)
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/fetchMGs/{MAG}-fetchMGs.done')
    params:
        fetchmgs = magpipe_config.dependencies["fetchmgs"],
        lib = magpipe_config.dependencies["fetchmgs_lib"],
        genes = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/gene_call/{MAG}-prokka.ffn',
        proteins = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/gene_call/{MAG}-prokka.faa',
        bestmgs = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/fetchMGs/{MAG}-bestMGs',
        allmgs = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/fetchMGs/{MAG}-allMGs'
    conda:
        magpipe_config.envs_dir + "/fetchmgs_env.yaml"
    threads:
        1
    resources:
        load = 4,
        mem = 8000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 60
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/fetchMGs/{MAG}-fetchMGs.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/fetchMGs/{MAG}-fetchMGs.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/fetchMGs/{MAG}-fetchMGs.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/fetchMGs/{MAG}-fetchMGs.qoutfile',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/fetchMGs/{MAG}-fetchMGs.qerrfile'
    shell:
        '''
        command="perl {params.fetchmgs} -m extraction -l {params.lib} -v -i -d {params.genes} -o {params.bestmgs} {params.proteins} -t {threads};
        perl {params.fetchmgs} -m extraction -l {params.lib} -d {params.genes} -o {params.allmgs} {params.proteins} -t {threads}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''

rule kegg_genome:
    input:
        ancient(dynamic('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/prokka/{MAG}-{domain}-prokka.done'))
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/kegg/{MAG}-kegg.done')
    params:
        magpipe_scripts = magpipe_config.scripts,
        faa = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/gene_call/{MAG}-prokka.faa',
        outprefix = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/kegg/{MAG}-kegg',
        kegg_db = "/science/paolil/databases/kegg/kegg_for_prokka-with_ko.dmnd",
        kegg_dict = "/science/paolil/databases/kegg/ko/"
        #kegg_db = "/cluster/scratch/paolil/db/kegg/kegg_for_prokka-with_ko.dmnd",
        #kegg_dict = "/cluster/scratch/paolil/db/kegg/ko/"
    conda:
        magpipe_config.envs_dir + "/kegg_env.yaml"
    threads:
        4
    resources:
        load = 2,
        mem = 3000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 60
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/kegg/{MAG}-kegg.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/kegg/{MAG}-kegg.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/kegg/{MAG}-kegg.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/kegg/{MAG}-kegg.qoutfile',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/kegg/{MAG}-kegg.qerrfile'
    shell:
        '''
        command="{params.magpipe_scripts}/annotate_kegg.py --query {params.faa} -o {params.outprefix} -d {params.kegg_db} -k {params.kegg_dict} -t {threads}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''

rule eggnog_genome:
    input:
        #ancient(magpipe_config.flag_gene_calling)
        ancient(dynamic('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/prokka/{MAG}-{domain}-prokka.done'))
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/eggnog/{MAG}-eggnog.done')
    params:
        magpipe_scripts = magpipe_config.scripts,
        faa = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/gene_call/{MAG}-prokka.faa',
        outprefix = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/eggnog/{MAG}-eggnog',
        eggnog_data = "/science/paolil/databases/eggnog/",
        eggnog_db = "/science/paolil/databases/eggnog/eggnog_proteins.dmnd"
        #eggnog_data = "/cluster/scratch/paolil/db/eggnog/",
        #eggnog_db = "/cluster/scratch/paolil/db/eggnog/eggnog_proteins.dmnd"
        #eggnog_data = "/cluster/work/biol/tmp/paolil/db/eggnog/",
        #eggnog_db = "/cluster/work/biol/tmp/paolil/db/eggnog/eggnog_proteins.dmnd"
    conda:
        magpipe_config.envs_dir + "/eggnog_env.yaml"
    threads:
        #8 
        16 # to finish the last jobs 
    resources:
        load = 4,
        mem = 3000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 180
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/eggnog/{MAG}-eggnog.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/eggnog/{MAG}-eggnog.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/eggnog/{MAG}-eggnog.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/eggnog/{MAG}-eggnog.qoutfile',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/eggnog/{MAG}-eggnog.qerrfile'
    shell:
        '''
        command="emapper.py -i {params.faa} --output {params.outprefix} -m diamond --data_dir {params.eggnog_data} --dmnd_db {params.eggnog_db} --query-cover 70 --subject-cover 70 --cpu {threads}";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''

rule antismash_genome:
    input:
        #ancient(magpipe_config.flag_gene_calling)
        ancient(dynamic('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/prokka/{MAG}-{domain}-prokka.done'))
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/antismash/{MAG}-antismash.done')
    params:
        antismash = magpipe_config.dependencies["antismash"],
        mag = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-genomes/{MAG}.fa',
        gff = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/prokka/{MAG}-*-prokka/{MAG}-prokka.gff',
        outdir = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/antismash/{MAG}-antismash'
    #conda:
    #    magpipe_config.envs_dir + "/antismash_env.yaml" # Currently using the module system
    threads:
        4 # Normally use 4, use 12 for the large ones
    resources:
        load = 4,
        mem = 6000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 180 # Was 45, 180 for last ones
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/antismash/{MAG}-antismash.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/antismash/{MAG}-antismash.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/antismash/{MAG}-antismash.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/antismash/{MAG}-antismash.qout',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/antismash/{MAG}-antismash.qerr'
    shell:
        '''
        command="{params.antismash} {params.mag} --output-dir {params.outdir} --genefinding-gff3 {params.gff} --genefinding-tool none --html-description '{wildcards.MAG}' --cb-general --cb-knownclusters --cb-subclusters --cpus {threads}";
        # FIXME Add process_antismash script
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''

rule traitar_genome:
    input:
        ancient(dynamic('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/prokka/{MAG}-{domain}-prokka.done'))
    output:
        outdir = directory('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/traitar/{MAG}-traitar'),
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/traitar/{MAG}-traitar.done')
    params:
        traitar_cmd = magpipe_config.dependencies["traitar"],
        mag_dir = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-genomes/',
        traitar_workdir = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/traitar/'
    conda:
        magpipe_config.envs_dir + "/traitar_env.yaml"
    threads:
        16 
    resources:
        load = 20,
        mem = 3000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 180 # Was 45, 180 for last ones
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/traitar/{MAG}-traitar.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/traitar/{MAG}-traitar.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/traitar/{MAG}-traitar.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/traitar/{MAG}-traitar.qout',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/traitar/{MAG}-traitar.qerr'
    shell:
        '''
        command="
if [ ! -d {params.traitar_workdir}/traitar_input ]; then
    mkdir {params.traitar_workdir}/traitar_input
    ln -s {params.mag_dir}/{wildcards.MAG}.fa {params.traitar_workdir}/traitar_input/{wildcards.MAG}.fna
fi
if [ ! -f {params.traitar_workdir}/traitar_sample_file.txt ]; then
    echo $'sample_file_name\tsample_name\tcategory' > {params.traitar_workdir}/traitar_sample_file.txt
    echo $'{wildcards.MAG}.fna\t{wildcards.MAG}\teremio_superproducer' >> {params.traitar_workdir}/traitar_sample_file.txt
fi
{params.traitar_cmd} phenotype {params.traitar_workdir}/traitar_input {params.traitar_workdir}/traitar_sample_file.txt from_nucleotides {output.outdir} -c {threads}
        ";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''

rule growthpred_genome:
    input:
        ancient(dynamic('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/prokka/{MAG}-{domain}-prokka.done'))
        #ancient(magpipe_config.flag_gene_calling)
    output:
        marker = touch('{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/growthpred/{MAG}-growthpred.done')
    params:
        growthpred_cmd = magpipe_config.dependencies["growthpred"],
        growthpred_workdir = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/growthpred/',
        faa = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/gene_call/{MAG}-prokka.faa',
        ffn = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/gene_call/{MAG}-prokka.ffn',
        outdir = '{MAG}-growthpred-output'
    conda:
        magpipe_config.envs_dir + "/growthpred_env.yaml"
    threads:
        1 
    resources:
        load = 1,
        mem = 8000,
        scratch = 500,
        time = lambda wildcards, attempt: attempt * 180 # Was 45, 180 for last ones
    benchmark:
        '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/growthpred/{MAG}-growthpred.benchmark'
    log:
        log = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/growthpred/{MAG}-growthpred.log',
        command = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/growthpred/{MAG}-growthpred.command',
        qoutfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/growthpred/{MAG}-growthpred.qout',
        qerrfile = '{path}' + '/' + magpipe_config.pooled_name + '-{binning_method}-cpl' + str(magpipe_config.cpl) + '_ctn' + str(magpipe_config.ctn) + '-annotations/{MAG}/growthpred/{MAG}-growthpred.qerr'
    shell:
        '''
        command="
cd {params.growthpred_workdir};
seqtk subseq {params.ffn} <(seqtk comp {params.faa} | cut -f1) > genes_for_growthpred.fna;
python {params.growthpred_cmd} -S -c 0 -o {params.outdir} -t -b -r -g genes_for_growthpred.fna;
        ";
        echo "$command" > {log.command};
        eval "$command" &> {log.log}
        '''

