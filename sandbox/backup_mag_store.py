#!/usr/bin/env python

import sys
from pathlib import Path

mag_store = Path("/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/MAGPIPE_MAGS_EAN/scratch/processed")
store_backup = Path("/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/MAGPIPE_MAGS_EAN/data/mag-store-backup")
method = "metabat2_c2k_e500-cpl50_ctn10" 
target_dir = "MAGs/INTERNAL"

if len(sys.argv) == 2:
    target_projects = [sys.argv[1]]
else:
    #target_projects = [f'{i.parent.name}/{i.name}' for i in mag_store.joinpath("MAGs").glob("*/*") if not i.name == "archive"]
    target_projects = [i.name for i in mag_store.joinpath(target_dir).glob("*") if not i.name == "archive"]

print(target_projects)

for project_name in target_projects: 
    project_path = mag_store.joinpath(f'{target_dir}/{project_name}')
    if not project_path.exists():
        raise Exception(f'Couldn\'t find {project_path}...')
    print(f'Processing project {project_name}...')
    project_backup = store_backup.joinpath(f'{target_dir}/{project_name}')
    project_backup.mkdir(parents = True, exist_ok = True)
    project_script = project_backup.joinpath(f'{project_name}-{method}-backup.sh') 
    print(f'Writing backup script: {project_script}')
    with open(project_script, 'w') as handle:
        handle.write('#!/bin/bash\nset -euo pipefail')
        # archive and export genomes
        handle.write(f'\n# Archive and backup genomes')
        genomes_path = project_path.joinpath(f'{project_name}-{method}-genomes')
        if not genomes_path.exists():
            raise Exception(f'Couldn\'t find {genomes_path}...')
        genomes_archive = project_path.joinpath(f'{project_name}-{method}-genomes.tar.gz')
        handle.write(f'\ntar -czf {genomes_archive} -C {genomes_path.parent} {genomes_path.name}')
        genomes_backup = project_backup.joinpath(f'{project_name}-{method}-genomes.tar.gz')
        handle.write(f'\nrsync -av {genomes_archive} {genomes_backup}')
        # export evaluate files
        handle.write(f'\n# Backup evaluation')
        for i in project_path.glob(f"{project_name}-{method}*eval*"):
            handle.write(f'\nrsync -av {project_path.joinpath(i)} {project_backup.joinpath(i.name)}')
        # if drep, export relevant drep files
        drep_folder = project_path.joinpath(f'{project_name}-{method}-drep_0.95', 'data_tables')
        if drep_folder.exists():
            handle.write("\n# Backup dRep")
            drep_backup = project_backup.joinpath(f'{project_name}-{method}-drep_0.95', 'data_tables')
            handle.write(f'\nmkdir -p {drep_backup}')
            for i in drep_folder.glob('*'):
                handle.write(f'\nrsync -av {i} {drep_backup.joinpath(i.name)}')
        # if gtdb, export relevant gtdb files
        gtdb_folder = project_path.joinpath(f'{project_name}-{method}-gtdbtk')
        if gtdb_folder.exists():
            handle.write("\n# Backup GTDB")
            gtdb_backup = project_backup.joinpath(f'{project_name}-{method}-gtdbtk')
            handle.write(f'\nmkdir {gtdb_backup}')
            gtdb_bacteria = "gtdbtk.bac120.summary.tsv"
            gtdb_archaea = "gtdbtk.ar122.summary.tsv"
            gtdb_log = "gtdbtk.log"
            handle.write(f'\nrsync -Lv {gtdb_folder.joinpath(gtdb_bacteria)} {gtdb_backup.joinpath(gtdb_bacteria)}')
            handle.write(f'\nrsync -Lv {gtdb_folder.joinpath(gtdb_archaea)} {gtdb_backup.joinpath(gtdb_archaea)}')
            handle.write(f'\nrsync -av {gtdb_folder.joinpath(gtdb_log)} {gtdb_backup.joinpath(gtdb_log)}')
        done_file = project_backup.joinpath(f'{project_name}-{method}-backup.done') 
        handle.write(f'\ntouch {done_file}\n')
