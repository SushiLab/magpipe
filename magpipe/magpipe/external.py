#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Functions associated with the processing of external genomes
"""

from pathlib import Path

def checkm_sized_genome_list(genome_list):
    """
    Yield lists of max 200 genomes for checkm to handle
    """
    for i in range(0, len(genome_list), 200):
        yield genome_list[i:i + 200]

def setup_ext_evaluation(dataset_dir, processed_genomes = 'processed_genomes'):
    """
    Takes an external dataset and processed_genomes folder as input
    Creates a checkm folder for that dataset with subfolders linking
    to max 200 genomes each
    Creates an anvio folder
    """
    dataset_path = Path(dataset_dir)
    dataset = dataset_path.name
    print('Preparing dataset {} for evaluation. Expects processed genomes in `{}`'.format(dataset, processed_genomes))
    checkm_path = dataset_path.joinpath('checkm')
    anvio_path = dataset_path.joinpath('anvio_db')
    if checkm_path.exists() or anvio_path.exists():
        print('Looks like this script has been ran already, skipping.')
        return
    anvio_path.mkdir()
    processed_path = dataset_path.joinpath(processed_genomes)
    genome_list = [i.name for i in processed_path.glob('*.fa')]
    chunk = 0
    for genome_chunk in checkm_sized_genome_list(genome_list):
        chunk += 1
        chunk_path = checkm_path.joinpath(dataset + '_chunk_' + str(chunk))
        chunk_path.mkdir(parents = True)
        for genome in genome_chunk:
            chunk_path.joinpath(genome).symlink_to(processed_path.joinpath(genome))
    print('Created {} chunks.'.format(chunk))

def get_ext_checkm_targets(dataset_dir):
    """
    For a given external dataset that has been processed gives the targets for snakemake
    """
    checkm_path = Path(dataset_dir).joinpath('checkm')
    if not checkm_path.exists():
        raise Exception('Did you setup the evaluation?')
    targets = []
    for chunk in checkm_path.glob('*_chunk_*[0-9]'):
        targets.append(str(chunk) + '-checkm_ext.done')
    return(targets)

def get_ext_anvi_evaluate_targets(dataset_dir):
    """
    For a given external dataset that has been processed gives the targets for snakemake
    """
    anvio_path = Path(dataset_dir).joinpath('anvio_db')
    genomes_path = Path(dataset_dir).joinpath('processed_genomes')
    if not anvio_path.exists():
        raise Exception('Did you setup the evaluation?')
    targets = []
    for genome in genomes_path.glob('*.fa'):
        targets.append(str(anvio_path.joinpath(genome.name.replace('.fa', ''))) + '-anvi_evaluate_ext.done')
    return(targets)

if __name__ == '__main__':
    setup_checkm_folder('delmont_mags')
    setup_checkm_folder('gorg')
    setup_checkm_folder('mardb')
