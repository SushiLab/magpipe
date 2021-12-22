#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

import re
import gzip
import argparse
from pathlib import Path
import Bio.SeqIO.FastaIO as FastaIO

import magpipe.log
import magpipe.pretty


def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser()
    # REQUIRED arguments:
    parser.add_argument('-m', '--mags', help = 'Path of the MAGs in need of renaming.', required = True, type = str)
    # return namespace
    return parser.parse_args()

def rename_mags(mags_directory, suffix = ".fa"):
    """
    Given the export MAGs directory, rename the MAGs and scaffolds according to:
    DATASET_SAMPLE_MAG_00000001
    DATASET_SAMPLE_MAG_00000001-scaffold_n
    """
    mags_path = Path(mags_directory)
    if not mags_path.exists():
        raise Exception("Looks like the mags directory doesn't exist...")
    mag_fastas = [i for i in mags_path.glob(f'*{suffix}')]
    if len(mag_fastas) >= 1:
        destination = mags_path.parent.joinpath(f"{mags_path.name}-renamed")
        if destination.exists():
            raise Exception(f"Yeah... something looks off. {destination} exists but it shouldn't.")
        destination.mkdir()
        genome_dict = destination.parent.joinpath(f"{destination.name}-genome_dict.tsv")
        scaffold_dict = destination.parent.joinpath(f"{destination.name}-scaffold_dict.tsv")
        with open(genome_dict, 'w') as genome_dict_file:
            genome_dict_file.write("old_mag\tnew_mag\n")
            with open(scaffold_dict, 'w') as scaff_dict_file:
                scaff_dict_file.write("old_scaffold\tnew_scaffold\n")
                for mag_fasta in mag_fastas:
                    mag_number = mag_fasta.name.split('.')[-2] 
                    sample = '-'.join(mag_fasta.name.split('-')[0:-1])
                    if "_METAG" in sample:
                        prefix = re.sub("_METAG$", "_MAG", sample)
                    else:
                        prefix = sample + "_MAG" 
                    blank = (8 - len(mag_number))*"0"
                    mag_name = f"{prefix}_{blank}{mag_number}" 
                    genome_dict_file.write(f"{'.'.join(mag_fasta.name.split('.')[0:-1])}\t{mag_name}\n")
                    with open(mag_fasta) as handle:
                        #with gzip.open(destination.joinpath(f"{mag_name}.fa.gz"), 'wt') as target:
                        with open(destination.joinpath(f"{mag_name}.fa"), 'w') as target:
                            counter = 0
                            for line in handle:
                                if line.startswith(">"):
                                    original_scaffold = line.strip().split('>')[-1]
                                    counter += 1
                                    new_scaffold = f"{mag_name}-scaffold_{counter}"
                                    scaff_dict_file.write(f"{original_scaffold}\t{new_scaffold}\n")
                                    target.write(f">{new_scaffold}\n")
                                else:
                                    target.write(f"{line}")
    else:
        print(f"No mags found in {mags_directory}, could be absolutely fine, just wanted to let you know so you could make sure.")


def main(args):
    """
    Rename a bunch of MAGs 
    """
    magpipe.pretty.print_single("Exporting and renaming the MAGs...")
    rename_mags(args.mags)

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
