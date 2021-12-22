#!/usr/bin/env python

# Import libraries -------------------------------------------------------------

import os
import sys
import pathlib
import argparse
import pandas as pd

# Define utilitary functions ---------------------------------------------------

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser(
    prog = 'prepare_anvio_collection.py',
    description = 'Takes a binning output and prepare a table anvio likes',
    formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # REQUIRED  arguments:
    parser.add_argument('-i', '--input',
    help = 'Path to bins directory.', 
    required = True,
    type = pathlib.Path)
    parser.add_argument('-o', '--output',
    help = 'Path and name of the output file.',
    required = True,
    type = pathlib.Path)

    return parser.parse_args()


def print_arguments(args):
    """
    Print arguments
    :param args
    """
    print("\nChecking arguments:")
    print("Input cluster file : {}".format(args.input))
    print("Output file : {}\n".format(args.output))


# Functions to identify the bins -----------------------------------------------

def prep_metabat2_clusters(input_directory, output_file, suffix = '*.fa'):
    """
    Read the bins fasta file to generate a table for anvio
    Replace '.' by '_' to make anvio happy
    """
    bins_path = pathlib.Path(input_directory)
    contigs_to_cluster = {}
    bins_fastas = [fasta for fasta in bins_path.glob(suffix)]
    for f in bins_fastas:
        cluster = 'BIN_' + f.name.rstrip(suffix).replace('.', '_')
        with open(f, "r") as handle:
            for line in handle:
                if line.startswith('>'):
                    contigs_to_cluster[line.lstrip('>').rstrip('\n')] = cluster

    df = pd.DataFrame.from_dict(contigs_to_cluster, orient = 'index') 
    if not pathlib.Path(output_file).parent.exists():
        pathlib.Path(output_file).parent.mkdir(parents = True)
    df.to_csv(output_file, sep = '\t', header = False)

# Run script -------------------------------------------------------------------

if __name__ == '__main__':
    args = get_arguments()
    print_arguments(args)
    prep_metabat2_clusters(args.input, args.output)
