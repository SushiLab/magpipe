#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Given a folder of fasta files, produces a map between scaffolds and their genome
"""

# Import libraries -------------------------------------------------------------

import re
import pandas
import pathlib
import argparse
import Bio.SeqIO.FastaIO as FastaIO
from collections import defaultdict

import magpipe.log

# Define utilitary functions ---------------------------------------------------

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser(description='Creates a big table with scaffolds and their bin',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # REQUIRED  arguments:
    parser.add_argument('-i', '--input', help='Path to the bins in fasta format.', required=True, type=pathlib.Path)
    parser.add_argument('-s', '--suffix', help='Suffix to add to the path to find the fasta files.', default='*.fa', type=str)
    parser.add_argument('-o', '--output', help='Path and name of the output file.', required=True, type=pathlib.Path)
    return parser.parse_args()

# Functions to identify the bins -----------------------------------------------

def export_bins_membership(input_directory, suffix, output_file):
    """
    Read the bins fasta file to generate a table of scaffolds and their bin
    """
    bins_path = pathlib.Path(input_directory)
    scaffolds_to_cluster = defaultdict()
    suffix = re.sub(r'^\*', '', suffix)
    bins_fastas = list(bins_path.glob("*" + suffix))
    for f in bins_fastas:
        cluster = re.sub(suffix + '$', '', f.name)
        with open(f, "r") as handle:
            for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                scaffolds_to_cluster[header] = cluster

    df = pandas.DataFrame.from_dict(scaffolds_to_cluster, orient='index')
    if not pathlib.Path(output_file).parent.exists():
        pathlib.Path(output_file).parent.mkdir(parents=True)
    df.to_csv(output_file, sep='\t', header=False)

# Run script -------------------------------------------------------------------

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    export_bins_membership(args.input, args.suffix, args.output)
