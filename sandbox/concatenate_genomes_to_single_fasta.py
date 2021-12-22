#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Given a folder of fasta files, generate a concatenate fasta file with headers like
>FILENAME-HEADER to replace >HEADER
"""

# Import libraries -------------------------------------------------------------

import re
import pathlib
import argparse
import Bio.SeqIO.FastaIO as FastaIO

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
    parser.add_argument('-i', '--input', help='Path to the genomes in fasta format.', required=True, type=pathlib.Path)
    parser.add_argument('-s', '--suffix', help='Suffix to add to the path to find the fasta files.', default='*.fa', type=str)
    parser.add_argument('-o', '--output', help='Path and name of the output file.', required=True, type=pathlib.Path)
    return parser.parse_args()

# Functions to identify the genomes -----------------------------------------------

def concat_genomes(input_directory, suffix, output_file):
    """
    Read the genomes fasta file to generate a table of scaffolds and their bin
    """
    genomes_path = pathlib.Path(input_directory)
    output_path = pathlib.Path(output_file)
    if not genomes_path.exists():
        raise Exception("Looks like there is something wrong with your input {}.".format(input_directory))
    if output_path.exists():
        raise Exception("Output {} already exists, please clean things up.".format(output_file))
    if not output_path.parent.exists():
        output_path.parent.mkdir(parents=True)
    suffix = re.sub(r'^\*', '', suffix)
    genomes_fastas = list(genomes_path.glob("*" + suffix))
    if not len(genomes_fastas) >= 1:
        raise Exception("Did not find any fasta file, please check input path and suffix.")
    with open(output_path, "w") as outfile:
        for f in genomes_fastas:
            genome = re.sub(suffix + '$', '', f.name)
            with open(f, "r") as handle:
                for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                    outfile.write(">" + genome + "-" + header + "\n")
                    outfile.write(sequence + "\n")

# Run script -------------------------------------------------------------------

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    concat_genomes(args.input, args.suffix, args.output)
