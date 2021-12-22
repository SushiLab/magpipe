#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
This script parses the annotations output to summarise all non CDS features in a nicely formatted
table
"""

import pandas
import argparse
from tqdm import tqdm
from pathlib import Path

import magpipe.log
import magpipe.pretty

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser()
    # REQUIRED  arguments:
    parser.add_argument('-a', '--annotations', help='Path to the annotation directory.', required=True, type=str)
    parser.add_argument('-o', '--output', help='Path and name of the output table.', required=True, type=str)
    return parser.parse_args()

def main(annotation_dir, output_table):
    """
    Main function of the pipeline
    Reads the list of genomes and process the prookka tables to make a nice table summary of non
    CDS features
    """
    annotation_path = Path(annotation_dir)
    if not annotation_path.exists():
        raise FileNotFoundError("Looks like there is something wrong with {}.".format(annotation_path))
    output_path = Path(output_table)
    if output_path.exists():
        raise Exception("Output {} already exists...".format(output_path))
    table_columns = ["genome", "gene", "feature", "length", "description"]
    output_path.write_text("\t".join(table_columns) + "\n")
    genome_count = 0
    for genome_path in tqdm(annotation_path.glob('*'), ncols=100):
        genome_count += 1
        genome = genome_path.name
        prokka_file = list(genome_path.joinpath("prokka").glob("{genome}-*-prokka/{genome}-prokka.tsv".format(genome=genome)))[0]
        #magpipe.pretty.print_pair("Processing file", prokka_file)
        prokka_table = pandas.read_csv(prokka_file, sep="\t")
        if not all([c in prokka_table.columns for c in ["locus_tag", "ftype", "length_bp", "gene", "EC_number", "COG", "product"]]):
            raise ValueError("Couldn't find the columns I'm looking for... did the names change?")
        if len(prokka_table.columns) > 7:
            raise ValueError("This table looks too big, I only expect 7 columns and got {}".format(len(prokka_table.columns)))
        prokka_table.insert(0, "genome", genome)
        prokka_table = prokka_table.loc[[i != "CDS" for i in prokka_table['ftype']]]
        prokka_table = prokka_table[["genome", "locus_tag", "ftype", "length_bp", "product"]]
        prokka_table.to_csv(output_path, sep="\t", index=False, header=False, mode='a')
    print("Finished processing {} genomes.".format(genome_count))

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args.annotations, args.output)
