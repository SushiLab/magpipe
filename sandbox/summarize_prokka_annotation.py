#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
This script parses the annotations output to summarise all prokka annotations 
table
"""

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
    Reads the list of genomes and process the prokka tables to generate a concatenated annotation file 
    """
    annotation_path = Path(annotation_dir)
    if not annotation_path.exists():
        raise FileNotFoundError(f"Looks like there is something wrong with {annotation_path}.")
    output_path = Path(output_table)
    if output_path.exists():
        raise Exception("Output {} already exists...".format(output_path))
    table_columns = ["locus_tag", "ftype", "length_bp", "gene", "EC_number", "COG", "product"]
    columns_string = "\t".join(table_columns)
    with open(output_path, 'w') as target:
        target.write(f"genome\t{columns_string}\n")
        genome_count = 0
        for genome_path in tqdm(annotation_path.glob('*'), ncols=100):
            genome_count += 1
            genome = genome_path.name
            prokka_file = list(genome_path.joinpath("prokka").glob("{genome}-*-prokka/{genome}-prokka.tsv".format(genome=genome)))[0]
            #magpipe.pretty.print_pair("Processing file", prokka__file)
            with open(prokka_file) as handle:
                for line in handle:
                    line = line.rstrip()
                    if line.startswith("locus_tag"):
                        if line != columns_string:
                            raise Exception(f"Looks like the headers are not what I expected... supposed to be {columns_string}")
                        continue
                    target.write(f"{genome}\t{line}\n")
    print("Finished processing {} genomes.".format(genome_count))

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args.annotations, args.output)
