#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
This script parses the annotations output to summarise all eggnog annotations 
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
    Reads the list of genomes and process the eggnog tables to generate a concatenated annotation file 
    """
    annotation_path = Path(annotation_dir)
    if not annotation_path.exists():
        raise FileNotFoundError(f"Looks like there is something wrong with {annotation_path}.")
    output_path = Path(output_table)
    if output_path.exists():
        raise Exception("Output {} already exists...".format(output_path))
    table_columns = ["query_name", "seed_eggNOG_ortholog", "seed_ortholog_evalue", "seed_ortholog_score", "best_tax_level", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "taxonomic scope", "eggNOG OGs", "best eggNOG OG", "COG Functional cat.", "eggNOG free text desc."]
    columns_string = "\t".join(table_columns)
    with open(output_path, 'w') as target:
        target.write(f"genome\t{columns_string}\n")
        genome_count = 0
        for genome_path in tqdm(annotation_path.glob('*'), ncols=100):
            genome_count += 1
            genome = genome_path.name
            eggnog_file = genome_path.joinpath(f"eggnog/{genome}-eggnog.emapper.annotations")
            #magpipe.pretty.print_pair("Processing file", eggnog__file)
            with open(eggnog_file) as handle:
                for line in handle:
                    line = line.rstrip()
                    if line.startswith("#"):
                        continue
                    target.write(f"{genome}\t{line}\n")
    print("Finished processing {} genomes.".format(genome_count))

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args.annotations, args.output)
