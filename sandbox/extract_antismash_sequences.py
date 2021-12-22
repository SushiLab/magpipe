#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
This is a convenience script to export anstismash biosynthetic regions into individual fasta files, e.g. for drep 
"""

import pandas
import argparse
from Bio import SeqIO
from tqdm import tqdm
from pathlib import Path

# Remove biopython warning as the header of the gbk files are too long and we get a crazy amount of warnings (although no consequences)
import warnings
from Bio import BiopythonParserWarning

import magpipe.log

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser()
    # REQUIRED  arguments:
    parser.add_argument('-a', '--annotations', help='Path to the annotation directory.', required=True, type=str)
    parser.add_argument('-t', '--table', help='Path to the summary antismash table, used to easily find the gbk files.', required=True, type=str)
    parser.add_argument('-o', '--output', help='Path and name of the output directory where the individual fasta will be stored.', required=True, type=str)
    return parser.parse_args()

def get_sequence_from_gbk(gbk_file):
    """
    A biopython based parser for antismash gbk output
    simply parse the file and return the sequence of the region
    """
    with open(gbk_file) as handle:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonParserWarning)
            record = list(SeqIO.parse(handle, "genbank"))
    if len(record) != 1:
        raise ValueError("Each file should be a single record... and for {} I've got {}...".format(gbk_file, len(record)))
    record = record[0]
    biosynth_seq = str(record.seq)
    return biosynth_seq

def main(annotation_dir, antismash_table, output_dir):
    """
    Main function of the pipeline
    Loops through the biosynthetic regions and write the corresponding fasta files
    """
    annotation_path = Path(annotation_dir)
    if not annotation_path.exists():
        raise FileNotFoundError("Looks like there is something wrong with {}.".format(annotation_path))
    output_path = Path(output_dir)
    if output_path.exists():
        raise Exception("Output {} already exists...".format(output_path))
    output_path.mkdir()
    antismash_table = pandas.read_csv(antismash_table, sep = "\t")
    expected_columns = ["genome", "scaffold", "region", "length", "contig edge", "products", "# candidate clusters", "candidate clusters type", "# protoclusters", "protoclusters products", "# CDS", "CDS list"]
    if not list(antismash_table.columns) == expected_columns:
        raise ValueError("The columns don't match what I was expecting, are you sure the antismash table is right ? It should have {}.".format(expected_columns))
    for i in tqdm(antismash_table.index, ncols=100):
        current = antismash_table.iloc[i]
        biosynth_name = "-".join([current.genome, current.region])
        biosynth_num = current.region.split("_")[-1]
        biosynth_num_formatted = "0" * (3 - len(biosynth_num)) + biosynth_num
        gbk_file = annotation_path.joinpath(current.genome, "antismash", current.genome + "-antismash", current.scaffold + ".region" + biosynth_num_formatted + ".gbk")
        if not gbk_file.exists():
            raise FileNotFoundError("Couldn't find file {}.".format(gbk_file))
        biosynth_seq = get_sequence_from_gbk(gbk_file)
        with open(output_path.joinpath(biosynth_name + ".fa"), "w") as fasta:
            fasta.write(">" + biosynth_name + "\n")
            fasta.write(biosynth_seq + "\n")

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args.annotations, args.table, args.output)
