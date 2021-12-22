#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Given a set of samples, generate atomized depth files. These are depth files for which we only have the profile of the sample in that sample, i.e. remove the co-abundances
"""

# Import libraries -------------------------------------------------------------

import re
import pandas
import pathlib
import argparse

import magpipe.log
import magpipe.pretty

# Define utilitary functions ---------------------------------------------------

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser(
        description='Read scaffold membership and approximate MAG abundances based on backmapping data.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # REQUIRED  arguments:
    parser.add_argument('-i', '--input', nargs='+', help='List of sample files to process. Each sample file corresponds to a dataset and should match the standard, e.g. TARA_OCEANS_prok.randN.samples.', required=True, type=str)
    parser.add_argument('-d', '--datasets', help='Path to the metagenomic datasets folder.', default="/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/MAGPIPE_MAGS_EAN/scratch/processed/metagenomes/", type=str)
    return parser.parse_args()

def check_file_or_path(file_or_path):
    """
    just a quick function to raise an error is a file or directory exists
    """
    if not isinstance(file_or_path, pathlib.Path):
        file_or_path = pathlib.Path(file_or_path)
    if not file_or_path.exists():
        raise FileNotFoundError("Looks like there is something wrong with {}.".format(file_or_path))

def atomize_depth_file(depth_file, sample):
    """
    read a depth file with pandas, remove columns with other samples, rewrite the necessary column and return table 
    """
    depth_table = pandas.read_csv(depth_file, sep="\t")
    sample_vs = sample + "_vs_" + sample + ".bam"
    depth_table = depth_table[["contigName", "contigName", "totalAvgDepth", sample_vs, sample_vs + "-var"]]
    depth_table["totalAvgDepth"] = depth_table[sample_vs]
    return depth_table

# Main function of the script --------------------------------------------------

def main(args):
    """
    The main function, loops through the sample files (and therefore datasets), reads the depth files and writes the depth files with only the sample 
    """
    datasets_path = pathlib.Path(args.datasets)
    check_file_or_path(datasets_path)
    # Looping through sample files
    for sample_file in args.input:
        sample_file_path = pathlib.Path(sample_file)
        check_file_or_path(sample_file)
        sample_list = sample_file_path.read_text().splitlines()
        dataset = sample_file_path.name.split(".")[0]
        dataset_path = datasets_path.joinpath(dataset)
        check_file_or_path(dataset_path)
        # we then want to guess the backmapped value of the normal depth file
        a_value = len(list(dataset_path.glob("*")))
        magpipe.pretty.print_pair("Processing " + dataset + " with 'a'", str(a_value), nl_before=1)
        for sample in sample_list:
            sample_path = dataset_path.joinpath(sample)
            check_file_or_path(sample_path)
            depth_path = sample_path.joinpath("depth_files")
            depth_file = list(depth_path.glob("{}*a{}*depth.t*".format(sample, a_value))) # Backcompatibility between _, . and - as well as txt and tsv
            if len(depth_file) != 1:
                raise ValueError("Expected to find one matching depth and found {}".format(depth_file))
            depth_file = depth_file[0]
            magpipe.pretty.print_pair(sample, depth_file)
            depth_atomized = atomize_depth_file(depth_file, sample)
            depth_atomized_file = re.sub("a" + str(a_value), "a1", str(depth_file))
            if pathlib.Path(depth_atomized_file).exists():
                raise Exception("{} already exists. Cowardly refusing to overwrite it.".format(depth_atomized_file))
            depth_atomized.to_csv(depth_atomized_file, index=False, sep="\t")
            pathlib.Path(".".join(depth_atomized_file.split(".")[0:-1]) + ".done").touch()
        
# Run script -------------------------------------------------------------------

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
