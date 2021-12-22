#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Combine the profiles in single tables, normalize by gene length if appropriate
"""

# Import libraries -------------------------------------------------------------

import argparse
from pathlib import Path
from itertools import zip_longest

import magpipe.log
import magpipe.pretty

# Define utilitary functions ---------------------------------------------------

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', help='Path to the profiling directory.', default='/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/MAGPIPE_MAGS_EAN/scratch/processed/quantification/go_microbiomics/go_microbiomics-integrated-cpl50_ctn10-genes-cds-w_eukarya/profile/', type=str)
    parser.add_argument('-o', '--output', help='Path and prefix of the output tables.', required=True, type=str)
    return parser.parse_args()

def process_multiple_lines(lines, column_to_summarize):
    """
    This functions takes the nth line of a number of zipped
    profiles. We want to extract the relevant columns
    and combine them in a single tsv file
    """
    lines = [l.strip().split("\t") for l in lines]
    gene = lines[0][0]
    gene_length = lines[0][1]
    #if not all([l[0] == gene for l in lines]) or not all([l[1] == gene_length for l in lines]):
    #    raise ValueError("Lines don't match.. {}".format(lines))
    values = [i[column_to_summarize] for i in lines]
    formatted_line = "\t".join([gene, gene_length] + values)
    return formatted_line

# Main function of the script --------------------------------------------------

def main(args):
    """
    This functions identifies all the profiles, opens them and produces the
    merged files for horizontal profiles, insert counts and length normalized 
    base counts
    """
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError("Looks like there is something wrong with {}.".format(input_path))
    # prepare outputs
    horizontal_out = Path(args.output  + "-horizontal.tsv")
    insert_out = Path(args.output  + "-insert.tsv")
    base_out = Path(args.output  + "-base.tsv")
    if any([horizontal_out.exists(), insert_out.exists(), base_out.exists()]):
        raise ValueError("Looks like some of your output already exists... please check {}, {}, {}.".format(horizontal_out, insert_out, base_out))
    # Find all the profiles
    file_name_list = list(input_path.glob("*/2*/*.profile")) # FIXME Might not be the most robust..
    samples = [".".join(s.name.split(".")[0:-2]) for s in file_name_list]
    if len(set(samples)) != len(samples):
        raise ValueError("Looks like there are more than one profile per sample, can't really handle that.")
    magpipe.pretty.print_pair("Number of samples", len(samples))
    # Try to open all of the files
    # Can't use 'with' as we have an unknown set of files
    # Open all of the files - trapping exceptions as we go
    try:
        open_files = [open(file_name) for file_name in file_name_list]
    except OSError:
        raise
    # Read all the lines at the same time
    # unpack the open_files list into zip
    # each open file is an iterator in it's own right
    try:
        with open(horizontal_out, 'w') as horizontal_target, open(insert_out, 'w') as insert_target, open(base_out, 'w') as base_target:
            for lines in zip_longest(*open_files):
                # lines will be a tuple of all of the
                # corresponding lines from the files
                if lines[0].startswith("#reference"):
                    horizontal_target.write("# horizontal coverage, no length normalisation\n")
                    insert_target.write("# insert count, no length normalisation\n")
                    base_target.write("# base count, length normalised\n")
                    horizontal_target.write("\t".join(["#reference", "length"] + samples) + "\n")
                    insert_target.write("\t".join(["#reference", "length"] + samples) + "\n")
                    base_target.write("\t".join(["#reference", "length"] + samples) + "\n")
                else:
                    horizontal_target.write(process_multiple_lines(lines, column_to_summarize=2) + "\n")
                    insert_target.write(process_multiple_lines(lines, column_to_summarize=6) + "\n")
                    base_target.write(process_multiple_lines(lines, column_to_summarize=3) + "\n")
    except:
        raise
    finally:
        # since we aren't using a context manager (with)
        # we have to ensure all files get closed.
        for open_file in open_files:
            open_file.close()

# Run script -------------------------------------------------------------------

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
