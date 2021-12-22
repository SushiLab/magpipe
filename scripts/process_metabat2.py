#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

import re
import gzip
import pandas
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
    parser.add_argument('-m', '--metabat2_dir', help = 'Path of the metabat2 output.', required = True, type = str)
    parser.add_argument('-a', '--assembly', help = 'Path to the assembly.', required = True, type = str)
    parser.add_argument('-c', '--seqtk_comp', help = 'Path of the seqtk comp tmp file.', required = True, type = str)
    parser.add_argument('-l', '--min_length', help = 'Min length of the scaffolds.', required = True, type = int)
    # return namespace
    return parser.parse_args()


def identify_missing_scaffolds(bins_membership, seqtk_comp, min_length):
    """
    reads the membership, identifies unbinned scaffolds and keep the ones above min length
    """
    bins_table = pandas.read_csv(bins_membership, sep = '\t', names = ["scaffold", "bin"])
    seqtk_table = pandas.read_csv(seqtk_comp, sep = '\t', header = None)
    table = bins_table.merge(seqtk_table, left_on = "scaffold", right_on = 0)
    table = table.loc[table.bin == 0]
    table = table.loc[table[1] >= min_length]
    return table.scaffold.to_list(), bins_table

def recover_missing_scaffolds(raw_bins_dir, seqtk_comp, assembly_file, min_length):
    """
    reads seqtk comp file as well as the binning results and recovers scaffolds over Xkbp that were left out
    """
    if not bins_path.exists():
        raise FileNotFoundError("Couldn't find the binning results {}".format(bins_path))
    if not seqtk_comp.exists():
        raise FileNotFoundError("Couldn't find the seqtk comp file".format(seqtk_comp))
    if not assembly_file.exists() or assembly_file.name.endswith(".gz"):
        raise ValueError("There was something wrong with {}. Currently expects a non gzipped file".format(assembly_file))
    bins_files = [f for f in bins_path.glob("*")]
    bins_fasta = [f for f in bins_files if f.name.endswith(".fa")]
    bins_membership = [f for f in bins_files if not f.name.endswith(".fa")]
    if len(bins_membership) != 1:
        raise ValueError("Something went wrong trying to identify the bin membership file in {}".format(bins_path))
    bins_membership = bins_membership[0]
    if len(bins_fasta) == 0:
        bins_number = 0
    else:
        bins_number = max([int(f.name.split(".")[-2]) for f in bins_fasta])
    if not bins_number == len(bins_fasta):
        raise ValueError("The number of bins doesn't match their name in {}".format(bins_path))
    magpipe.pretty.print_pair("Number of bins initially", len(bins_fasta))
    scaffolds_to_recover, membership_table = identify_missing_scaffolds(bins_membership, seqtk_comp, args.min_length)
    bin_counter = bins_number
    # loop through the assembly 
    with open(assembly_file, "rt") as handle:
        for (header, sequence) in FastaIO.SimpleFastaParser(handle):
            if header in scaffolds_to_recover:
                bin_counter += 1
                dest = bins_path.joinpath(".".join([bins_membershhip.name, bin_counter, "fa"]))
                membership_table.loc[membership_table.scaffold == header, "bin"] = "recovered_{}".format(bin_counter)
                with open(dest, "w") as target:
                   target.write(">" + header + "\n") 
                   target.write(sequence + "\n") 
            if len(sequence) < args.min_length:
                break
    membership_table.to_csv(bins_path.joinpath(bins_membership.name + ".fixed"), sep = '\t', index = False)
    magpipe.pretty.print_pair("Number of bins recovered", bin_counter - bins_number)

def rename_bins(bins_directory, suffix = ".fa"):
    """
    Given the metabat output bins directory, rename the BINs and scaffolds according to:
    DATASET_SAMPLE_BIN_00000001
    DATASET_SAMPLE_BIN_00000001-scaffold_n
    """
    bins_path = Path(bins_directory)
    if not bins_path.exists():
        raise Exception("Looks like the bins directory doesn't exist...")
    bin_fastas = [i for i in bins_path.glob(f'*{suffix}')]
    if len(bin_fastas) >= 1:
        destination = bins_path.parent.joinpath("bins")
        if destination.exists():
            raise Exception(f"Yeah... something looks off. {destination} exists but it shouldn't.")
        destination.mkdir()
        for bin_fasta in bin_fastas:
            bin_number = bin_fasta.name.split('.')[-2] 
            sample = '-'.join(bin_fasta.name.split('-')[0:-1])
            prefix= re.sub("_METAG$", "_BIN", sample)
            blank = (8 - len(bin_number))*"0"
            bin_name = f"{prefix}_{blank}{bin_number}" 
            scaffold_dict = destination.joinpath(f"{sample}-scaffold_dict.tsv")
            with open(bin_fasta) as handle:
                with gzip.open(destination.joinpath(f"{bin_name}.fa.gz"), 'wt') as target:
                    with open(scaffold_dict, 'w') as dict_file:
                        counter = 0
                        for line in handle:
                            if line.startswith(">"):
                                original_scaffold = line.strip().split('>')[-1]
                                counter += 1
                                new_scaffold = f"{bin_name}-scaffold_{counter}"
                                dict_file.write(f"{original_scaffold}\t{new_scaffold}\n")
                                target.write(f">{new_scaffold}\n")
                            else:
                                target.write(f"{line}")
    else:
        print(f"No bins found in {bins_directory}, could be absolutely fine, just wanted to let you know so you could make sure.")


def main(args):
    """
    Recover missed bins (older version of metabat) and rename bins
    """
    raw_bins_path = Path(args.metabat2_dir)
    if not raw_bins_path.name == "raw_bins":
        raise Exception("Expects metabat2 output to be 'raw_bins' so the processed one is 'bins'.")
    seqtk_comp = Path(args.seqtk_comp)
    assembly_file = Path(args.assembly)
    recover_missing_scaffolds(raw_bins_path, seqtk_comp, assembly_file, args.min_length)
    magpipe.pretty.print("Exporting and renaming the resuts...")
    rename_bins(args.metabat2_dir)

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
