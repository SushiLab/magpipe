#!/usr/bin/env python

# Imports ----------------------------------------------------------------------

import argparse
from pathlib import Path

import pandas
from tqdm import tqdm
import Bio.SeqIO.FastaIO as FastaIO

import magpipe.log
import magpipe.pretty

# Define utilitary functions ---------------------------------------------------

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser()
    # REQUIRED arguments:
    parser.add_argument('-a', '--annotations', help='Path to the annotations directory, which contains FetchMGs/{genomes}-bestMGs and -allMGs.', required=True, type=str)
    parser.add_argument('-o', '--output', help='Path and name of the output directory.', required=True, type=str)
    # OPTIONAL arguments:
    #parser.add_argument('-g', '--genomes', help='A file with the list of genomes to process.', required=False, default='/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/MAGPIPE_MAGS_EAN/scratch/processed/integrated/go_microbiomics/go_microbiomics-integrated-cpl50_ctn10-dictionaries/test-genomes-specI', type=str)
    parser.add_argument('-g', '--genomes', help='A file with the list of genomes to process.', required=False, default='/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/MAGPIPE_MAGS_EAN/scratch/processed/integrated/go_microbiomics/go_microbiomics-integrated-cpl50_ctn10-dictionaries/go_microbiomics-integrated-cpl50_ctn10-prokarya.txt', type=str)
    # Return namespace
    return parser.parse_args()

def load_file_as_list(input_file): 
    """
    Takes a file file and lists all the entries in it
    """
    magpipe.pretty.print_pair("Loading file", input_file)
    file_list = []
    file_path = Path(input_file) 
    if not file_path.exists():
        raise ValueError('Well then, input file {} doesn\'t exist.'.format(file_path))
    with open(file_path) as handle:
        for line in handle: 
            file_list.append(line.strip())
    return file_list

def main(args):
    """
    Main function to summarize marker genes results 
    """
    output = Path(args.output)
    if output.exists():
        raise ValueError('Seems the script has already been ran... cleanup and start again!')
    genomes_list = load_file_as_list(args.genomes)
    res_dict = {}
    mgs = None
    for genome in tqdm(genomes_list, ncols=100):
        allmgspath = Path(args.annotations).joinpath(genome, "fetchMGs", genome + "-allMGs") 
        bestmgspath = Path(args.annotations).joinpath(genome, "fetchMGs", genome + "-bestMGs")
        if not allmgspath.exists() or not bestmgspath.exists():
            raise FileNotFoundError("Couldn't find {} or {}...".format(allmgspath, bestmgspath))
        if mgs is None:
            mgs = [f.name.replace(".fna", "") for f in allmgspath.glob("*.fna")]
            if len(mgs) != 40:
                raise ValueError("Well, I was expecting 40 MGs and found {} when processing {}.".format(len(mgs), genome))
        res_dict[genome] = {}
        for mg in mgs:
            allmg_fna = allmgspath.joinpath(mg + ".fna")
            bestmg_fna = bestmgspath.joinpath(mg + ".fna")
            with open(allmg_fna) as allmg_in, open(bestmg_fna) as bestmg_in:
                alllist = [{'h': all_h, 's': all_s} for (all_h, all_s) in FastaIO.SimpleFastaParser(allmg_in)]
                bestlist = [{'h': best_h, 's': best_s} for (best_h, best_s) in FastaIO.SimpleFastaParser(bestmg_in)]
                res_dict[genome][mg + "-all"] = len(alllist)
                res_dict[genome][mg + "-best"] = len(bestlist)
    # convert to df and write 
    res_df = pandas.DataFrame.from_dict(res_dict, orient = 'index')
    res_df.index.rename("genome", inplace = True)
    res_df.to_csv(output, sep = "\t")

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
