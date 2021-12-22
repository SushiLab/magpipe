#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Give scaffold membership and a dataset this script can approximate the abundance of a mag across the
dataset's samples
"""

# Import libraries -------------------------------------------------------------

import re
import json
import pandas
import pathlib
import argparse
import itertools

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
    parser.add_argument('-i', '--input', help='Path to the scaffold membership table.', required=True, type=str)
    parser.add_argument('-o', '--output', help='Path and prefix of the output table.', required=True, type=str)
    parser.add_argument('-f', '--scaffold_mapping', help='Path to the scaffold mapping files if scaffold were renamed (can be a list). Expects old{tab}new.', nargs='+', default=None, type=str)
    parser.add_argument('-l', '--sample_mapping', help='Path to the sample mapping files if scaffold were renamed (list). Expects old{tab}new.', nargs='+', default=None, type=str)
    parser.add_argument('-a', '--dataset_abbrev', help='Dictionary of the dataset abreviations. (str -> list)', default='{"GEOTRACES": ["BATS", "BGEO", "HOTS"], "MALASPINA": ["MALA"], "TARA_OCEANS_prok": ["TARA"], "TARA_OCEANS_vir": ["TARA"]}', type=json.loads)
    parser.add_argument('-d', '--datasets', help='Path to the metagenomic datasets folder.', default="/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/MAGPIPE_MAGS_EAN/scratch/processed/metagenomes/", type=str)
    parser.add_argument('-s', '--suffix', help='Suffix to add to the path to find the depth files.', default='.depth.txt', type=str)
    parser.add_argument('-x', '--index', help='Minimal vertical coverage for the sample to be included in the diff coverage index.', default=1, type=int)
    return parser.parse_args()

def load_sample_mapping(sample_mapping_list):
    """
    Load the sample mapping dictionary
    The funtion takes a list of tab-separated files (with header) as input and returns a dictionary with the old names (first column)
    as keys and the new names as values (second column).
    """
    sample_mapping = {}
    for sample_mapping_file in sample_mapping_list:
        magpipe.pretty.print_pair("Loading", sample_mapping_file)
        sample_mapping_path = pathlib.Path(sample_mapping_file)
        if not sample_mapping_path.exists():
            raise FileNotFoundError("Well, couldn't find {}.".format(sample_mapping_path))
        mapping_df = pandas.read_csv(sample_mapping_path, sep='\t')
        mapping_df.columns = ["old", "new"]
        # Make sure everything is alright, use METAG do desambiguate METAG and METAT
        mapping_df.old = [s + "_METAG" if s.startswith("TARA_") and not s.endswith("_METAG") else s for s in mapping_df.old]
        mapping_df.new = [s if s.endswith("_METAG") else s + "_METAG" for s in mapping_df.new]
        mapping_df.set_index("old", inplace=True)
        mapping_dict = mapping_df["new"].to_dict() # Needs to be from Series to dict
        sample_mapping.update(mapping_dict)
    return sample_mapping

def load_scaffold_mapping(scaffold_mapping_list):
    """
    Load the scaffold mapping dictionary
    The funtion takes a list of tab-separated files (w/o header) as input and returns a dictionary with the new names (second column)
    as keys and the old names as values (first column).
    """
    scaffold_mapping = {}
    for scaffold_mapping_file in scaffold_mapping_list:
        magpipe.pretty.print_pair("loading", scaffold_mapping_file)
        scaffold_mapping_path = pathlib.Path(scaffold_mapping_file)
        if not scaffold_mapping_path.exists():
            raise FileNotFoundError("Well, couldn't find {}.".format(scaffold_mapping_path))
        mapping_df = pandas.read_csv(scaffold_mapping_path, sep='\t', header=None)
        mapping_df.set_index(1, inplace=True)
        scaffold_mapping.update(mapping_df[0].to_dict())
    return scaffold_mapping

def select_columns(sample, columns):
    """
    Filter the column of a depth file based on the sample currently being processed
    They look like:
    [contigName, length, totalAbundance, sample_vs_otherSample, sample_vs_otherSample-var ...]
    And we only want to keep sample_vs_otherSample
    """
    res = list()
    for col in columns:
        if sample in col and '-var' not in col:
            res.append(col)
    return res

def process_dataset(dataset_path, mags, suffix, membership, scaffold_map=None, sample_map=None, index=1):
    """
    For a given dataset, identify all the depth files and process them to approximate the
    abundance of the mags in the list. For this it uses the scaffold membership and takes the median
    values of all the scaffolds.
    It needs for that a couple of maps in case there was renaming going on:
    scaffold_map, a dict from the new names (genomes) to the old one (depth files)
    samples_map, a dict from the old names (the ones in the dataset path) and the new one (the ones in the genomes)
    """
    dataset = dataset_path.name
    output_path = pathlib.Path("{prefix}-approx_abundance-{dataset}.tsv".format(prefix=args.output, dataset=dataset.lower()))
    if output_path.exists():
        raise Exception("Output {} already exists...".format(output_path))
    magpipe.pretty.print_pair("Datatet {} located in".format(dataset), dataset_path, nl_before=1)
    samples = [s.name for s in dataset_path.glob('*')]
    abundances = []
    for s in samples:
        # Before going, we can only iterate over the mags derived from that sample so let's filter that
        if sample_map:
            filtered_mags = [m for m in mags if sample_map[s] in m]
        else:
            filtered_mags = [m for m in mags if s in m]
        if sample_map:
            magpipe.pretty.print_pair('Processing {} genomes from sample'.format(len(filtered_mags)), "{old} -> {new}".format(old=s, new=sample_map[s]))
        else:
            magpipe.pretty.print_pair('Processing {} genomes from sample'.format(len(filtered_mags)), s)
        if len(filtered_mags) > 0:
            depth_file = list(dataset_path.joinpath(s, 'depth_files').glob(s + "*" + suffix))
            if len(depth_file) != 1:
                raise ValueError("It seems there was something wrong with finding the depth file for sample {}. Found: {}".format(s, depth_file))
            depth_table = pandas.read_csv(depth_file[0], sep='\t')
            depth_table.set_index('contigName', inplace=True)
            # Now let's use the table to get the mags abundances
            for mag in filtered_mags:
                # then we want to subset the table for only the relevant scaffolds
                mag_membership = membership.loc[mag]
                if isinstance(mag_membership, pandas.DataFrame):
                    scaffolds = mag_membership["scaffold"].to_list()
                elif isinstance(mag_membership, pandas.Series): # Needs a special case if the genome has only 1 scaffold, which returns a Series and then scaffold as a string 
                    scaffolds = [mag_membership["scaffold"]]
                else:
                    raise Exception("Something went wrong with mag {}. Expected to get a pandas DataFrame or Series, got {}.".format(mag, type(mag_membership))) 
                if scaffold_map:
                    scaffolds = [scaffold_map[s] for s in scaffolds]
                subset_table = depth_table[select_columns(s, depth_table.columns)]
                if sample_map:
                    subset_table.columns = [sample_map[col.split('_vs_')[0]] for col in subset_table.columns]
                else:
                    subset_table.columns = [col.split('_vs_')[0] for col in subset_table.columns]
                mag_abundance = subset_table.loc[scaffolds].median(axis=0).to_dict()
                mag_abundance['genome'] = mag
                abundances.append(mag_abundance)
    df = pandas.DataFrame(abundances)
    df.set_index('genome', inplace=True)
    df.round(4).to_csv(output_path, sep='\t')
    # Count the number of samples where approx abundance is >= index
    df['diff_coverage_index'] = [int(sum([a >= index for a in df.loc[i]])) for i in df.index]
    df = pandas.DataFrame(df['diff_coverage_index']) # Otherwise a 1 col selection gives a series
    return df

# Main function of the script --------------------------------------------------

def main(args):
    """
    The main function, loops through the mags and the datasets to extract the relevant information
    """
    datasets_path = pathlib.Path(args.datasets)
    if not datasets_path.exists():
        raise FileNotFoundError("Looks like there is something wrong with {}.".format(datasets_path))
    # Load the mapping files
    sample_map = load_sample_mapping(args.sample_mapping)
    scaffold_map = load_scaffold_mapping(args.scaffold_mapping)
    # Start by loading the scaffold membership
    membership = pandas.read_csv(args.input, sep='\t', header=None, names=["scaffold", "genome"])
    # Set genomes as index
    membership.set_index('genome', drop=False, inplace=True)
    # Remove external genomes (i.e. genomes not reconstructed from the relevant datasets)
    if not all([isinstance(v, list) for v in args.dataset_abbrev.values()]):
        raise TypeError("Expected the dataset abbreviations to contain lists... got {}".format(args.dataset_abbrev))
    all_abbrevs = list(itertools.chain.from_iterable(args.dataset_abbrev.values()))
    membership = membership.loc[[genome.split('_')[0] in all_abbrevs for genome in membership.index]]
    # Quick fix to remove external MAGs (the diff is that they should have numbers in the MAG unique ID)
    membership = membership.loc[[not bool(re.search('[0-9]', genome.split('_')[-1])) for genome in membership.index]]
    # now we can get the list of mags
    mags = list(set(membership.index))
    magpipe.pretty.print_pair("Numer of MAGs to process:", len(mags))
    # Init df to store the mag + differential coverage index information
    res = pandas.DataFrame()
    # Loop over the datasets
    for dataset in args.dataset_abbrev.keys():
        magpipe.pretty.print_pair("Processing dataset", dataset, nl_before=1)
        df = process_dataset(dataset_path=datasets_path.joinpath(dataset),
                             mags=[m for m in mags if m.split('_')[0] in args.dataset_abbrev[dataset]],
                             suffix=args.suffix,
                             membership=membership,
                             scaffold_map=scaffold_map,
                             sample_map=sample_map,
                             index=args.index)
        res = res.append(df)
    output_diff_index = pathlib.Path("{prefix}-diff_cov_index.tsv".format(prefix=args.output))
    res.to_csv(output_diff_index, sep='\t')

# Run script -------------------------------------------------------------------

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
