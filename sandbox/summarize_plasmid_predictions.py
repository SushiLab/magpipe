#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Summarizes the plasmid/viral predictions of the scaffolds
## FIXME: Initialize columns when results table may be empty or it breaks dowstream with key error
i.e. ccontig, plasmidfinder and virsorter?
"""

import json
import pathlib
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm

import magpipe.log
import magpipe.pretty

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser()
    # REQUIRED arguments:
    parser.add_argument('-r', '--results_path', help = 'Path of the results directory.', required = True, type = pathlib.Path)
    parser.add_argument('-d', '--datasets', nargs = '+', help = 'List of the datasets to process.', required = True, type = str)
    parser.add_argument('-o', '--output', help = 'Name of the summary plasmid prediction table.', required = True, type = str)
    # Optional arguments
    parser.add_argument('-s', '--samples', help = 'List of samples to override the automatic lookup.', required = False, type = pathlib.Path, default = None)
    return parser.parse_args()

def read_ccontig(table_path):
    """
    Read the output of ccontig for a given sample
    and return it in a nicely formatted table
    """
    res = pd.read_csv(table_path, sep = "\t", names = ["scaffold_name", "ccontig_score"])
    return(res)

def read_cbar(table_path):
    """
    Read the output of cbar for a given sample
    and return it in a nicely formatted table
    """
    res = pd.read_csv(table_path, sep = "\t")
    res = res.rename(columns = {"#SeqID": "scaffold_name", "Prediction": "cbar_prediction"})
    res = res.drop(columns = ["Length"])
    return(res)

def read_plasflow(table_path):
    """
    Read the output of plasflow for a given sample
    and return it in a nicely formatted table
    """
    res = pd.read_csv(table_path, sep = "\t", index_col = 0)
    res = res.drop(res.columns.difference(['contig_name', 'contig_length', 'label']), 1)
    res = res.rename(columns = {"contig_name": "scaffold_name", "contig_length": "scaffold_length"})
    #split existing label column in two new ones
    new_labels = res["label"].str.split(".", n = 1, expand = True)
    res["plasflow_prediction"] = new_labels[0]
    res["plasflow_taxonomy"] = new_labels[1]
    res = res.drop(columns = ["label"])
    return(res)

def read_plasmidfinder(table_path):
    """
    Read the output of read_plasmidfinder for a given sample
    and return it in a nicely formatted table
    """
    #FIXME failsafe if duplicated hits?
    res = pd.read_csv(table_path, sep = "\t")
    res = res.drop(columns = ["Position in contig", "Note", "Accession number"])
    res = res.rename(columns = {"Database": "plasmidfinder_database", "Plasmid": "plasmidfinder_plasmid",
    "Identity": "plasmidfinder_identity", "Query / Template length": "plasmidfinder_coverage", "Contig": "scaffold_name"})
    if len(res.index) > 0:
        # streamlining the scaffold names
        res['scaffold_name'] = res['scaffold_name'].str.split(' ').str[0]
    return(res)

def read_virsorter(table_path):
    """
    Read the output of virsorter for a given sample
    and return it in a nicely formatted table
    """
    res = pd.read_csv(table_path, sep = ",", skiprows = 1, engine = "python")
    res = res.drop(res.columns.difference(['## Contig_id', 'Category', 'Fragment', 'Nb genes', 'Nb phage hallmark genes']), 1)
    res = res.rename(columns = {"## Contig_id": "scaffold_name", "Category": "virsorter_score",
    "Fragment": "virsorter_fragment", "Nb genes": "virsorter_fragment_n_genes", "Nb phage hallmark genes": "virsorter_n_hallmark_genes"})
    # Drop rows starting with '## Contig_id'
    res = res[~res.scaffold_name.str.contains("## Contig_id")]
    category_indexes = list(res.index[res.scaffold_name.str.startswith('##')])
    res.loc[0:category_indexes[2], "virsorter_category"] = "phage"
    res.loc[category_indexes[2]:res.index[-1], "virsorter_category"] = "prophage"
    res = res[~res.scaffold_name.str.startswith("##")]
    # streamline names
    res['scaffold_name'] = res['scaffold_name'].str.split('VIRSorter_').str[1]
    res['scaffold_name'] = res['scaffold_name'].str.split('_type=').str[0]
    return(res)

def read_deepvirfinder(table_path):
    """
    Read the output of deepvirfinder for a given sample
    and return it in a nicely formatted table
    """
    res = pd.read_csv(table_path, sep = "\t")
    res = res.rename(columns = {"name": "scaffold_name", "score": "deepvirfinder_score", "pvalue": "deepvirfinder_pvalue"})
    res['scaffold_name'] = res['scaffold_name'].str.split(' type=').str[0]
    res = res.drop(columns = ["len"])
    return(res)

def add_eukrep(table, eukrep_output):
    """
    Takes a table with the results of plasmid predictions
    and adds the eukrep information (i.e. euk or prok) to it
    """
    with open(str(eukrep_output) + '.euk') as f:
        euks = f.read().splitlines()
    with open(str(eukrep_output) + '.prok') as f:
        proks = f.read().splitlines()
    table = table.set_index('scaffold_name', drop = False)
    table.loc[euks, "eukrep"] = "Eukarya"
    table.loc[proks, "eukrep"] = "Prokarya"
    table = table.reset_index(drop = True)
    return(table)

def add_classification(table):
    """
    Add a classification column to the plasmid prediction table
    We consider plasmidfinder prediction highly specific, so sufficient for the prediction
    We additionally consider that if it's predicted as plasmid by both cbar and
    plasflow without hits for viruses (FIXME and eukaryotes...) we predict a plasmid
    i.e. congruent evidence without contradicting evidence
    Note: to test for virsorter_score, keep in mind that due to the input the type is str
    """
    plasmids = (table.loc[:, 'plasmidfinder_identity'] >= 80) | \
    ((table.loc[:, 'cbar_prediction'] == "Plasmid") & (table.loc[:, 'plasflow_prediction'] == "plasmid") & \
    (table.loc[:, 'deepvirfinder_pvalue'] > 0.05) & (table.loc[:, 'virsorter_score'] != "1") & (table.loc[:, 'virsorter_score'] != "2"))
    """
    We create a more lenient category for putative plasmids
    if there is some evidence that this is a plasmid, and not contradicting information
    Note that a virsorter score of 3 often predicts plasmids
    i.e. some evdence and no strong contradiction
    """
    putative_plasmids = (table.loc[:, 'plasmidfinder_identity'] >= 50) | (((table.loc[:, 'cbar_prediction'] == "Plasmid") | (table.loc[:, 'plasflow_prediction'] == "plasmid") | (table.loc[:, 'virsorter_score'] == 3)) & \
        (table.loc[:, 'plasflow_prediction'] != "chromosome") & \
        (table.loc[:, 'deepvirfinder_pvalue'] > 0.01) & \
        (table.loc[:, 'virsorter_score'] != "1") & (table.loc[:, 'virsorter_score'] != "2"))
    """
    We also add a case for any viral signal
    """
    viral_signal = (table.loc[:, 'deepvirfinder_pvalue'] <= 0.05) | (table.loc[:, 'virsorter_score'] == "1") | (table.loc[:, 'virsorter_score'] == "2") | (table.loc[:, 'virsorter_score'] == "3")
    # Add that to the table, note that the order matters
    table.loc[viral_signal, 'plasmid_prediction'] = 'Viral_Signal'
    table.loc[putative_plasmids, 'plasmid_prediction'] = 'Putative_Plasmid'
    table.loc[plasmids, 'plasmid_prediction'] = 'Plasmid'
    # Correct wrong prediction in case eukrep doesn't predict a prokaryotic sequence
    # i.e. pedicts Eukarya or tie (None)
    table.loc[(table.loc[:, 'eukrep'] != "Prokarya"), 'plasmid_prediction'] = None
    # write result to file
    return table

def main(args):
    """
    Loop through datasets and sample to merge and concatenate all the results
    Finalise the plasmid prediction
    """
    output = pathlib.Path(args.output)
    if output.exists():
        raise ValueError("Output ({}) already exists...".format(output))
    ordered_columns = ['scaffold_name', 'scaffold_length', 'ccontig_score', "eukrep",
                       'plasflow_prediction', 'plasflow_taxonomy', 'cbar_prediction',
                       'plasmidfinder_database', 'plasmidfinder_plasmid', 'plasmidfinder_identity', 'plasmidfinder_coverage',
                       'deepvirfinder_score', 'deepvirfinder_pvalue',
                       'virsorter_score', 'virsorter_category', 'virsorter_fragment', 'virsorter_fragment_n_genes', 'virsorter_n_hallmark_genes', 
                       'plasmid_prediction']
    res_df = pd.DataFrame(columns = ordered_columns)
    res_df.to_csv(output, sep = '\t', index = False)
    # Loop over datasets
    for dataset in args.datasets:
        dataset_path = args.results_path.joinpath(dataset)
        if args.samples:
            magpipe.pretty.print_pair("Looking for samples given by", args.samples)
            with open(str(args.samples)) as f:
                samples = f.read().splitlines()
        else:
            magpipe.pretty.print_pair("Looking for samples in", dataset_path, nl_before=1)
            samples = [sample.name for sample in dataset_path.glob('*')]
        # Loop over the samples of the dataset
        magpipe.pretty.print_pair("Number of samples to loop through", len(samples))
        for sample in tqdm(samples, ncols = 100):
            ccontig_df = read_ccontig(args.results_path.joinpath(dataset, sample, 'ccontig', sample + '-ccontig.tsv'))
            cbar_df = read_cbar(args.results_path.joinpath(dataset, sample, 'cbar', sample + '-cbar.txt'))
            plasflow_df = read_plasflow(args.results_path.joinpath(dataset, sample, 'plasflow', sample + '-plasflow.tsv'))
            plasmidfinder_df = read_plasmidfinder(args.results_path.joinpath(dataset, sample, 'plasmidfinder', sample + '-plasmidfinder', 'results_tab.tsv'))
            dvf_files = list(args.results_path.joinpath(dataset, sample, 'deepvirfinder', sample + '-deepvirfinder/').glob('*dvfpred.txt'))
            if len(dvf_files) > 1:
                raise ValueError("looks like there was more than one result file here... {}".format(sample))
            deepvirfinder_df = read_deepvirfinder(dvf_files[0])
            virsorter_df = read_virsorter(args.results_path.joinpath(dataset, sample, 'virsorter', sample + '-virsorter', 'VIRSorter_global-phage-signal.csv'))
            # merging all the tools
            # Start with plasflow but reorder
            sample_predictions = pd.merge(plasflow_df, ccontig_df, how = "outer", on = "scaffold_name")
            sample_predictions = pd.merge(sample_predictions, cbar_df, how = "outer", on = "scaffold_name")
            sample_predictions = pd.merge(sample_predictions, plasmidfinder_df, how = "outer", on = "scaffold_name")
            sample_predictions = pd.merge(sample_predictions, deepvirfinder_df, how = "outer", on = "scaffold_name")
            sample_predictions = pd.merge(sample_predictions, virsorter_df, how = "outer", on = "scaffold_name")
            # Add eukrep results
            sample_predictions = add_eukrep(sample_predictions, args.results_path.joinpath(dataset, sample, 'eukrep', sample))
            # Append the sample-specific table to the overall summary
            # Reorder the columns
            sample_predictions = add_classification(sample_predictions)
            sample_predictions = sample_predictions[ordered_columns]
            #res_df = res_df.append(sample_predictions, ignore_index = True)
            sample_predictions.to_csv(output, sep = "\t", index = False, header = False, mode = 'a')

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
