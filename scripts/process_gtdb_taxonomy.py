#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

import sys
import pandas
import argparse
from pathlib import Path

import magpipe.log
import magpipe.pretty

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser()
    # REQUIRED arguments:
    parser.add_argument('-e', '--evalutation_table', help = 'Path of the evaluation table.', required = True, type = str)
    parser.add_argument('-g', '--gtdbtk_output', help = 'Path to the gtdbtk output, expects gtdb.ar122.summary.tsv etc in there.', required = True, type = str)
    parser.add_argument('-o', '--output', help = 'Path and prefix of the output files.', required = True, type = str)
    parser.add_argument('-m', '--completeness', help = 'Completeness value for filtering.', default = 30, type = int)
    parser.add_argument('-n', '--contamination', help = 'Contamination value for filtering.', default = 10, type = int)
    # return namespace
    return parser.parse_args()

def main(args):
    """
    Identifies Eukarya genomes based on Anvio and logs it so you can make sure it works
    Creates a list of mags per domain (Eukarya, Bacteria, Archaea).
    As well as two convenient lists All and Prokarya
    """
    evaluation_table = pandas.read_csv(args.evalutation_table, sep = '\t')
    # Check if the Anvio Domain Confidence table exists for retro-active compatibility (older anvio version)
    if not 'Anvio Domain Confidence' in evaluation_table.columns:
        raise ValueError("Couldn't find 'Anvio Domain Confidence' in the coluns of {}.".format(args.evaluation_table))
    # Identifies eukaryotes and simplify the table
    candidate_euks = evaluation_table.loc[evaluation_table["Anvio Domain"] == 'EUKARYA']
    sub_columns = ['Bin Id', 'Anvio Domain', 'Anvio Domain Confidence', 'CheckM Marker lineage', 'CheckM Completeness', 'CheckM Contamination', 'Anvio Completion', 'Anvio Redundancy', 'Anvio # Scaffolds', 'Anvio Length']
    candidate_euks = candidate_euks[sub_columns]
    candidate_euks.columns = ['Bin Id', 'Anvio Domain', 'Anvio Confidence', 'CheckM lineage', 'CheckM Compl.', 'CheckM Cont.', 'Anvio Compl.', 'Anvio Red.', '# Scaffolds', 'Length (bp)']
    # Print some information to evaluate the euk identification
    magpipe.pretty.print_pair('Printing Eukaryotes', 'All', nl_before =1)
    magpipe.pretty.print_table(candidate_euks)
    magpipe.pretty.print_pair('Printing Eukaryotes', 'Low confidence (Confidence >= 0.4)', nl_before =1)
    magpipe.pretty.print_table(candidate_euks[candidate_euks['Anvio Confidence'] <= 0.4])
    magpipe.pretty.print_pair('Printing Eukaryotes', 'Medium confidence (0.4 <= Confidence <= 0.6)', nl_before =1)
    magpipe.pretty.print_table(candidate_euks[(candidate_euks['Anvio Confidence'] <= 0.6) & (candidate_euks['Anvio Confidence'] >= 0.5)])
    magpipe.pretty.print_pair('Printing Eukaryotes', 'High confidence (Confidence >= 0.7)', nl_before =1)
    magpipe.pretty.print_table(candidate_euks[candidate_euks['Anvio Confidence'] >= 0.7])
    magpipe.pretty.print_single("From that, it appears legitimate to consider the list of eukaryotes to be EUK + Anvio pass, which removes the false positive (i.e. ChecM pass, Anvio not pass but EUK)", nl_before=1, nl_after=1)

    # Getting the eukarya list
    eukarya_table = evaluation_table.loc[(evaluation_table["Anvio Domain"] == 'EUKARYA') & (evaluation_table["Anvio Completion"] >= args.completeness) & (evaluation_table["Anvio Redundancy"] <= args.contamination)]
    eukarya_list = list(eukarya_table["Bin Id"])

    # Getting the complete list
    all_list = list(evaluation_table["Bin Id"])

    # Use gtdbtk output for the lists of other domains
    # Using classify instead of summary means you can manually run the script before gtdbtk finishes
    prokarya_table = evaluation_table.loc[~evaluation_table["Bin Id"].isin(eukarya_list)]
    prokarya_list = list(prokarya_table["Bin Id"])
    archaea = pandas.read_csv(args.gtdbtk_output + '/classify/intermediate_results/gtdbtk.ar122.classification_pplacer.tsv', sep = '\t', header = None)
    archaea_table = prokarya_table.loc[prokarya_table["Bin Id"].isin(list(archaea[0]))]
    archaea_list = list(archaea_table["Bin Id"])
    bacteria = pandas.read_csv(args.gtdbtk_output + '/classify/intermediate_results/gtdbtk.bac120.classification_pplacer.tsv', sep = '\t', header = None)
    bacteria_table = prokarya_table.loc[prokarya_table["Bin Id"].isin(list(bacteria[0]))]
    bacteria_list = list(bacteria_table["Bin Id"])

    # Some safety check
    test_list_lengths = (len(all_list) == len(eukarya_list) + len(archaea_list) + len(bacteria_list))
    test_prok_lengths = (len(prokarya_list) == len(archaea_list) + len(bacteria_list))
    test_table_list = (len(evaluation_table.index) == len(all_list))
    if not all([test_list_lengths, test_prok_lengths, test_table_list]):
        magpipe.pretty.print_single("Some MAGs are missing, please check that it makes sense based on the gtdbtk log...", nl_before=1, nl_after=1, color="yellow")

    # Some more log
    magpipe.pretty.print_pair("Number of Archaea", len(archaea_list))
    magpipe.pretty.print_pair("Number of Bacteria", len(bacteria_list))
    magpipe.pretty.print_pair("Number of Eukarya", len(eukarya_list))

    # Write the files
    magpipe.pretty.print_single('Writing the files...')
    Path(args.output).parent.mkdir(exist_ok=True)
    with open(args.output + '-all.txt', 'w') as f:
        f.write("\n".join(all_list) + "\n")
    with open(args.output + '-archaea.txt', 'w') as f:
        f.write("\n".join(archaea_list) + "\n")
    with open(args.output + '-bacteria.txt', 'w') as f:
        f.write("\n".join(bacteria_list) + "\n")
    with open(args.output + '-eukarya.txt', 'w') as f:
        f.write("\n".join(eukarya_list) + "\n")
    with open(args.output + '-prokarya.txt', 'w') as f:
        f.write("\n".join(prokarya_list) + "\n")

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
