#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long


import argparse
import pandas as pd

import magpipe.log
import magpipe.evaluation

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser()
    # REQUIRED arguments:
    parser.add_argument('-r', '--results_dir', help='Path of the results directory.', required=True, type=str)
    parser.add_argument('-o', '--output_dir', help='Path of the directory where to write the output files. Base folder name is used as prefix.', required=True, type=str)
    parser.add_argument('-d', '--datasets', nargs='+', help='List of the dataset to process.', required=True, type=str)
    parser.add_argument('-m', '--method', help='Method used for the binning step (you need to pick just one).', required=True, type=str)
    # OPTIONAL arguments:
    parser.add_argument('--filter', dest='filter_table', action='store_true', help='Filter the table based on completeness and contamination (default doesn\'t).')
    parser.add_argument('--completeness', help="Minimum completeness percentage of a bin to be exported", type=int, default=50)
    parser.add_argument('--contamination', help="Maximum contamination score of a bin to be exported", type=int, default=10)
    parser.add_argument('--anvio', dest='anvio', action='store_true', help='Summarise anvi summaries (default).')
    parser.add_argument('--no_anvio', dest='anvio', action='store_false', help='Don\'t summarise anvio summaries.')
    parser.add_argument('--export_sequences', dest='export_sequences', action='store_true', help='Export the sequences (default doesn\'t).')
    parser.add_argument('--table_for_drep', dest='table_for_drep', action='store_true', help='Generate the table for drep (default doesn\'t).')
    parser.add_argument('--datasets_table', help="Table with datasets details.", type=str, default='/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/DEV/MAGPIPE_DEV/code/resources/datasets.tsv')
    parser.add_argument('--methods_table', help="Table with methods details.", type=str, default='/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/DEV/MAGPIPE_DEV/code/resources/implemented_methods.tsv')
    parser.set_defaults(filter_table=False, anvio=True, export_sequences=False, table_for_drep=False)
    # parse and return the arguments
    return parser.parse_args()

def main(args):
    """
    Read the checkm and anvio summary tables for the listed datasets
    merge and export the result, including the bins passing filtering
    write a formatted version of the checkm table for drep
    """
    methods_table = pd.read_csv(args.methods_table, sep='\t').set_index('method')
    methods_table.fillna("", inplace=True)
    methods_table = methods_table.loc[args.method].to_dict()
    datasets_table = pd.read_csv(args.datasets_table, sep='\t').set_index('datasets')
    magpipe.evaluation.wrapper_export_evaluation(output_dir=args.output_dir,
                                                 results_dir=args.results_dir,
                                                 datasets=args.datasets,
                                                 datasets_table=datasets_table,
                                                 methods_table=methods_table,
                                                 filter_table=args.filter_table,
                                                 completeness=args.completeness,
                                                 contamination=args.contamination,
                                                 anvio=args.anvio,
                                                 export_sequences=args.export_sequences,
                                                 table_for_drep=args.table_for_drep)

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
