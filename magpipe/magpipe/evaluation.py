# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Functions to manipulate and summarise the output of magpipe evaluation rules
"""

import shutil
import pandas
from pathlib import Path

import magpipe.drep
import magpipe.pretty
import magpipe.readers
import magpipe.files_and_paths

def wrapper_export_evaluation(output_dir, results_dir, datasets, datasets_table, methods_table, filter_table=True, completeness=50, contamination=10, anvio=True, export_sequences=True, table_for_drep=True):
    """
    This functions is directly called by the magpipe script and is wrapped around the different
    functons of the magpipe.evaluation module.
    The goal is to summarise the binning efforts of the phase 1 and move to phase 2.
    Effectvely, it's about going through all the metagenomic samples of the provided datasets,
    identify the MAGs (passing filtering) and export their sequences as well as create two summary
    tables, one with all the information and the other formatted for drep
    """
    # Some inital testing
    if not isinstance(datasets_table, pandas.DataFrame):
        raise TypeError("Expects datasets to be provided as a pandas.DataFrame and got {}.".format(type(datasets_table)))
    if not isinstance(methods_table, dict):
        raise Exception("If you're planning on filtering the bins, exporting the MAGs and moving on with analyses, I would have expected you to settle on a given binning method. Yet I got {} instead of a dictionary.".format(type(methods_table)))
    if not filter_table:
        completeness = "NaN"
        contamination = "NaN"
    # format input/output
    results_path = Path(results_dir)
    output_path = Path(output_dir)
    output_prefix = "{MAGs_name}-{method}{suffix}-cpl{cpl}_ctn{ctn}".format(MAGs_name=output_path.name,
                                                                            method=methods_table['tool'],
                                                                            suffix=methods_table['suffix'],
                                                                            cpl=completeness,
                                                                            ctn=contamination)
    output_evaluation_table = output_path.joinpath(output_prefix + "-evaluate_summary.tsv")
    output_drep_table = output_path.joinpath(output_prefix + "-evaluate_for_drep.csv")
    output_bins_sequences = output_path.joinpath(output_prefix + "-genomes")
    # Checking input/output
    if output_evaluation_table.exists():
        raise Exception('The targeted output file {} exists, cannot handle that so just exiting.'.format(output_evaluation_table))
    if output_drep_table.exists():
        raise Exception('The targeted output file {} exists, cannot handle that so just exiting.'.format(output_drep_table))
    if output_bins_sequences.exists():
        raise Exception("The target folder to copy the MAGs sequences {} already exists, probably been copied already. Please check.".format(output_bins_sequences))
    if not results_path.is_dir():
        raise FileNotFoundError('The prodived results directory {} is not a directory'.format(results_dir))
    if not isinstance(datasets, list):
        raise TypeError("Datasets must be given as a list, even if only a single one.  You gave a {}.".format(type(datasets)))
    # Now we are ready to loop through datasets and get things going
    for dataset in datasets:
        magpipe.pretty.print_pair("Processing dataset", dataset, nl_before=1)
        method = "{tool}_a{depth}{suffix}".format(tool=methods_table['tool'],
                                                  depth=int(datasets_table.loc[dataset, 'diffcov']),
                                                  suffix=methods_table['suffix'])
        dataset_export_evaluation(dataset_path=results_path.joinpath(dataset),
                                  method=method,
                                  output_evaluation_table=output_evaluation_table,
                                  output_drep_table=output_drep_table,
                                  output_bins_sequences=output_bins_sequences,
                                  filter_table=filter_table,
                                  completeness=completeness,
                                  contamination=contamination,
                                  anvio=anvio,
                                  export_sequences=export_sequences,
                                  table_for_drep=table_for_drep)

def dataset_export_evaluation(dataset_path, method, output_evaluation_table, output_drep_table=None, output_bins_sequences=None, filter_table=True, completeness=50, contamination=10, anvio=True, export_sequences=True, table_for_drep=True):
    """
    Summarises the evaluation of a binning effort for a given dataset
    1 - Recursively go through the metagenomic samples
    2 - Identify checkm summaries, anvio summaries and checkm analyses as specified
    3 - Merge the three tables
    4 - Update completeness and contamination as the mean of both checkm and anvio if anvio avail
    5 - Optionally filter the table based on completeness and contamination
    6 - Append the dataframe to the target file
    7 - Optionally append the drep format dataframe to the corresponding file
    8 - Optionally export the fasta files to the target path
    """
    # quick check
    if not all([isinstance(anvio, bool), isinstance(anvio, bool), isinstance(table_for_drep, bool), isinstance(filter_table, bool), isinstance(export_sequences, bool)]):
        raise TypeError("Needs a boolean to know whether to look for anvio results, generate drep table, filter or export bins. Got {}.".format([type(anvio), type(table_for_drep), type(filter_table), type(export_sequences)]))
    # Iterate over the samples
    n_samples = 0
    for sample_dir in dataset_path.iterdir():
        n_samples += 1
        sample_path = Path(sample_dir)
        sample = sample_path.name
        checkm_summary = sample_path.joinpath(method, "{sample}-{method}-checkm.summary".format(sample=sample, method=method))
        checkm_analysis = sample_path.joinpath(method, 'checkm_lineage', 'storage', 'bin_stats.analyze.tsv')
        if not checkm_summary.exists() or not checkm_analysis.exists():
            magpipe.pretty.print_single("Couldn't find one of the CheckM files for sample {}. Please check {} and {}. I will skip this sample for now, assuming there was no bins".format(sample, checkm_summary, checkm_analysis), color = "yellow") # Need a warning as some samples may not have them (if no bins)
            continue
        table_checkm_summary = magpipe.readers.read_checkm_summary(checkm_summary)
        table_checkm_analysis = magpipe.readers.read_checkm_analysis(checkm_analysis)
        table_summary = table_checkm_summary.merge(table_checkm_analysis, left_on='Bin Id', right_index=True)
        if anvio:
            anvio_summary = sample_path.joinpath(method, "{sample}-{method}-anvi_evaluate.summary".format(sample=sample, method=method))
            if not anvio_summary.exists():
                raise FileNotFoundError("Couldn't find anvio summary for sample {}. Please check {}.".format(sample, anvio_summary))
            table_anvio_summary = magpipe.readers.read_anvio_summary(anvio_summary, sample=sample, method=method)
            table_summary = table_summary.merge(table_anvio_summary, left_on='Bin Id', right_on='Bin Id')
            # Some sanity check
            check_size = all([table_summary.loc[index, 'CheckM Genome size'] == table_summary.loc[index, 'Anvio Length'] for index in table_summary.index])
            check_nscaff = all([table_summary.loc[index, 'CheckM # scaffolds'] == table_summary.loc[index, 'Anvio # Scaffolds'] for index in table_summary.index])
            if not check_size and check_nscaff:
                raise Exception('CheckM and Anvio genome size and/or number of scaffolds don\'t match... shouldn\'t happen.')
            # Generate mean values, convenient and needed for drep
            mean_contamination = table_summary[['CheckM Contamination', 'Anvio Redundancy']].mean(axis=1).round(2)
            mean_completeness = table_summary[['CheckM Completeness', 'Anvio Completion']].mean(axis=1).round(2)
            table_summary.insert(loc=1, column='Mean Contamination', value=mean_contamination)
            table_summary.insert(loc=1, column='Mean Completeness', value=mean_completeness)
        # Add the location
        table_summary.loc[:, "Folder Location"] = str(sample_path.joinpath(method))
        # If asked for, filter the table
        if filter_table:
            table_summary = filter_bins(table_summary, completeness=completeness, contamination=contamination, anvio=anvio)
        # Append the sample table to file
        with open(str(output_evaluation_table), 'a') as f:
            table_summary.to_csv(f, header=(f.tell() == 0), sep='\t', mode='a', index=False)
        # Generate the drep table if required
        if table_for_drep:
            table_summary_for_drep = magpipe.drep.format_evaluate_for_drep(table_summary)
            with open(str(output_drep_table), 'a') as f:
                table_summary_for_drep.to_csv(f, header=(f.tell() == 0), mode='a', index=False)
        # Export the fasta files if necessary
        if export_sequences:
            export_bins(table_summary, output_bins_sequences)
    magpipe.pretty.print_pair("Number of samples processed", n_samples)

def filter_bins(table, completeness=50, contamination=10, anvio=False):
    """
    Filter bins based on CheckM Completeness and Contamination
    if Anvio is provided, will filter on an either or basis
    """
    # Note use | and & for compatiblity with pandas series ('and' and 'or' don't work)
    if anvio:
        pass_checkm = (table.loc[:, 'CheckM Completeness'] >= completeness) & (table.loc[:, 'CheckM Contamination'] <= contamination)
        pass_anvio = (table.loc[:, 'Anvio Completion'] >= completeness) & (table.loc[:, 'Anvio Redundancy'] <= contamination)
        pass_quality = table.index[pass_checkm | pass_anvio]
    else:
        pass_quality = (table.loc[:, 'CheckM Completeness'] >= completeness) & (table.loc[:, 'CheckM Contamination'] <= contamination)
    table = table.loc[pass_quality]
    return table

def export_bins(table, output_path):
    """
    Takes the evaluate summary table as input (from summarise_dataset_evaluation)
    Extracts the positions of the fasta file from Sample Folder and Bin Id
    Export the '.fa' bins from the filtered table into the output_path
    """
    # Prepare the output directory
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    # Prepare the origins:
    paths_list = table['Folder Location'].tolist()
    bins_list = table['Bin Id'].tolist()
    origins = [Path(p).joinpath('bins/' + b + '.fa') for p, b in zip(paths_list, bins_list)]
    destinations = [output_path.joinpath(b + '.fa') for b in bins_list]
    if any([ori.name != des.name for ori, des in zip(origins, destinations)]):
        raise Exception("Origins bin names don't match destinations bin names")
    for ori, des in zip(origins, destinations):
        shutil.copy(str(ori), str(des))
