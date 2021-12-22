# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Functions to read the output of programs used in the magpipe snakemake pipeline
when they don't have a dedicated submodule
"""

import re
import ast
import pandas
from pathlib import Path

def read_checkm_summary(checkm_file):
    """
    Function to read the summary output of CheckM
    Try to read it as tsv, if it doesn't match, read as space delimited.
    This enables compatibility between versions of the rule that specified the
    --tab_table option or not
    return standardized pandas dataframe.
    """
    # Define the general expectations, may need to be updated if checkm changes...
    expected_names = ["Bin Id", "Marker lineage", "# genomes", "# markers", "# marker sets", "0", "1", "2", "3", "4", "5+", "Completeness", "Contamination", "Strain heterogeneity"]
    wanted_names = ["Bin Id", "Marker lineage", "lineage id", "# genomes", "# markers", "# marker sets", "0", "1", "2", "3", "4", "5+", "Completeness", "Contamination", "Strain heterogeneity"]
    # try reading a tsv
    table = pandas.read_csv(checkm_file, sep='\t')
    if list(table.columns) == expected_names:
        # Split "Marker lineage" in two columns
        splitted = splitted = table["Marker lineage"].str.split(' ', n=1, expand=True)
        table.loc[:, 'Marker lineage'] = splitted[0]
        table.insert(loc=2, column='lineage id', value=splitted[1])
    else: # Read as space delim...
        table = pandas.read_csv(checkm_file, delim_whitespace=True, skiprows=3, header=None)
        # Use the internal definition of the columns names...
        table.columns = wanted_names
        nan_row = table.index[table.iloc[:, 1:].isnull().all(1)] # find rows where all the non first column values are nans
        if len(nan_row) == 1 and nan_row[0] == (table.shape[0] - 1):
            table = table.drop(nan_row)
        else:
            raise Exception("This file {} looks bad.\nNot just one empty row or empty row is not last row... shoudln\'t happen.".format(checkm_file))
    # Add 'CheckM ' in front of columns beside 'Bin Id' for desambiguation
    table.columns = ["CheckM {}".format(name) if name != "Bin Id" else name for name in table.columns]
    return table

def read_checkm_analysis(analysis_file):
    """
    Function to read the analysis output of checkm.
    Reads the checkm bin_analysis file as a dataframe.
    Initially reads it as a dataframe and expands the dictionary into the columns.
    """
    init_table = pandas.read_csv(analysis_file, sep='\t', header=None)
    table = pandas.DataFrame()
    for i in init_table.index:
        row = init_table.loc[i]
        index = [row[0]]
        init_dict = ast.literal_eval(row[1])
        df = pandas.DataFrame(init_dict, index=index)
        table = table.append(df)
    # Add 'CheckM ' in front of columns beside 'Bin Id' for desambiguation
    table.columns = ["CheckM {}".format(name) if name != "Bin Id" else name for name in table.columns]
    return table

def read_anvio_summary(anvio_file, sample, method, external=False, ext_pattern='.contigs'):
    """
    Reads the anvio summary file
    test for the two anvio output versions and returns 7 colunns
    Converts the Bin Id back to orginal name
    Handles external cases
    """
    col_rename = ['Anvio Bin Name', 'Anvio Domain', 'Anvio Domain Confidence', 'Anvio Completion', 'Anvio Redundancy', 'Anvio # Scaffolds', 'Anvio Length']
    # Read table
    table = pandas.read_csv(anvio_file, sep='\t')
    if len(table.columns) == 7:
        table.columns = col_rename
    elif len(table.columns) == 6:
        new_domain = table["domain"].str.split(" ", n=1, expand=True)
        table.loc[:, "domain"] = new_domain[0]
        table.insert(loc=2, column='confidence', value=[float(c.lstrip(r'\(').rstrip(r'\)')) for c in list(new_domain[1])]) # Use r'string' with backslashes
        table.columns = col_rename
    else:
        raise Exception("Unexpected column behavior in file {}.\nExpected 6 or 7 columns and got {}.".format(anvio_file, len(table.columns)))
    # Create the matching Bin Id
    # Prepare method strings to rename method + bin number
    # Also remove the first 'BIN_' artificially added (anvio wants weird things sometimes)
    if not external:
        anvio_method = method.strip('_').strip('.') + '_'
        anvio_method = anvio_method.replace('.', '_')
        method = method.strip('_').strip('.') + '.'
        sample_fixed = re.sub("[^0-9a-zA-Z]+", "_", sample)
        if not all(anvio_method in bin_name for bin_name in list(table['Anvio Bin Name'])):
            raise Exception('Anvio transformed method provided not in all anvio bin names... which it should.')
        table.insert(loc=0, column='Bin Id', value=[bin.replace(sample_fixed, sample).replace(anvio_method, method).replace('BIN_', '', 1) for bin in list(table['Anvio Bin Name'])])
    elif external:
        table.insert(loc=0, column='Bin Id', value=[bin.replace(ext_pattern, '') for bin in list(table['Anvio Bin Name'])])
    return table

def read_append_tables(list_of_files, read_type='checkm_summary', method='method', external=False, ext_pattern='.contigs'):
    """
    Takes a list of files and a function to read them with pandas
    and merge all the tables to return a single one
    Supports different type of files, listed in 'available'
    Generally deprecated because pandas append is too slow, faster to append to file
    But used in specific cases, e.g. external genomes
    """
    available = ['checkm_summary', 'checkm_analysis', 'anvio_summary']
    if read_type not in available:
        raise Exception("Read type must be in {}, and you've given {}.".format(available, read_type))
    res = pandas.DataFrame()
    for table_file in list_of_files:
        # FIXME Maybe improve by defining read_fun outside loop, using partial for anvio's method
        if read_type == 'checkm_summary':
            table = read_checkm_summary(table_file)
        if read_type == 'checkm_analysis':
            table = read_checkm_analysis(table_file)
        if read_type == 'anvio_summary':
            table = read_anvio_summary(table_file, method, external, ext_pattern)
        init_len = len(res.index)
        table_len = len(table.index)
        res = res.append(table)
        if len(res.index) != init_len + table_len:
            raise Exception("Summary table as fewer or more rows than the two compiled tables.\nWas processing {} when this happened.".format(table_file))
    return res
