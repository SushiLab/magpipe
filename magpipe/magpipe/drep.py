# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Functions to prepare the magpipe evaluation summary for drep
"""

import gzip
import shutil
import pandas
from pathlib import Path

def format_evaluate_for_drep(input_table):
    """
    Takes an evaluation table as input and formats it for drep
    This means reducing the number of columns to:
    bin id, completeness, contamination, strain heterogeneity
    """
    if not isinstance(input_table, pandas.DataFrame):
        raise TypeError("Expected and pandas dataframe and got {}.".format(type(input_table)))
    if  all([c in input_table.columns for c in ["Mean Completeness", "Mean Contamination", "CheckM Strain heterogeneity"]]):
        formatted_table = input_table[['Bin Id', 'Mean Completeness', 'Mean Contamination', 'CheckM Strain heterogeneity']]
    elif all([c in input_table.columns for c in ["CheckM Completeness", "CheckM Contamination", "CheckM Strain heterogeneity"]]):
        formatted_table = input_table[['Bin Id', 'CheckM Completeness', 'CheckM Contamination', 'CheckM Strain heterogeneity']]
    else:
        raise ValueError("Expected Mean/CheckM Completeness, Contamination and Strain heterogeneity columns. Got {}.".format(input_table.columns))
    formatted_table.columns = ['genome', 'completeness', 'contamination', 'strain_heterogeneity']
    formatted_table.genome += ".fa"
    return formatted_table

def export_drep_data_tables(origin_path, dest_path, sublevel='data_tables', files=['genomeInformation.csv', 'Cdb.csv', 'Wdb.csv', 'Sdb.csv', 'Ndb.csv']):
    """
    A convenience function to export the drep data tables of interest,
    compressing them on the way
    It takes the drep output as an origin path and looks for the files in the
    sublevel folder
    """
    list_of_drep_tables = ["Bdb.csv", "Cdb.csv", "genomeInfo.csv", "genomeInformation.csv", "Mdb.csv", "Ndb.csv", "Sdb.csv", "Wdb.csv", "Widb.csv"]

    origin_path = Path(origin_path)
    dest_path = Path(dest_path)

    if not dest_path.is_dir():
        dest_path.mkdir(parents=True)

    for f in files:
        o = origin_path.joinpath(sublevel, f)
        if not o.is_file():
            raise FileNotFoundError("File {} doesn't seem to exist.\nIs it a standard drep table {}?".format(o, list_of_drep_tables))
        d = dest_path.joinpath(f + '.gz')
        with open(o, 'rb') as f_in, gzip.open(d, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
