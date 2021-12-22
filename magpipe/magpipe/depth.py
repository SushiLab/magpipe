# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Functions to read, split and manipulate depth files (the output of metabat2
jgi script to extract abundances from bamfiles)
"""

import re
import pandas
from pathlib import Path

def read_depth_file(file_name):
    """
    TODO: implement
    """
    table = pd.read_csv(file_name, sep = '\t')
