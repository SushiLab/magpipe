# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Functions related to files and directory manipulations
"""

import re
from pathlib import Path

def find_suffixed_files(path, file_suffix, sublevels=1, iter_subfolders=True):
    """
    This function takes a path as input and looks for files with a given suffix
    it should be able to look through a specified amound of sublevels and/or
    iterate through all the subfolders in the path. If both are specified, it
    will first iterate through the subfolders and look then add the sublevels
    Returns a list of the files found
    """
    # Some safety check
    path = Path(path)
    if not path.is_dir():
        raise FileNotFoundError('Cannot find the target directory: {}.'.format(path))
    if not isinstance(sublevels, int) and not isinstance(sublevels, str):
        raise TypeError("Argument 'sublevels' must be a integer or a string, you gave {}.".format(type(sublevels)))
    if not isinstance(iter_subfolders, bool):
        raise TypeError("Argument 'iter_subfolders' must be a boolean, you gave {}.".format(type(iter_subfolders)))
    # define the pattern to glob
    if isinstance(sublevels, int):
        pattern = "{}/*{}".format("/".join(['*'] * sublevels), file_suffix)
    elif isinstance(sublevels, str):
        pattern = "{}/*{}".format(sublevels, file_suffix)
    # remove doubled (or more) chars ('//' and '**') as well as leading '/'
    pattern = re.sub('/+', '/', pattern)
    pattern = re.sub(r'\*+', '*', pattern) # use r'string' when escaping characters for regex
    pattern = re.sub('^/', '', pattern)
    # Look for the files
    files_list = []
    if iter_subfolders:
        for subfolder in path.iterdir(): # go through subfolders first if specified
            for res_file in path.joinpath(subfolder.parts[-1]).glob(pattern):
                files_list.append(str(res_file))
    else:
        for res_file in path.glob(pattern):
            files_list.append(str(res_file))
    return files_list
