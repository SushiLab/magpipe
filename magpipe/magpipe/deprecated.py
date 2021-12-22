# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Functions that were used at some point but lost relevance
"""

from pathlib import Path

def get_repository(repo_name='Ocean_MAGs'):
    """
    This functions finds the folder in the path of the file to infer the path of the overall repo
    The idea was to directly infer the output directory (repo/results)
    """
    file_path = Path(__file__).resolve()
    if repo_name in file_path.parts:
        repo_index = file_path.parts.index(repo_name)
    else:
        raise Exception('Directory of interest not in path...')
    repo_path = file_path.parents[(len(file_path.parts) - 1) - (repo_index + 1)]
    if repo_path.parts[-1] is not repo_name:
        raise NotImplementedError('Oops... wrong repo path, I made a mistake somewhere...')
    return repo_path

