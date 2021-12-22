# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Wrapper function to run a commandline
"""

import subprocess

import magpipe.pretty

def run_command(command: str):
    """
    Simple wrapper to run a command using check_call and catch exceptions
    :param command:
    :return:
    """
    magpipe.pretty.print_pair('Runnning command', command, nl_after = 1)
    try:
        returncode = subprocess.check_call(command, shell = True)
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError('Command {} failed with message:\t{}'.format(e.cmd, e.stderr))
    return returncode
