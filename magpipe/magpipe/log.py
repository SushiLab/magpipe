# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Logging functions for the magpipe module
TODO: implement the logging
"""

import sys
import shutil
import logging
import argparse
from pathlib import Path

import magpipe.pretty

def print_arguments(args):
    """
    print arguments for logging and debugging
    """
    if not type(args) is argparse.Namespace:
        raise TypeError("Expects an argparse.Namespace and got {}...".format(type(args)))
    magpipe.pretty.print_header("Printing arguments", nl_after = 0)
    for k, v in args.__dict__.items():
        magpipe.pretty.print_pair(k, v)
    print()

def startup():
    logging.basicConfig(format = '%(asctime)s %(levelname)s: %(message)s', level = logging.INFO)
    logging.info('Starting magpipe')

def shutdown(status = 0):
    logging.info('Finishing magpipe with status:\t{}'.format(status))
    sys.exit(status)
