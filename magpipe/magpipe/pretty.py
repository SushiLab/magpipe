#!/usr/bin/bash python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Functions to print prettier things
#FIXME do not write colors when writing to file?
"""

import sys
from tabulate import tabulate

def prettify_text(text, color = 'white', type = 'normal'):
    """
    Takes a text string a return it prettier (colored, bold)
    We want to return something like '\033[1;37mTEXT\033[0m'
    Inspired by the Anvio ttycolors
    Only formats the text if tty is detected
    """
    if sys.stdout.isatty():
        escape = '\033'
        reset = '\033[0m'
        colors_dict = {'gray'     : '30',
                       'red'      : '31',
                       'green'    : '32',
                       'yellow'   : '33',
                       'blue'     : '34',
                       'magenta'  : '35',
                       'cyan'     : '36',
                       'white'    : '37',
                       'crimson'  : '38'}
        type_dict = {'normal'     : '0',
                     'bold'       : '1',
                     'faint'      : '2',
                     'italic'     : '3',
                     'underlined' : '4'}
        prettify = escape + '[{}m'.format(';'.join([type_dict[type], colors_dict[color]]))
        pretty_text = "{prefix}{text}{reset}".format(prefix = prettify, text = text, reset = reset)
        return(pretty_text)
    else:
        return(text)


def print_pair(key, value, nl_before = 0, nl_after = 0, k_color = 'green', k_type = 'normal', v_color = 'blue', v_type = 'normal'):
    """
    Nice print of a value pair in Anvio style
    """
    n_chars = 50
    if len(key) > n_chars - 5:
        raise NotImplementedError('Currently did not plan the key to be that long...')
    pretty_key = prettify_text(key, color = k_color, type = k_type)
    pretty_value = prettify_text(value, color = v_color, type = v_type)
    middle = ' {}: '.format('.' * (n_chars - len(key) - 3))
    pretty_pair = "{before}{key}{middle}{value}{after}".format(before = ('\n' * nl_before), key = pretty_key, middle = middle, value = pretty_value, after = ('\n' * nl_after))
    print(pretty_pair)


def print_single(text, nl_before = 0, nl_after = 0, color = 'white', type = 'normal'):
    """
    Print a nice single line
    """
    pretty_text = ('\n' * nl_before) + prettify_text(text, color = color, type = type) + ('\n' * nl_after)
    print(pretty_text)


def print_header(text, nl_before = 1, nl_after = 1, color = 'blue', type = 'bold', line_length = None):
    """
    Print a nice header line with double dashed line below
    """
    if not line_length:
        line_length = len(text)
    header_underlined = text + '\n' + ('=' * line_length)
    pretty_header = prettify_text(header_underlined, color = color, type = type)
    print(('\n' * nl_before) + pretty_header + ('\n' * nl_after))


def print_table(df, nl_before = 1, nl_after = 1, color = 'blue', type = 'normal'):
    """
    Print a pandas df nicely
    """
    pretty_table = prettify_text(tabulate(df, headers = df.columns, showindex = False), color = color, type = type)
    print(('\n' * nl_before) + pretty_table + ('\n' * nl_after))


def melt_table(df, k_color = 'green', k_type = 'normal', v_color = 'blue', v_type = 'normal'):
    """
    Melt a pandas df into pairs to print large strings nicely
    """
    for i in df.index:
        sub = df.loc[i]
        print_single('Dataset #{}'.format(list(df.index).index(i)), color = 'green', type = 'underlined')
        for key in df.columns:
            print_pair(key, str(sub[key]), k_color = k_color, k_type = k_type, v_color = v_color, v_type = v_type)
