#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

import re
import argparse
from pathlib import Path
from collections import defaultdict, Counter, OrderedDict

import magpipe.log
import magpipe.pretty

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser(description='Go through a prokka output folder rename the genes id so they are consistent with the rest of the project for relevant files and create the gene_call folder.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # REQUIRED  arguments:
    parser.add_argument('-p', '--prokka_path', help='path to the prokka output.', required=True, type=Path)
    # OPTIONAL arguments:
    parser.add_argument('--original_suffix', help="prokka suffix to be removed, if appropriate.", type=str, default="-original")
    # Return namespace
    return parser.parse_args()

# From https://stackoverflow.com/a/23747652
class OrderedCounter(Counter, OrderedDict):
    """
    Create a OrderedCounter class
    A Counter that orders elements based on the appearance order
    """
    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, OrderedDict(self))

    def __reduce__(self):
        return self.__class__, (OrderedDict(self),)

def file_checkup(args, file_format):
    """
    Some safety checks on a prokka output file
    """
    files = list(args.prokka_path.glob("*{}.{}".format(args.original_suffix, file_format)))
    if len(files) > 1:
        raise ValueError("Looking for {format} than one file matching *{suffix}.{format} in the output. That shouldn't be.".format(suffix=args.original_suffix, format=file_format))
    file = files[0]
    if file.stat().st_size == 0:
        raise ValueError("Found {} file {} but it's empty...".format(file_format, file))
    return file

def extract_locus_tag(args):
    """
    Function that extracts the locus_tag from the prokka log and returns it
    """
    log_file = file_checkup(args, 'log')
    # Extracting the the locus tag from log file
    with open(log_file) as log:
        for line in log:
            if "Setting --locustag" in line.rstrip():
                locus_tag = re.sub('.*Setting --locustag | from MD5.*', '', line.rstrip())
                break
    if not bool(re.match('^[A-Z]{8}$', locus_tag)):
        raise ValueError("Locus tag expected to be 8 captial letters... and got {}.".format(locus_tag))
    return locus_tag

def dict_from_gff(gff_file):
    """
    Takes the gff file from a prokka output as input
    returns a dictionary between the prokkka gene id and scaffold id
    This will include rRNA, tRNA, tmRNA and misc_RNA but
    not repeat_regions as they don't have an ID
    """
    gff_dict = defaultdict()
    with open(gff_file) as f:
        for line in f:
            if 'ID=' in line:
                prokka_gene_id = re.sub('.*ID=|;.*', '', line.rstrip())
                scaffold_id = re.sub(r'\t.*', '', line.strip())
                gff_dict.update({prokka_gene_id : scaffold_id})
    return gff_dict

def correct_gene_names(gff_dict):
    """
    From a prokka gene id to scaffold dictionary
    return a prokka gene id to corrected gene id
    """
    scaffolds_counter = OrderedCounter(gff_dict.values())
    scaffolds_gene_numbers = []
    for n in scaffolds_counter.values():
        scaffolds_gene_numbers.extend(list(range(1, n + 1)))
    correct_gene_ids = [scaffold + '-gene_' + str(number) for scaffold, number in zip(gff_dict.values(), scaffolds_gene_numbers)]
    corrected_dict = defaultdict()
    for prokka_gene_id, corrected_gene_id in zip(gff_dict.keys(), correct_gene_ids):
        corrected_dict.update({prokka_gene_id : corrected_gene_id})
    if len(corrected_dict) != len(set(gff_dict.keys())):
        raise ValueError("Something went wrong... We're missing some genes. Corrected {} vs orginal {} # genes.".format(len(corrected_dict), len(set(gff_dict.keys()))))
    return corrected_dict

def write_file_with_correct_genes_id(args, file_format):
    """
    Given a file format, create a copy without the suffix and
    replace the prokka genes id by the correct gene ids
    """
    input_file = file_checkup(args, file_format)
    magpipe.pretty.print_pair("Fixing file", input_file)
    output_file = Path(re.sub("{suffix}.{format}$".format(suffix=args.original_suffix, format=file_format),
                              ".{format}".format(format=file_format),
                              str(input_file)))
    if file_format in ["ffn", "faa"]:
        magpipe.pretty.print_single("Moving to gene_call folder.", color='yellow')
        gene_call_path = args.prokka_path.parent.parent.joinpath("gene_call") # need to goo up twice
        gene_call_path.mkdir(exist_ok=True)
        output_file = gene_call_path.joinpath(output_file.name)

    output_file.touch()
    with open(input_file) as i, open(output_file, 'w') as o:
        regex = args.locus_tag + '_[0-9]{5}'
        for line in i:
            if bool(re.search(regex, line)):
                keys = re.findall(regex, line)
                if len(set(keys)) > 1:
                    raise ValueError("More than one locus-tag in a line of file {}... did not expect that".format(input_file))
                key = keys[0]
                if file_format in ["ffn", "faa"]:
                    newline = ">{}\n".format(args.corrected_dict[key])
                else:
                    newline = line.replace(key, args.corrected_dict[key])
                o.write(newline)
            else:
                o.write(line)

def main(args):
    """
    Main function, creates corrected copies of relevant files
    Files affected
    .faa
    .ffn
    .gff
    .tsv
    Note that faa and ffn are moved to a gene_call folder
    """
    args.locus_tag = extract_locus_tag(args)
    args.gff_file = file_checkup(args, 'gff')
    args.gff_dict = dict_from_gff(args.gff_file)
    args.corrected_dict = correct_gene_names(args.gff_dict)
    for f in ['ffn', 'faa', 'gff', 'tsv']:
        write_file_with_correct_genes_id(args, f)

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
