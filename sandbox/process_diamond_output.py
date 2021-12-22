#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
This script is a diamond processing script based on the NCBI PGAP criteria
It uses diamond to blastp a set of protein query against a database database
Hits are first filtered based on id >= 25%, qcov and scov >= 70%
Then, we compute a normalised bitscore (bitscore normalised by the maximal
bitsore of the sseq)
Only hits with normalised bitscore >= 0.5 are kept
"""

import re
import math
import shutil
import pandas
import argparse
from pathlib import Path
from collections import defaultdict

import magpipe.log
import magpipe.pretty
import magpipe.commandline

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser()
    # REQUIRED arguments:
    parser.add_argument("-d", "--diamond_output", help="Path of the diamond output, needs the following format: qseqid sseqid qlen slen length pident qcovhsp evalue bitscore full_qseq/fullsseq", required=True, type=str)
    parser.add_argument("-o", "--filtered_output", help="Specify the output file name.", required=True, type=str)
    # parse and return the arguments
    return parser.parse_args()

def get_max_bitscore(seq):
    """
    This function takes a sequence as input and
    returns a maximal theoretical bitscore
    Basically, this is the sum of the diagonal of
    BLOSUM62 scores for each AA of the sequence
    Which is converted into bit-score using
    Lambda = 0.267 and  K = 0.041 for BLOSUM62
    see diamond log for values and
    https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
    for the equation
    ===== Ambiguous amino acids:
    B = D or N
    Z = E or Q
    J = I or L
    ===== Missing values and stop codons:
    Unknown = X
    STOP = *
    ===== Potential problems:
    Selenocysteine, Sec, U
    Pyrrolysine, Pyl, O
    """
    blosum62_diag = {"A": 4,
                     "R": 5,
                     "N": 6,
                     "D": 6,
                     "C": 9,
                     "Q": 5,
                     "E": 5,
                     "G": 6,
                     "H": 8,
                     "I": 4,
                     "L": 4,
                     "K": 5,
                     "M": 5,
                     "F": 6,
                     "P": 7,
                     "S": 4,
                     "T": 5,
                     "W": 11,
                     "Y": 7,
                     "V": 4,
                     # Ambiguous
                     "B": 4,
                     "Z": 4,
                     "J": 4,
                     # Unknown
                     "X": -1,
                     # Stop
                     "*": 1}
    Lambda = 0.267
    K = 0.041
    scores = [blosum62_diag[AA] for AA in seq]
    max_score = sum(scores)
    max_bitscore = (Lambda*max_score - math.log(K)) / math.log(2)
    return round(max_bitscore, 1)

def process_diamond_output(diamond_output, diamond_filtered):
    """
    Read a diamond output (.tsv) of the following format:
    qseqid sseqid qlen slen length pident qcovhsp evalue bitscore full_sseq
    Replaces the column that contains the sequences with a normalised bitscore
    which is the bitscore divided by max theoretical bitscore of the sequence
    It also applies the filter: norm_bitscore >= 0.5
    """
    magpipe.pretty.print_single("Processing the diamond output...", nl_after=1, color='blue')
    with open(diamond_output) as handle, open(diamond_filtered, "w") as target: # Use the filtered out for debugging
    #with open(diamond_temp) as handle, open(diamond_processed, "w") as target, open(diamond_processed.replace('.tsv', '-non_sig.tsv'), "w") as trash:
        header = ["qseqid", "sseqid", "qlen", "slen", "length", "pident", "qcovhsp", "evalue", "bitscore", "max_bitscore", "norm_bitscore"]
        target.write("\t".join(header) + "\n")
        for line in handle:
            line = line.rstrip().split("\t")
            seq = line[-1]
            line[-1] = get_max_bitscore(seq)
            norm_bitscore = round(float(line[-2]) / line[-1], 3)
            norm_bitscore = min(norm_bitscore, 1)
            line.append(norm_bitscore)
            line = [str(i) for i in line]
            if norm_bitscore >= 0.5:
                target.write("\t".join(line) + "\n")
            #else:
                #trash.write("\t".join(line) + "\n")

def main(args):
    """
    The main function that patches everything together
    """
    process_diamond_output(args.diamond_output, args.filtered_output)

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
