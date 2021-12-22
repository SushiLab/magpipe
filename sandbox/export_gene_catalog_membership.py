#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Combine the profiles in single tables, normalize by gene length if appropriate
"""

# Import libraries -------------------------------------------------------------

import argparse
from pathlib import Path

import magpipe.log
import magpipe.pretty

# Define utilitary functions ---------------------------------------------------

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--dedupe_clstr', help='Path to the dedupe clustering file (optional).', default="skip", type=str)
    parser.add_argument('-c', '--cdhit_clstr', help='Path to the cdhit clustering file.', default='/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/MAGPIPE_MAGS_EAN/scratch/processed/gene_catalog/go_microbiomics/go_microbiomics-integrated-cpl50_ctn10-gene_catalog/3cdhit/tmp/go_microbiomics-integrated-cpl50_ctn10-genes-cds-w_eukarya.lengthsorted.dedupe100.cdhit9590.fasta.clstr', type=str)
    parser.add_argument('-o', '--output', help='Path the output dictionary.', required=True, type=str)
    return parser.parse_args()

def read_clstr_as_dict(clstr):
    """
    Takes a clstr file as input and returns dict of repr, list of members
    >Cluster 1
    0       5796nt, >TARA_SAMEA4397330_METAG-scaffold_20852-gene_1... at 245:5750:29172:34679/+/99.60%
    1       3066nt, >TARA_SAMEA4397174_METAG-scaffold_8551-gene_6... at 1:3062:14842:17903/+/100.00%
    2       4458nt, >TARA_SAMEA4397174_METAG-scaffold_16015-gene_1... at 327:4458:4689:8821/+/99.71%
    3       51870nt, >TARA_SAMEA4397118_METAG-scaffold_132-gene_52... *
    4       2559nt, >TARA_SAMEA4397330_METAG-scaffold_56906-gene_1... at 1:2559:12703:15262/+/99.57%
    5       2331nt, >TARA_SAMEA4397174_METAG-scaffold_16015-gene_2... at 1:2331:8830:11162/+/99.49%
    6       2094nt, >TARA_SAMEA4397330_METAG-scaffold_60567-gene_1... at 1:2094:16783:18895/+/97.11%
    7       1965nt, >TARA_SAMEA4397174_METAG-scaffold_8551-gene_5... at 1:1965:17923:19898/+/98.94%
    8       1911nt, >TARA_SAMEA4397330_METAG-scaffold_55717-gene_3... at 1:1747:47650:49396/+/99.94%
    9       1710nt, >TARA_SAMEA4397174_METAG-scaffold_20078-gene_1... at 1:1710:50161:51870/+/99.94%
    10      630nt, >TARA_SAMEA4397330_METAG-scaffold_55717-gene_1... at 1:627:49765:50394/+/99.37%
    11      402nt, >TARA_SAMEA4397174_METAG-scaffold_58077-gene_2... at 1:402:29065:29470/+/96.55%
    12      327nt, >TARA_SAMEA4397330_METAG-scaffold_55717-gene_2... at 1:327:49474:49793/+/96.94%

    with the star indicating the repr.
    """
    magpipe.pretty.print_pair("Reading clstr file", clstr)
    res = {}
    rep = None
    members = []
    with open(clstr) as handle:
        for line in handle:
            if line.startswith(">"):
                if rep is not None:  
                    for m in members:
                        res[m] = rep
                    rep = None
                    members = []
            else:
                seq = line.split(">")[1].split("...")[0]
                if "*" in line:
                    rep = seq
                members.append(seq)
        # catching the last cluster
        for m in members:
            res[m] = rep
    magpipe.pretty.print_pair("Number of sequences with membership", len(res), nl_after=1)
    return res

# Main function of the script --------------------------------------------------

def main(args):
    """
    Reads the 100% derep and cdhit clstr files and merge them to get the full mapping
    """
    output_path = Path(args.output)
    if output_path.exists():
        raise ValueError("Looks like some of your output already exists... please check {}.".format(output_path))
    if args.dedupe_clstr == "skip":
        cdhit_path = Path(args.cdhit_clstr)
        if not cdhit_path.exists(): 
            raise ValueError("Can't really find the expected files...")
        cdhit_dict = read_clstr_as_dict(cdhit_path) 
        magpipe.pretty.print_pair("Writing output to", output_path)
        with open(output_path, 'w') as target:
            target.write("gene\trepresentative\n")
            for k in cdhit_dict.keys():
                target.write(k + "\t" + cdhit_dict[k] + "\n")
    else:
        derep_path = Path(args.dedupe_clstr) 
        cdhit_path = Path(args.cdhit_clstr)
        if not derep_path.exists() or not cdhit_path.exists(): 
            raise ValueError("Can't really find the expected files...")
        derep_dict = read_clstr_as_dict(derep_path) 
        cdhit_dict = read_clstr_as_dict(cdhit_path) 
        magpipe.pretty.print_pair("Writing output to", output_path)
        with open(output_path, 'w') as target:
            target.write("gene\trepresentative\n")
            for k in derep_dict.keys():
                target.write(k + "\t" + cdhit_dict[derep_dict[k]] + "\n")

# Run script -------------------------------------------------------------------

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
