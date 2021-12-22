#!/usr/bin/env python

# Import libraries -------------------------------------------------------------

import os
import sys
import gzip
import pathlib
import argparse
import Bio.SeqIO.FastaIO as FastaIO
from collections import defaultdict

#cluster_file = '/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/scratch/paolil/Ocean_MAGs/results/EMOSE_MAGs/ETHSEQ0000004025/vamb_bins_a26_m2k/vamb_out/clusters.tsv'
#outdir = '/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/scratch/paolil/Ocean_MAGs/results/EMOSE_MAGs/ETHSEQ0000004025/vamb_bins_a26_m2k/'
#assembly = '/nfs/cds-shini.ethz.ch/exports/biol_micro_cds_gr_sunagawa/SequenceStorage/analysis/ETHSEQ0000004025/assembly/spades/ETHSEQ0000004025_scaffolds.min1000.fasta.gz'
#prefix = 'ETHSEQ0000004025_a26_m2k'

# Define utilitary functions ---------------------------------------------------

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser(
    prog = 'write_vamb_bins.py',
    description = 'Use vamb cluster output and write bins as fasta files',
    formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    # Not Shown:
    parser.add_argument('--version',
    action = 'version',
    version = '%(prog)s 0.1',
    help = argparse.SUPPRESS)

    parser.add_argument('--debug',
    action = 'store_true',
    help = argparse.SUPPRESS)

    # REQUIRED  arguments:
    parser.add_argument('-c', '--cluster_file',
    help = 'Path to the cluster.tsv file from VAMB output.',
    required = True,
    type = pathlib.Path)
    parser.add_argument('-o', '--outdir',
    help = 'Output directory to write the bins in fasta format.',
    required = True,
    type = pathlib.Path)
    parser.add_argument('-a', '--assembly',
    help = 'Assembly file.',
    required = True,
    type = pathlib.Path)
    parser.add_argument('-p', '--prefix',
    help = 'Prefix for the bins names.',
    required = True,
    type = str)

    # OPTIONAL arguments:
    parser.add_argument('-m', '--min_size',
    help = "Minimum size for a bin to be written (in base pair).",
    type = int,
    default = 200000)

    return parser.parse_args()


def print_arguments(args):
    """
    Print arguments
    :param args
    """
    print("\nChecking arguments:")
    print("Cluster file (VAMB output) : {}".format(args.cluster_file))
    print("Output directory : {}".format(args.outdir))
    print("Prefix of the bins : {}".format(args.prefix))
    print("Assembly file : {}".format(args.assembly))
    print("Minimum size of the written bins : {}".format(args.min_size))


def remove_unzipped_assembly(outdir, prefix):
    unzassembly = outdir.parent.joinpath(prefix.split('_vamb')[0] + '.assembly.fa')
    print("\nLooking for unzipped assembly in output's parent directory...")
    if unzassembly.is_file():
        print("removing '{}'...".format(unzassembly))
        os.remove(unzassembly)
    else:
        print("Nothing found matching '{}'...".format(unzassembly))


def check_arguments(args):
    """
    Check arguments, i.e. wether directories and file exist
    :param args
    """
    if not args.cluster_file.is_file():
        sys.exit("\nERROR: '{}' is not a file.".format(args.cluster_file))
    if not args.assembly.is_file():
        sys.exit("\nERROR: '{}' is not a file.".format(args.assembly))
    if args.outdir.parent != args.cluster_file.parents[1]:
        sys.exit("\nERROR: outdir ('{}') is expected to be a subfolder of 'sample/workflow/'.".format(args.outdir))
    if not args.outdir.is_dir():
        print("\n'{}' doesn't exist... creating it.".format(args.outdir))
        args.outdir.mkdir()
    if any(files.endswith('.fa') for files in os.listdir(args.outdir)):
        sys.exit("\nERROR: File(s) ending in '.fa' in outdir... check if bins already exist.")
    print('\n... looks good, let\'s go on.\n')


# Functions to write the VAMB bins ---------------------------------------------

def read_vamb_clusters(cluster_file):
    """
    Read the vamb clusters
    """
    clusters = defaultdict(lambda: defaultdict(list))
    with open(cluster_file, "r") as handle:
        for line in handle:
            cluster_name = line.split()[0]
            contig_name = line.split()[1]
            contig_length = int(line.split()[3].strip('length='))
            clusters[cluster_name]['contig_names'].append(contig_name)
            clusters[cluster_name]['contig_lengths'].append(contig_length)
    return(clusters)


def filter_clusters(clusters, threshold = 200000):
    contig_to_cluster = defaultdict()
    for cluster in clusters.keys():
        if sum(clusters[cluster]['contig_lengths']) >= threshold:
            for contig in clusters[cluster]['contig_names']:
                contig_to_cluster[contig] = cluster
    return(contig_to_cluster)


def write_vamb_bins(contig_to_cluster, assembly, outdir, prefix):
    with gzip.open(assembly, "rt") as handle:
        for (header, sequence) in FastaIO.SimpleFastaParser(handle):
            if contig_to_cluster.get(header.split()[0]) is not None:
                contig = header.split()[0]
                cluster = contig_to_cluster[contig]
                outfile = outdir.joinpath(prefix + '.' + cluster.strip('cluster_') + '.fa')
                with open(outfile, 'a') as out:
                    out.write(''.join(['>', contig, '\n', sequence, '\n']))


# Run script -------------------------------------------------------------------

if __name__ == '__main__':
    args = get_arguments()
    print_arguments(args)
    remove_unzipped_assembly(args.outdir, args.prefix)
    check_arguments(args)
    clusters = read_vamb_clusters(args.cluster_file)
    contig_to_cluster = filter_clusters(clusters)
    write_vamb_bins(contig_to_cluster, args.assembly, args.outdir, args.prefix)
