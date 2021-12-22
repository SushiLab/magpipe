#!/usr/bin/env python

# Imports ----------------------------------------------------------------------

import argparse
from pathlib import Path

import pandas
from tqdm import tqdm
import Bio.SeqIO.FastaIO as FastaIO

import magpipe.log
import magpipe.pretty

# Define utilitary functions ---------------------------------------------------

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser()
    # REQUIRED arguments:
    parser.add_argument('-a', '--annotations', help='Path to the annotations directory, which contains FetchMGs/{genomes}-bestMGs and -allMGs.', required=True, type=str)
    parser.add_argument('-o', '--output', help='Path and name of the output directory.', required=True, type=str)
    # OPTIONAL arguments:
    #parser.add_argument('-g', '--genomes', help='A file with the list of genomes to process.', required=False, default='/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/MAGPIPE_MAGS_EAN/scratch/processed/integrated/go_microbiomics/go_microbiomics-integrated-cpl50_ctn10-dictionaries/test-genomes-specI', type=str)
    parser.add_argument('-g', '--genomes', help='A file with the list of genomes to process.', required=False, default='/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/MAGPIPE_MAGS_EAN/scratch/processed/integrated/go_microbiomics/go_microbiomics-integrated-cpl50_ctn10-dictionaries/go_microbiomics-integrated-cpl50_ctn10-prokarya.txt', type=str)
    parser.add_argument('-t', '--taxo', help='Path to the gtdbtk taxonomy of the genomes. Just path/prefix (usually gtdbtk).', required=False, default='/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/EAN/MAGPIPE_MAGS_EAN/scratch/processed/integrated/go_microbiomics/go_microbiomics-integrated-cpl50_ctn10-gtdbtk/gtdbtk', type=str)
    # Return namespace
    return parser.parse_args()

def load_file_as_list(input_file): 
    """
    Takes a file file and lists all the entries in it
    """
    magpipe.pretty.print_pair("Loading file", input_file)
    file_list = []
    file_path = Path(input_file) 
    if not file_path.exists():
        raise ValueError('Well then, input file {} doesn\'t exist.'.format(file_path))
    with open(file_path) as handle:
        for line in handle: 
            file_list.append(line.strip())
    return file_list

def process_marker_genes(annotations, genome, fetch_type, output):
    """
    Takes all the marke genes (mgs) list for a genome, renames them and move them over
    refg:
    Adds the Marker Genes of a reference genome to thee combined marker genes file
    For that we need to read the bestMGs of that genome, rename them and append them
    to the correct input file for specI
    else:
    Adds the marker genes of a MAG or a SAG to the combined marker genes
    For that we read the marker genes for that genome, ignore the duplicated ones
    rename and write to file for specI
    """
    output = Path(output)
    mgspath = Path(annotations).joinpath(genome, "fetchMGs", genome + "-" + fetch_type.lstrip('-'))
    if not mgspath.exists():
        raise FileNotFoundError("Couldn't find {}...".format(mgspath))
    mgs = [f.name.replace(".fna", "") for f in mgspath.glob("*.fna")]
    if len(mgs) != 40:
        raise ValueError("Well, I was expecting 40 MGs and found {} when processing {}.".format(len(mgs), genome))
    genome_id = genome.replace(".", ":") + "." + genome.replace(".", ":")
    for mg in mgs:
        mg_fna = mgspath.joinpath(mg + ".fna")
        mg_faa = mgspath.joinpath(mg + ".faa")
        if mg_fna.stat().st_size == 0 and mg_faa.stat().st_size == 0:
            continue # skip empty files, i.e. the mg wasn't found
        with open(mg_fna) as fna_in, open(mg_faa) as faa_in:
            seqlist = [{'fna_h': fna_h, 'fna_s': fna_s, 'faa_h': faa_h, 'faa_s': faa_s} for (fna_h, fna_s), (faa_h, faa_s) in zip(FastaIO.SimpleFastaParser(fna_in), FastaIO.SimpleFastaParser(faa_in))]
            if len(seqlist) > 1:
                continue # skipping duplicated genes as we can't tell which is the right one and they could mess up specI
            seqdict = seqlist[0]
            if seqdict['fna_h'] != seqdict['faa_h']:
                raise ValueError("Mmmh, I've got gene names that don't match between fna and faa here... See {}, {} gave {} and {}.".format(genome, mg, seqdict['fna_h'], seqdict['faa_h']))
            header = genome_id + "." + seqdict['fna_h'].replace(".", ":")
            with open(output.joinpath(mg + ".fna"), "a") as fna_out, open(output.joinpath(mg + ".faa"), "a") as faa_out:
                fna_out.write(">" + header + "\n")
                fna_out.write(seqdict['fna_s'] + "\n")
                faa_out.write(">" + header + "\n")
                faa_out.write(seqdict['faa_s'] + "\n")
    return genome_id

def main(args):
    """
    Main function to process the marker genes for checkm 
    """
    output = Path(args.output)
    if output.exists():
        raise ValueError('Seems the script has already been ran... cleanup and start again!')
    mgs_out = output.joinpath("MGs")
    mgs_out.mkdir(parents=True)
    genomes_list = load_file_as_list(args.genomes)
    ids_list = []
    for genome in tqdm(genomes_list, ncols=100):
        if "REFG" in genome:
            genome_id = process_marker_genes(args.annotations, genome, "bestMGs", mgs_out)
        else:
            genome_id = process_marker_genes(args.annotations, genome, "allMGs", mgs_out)
        ids_list.append(genome_id)
    # Write the id list
    with open(output.joinpath("genome_ids_for_specI.txt"), "w") as handle:
        handle.write("\n".join(ids_list))
    # get the taxo for specI
    bac_table = pandas.read_csv(args.taxo + ".bac120.summary.tsv", sep="\t")
    bac_table = bac_table[["user_genome", "classification"]]
    arc_table = pandas.read_csv(args.taxo + ".ar122.summary.tsv", sep="\t")
    arc_table = arc_table[["user_genome", "classification"]]
    taxo = bac_table.append(arc_table)
    taxo.user_genome = [i.replace('.', ':') for i in taxo.user_genome]
    taxo.to_csv(output.joinpath("genome_taxo_for_specI.tsv"), index=False, sep="\t")

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
