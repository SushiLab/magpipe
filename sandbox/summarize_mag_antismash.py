#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
This script can be used to nicely summarize the  biosynthetic regions found by antismash in all the
mags
"""

import gc
import pandas
import argparse
from Bio import SeqIO
from tqdm import tqdm
from pathlib import Path

# Remove biopython warning as the header of the gbk files are too long and we get a crazy amount of warnings (although no consequences)
import warnings
from Bio import BiopythonParserWarning

import magpipe.log

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser()
    # REQUIRED  arguments:
    parser.add_argument('-a', '--annotations', help='Path to the annotation directory.', required=True, type=str)
    parser.add_argument('-o', '--output', help='Path and name of the output table.', required=True, type=str)
    return parser.parse_args()

def read_gbk_biopython(gbk_file):
    """
    A biopython based parser for antismash gbk output
    Designed for Antismash 5.1
    Each gbk file represents a biosynthetic region
    that can contain several protocluster that are
    organised in candidate clusters
    """
    # scaffold, region, length, contig edge, poducts, # cand cluster, cand cluster types, # protocluster, protocluster products, # genes, CDS list
    biosynt_cluster = {}
    #print("Processing file {} ...".format(gbk_file))
    with open(gbk_file) as handle:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonParserWarning)
            record = list(SeqIO.parse(handle, "genbank"))
    if len(record) != 1:
        raise ValueError("Each file should be a single record... and for {} I've got {}...".format(gbk_file, len(record)))
    record = record[0]
    region = [r for r in record.features if r.type == "region"]
    if len(region) != 1:
        raise ValueError("Each file should be a single record... and for {} I've got {}...".format(gbk_file, len(region)))
    region = region[0]
    biosynt_cluster["scaffold"] = record.id
    region_number = region.qualifiers["region_number"]
    if len(region_number) != 1:
        raise ValueError("Expected a length 1 list as region id and got a length of {}.".format(len(region_number)))
    region_number = region_number[0]
    biosynt_cluster["region"] = "{}-biosynth_{}".format(record.id, region_number)
    biosynt_cluster["length"] = len(record.seq)
    biosynt_cluster["contig edge"] = ";".join(region.qualifiers["contig_edge"])
    biosynt_cluster["products"] = ";".join(region.qualifiers["product"])
    biosynt_cluster["# candidate clusters"] = len(region.qualifiers["candidate_cluster_numbers"])
    cand_clusters_kind = [";".join(f.qualifiers["kind"]) for f in record.features if f.type == "cand_cluster"]
    if not len(cand_clusters_kind) == biosynt_cluster["# candidate clusters"]:
        raise ValueError("Mismatch in # cand clusters, region says {} but there are {} in the gbk file {}.".format(len(cand_clusters_kind), biosynt_cluster["# candidate clusters"], gbk_file))
    biosynt_cluster["candidate clusters type"] = ";".join(cand_clusters_kind)
    protoclusters_product = [";".join(f.qualifiers["product"]) for f in record.features if f.type == "protocluster"]
    biosynt_cluster["# protoclusters"] = len(protoclusters_product)
    biosynt_cluster["protoclusters products"] = ";".join(protoclusters_product)
    try:
        genes = [";".join(f.qualifiers["ID"]) for f in record.features if f.type == "CDS" and "ID" in f.qualifiers.keys()]
        biosynthetic = [";".join(f.qualifiers['ID']) for f in record.features if f.type == "CDS" and "ID" in f.qualifiers.keys() and "gene_kind" in f.qualifiers.keys() and "biosynthetic" in f.qualifiers["gene_kind"]]
        biosynthetic_additional = [";".join(f.qualifiers['ID']) for f in record.features if f.type == "CDS" and "ID" in f.qualifiers.keys() and "gene_kind" in f.qualifiers.keys() and "biosynthetic-additional" in f.qualifiers["gene_kind"]]
        regulatory = [";".join(f.qualifiers['ID']) for f in record.features if f.type == "CDS" and "ID" in f.qualifiers.keys() and "gene_kind" in f.qualifiers.keys() and "regulatory" in f.qualifiers["gene_kind"]]
        resistance = [";".join(f.qualifiers['ID']) for f in record.features if f.type == "CDS" and "ID" in f.qualifiers.keys() and "gene_kind" in f.qualifiers.keys() and "resistance" in f.qualifiers["gene_kind"]]
        transport = [";".join(f.qualifiers['ID']) for f in record.features if f.type == "CDS" and "ID" in f.qualifiers.keys() and "gene_kind" in f.qualifiers.keys() and "transport" in f.qualifiers["gene_kind"]]
    except:
        raise Exception("There was an error processing a CDS sequence from {}.".format(gbk_file))
    biosynt_cluster["# CDS"] = len(genes)
    biosynt_cluster["CDS list"] = ";".join(genes)
    biosynt_cluster["CDS biosynthetic"] = ";".join(biosynthetic)
    biosynt_cluster["CDS biosynthetic additional"] = ";".join(biosynthetic_additional)
    biosynt_cluster["CDS regulatory"] = ";".join(regulatory)
    biosynt_cluster["CDS resistance"] = ";".join(resistance)
    biosynt_cluster["CDS transport"] = ";".join(transport)
    return biosynt_cluster

def main(annotation_dir, output_table):
    """
    Main function of the pipeline
    Reads the list of genomes and process the antismash output to make a nice table summary
    """
    annotation_path = Path(annotation_dir)
    if not annotation_path.exists():
        raise FileNotFoundError("Looks like there is something wrong with {}.".format(annotation_path))
    output_path = Path(output_table)
    if output_path.exists():
        raise Exception("Output {} already exists...".format(output_path))
    table_columns = ["genome", "scaffold", "region", "length", "contig edge", "products", "# candidate clusters", "candidate clusters type", "# protoclusters", "protoclusters products", "# CDS", "CDS list", "CDS biosynthetic", "CDS biosynthetic additional", "CDS regulatory", "CDS resistance", "CDS transport"]
    output_path.write_text("\t".join(table_columns) + "\n")
    res_list = []
    genome_count = 0
    for genome_path in tqdm(annotation_path.glob('*'), ncols=100):
        genome_count += 1
        genome = genome_path.name
        #print("Processing genome {}...".format(genome))
        antismash_path = genome_path.joinpath('antismash', genome + '-antismash')
        # find the gbk files
        gbk_files = [f.name for f in antismash_path.glob('*region*.gbk')]
        for gbk in gbk_files:
            gbk_as_dict = read_gbk_biopython(antismash_path.joinpath(gbk))
            gbk_as_dict["genome"] = genome
            res_list.append(gbk_as_dict)
        if len(res_list) > 5000:
            print("Finished processing {} genomes. Now writing to file.".format(genome_count))
            res_df = pandas.DataFrame(res_list)
            res_df = res_df[table_columns]
            res_df.to_csv(output_table, sep="\t", header=False, index=False, mode='a')
            res_df = pandas.DataFrame()
            res_list = []
            gc.collect()

    print("Finished processing {} genomes. Now writing to file.".format(genome_count))
    res_df = pandas.DataFrame(res_list)
    res_df = res_df[table_columns]
    res_df.to_csv(output_table, sep="\t", header=False, index=False, mode='a')

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args.annotations, args.output)
