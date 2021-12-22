#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
This script is a kegg annotator based on the NCBI PGAP criteria
It uses diamond to blastp a set of protein query against the KEGG database
Hits are first filtered based on id >= 25%, qcov and scov >= 70%
Then, we compute a normalised bitscore (bitscore normalised by the maximal
bitsore of the sseq)
Only hits with normalised bitscore >= 0.5 are kept
Then, we attribute the set of KOs within the 0.1 range of the best normalised
bitscore to the query
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
    parser.add_argument("-q", "--query", help="Query file in a .faa format", required=True, type=str)
    parser.add_argument("-o", "--output_prefix", help="Path and prefix of the output files", required=True, type=str)
    # OPTIONAL arguments:
    parser.add_argument("-d", "--database", help="Path to the (processed) kegg database", required=False, default='/science/paolil/databases/kegg/kegg_for_prokka-with_ko.dmnd', type=str)
    parser.add_argument("-k", "--ko_dictionaries", help="Path to the ko directory containing the mapping files.", required=False, default='/nfs/cds/Databases/KEGG/February_2020/ftp.bioinformatics.jp/kegg/genes/ko/', type=str)
    parser.add_argument("-t", "--threads", help="Number of threads to use.", required=False, default=8, type=int)
    parser.add_argument("-m", "--mode", help="Specify if diamond should run in default, 'sensitive' or 'more-sensitive' mode", required=False, default=None, type=str)
    # parse and return the arguments
    return parser.parse_args()

def check_diamond_db(args):
    """
    This function just checks that the diamond db is in the right format
    and creates a new one if necessary
    returns the name of the database
    """
    if shutil.which('diamond') is None:
        raise Exception("Couldn't find diamond in your path...")
    magpipe.pretty.print_pair("Diamond executable", shutil.which('diamond'))
    if args.database.endswith(".dmnd"):
        magpipe.pretty.print_single("The database {} ends with '.dmnd' so I assume it's all good to go.".format(args.database), nl_after=1, color='blue')
        db_name = args.database
    else:
        db_prefix = "".join(args.database.split(".")[0:-1])
        db_name = db_prefix + '.dmnd'
        command = "diamond makedb --in {} --db {}".format(args.database, db_prefix)
        returncode = magpipe.commandline.run_command(command)
        if returncode != 0:
            raise Exception("Diamond commandline {} finished with returncode {}...".format(command, returncode))
        if not Path(db_name).exists():
            raise FileNotFoundError("The diamond database file {} doesn't exist... something went wrong.".format(db_name))
    return db_name

def run_diamond(args):
    """
    This functions sets up the diamond commandline and runs it
    """
    # Note that we can set the --top 20, which means it will only report alignments with a bitscore at most 20% lower than the top alignment
    # This speeds things up without any cost as we only look at alignements with a normalised bitscore of normm > = 0.5 and norm >= (max(norm) - 0.1)
    fixed_params = "--id 25 --query-cover 70 --subject-cover 70 --top 20 --outfmt 6 qseqid sseqid qlen slen length pident qcovhsp evalue bitscore full_sseq"
    if args.mode == 'sensitive':
        fixed_params = "--sensitive " + fixed_params
    elif args.mode == 'more-sensitive':
        fixed_params = "--more-sensitive " + fixed_params
    temp_output = args.output_prefix + "-diamond-temp.tsv"
    command = "diamond blastp -q {} -d {} -o {} --threads {} {} &> {}.log".format(args.query, args.database, temp_output, args.threads, fixed_params, args.output_prefix)
    returncode = magpipe.commandline.run_command(command)
    if returncode != 0:
        raise Exception("Diamond commandline {} finished with returncode {}...".format(command, returncode))
    if not Path(temp_output).exists():
        raise FileNotFoundError("The diamond output file {} doesn't exist... something went wrong.".format(temp_output))
    return temp_output

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

def process_diamond_output(diamond_temp, diamond_processed):
    """
    Read a diamond output (.tsv) of the following format:
    qseqid sseqid qlen slen length pident qcovhsp evalue bitscore full_sseq
    Replaces the column that contains the sequences with a normalised bitscore
    which is the bitscore divided by max theoretical bitscore of the sequence
    It also applies the filter: norm_bitscore >= 0.5
    """
    magpipe.pretty.print_single("Processing the diamond output...", nl_after=1, color='blue')
    with open(diamond_temp) as handle, open(diamond_processed, "w") as target: # Use the filtered out for debugging
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

def read_ko_mapping(path, mapping_file):
    """
    This functions takes a file from the ko mappings
    and return a dictionary
    currently work for genes, cog and ec
    """
    path = Path(path)
    dict = defaultdict(lambda: None)
    with open(path.joinpath(mapping_file)) as handle:
        for line in handle:
            line = line.strip().split('\t')
            ko = re.sub("^ko:", "", line[0])
            value = re.sub("^cog:|^ec:", "", line[1])
            if mapping_file == "ko_genes.list": # reverse the dictionary in this specific case
                value, ko = (ko, value)
            if ko in dict.keys():
                if isinstance(dict[ko], list):
                    dict[ko].append(value)
                else:
                    raise TypeError("{} is not a list...".format(dict[ko]))
            else:
                dict.update({ko: [value]})
    return dict

def read_ko_descriptions(path):
    """
    This functions reads the ko description file and
    returns a mapping between the KO and their description
    """
    path = Path(path)
    dict = defaultdict(lambda: None)
    with open(path.joinpath("ko")) as handle:
        for line in handle:
            # The firt filterings is because we want to keep lines that start
            # by spaces and then  have the KO number with the description
            if line[0] != " ":
                continue
            line = line.strip()
            if line[0] != "K" or not line[1].isdigit():
                continue
            line = line.split("  ")
            if len(line) != 2:
                continue
            ko = line[0]
            descr = line[1]
            if ko in dict.keys():
                if isinstance(dict[ko], list):
                    if not descr in dict[ko]:
                        dict[ko].append(descr)
                else:
                    raise TypeError("{} is giving us some trouble...\nThe line looks like that:\n{}\nAnd the dict value is {}: {}.".format(ko, line, dict[ko], type(dict[ko])))
            else:
                dict.update({ko: [descr]})
    return dict

def match_ko(ko_mapping, ko_list, sep = ","):
    """
    This function takes a dictionary mapping KOs to external IDs,
    e.g. KO to COG
    along with a list of KOs (string), which can be of the form "KO1,KO2"
    and returns the transformed list, with e.g. "COG1,COG2"
    """
    res = []
    for kos in ko_list:
        kos = kos.split(",")
        others = []
        for ko in kos:
            if ko_mapping[ko] is not None:
                for other in ko_mapping[ko]:
                    others.append(other)
        res.append(sep.join(others))
    return res

def get_kegg_annotations(args, processed_diamond_output):
    """
    This functions translate the processed diamond results into an annotation table
    """
    magpipe.pretty.print_single("Retrieving KO information...", nl_after=1, color='blue')
    # let's first read mappings from KEGG
    ko_dict = read_ko_mapping(args.ko_dictionaries, "ko_genes.list")
    ko_cog = read_ko_mapping(args.ko_dictionaries, "ko_cog.list")
    ko_ec = read_ko_mapping(args.ko_dictionaries, "ko_enzyme.list")
    ko_descr = read_ko_descriptions(args.ko_dictionaries)
    magpipe.pretty.print_single("Formatting annotations...", nl_after=1, color='blue')
    diamond_df = pandas.read_csv(processed_diamond_output, sep="\t")
    diamond_df = diamond_df.set_index("qseqid", drop=False)
    diamond_df["KO"] = [",".join(ko_dict[i]) for i in diamond_df.sseqid]
    diamond_df["COG"] = match_ko(ko_cog, diamond_df.KO)
    diamond_df["EC"] = match_ko(ko_ec, diamond_df.KO)
    diamond_df["Description"] = match_ko(ko_descr, diamond_df.KO, sep=" / ")
    annotation_dict = defaultdict(lambda: None)
    for gene in set(diamond_df.index):
        gene_df = diamond_df.loc[gene]
        if isinstance(gene_df, pandas.DataFrame):
            gene_df = gene_df[gene_df.norm_bitscore >= (max(gene_df.norm_bitscore) - 0.1)]
            KOs = ','.join([i for i in set(gene_df.KO) if i != ""]) # We want to remove empty strings (from None) that make things ugly
            COGs = ','.join([i for i in set(gene_df.COG) if i != ""])
            ECs = ','.join([i for i in set(gene_df.EC) if i != ""])
            Descriptions = ' / '.join([i for i in set(gene_df.Description) if i != ""])
            annotation_dict.update({gene: [KOs, COGs, ECs, Descriptions]})
        elif isinstance(gene_df, pandas.Series):
            KO = gene_df.KO
            COG = gene_df.COG
            EC = gene_df.EC
            Description = gene_df.Description
            annotation_dict.update({gene: [KO, COG, EC, Description]})
        else:
            raise TypeError("Unexpected type {}...".format(type(gene_df)))
    annotation_df = pandas.DataFrame.from_dict(annotation_dict, orient="index", columns=["KO", "COG", "EC", "Description"])
    annotation_df.index.name = "gene"
    annotation_df.to_csv(args.output_prefix + "-annotations.tsv", sep="\t")

def main(args):
    """
    The main function that patches everything together
    """
    args.database = check_diamond_db(args)
    diamond_temp = run_diamond(args)
    process_diamond_output(diamond_temp, args.output_prefix + '-diamond.tsv')
    if Path(diamond_temp).is_file():
        Path(diamond_temp).unlink()
    get_kegg_annotations(args, args.output_prefix + '-diamond.tsv')

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
