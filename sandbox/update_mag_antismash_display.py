#!/usr/bin/env python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
The role of this script is to come after the first round of functional
annotations of the genomes, i.e. (1) prokka, (2) KEGG and (3) eggNOG and (4) antiSMASH
The goal is to integrate annotations from prokka, KEGG and eggNOG in a standardized format to the
antiSMASH display.
"""

import re
import shutil
import pandas
import argparse
from pathlib import Path
from collections import defaultdict

import magpipe.pretty
import magpipe.log

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser()
    # REQUIRED  arguments:
    parser.add_argument('-a', '--annotations', help='Path to the annotation directory.', required=True, type=str)
    return parser.parse_args()

def format_annotations_for_html(tool, dict, descr,
                                tab="&nbsp;&nbsp;",
                                base_url_ko="https://www.genome.jp/dbget-bin/www_bget?ko:",
                                base_url_cog="https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=",
                                base_url_ec="https://enzyme.expasy.org/EC/"): #FIXME add a EC website?
    """
    The goal of this function is to format the processed annotations for the different sources
    Because the different sources have different attributes, read the keys and go from there
    """
    annotation_list = ["<b>{}:</b>".format(tool)]
    annotation_list.append("{tab}<u>Descr:</u> {descr}".format(tab=tab, descr=descr))
    for k in dict.keys():
        if k not in ["KO", "COG", "EC"]:
            raise ValueError("Got an expected keys in the annotation dictinoary: {}. Should be in {}.".format(k, ["KO", "COG", "EC"]))
        if k == "KO" and dict[k] != "nan":
            ko_keys = re.findall("K[0-9]{5}", dict[k])
            ko_urls = [r"<a href=\"{url}{KO}\" target=\"_blank\">{KO}</a>".format(url=base_url_ko, KO=i) for i in ko_keys]
            dict[k] = ",".join(ko_urls)
        if k == "COG" and dict[k] != "nan":
            cog_keys = re.findall("COG[0-9]{4}", dict[k])
            cog_urls = [r"<a href=\"{url}{COG}\" target=\"_blank\">{COG}</a>".format(url=base_url_cog, COG=i) for i in cog_keys]
            dict[k] = ",".join(cog_urls)
        if k == "EC" and dict[k] != "nan":
            ec_keys = dict[k].split(',')
            ec_urls = [r"<a href=\"{url}{EC}\" target=\"_blank\">{EC}</a>".format(url=base_url_ec, EC=i) for i in ec_keys]
            dict[k] = ",".join(ec_urls)
        annotation_list.append("{tab}<u>{key}:</u> {value}".format(tab=tab, key=k, value=dict[k]))
    return "<br>".join(annotation_list)

def load_kegg_annotations(kegg_annotations):
    """
    The goal of this function is to process the kegg annotation table and
    return a standard annotation in a dictionary with locus tags as keys
    format: "tool: KO:KOXXX<br> COG:COGXXX<br> EC:X.X.XX<br> description:text"
    """
    magpipe.pretty.print_pair("Loading annotations", "KEGG")
    kegg_table = pandas.read_csv(kegg_annotations, sep = "\t")
    if not all([c in kegg_table.columns for c in ["gene", "KO", "COG", "EC", "Description"]]):
        raise ValueError("Couldn't find the columns I'm looking for... did the names change?")
    if len(kegg_table.columns) > 5:
        raise ValueError("This table looks too big, I only expect 5 columns (gene, ko, cog, ec and descr) and got {}".format(len(table.columns)))
    kegg_table = kegg_table.set_index("gene")
    kegg_annotations_dict = defaultdict(lambda: "<b>KEGG:</b> nan")
    for i in kegg_table.index:
        kegg_annotation = format_annotations_for_html(tool="KEGG",
                                                      dict={"KO": str(kegg_table.loc[i, "KO"])},
                                                      #KO=str(kegg_table.loc[i, "KO"]),
                                                      #COG=str(kegg_table.loc[i, "COG"]),
                                                      #EC=str(kegg_table.loc[i, "EC"]),
                                                      descr=str(kegg_table.loc[i, "Description"]))
        kegg_annotations_dict.update({i: kegg_annotation})
    return(kegg_annotations_dict)

def load_eggnog_annotations(eggnog_annotations):
    """
    The goal of this function is to process the eggnog annotation table and
    return a standard annotation in a dictionary with locus tags as keys
    format: "tool: KO:KOXXX, COG:COGXXX, EC:X.X.XX, description:text"
    """
    magpipe.pretty.print_pair("Loading annotations", "EggNOG")
    eggnog_table = pandas.read_csv(eggnog_annotations, sep = "\t", header = 3, skipfooter = 3, engine = "python")
    if not all([c in eggnog_table.columns for c in ["#query_name", "KEGG_ko", "eggNOG OGs", "Preferred_name", "eggNOG free text desc."]]):
        raise ValueError("Couldn't find the columns I'm looking for... did the names change?")
    eggnog_table = eggnog_table.set_index("#query_name")
    eggnog_annotations_dict = defaultdict(lambda: "<b>eggNOG:</b> nan")
    for i in eggnog_table.index:
        KO = str(eggnog_table.loc[i, "KEGG_ko"]).replace('ko:', '') # str() cause there are some nan
        COG = [cog.split("@")[0] for cog in eggnog_table.loc[i, "eggNOG OGs"].split(",") if "COG" in cog]
        EC = str(eggnog_table.loc[i, "EC"])
        name = str(eggnog_table.loc[i, "Preferred_name"])
        text = str(eggnog_table.loc[i, "eggNOG free text desc."])
        if len(COG) == 0:
            COG = "nan"
        else:
            COG = ",".join(set(COG)) # Note that we lose the order here
        if any([KO != "nan", COG != "nan", EC != "nan", name != "nan", text != "nan"]):
            if name != "nan":
                Description = "{name}; {text}".format(name=name, text=text)
            else:
                Description = text
            eggnog_annotation = format_annotations_for_html(tool="EggNOG",
                                                            dict={"COG": COG},
                                                            #KO=KO,
                                                            #COG=COG,
                                                            #EC=EC,
                                                            descr=Description)
            eggnog_annotations_dict.update({i: eggnog_annotation})
    return(eggnog_annotations_dict)

def load_prokka_annotations(prokka_annotations):
    """
    The goal of this function is to process the prokka annotation table and
    return a standard annotation in a dictionary with locus tags as keys
    format: "tool: KO:KOXXX, COG:COGXXX, EC:X.X.XX, description:text"
    """
    magpipe.pretty.print_pair("Loading annotations", "Prokka")
    prokka_table = pandas.read_csv(prokka_annotations, sep = "\t")
    if not all([c in prokka_table.columns for c in ["locus_tag", "gene", "EC_number", "COG", "product"]]):
        raise ValueError("Couldn't find the columns I'm looking for... did the names change?")
    if len(prokka_table.columns) > 7:
        raise ValueError("This table looks too big, I only expect 7 columns and got {}".format(len(table.columns)))
    prokka_table = prokka_table.dropna(subset=['locus_tag']) # remove features without oocus tags, e.g. CRISPR regions
    prokka_table = prokka_table.set_index("locus_tag")
    prokka_annotations_dict = defaultdict(lambda: "<b>Prokka:</b> nan")
    for i in prokka_table.index:
        KO = "nan" # FIXME Add the KO mapping?
        COG = str(prokka_table.loc[i, "COG"]) # str() cause there are some nan
        EC = str(prokka_table.loc[i, "EC_number"]) # str() cause there are some nan
        gene = str(prokka_table.loc[i, "gene"])
        text = str(prokka_table.loc[i, "product"])
        Description = "{name}; {text}".format(name=gene, text=text)
        if any([KO != "nan", COG != "nan", EC != "nan", gene != "nan", text != "hypothetical protein"]):
            prokka_annotation = format_annotations_for_html(tool="prokka",
                                                            dict={"COG": COG,
                                                                  "EC": EC},
                                                            #KO=KO,
                                                            #COG=COG,
                                                            #EC=EC,
                                                            descr=Description)
            prokka_annotations_dict.update({i: prokka_annotation})
    return(prokka_annotations_dict)

def update_antismash_index(html_index, genome, url="https://sunagawalab.ethz.ch/share/MARINE_METABOLITES/db", version="1.0"):
    """
    A quick function to replace the page description by a link back to the genome's home folder
    """
    pattern="<div class=\"custom-description\">{}</div>".format(genome)
    replacement="<div class=\"custom-description\"><a href=\"{url}/{version}/{genome}/\" target=\"_blank\">{genome}</a></div>".format(url=url, version=version, genome=genome)
    index_path = Path(html_index)
    magpipe.pretty.print_pair("Updating antismash index", index_path)
    if not index_path.exists():
        raise FileNotFoundError("Well, seems that {} doesn't exsit.".format(index_path))
    index_content = index_path.read_text().splitlines()
    with open(index_path, 'w') as handle:
        for line in index_content:
            if pattern in line:
                handle.write(re.sub(pattern, replacement, line) + "\n")
            else:
                handle.write(line + "\n")

def update_antismash_display(js_file, prokka_dict, kegg_dict, eggnog_dict): # (annotation_dictionary, js_file)
    """
    Because everything else did not work...
    Let's just edit the regions.js file and add the annotation that we want:
    """
    js_file = str(js_file)
    magpipe.pretty.print_pair("Updating antismash display", js_file)
    formatted_file = Path(js_file)
    backup_file = Path(re.sub(".js$", "-original.js", js_file))
    if backup_file.exists():
        raise ValueError("Backup file {} already exists".format(backup_file))
    shutil.move(formatted_file, backup_file)
    with open(backup_file) as input, open(formatted_file, "w") as target:
        flag = False
        for line in input:
            line = line.strip("\n")
            if line.startswith("var all_regions = {"): # We're just changing that variable
                flag = True
            if line.startswith("var details_data = {"):
                flag = False
            if "Locus tag:" in line:
                # get line in a nice format, need to add "                " at the beginning again though
                line = [i.replace("\\n", "") for i in line.split("<br>")]
                indent = line[0].split("\"")[0]
                line = [i.strip() for i in line if i.strip() != '']
                if not line[1].startswith("Locus tag:"): # checking that there are no annotatinos already
                    del line[1]
                if not all([line[1].startswith("Locus tag:"), line[2].startswith("Protein ID:"), line[3].startswith("Gene:")]):
                    raise ValueError("something seems wrong processing line:\n{}".format(line))
                del line[1:4] # Removing these 3 elements we don't want
                # removing the formating around the gene id and add newline before annotations
                line[0] = re.sub("<span class=.\"serif.\">", "", line[0])
                line[0] = re.sub("</span></strong>", "</strong><br>", line[0])
                # Extract and cleanup locus tag
                locus_tag = line[0].split("strong")[1]
                locus_tag = re.sub(".*>", "", locus_tag)
                locus_tag = re.sub("<.*", "", locus_tag)
                # Adding the annotation
                line[0] = "{indent}{init}<br>{keg}<br>{prk}<br>{egg}<br>".format(indent=indent,
                                                                                 init=line[0],
                                                                                 keg=kegg_dict[locus_tag],
                                                                                 prk=prokka_dict[locus_tag],
                                                                                 egg=eggnog_dict[locus_tag])
                # Adding a newline after the location
                line[1] = line[1] + "<br>"
                # Tried to add line break before blasts, did not work
                #newline = re.sub("<div class=\"focus-urls\">", "<div class=\"focus-urls\"> &nbsp;<br>", "<br>".join(line))
                target.write("<br>".join(line) + "\n")
            else:
                target.write(line + "\n")

def main(args):
    """
    This is the main function of the script
    It will loop through the genomes found in the annotation folder and update the relevant antiSMASH
    files
    """
    annotation_path = Path(args.annotations)
    if not annotation_path.exists():
        raise FileNotFoundError("Looks like there is something wrong with {}.".format(annotation_path))
    genome_count = 0
    for genome_path in annotation_path.glob('*'):
        genome = genome_path.name
        # get the dictionaries
        prokka_file = [f for f in genome_path.joinpath("prokka").glob("{genome}-*-prokka/{genome}-prokka.tsv".format(genome=genome))][0]
        prokka_dict = load_prokka_annotations(prokka_file)
        kegg_dict = load_kegg_annotations(genome_path.joinpath("kegg", "{}-kegg-annotations.tsv".format(genome)))
        eggnog_dict = load_eggnog_annotations(genome_path.joinpath("eggnog", "{}-eggnog.emapper.annotations".format(genome)))
        # update the relevant files
        update_antismash_index(html_index=genome_path.joinpath("antismash", "{}-antismash".format(genome), "index.html"),
                               genome=genome)
        update_antismash_display(js_file=genome_path.joinpath("antismash", "{}-antismash".format(genome), "regions.js"),
                                 prokka_dict=prokka_dict,
                                 kegg_dict=kegg_dict,
                                 eggnog_dict=eggnog_dict)

if __name__ == '__main__':
    args = get_arguments()
    magpipe.log.print_arguments(args)
    main(args)
