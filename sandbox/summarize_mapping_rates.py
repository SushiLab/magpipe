#!/usr/bin/env python

import argparse
from pathlib import Path

# we need to glob all the files endind in '.mappingrates', extract the input reads and output reads and keep sample in memory.

def get_arguments():
    """
    Get commandline arguments and return namespace
    """
    # Initialize Parser
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', help='Path to the mapping rates folder. Expects samples folder with *.mappingrates files.', required=True, type=str)
    parser.add_argument('-o', '--output', help='Path to the ouptut table.', required=True, type=str)
    return parser.parse_args()

def main(args):
    """
    glob the given directory
    """
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError("Couldn't find the error file...")
    res = []
    for f in input_path.glob("*/1mappingrates/*/*.mappingrates"):
        print("Processing file {}".format(f))
        sample = f.parent.parent.parent.name
        bamname = f.name
        with open(f) as handle:
            for line in handle:
                if line.startswith("Input Reads"):
                    input_reads = line.strip().split('\t')[-1]
                if line.startswith("Output Reads"):
                    output_reads = line.strip().split('\t')[-1]
            res.append("\t".join([sample, bamname, input_reads, output_reads]))
    with open(args.output, "w") as target:
        target.write("\t".join(["Sample", "File", "Input_Reads", "Output_Reads"]) + "\n")
        target.write("\n".join(res) + "\n")

if __name__ == '__main__':
    args = get_arguments()
    main(args) 
