#!/usr/bin/env python
# -*- coding: utf-8

"""
name: magpipe, phases 1, 2 and 3
description: Snakemake pipeline to build MAGs, sample-specific processing
author: Lucas Paoli (paolil@ethz.ch)
dependencies:
    - depth
    - binning
    - evaluate
    - export
    - annotate
    - dereplicate
"""

import re
from pathlib import Path
from magpipe import pretty
from magpipe import external
from magpipe.setup_snake import setup_magpipe 

# Setup ========================================================================
# Set things up and prepare variables used throughout the pipeline
# config object coming from the --configfile flag when running the Snakefile

# Updating workdir and getting ready
workdir: config["snake_workdir"]
magpipe_config = setup_magpipe(config, workflow)

# Define the targets ===========================================================

TARGETS = magpipe_config.get_targets()

# One rule to rule them all ====================================================

#pretty.print_single("Starting snakemake 🐍",  nl_after = 1, color = 'green')
pretty.print_single("Starting snakemake",  nl_after = 1, color = 'green')
rule all:
    input: TARGETS

# Include the required rules ===================================================

include: str(magpipe_config.magpipe_path) + "/rules/depth.smk"
include: str(magpipe_config.magpipe_path) + "/rules/export.smk"
include: str(magpipe_config.magpipe_path) + "/rules/binning.smk"
include: str(magpipe_config.magpipe_path) + "/rules/evaluate.smk"
include: str(magpipe_config.magpipe_path) + "/rules/dereplicate.smk"
include: str(magpipe_config.magpipe_path) + "/rules/annotate-genomes.smk"
include: str(magpipe_config.magpipe_path) + "/rules/annotate-metagenomes.smk"
include: str(magpipe_config.magpipe_path) + "/rules/identify-elements-genomes.smk"
include: str(magpipe_config.magpipe_path) + "/rules/identify-elements-metagenomes.smk"
