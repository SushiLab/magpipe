#!/usr/bin/bash python
# -*- coding: utf-8
# pylint: disable=line-too-long

"""
Functions needed to set up the pipeline
"""

import yaml
import pandas
import subprocess
from pathlib import Path

from magpipe import pretty
from magpipe import external as magext

class setup_magpipe:
    """
    A handy class that sets up the pipeline and stores what's needed to get goin'
    """
    def __init__(self, config, workflow):
        """
        This functions takes the snakemake config and workflow objects
        config is based on the configfile used to specify the pipeline's parameters
        The functions adds all the needed variables for the pipeline to run smoothly.
        """
        # Just setting up things
        pretty.print_header("MAGPIPE", color = 'green', type = 'normal')
        pretty.print_single("Starting by purging the module system.", nl_after = 0, color = 'green')
        subprocess.run(["ml", "purge"], shell = True)

        # Stuff regarding snakemake itself
        self.magpipe_path = Path(workflow.snakefile).resolve().parent.parent
        self.config = config
        if config['envs'] is not None:
            self.envs_dir = config['envs'].rstrip('/')
        else:
            self.envs_dir = str(self.magpipe_path.joinpath('envs'))
        pretty.print_pair("Running snakemake file", Path(workflow.snakefile).name, nl_before = 1)
        pretty.print_pair("Located in folder", Path(workflow.snakefile).resolve().parent)
        pretty.print_pair("Using working directory", workflow._workdir) 
        pretty.print_pair("With magpipe path", self.magpipe_path)
        pretty.print_pair("Using environments", self.envs_dir)
        pretty.print_pair("Running in debug mode", config["debug"], nl_after = 1)

        # Prepare dependencies and custom scripts
        self.scripts = str(self.magpipe_path.joinpath('scripts'))
        with open(self.magpipe_path.joinpath('resources/dependencies.yaml')) as file:
            self.dependencies = yaml.full_load(file)

        # Get input ready
        if "datasets" in config.keys() and config["datasets"] is not None:
            # Checking if there are ctually some metagenomic datasets  to process
            self.datasets = sorted(config['datasets']) # sorted list of datasets to be processed
            if "dataset_info" in config.keys() and config["dataset_info"] is not None: # option to override default table
                self.dataset_info = config["dataset_info"]
                self.input = pandas.read_csv(config["dataset_info"], sep="\t")
            else:
                self.dataset_info = str(self.magpipe_path.joinpath('resources/datasets.tsv'))
                self.input = pandas.read_csv(self.dataset_info, sep="\t")
            self.input = self.input.set_index("datasets", drop = False)
            self.input = self.input.loc[self.datasets]
            # format depth suffixes:
            self.input["depth_sfx"] = ["_a{}".format(int(i)) for i in list(self.input["diffcov"])]
            # use str(Path()) to remove trailing '/' and uniformise paths:
            self.input.loc[:, 'assemblies'] = [str(Path(i)) for i in list(self.input.loc[:, 'assemblies'])]
            self.input.loc[:, 'backmapping'] = [str(Path(i)) for i in list(self.input.loc[:, 'backmapping'])]
            # FIXME: Add schemas for validation
            # e.g. validate(self.input, schema = "schemas/samples.schema.yaml")
            pretty.print_single('The pipeline is moving forward using this input:', color = 'green')
            pretty.melt_table(self.input)
        else:
            pretty.print_single('No metagenomic datasets to process. That can make sense if you only want to work with external genomes. Just checking with you.', color = 'yellow')
            self.dataset_info = "EXTERNAL_DUMMY" #FIXME dummy to prevent unassigned variable in snakemake..

        # Prepare Output
        self.output_path = Path(config['output'])
        self.metag_path = self.output_path.joinpath('metagenomes')
        self.metag_dir = str(self.metag_path)
        self.fast_path = Path(config['fast_drive'])
        self.fast_dir = str(self.fast_path)
        self.mags_path = self.output_path.joinpath('MAGs', config['pooled_name'])
        self.mags_dir = str(self.mags_path)
        self.pooled_name = config['pooled_name'].split("/")[-1] # account for the possibility of 'folder/pooled_name'
        self.integrate_path = self.output_path.joinpath('integrated', config['pooled_name'])
        self.integrate_dir = str(self.integrate_path)
        pretty.print_single('The output will be structured as follows:', color = 'green', nl_before = 1)
        pretty.print_pair("The main output folder is", str(self.output_path))
        pretty.print_pair("Each metagenomic dataset will be processed in", self.metag_dir)
        pretty.print_pair("MAGs will be pooled in", self.mags_dir)
        pretty.print_pair("Integrations (ext genomes + MAGs) will be in", self.integrate_dir)
        pretty.print_pair("The temp output / anvio databases will be in", self.fast_dir)

        # Specify a binning method
        self.method_info = str(self.magpipe_path.joinpath('resources/implemented_methods.tsv'))
        self.methods = pandas.read_csv(self.method_info, sep="\t")
        self.methods = self.methods.set_index("method", drop = False).fillna("")
        if len(config['binning']) != 1:
            raise NotImplementedError("Can only run one and only one binning method.")
        self.methods = self.methods.loc[config['binning']]
        if not type(self.methods) is pandas.core.frame.DataFrame:
            raise TypeError("Needs methods to be a dataframe type. It is a {}.".format(type(self.methods)))
        pretty.print_single('The pipeline is using the following binning method(s):', color = 'green', nl_before = 1)
        pretty.print_table(self.methods, nl_before = 0, nl_after = 1)
        #pretty.print_single('DEBUG: magpipe_config.methods = {}.'.format(type(self.methods)), color = 'red', nl_after = 1)
        #pretty.print_single('DEBUG: magpipe_config.methods.suffix = ' + self.methods.suffix[0], color = 'red', nl_after = 1)

        # MAGs quality filters
        self.cpl = int(config['filter']['completeness'])
        self.ctn = int(config['filter']['contamination'])
        pretty.print_single('Filtering settings:', color = 'green')
        pretty.print_pair("Completeness", '{}%'.format(self.cpl))
        pretty.print_pair("Contamination", '{}%'.format(self.ctn), nl_after = 1)

    def get_targets(self):
        """
        This functions generates and return the targets
        Additionally stores convenient flags, could be replaced by checkpoints
        The pipeline is divided in 3 phases
        Phase 1 is for each metagenome
            Phase 1 is divided in two subparts
            Phase 1-A is at the scaffold-level
            Phase 1-B is at the bin-level
        Phase 2 is for the analyses on the pooled MAGs
        Phase 3 is for analyses on each MAG individually
        """
        #FIXME also can probably make a better use of functions and expand()..
        # Init lists, FIXME maybe generators wouuld make things smoother
        targets = []
        self.flag_export_mags = []
        self.flag_integrate_external = []
        self.flag_gene_calling = []

        # Some convenient stuff
        self.targets = self.config["targets"] # These are the targets specified in config, used to generate the snakemake targets

        num_phases = len(self.targets)
        if num_phases != 1:
            raise NotImplementedError("The pipeline can only run one phase at a time. You gave {} phases... Please comment unused phases.".format(num_phases))

        # Phase 1, metagenome specific
        if "phase_1" in self.targets.keys() and self.targets["phase_1"] is not None:
            pretty.print_pair("Running magpipe", 'Phase 1', nl_after = 1)
            # Loop over datasets and binning methods
            for d in self.datasets:
                # load samples (check for debug state)
                if self.config['debug']: # Load only a small sample
                    magpipe.pretty.print_single("Running in debug mode!")
                    samples = [str(self.input.loc[d, 'debug'])]
                else:
                    sample_file = Path(self.input.loc[d, 'samples'])
                    samples = set(sample_file.read_text().splitlines())
                # loop over samples
                for s in samples:
                    # add phase 1-A targets
                    if "A" in self.targets["phase_1"].keys() and self.targets["phase_1"]["A"] is not None:
                        for t in self.targets['phase_1']["A"]:
                            if t == "depth":
                                if not "depth" in self.input.loc[d, "backmapping"]:
                                    targets.append(self.metag_path.joinpath(d, s, t, "{s}-{t}.done".format(s = s, t = t)))
                            else:
                                targets.append(self.metag_path.joinpath(d, s, t, "{s}-{t}.done".format(s = s, t = t)))
                    # add phase 1-B targets
                    if "B" in self.targets['phase_1'].keys() and self.targets['phase_1']["B"] is not None:
                        # loop over binning methods
                        for b in self.config["binning"]:
                            # write down the method's name
                            m = self.methods.loc[b, "tool"] + self.input.loc[d, "depth_sfx"] + self.methods.loc[b, "suffix"]
                            for t in self.targets['phase_1']["B"]:
                                targets.append(self.metag_path.joinpath(d, s, m, "{s}-{m}-{t}.done".format(s = s, m = m, t = t)))

        # Optionally adding some genomes
        # FIXME combine external genomes and exported MAGs for bulk processing
        if self.config["external_genomes"] is not None:
            pretty.print_single('Running magpipe with external genomes:', color = 'green')
            for external_dataset in self.config["external_genomes"]:
                pretty.print_pair("Processing external dataset", external_dataset)
                magext.setup_ext_evaluation(external_dataset)
                targets.extend(magext.get_ext_checkm_targets(external_dataset))
                self.flag_integrate_external.extend(magext.get_ext_checkm_targets(external_dataset))
                targets.extend(magext.get_ext_anvi_evaluate_targets(external_dataset))
                self.flag_integrate_external.extend(magext.get_ext_anvi_evaluate_targets(external_dataset))
            pretty.print_single('Done.', color = 'green', nl_after = 1)

        # Phase 2, processing genomes in bulk
        if "phase_2" in self.targets.keys() and self.targets["phase_2"] is not None:
            pretty.print_pair("Running magpipe", 'Phase 2', nl_after = 1)
            # Get input
            # Loop over datasets and binning methods
            for d in self.datasets:
                # load samples (check for debug state)
                if self.config['debug']: # Load only a small sample
                    pretty.print_single("Running in debug mode!", color = "yellow", nl_after = 1)
                    samples = [str(self.input.loc[d, 'debug'])]
                else:
                    sample_file = Path(self.input.loc[d, 'samples'])
                    samples = set(sample_file.read_text().splitlines())
                # loop over samples
                for s in samples:
                    # loop over binning methods
                    for b in self.config["binning"]:
                        # write down the method's name
                        m = self.methods.loc[b, "tool"] + self.input.loc[d, "depth_sfx"] + self.methods.loc[b, "suffix"]
                        # add the phase 2 input
                        self.flag_export_mags.append(self.metag_path.joinpath(d, s, m, "{s}-{m}-checkm.done".format(s = s, m = m)))
                        #FIXME If anvio fails because there is no bins, you get stuck / currently manually touch the markers that failed
                        self.flag_export_mags.append(self.metag_path.joinpath(d, s, m, "{s}-{m}-anvi_evaluate.done".format(s = s, m = m)))
            # Loop through binning methods
            for b in self.config["binning"]:
                m = self.methods.loc[b, "tool"] + self.methods.loc[b, "suffix"]
                if self.config["external_genomes"] is not None:
                    pretty.print_single('Assuming external genomes should be included here.', color = 'green')
                    continue
                # This one is always here if phase_2 because... otherwise there is no phase_2
                targets.append(self.mags_path.joinpath("{mags}-{m}-cpl{cpl}_ctn{ctn}-export_evaluation.done".format(mags = self.pooled_name, m = m, cpl = self.cpl, ctn = self.ctn)))
                # Add the phase 2 targets
                for t in self.targets["phase_2"]:
                    if t == 'drep':
                        t = "{t}_{ani}".format(t = t, ani = self.config["drep_ani"])
                    targets.append(self.mags_path.joinpath("{mags}-{m}-cpl{cpl}_ctn{ctn}-{t}.done".format(mags = self.pooled_name, m = m, cpl = self.cpl, ctn = self.ctn, t = t)))
            if self.config["external_genomes"] is not None:
                if len(self.config["binning"]) != 1:
                    #  Need to break if several methods as interate_external cannot deal with that
                    raise ValueError("You want me to integrate external genomes at phase_2 but gave me {} binning methods. I can't really deal with that, you need to make up your mind on a signle one.".format(len(self.config["binning"])))
                m = "integrated"
                # This one is always here if phase_2 because... otherwise there is no phase_2 FIXME paused for now
                #targets.append(self.integrate_path.joinpath("{mags}-{m}-cpl{cpl}_ctn{ctn}-export_evaluation.done".format(mags = self.pooled_name, m = m, cpl = self.cpl, ctn = self.ctn)))
                # This one is always here if you want to include your external genomes.
                #targets.append(self.integrate_path.joinpath("{mags}-{m}-cpl{cpl}_ctn{ctn}-integrate_external.done".format(mags = self.pooled_name, m = m, cpl = self.cpl, ctn = self.ctn)))
                # Add the phase 2 targets
                for t in self.targets["phase_2"]:
                    if t == 'drep':
                        t = "{t}_{ani}".format(t = t, ani = self.config["drep_ani"])
                    targets.append(self.integrate_path.joinpath("{mags}-{m}-cpl{cpl}_ctn{ctn}-{t}.done".format(mags = self.pooled_name, m = m, cpl = self.cpl, ctn = self.ctn, t = t)))

        # Phase 3, genome specific
        if "phase_3" in self.targets.keys() and self.targets["phase_3"] is not None:
            pretty.print_pair("Running magpipe", 'Phase 3', nl_after = 1)
            if len(self.config["binning"]) != 1:
                raise ValueError("You gave {} binning methods for phase_3. I can't really see why... so you need to make up your mind on a single one. Run this several times if you really need the different binning methods.".format(len(self.config["binning"])))
            b = self.config["binning"][0]
            if self.config["external_genomes"] is not None:
                pretty.print_single('Assuming external genomes should be included here.', color = 'green')
                m = "integrated"
                phase_3_path = self.integrate_path
            else:
                m = self.methods.loc[b, "tool"] + self.methods.loc[b, "suffix"]
                phase_3_path = self.mags_path
            # Update the output path
            phase_3_path = phase_3_path.joinpath("{mags}-{m}-cpl{cpl}_ctn{ctn}-annotations".format(mags = self.pooled_name, m = m, cpl = self.cpl, ctn = self.ctn))
            # Read the lists of bacterial and archaeal MAGs
            def _get_targets_from_taxo_files(self, phase_3_path, type):
                """
                Convenience function that differentiates archaea and bacteria
                Needed for prokka
                type = 'Bacteria' or 'Archaea'
                """
                targets = []
                file_ext_dict = {"Bacteria" : "bacteria", "Archaea" : "archaea"}
                mags_file = Path(phase_3_path.parent.joinpath("{mags}-{m}-cpl{cpl}_ctn{ctn}-dictionaries/{mags}-{m}-cpl{cpl}_ctn{ctn}-{type}.txt".format(mags=self.pooled_name, m=m, cpl=self.cpl, ctn=self.ctn, type=file_ext_dict[type])))
                # Just checking and returning approriate info if file missing or empty
                if not mags_file.exists():
                    pretty.print_single("Couldn't file for {}, ({}). Assuming there are none and moving on. Please check if that seems wrong".format(type, mags_file), color="yellow")
                    return
                mags = set(mags_file.read_text().splitlines())
                if len(mags) == 0:
                    pretty.print_single("File {} seems empty. Assuming there are no {} and moving on. Please check if that seems wrong".format(mags_file, type), color="yellow")
                    return
                # Loop through the MAGs
                for mag in mags:
                    for t in self.targets["phase_3"]:
                        target_path = phase_3_path.joinpath(mag, t)
                        if t == 'prokka':
                            mag_prokka_target = target_path.joinpath("{mag}-{type}-{t}.done".format(mag=mag, type=type, t=t))
                            targets.append(mag_prokka_target)
                            self.flag_gene_calling.append(mag_prokka_target)
                        else:
                            targets.append(target_path.joinpath("{mag}-{t}.done".format(mag=mag, t=t)))
                # Just a helper here #FIXME did I actually end up using that? if not remove
                if type  == 'Bacteria':
                    self.bacterial_mags = mags
                if type == 'Archaea':
                    self.archaeal_mags = mags
                # Return the targets
                return targets
            # Using the function for both Archaea and Bacteria
            phase_3_archaeal_targets = _get_targets_from_taxo_files(self, phase_3_path, "Archaea") 
            if phase_3_archaeal_targets is not None:
                targets.extend(phase_3_archaeal_targets)
            phase_3_bacterial_targets = _get_targets_from_taxo_files(self, phase_3_path, "Bacteria")
            if phase_3_bacterial_targets is not None:
                targets.extend(phase_3_bacterial_targets)

        # Exit if something went wrong
        if len(targets) == 0:
            raise ValueError("Didn't find any target...")

        # All set, just return the targets, which should be used for rule all
        return(targets)
