# magpipe 

*Like a bagpipe, but for MAGs.*

# Code associated with the paper 'Biosynthetic potential of the global ocean microbiome' and the Ocean Microbiomics Database.

Link to the [preprint](https://www.biorxiv.org/content/10.1101/2021.03.24.436479v1) and the [paper](). The Ocean Microbiomics Database is available here: [https://www.microbiomics.io/ocean/](https://www.microbiomics.io/ocean/).

This repo contains the code behing the snakemake pipeline used for the generation of the MAGs, their analyses and the R code used for the figures. 

## Structure

- [**configs**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/configs): config files used by the launchers to run the pipeline.
- [**envs**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/envs): conda environments for the rules.
- [**figures**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/figures): R code behind the figures.
- [**launchers**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/launchers): bash scripts to run the pipeline (on sge etc).
- [**magpipe**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/magpipe): python module for the pipeline.
- [**resources**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/resources): some fixed input information on methods and datasets.
- [**rules**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/rules): the rules of the snakemake pipeline. 
- [**sandbox**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/sandbox): scripts used for postprocessing and analyses.
- [**scripts**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/scripts): ad-hoc python scripts used by the pipeline.
- [**snakes**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/snakes): snakemake files.

You can find additional documentation of the different metagenomic analyses steps here: [https://methods-in-microbiomics.readthedocs.io/en/latest/](https://methods-in-microbiomics.readthedocs.io/en/latest/). 

# Running the pipeline:

## Installation

We highly recommend using conda to create a dedicated environment, as follows:

```
conda create -n mapgipe
conda activate magpipe
conda install snakemake
```

For the pipeline to work, you need to install the python module as follows:
 
```
git clone git@github.com:SushiLab/magpipe.git
cd magpipe/magpipe
python -m pip install -r requirements.txt -e .
```

## Configuration

You now need to setup a few things:

- In `resources`, you'll need to update the dataset table and add the name and path to a given metagenomic dataset, the corresponding assemblies and depth files.
- In config, you'll need to specify the paths for the output, snakemake workdir, and fast drive (if appropriate).
- Once configured, you can use the snakemake command line to start the snakemake pipeline. 
