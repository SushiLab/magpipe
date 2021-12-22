# magpipe 

*Like a bagpipe, but for MAGs.*

This repo contains the code behing the snakemake pipeline used for the generation of the MAGs used in this [paper](https://www.biorxiv.org/content/10.1101/2021.03.24.436479v1).

## Structure
a
- [**magpipe**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/magpipe): python module for the pipeline.
- [**snakes**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/snakes): snakemake files.
- [**rules**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/rules): the rules of the snakemake pipeline. 
- [**configs**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/configs): config files used by the launchers to run the pipeline.
- [**scripts**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/scripts): ad-hoc python scripts and the magpipe python module.
- [**launchers**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/launchers): bash scripts to run the pipeline (on sge etc).
- [**resources**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/resources): some fixed input information on methods and datasets.
- [**envs**](https://github.com/SushiLab/MAGPIPE_DEV/tree/master/code/envs): conda environments for the rules.

## Installation

We highly recommend using conda

```
conda create -n mapgipe
conda activate magpipe
conda install snakemake
```

For the pipeline to work, you need to install the python module as follows
 
```
git clone git@github.com:SushiLab/magpipe.git
cd magpipe
python -m pip install -r requirements.txt -e .
```

## Configuration

You now need to setup a few things:

- In `resources`, you'll need to update the dataset table and add the name and path to a given metagenomic dataset, the corresponding assemblies and depth files.
- In config, you'll need to specify the paths for the output, snakemake workdir, and fast drive (if appropriate).
- Once configured, you can use the snakemake command line to start the snakemake pipeline. 
