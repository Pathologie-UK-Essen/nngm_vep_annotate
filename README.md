# nNGM variant file annotation

This workflow annotates vcf files as definied by the nNGM-project.


## Requirements

* [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge)
* on Apple Silicon devices [Rosetta2](https://support.apple.com/en-us/HT211861) needs to be installed
## Setup

### Clone this repository git clone
`https://github.com/FelixMoelder/nngm_vep_annotate`

### [Optional] Change default input- and output-directories

By default vcf-files need to be placed in a subdirectory called `input` and will be written to an `output` folder.
An alternative input- and output-path can be set in `config.yaml`. Paths can be relative and absolute.

### Make the script executable
`chmod +x init_cronjob.sh`

### Setup a cronjob by calling

To initialize a conda environment and setup a cronjob executeing the workflow once an hour call:

`sh ./init_cronjob.sh`
