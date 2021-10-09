#!/bin/bash

if ! command -v conda &>/dev/null; then
    printf "conda has not been found. Please install mambaforge via https://github.com/conda-forge/miniforge\n"
    exit
fi

if ! command -v mamba &>/dev/null; then
    printf "mamba has not been found. Please install mamba to your base environment by calling 'conda install -n base -c conda-forge mamba'\n"
    exit
fi
#Create conda environment if it does not exist
if ! conda info --envs | grep vep_annotation >/dev/null; then
    mamba create -n vep_annotation snakemake -y
fi

#make cronjob if not existing

if [ $(crontab -l | grep snakefile_vep | wc -l) -eq 0 ]; then
    crontab -l >update_cronjob
    #if Apple Silicon cpu set conda channels to x86 architecture
    if [ $(uname) = "Darwin" ] && [ $(uname -m) = "arm64" ]; then
        CONDA_SUBDIR="CONDA_SUBDIR=osx-64"
    fi
    echo "47 * * * * (eval "$(conda shell.bash activate vep_annotation)" && cd $(pwd) && $CONDA_SUBDIR snakemake -s snakefile_vep.smk --use-conda --conda-frontend mamba -j 1 -p) >> $(pwd)/workflow.log 2>&1" >>updated_cronjob
    crontab updated_cronjob
    rm updated_cronjob
fi

echo "Workflow successfully initialized. If anything failes please file an issue on https://github.com/FelixMoelder/nngm_vep_annotate."
