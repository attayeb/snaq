# Snaq Pipeline

Snaq is a snakemake pipeline for Microbiome data analsysis using QIIME2

## Citation:

## Installation

This pipeline works natively in Linux, Mac and windows on WSL (basically ubuntu). It also can run using docker container system in windows.

### Linux and Mac and windows using [WSL](https://docs.microsoft.com/en-us/windows/wsl/install)
* [Install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

* Install mamba:
```bash
conda install -n base -c conda-forge mamba
```
* [Install snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html):
```bash
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```
clone this repository or download zipped release file:
```bash
$ To be determined soon
```

### Windows with docker:
* [Install docker desktop in windows](https://docs.docker.com/desktop/windows/install/)
* Using command  or windows [PowerShell](https://en.wikipedia.org/wiki/PowerShell) clone docker image for snakemake by sending this command:
```
docker pull snakemake/snakemake
```
* clone this repository or copy the zipped release file and extract it in your favourite folder.

