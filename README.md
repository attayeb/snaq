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
* Clone this repository or copy the zipped release file and extract it in your favourite folder.

## How to do the analysis:

* Create a new folder inside data folder use only letters in capital to name it, no spaces. Eg: AB, CONTROL, COHORTONE etc.


* Copy your paired-end fastq files to the folder, check the identifier of the R1 and R2 is it is \_R1\_ and \_R2\_ or _1 and _2 then Snaq will understand and differentate the R1 and R2 files. If any other identifiers are used, please prepare a manifest file and copy it to results/COHORT/ folder.\
Snaq will follow that manifest file if you provide it. Keep a copy of that manifest file somewhere outside the pipeline folder, because it could be overwritten by mistake.
* You need to send the snakemake command with basically two needed parameters ```--cores 10 --use-conda```, These two parameters are essential to run the analysis.
* After applied these two parameters you write the target. For example for an artifact of AB the target should be ```results/AB/AB.qza```. Snakemake will understand to import the data set saved in ```data/AB``` to ```AB.qza```; That will be done in two steps, first a manifest file is created, ```result/AB/AB_manifest.qza``` and then the real ```results/AB/AB.qza``` will be created.
