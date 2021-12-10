# Snaq Pipeline

Snaq is a snakemake pipeline for Microbiome data analsysis using QIIME2

This pipeline works natively in Linux, Mac and windows on WSL (basically ubuntu). It also can run using docker container system in windows.

## Citation:
TBD



## Installation



### Linux and Mac and windows using [WSL](https://docs.microsoft.com/en-us/windows/wsl/install)

All these steps should be executed in the terminal

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
download the latest release file of this repository:
```bash
$ wget https://github.com/attayeb/snaq/archive/refs/tags/testing.zip
$ unzip testing.zip
```

Test the integrity of the pipeline by executing:
```
cd snaq-testing2
snakemake -lt
```
### Windows with docker:
* [Install docker desktop in windows](https://docs.docker.com/desktop/windows/install/)
* Using command  or windows [PowerShell](https://en.wikipedia.org/wiki/PowerShell) clone docker image for snakemake by sending this command:
```
docker pull snakemake/snakemake
```
* Download the latest release [Source code (zip)](https://github.com/attayeb/snaq/archive/refs/tags/testing2.zip) from github repository, and extract it.

* check the integrity of the pipeline by sending this command:
```
docker run -it snakemake/snakemake snakemake -lt
```

## How to do the analysis:

* Create a new folder inside data folder use only letters in capital to name it, no spaces. Eg: AB, CONTROL, COHORTONE etc.


* Copy your paired-end fastq files to the folder, check the identifier of the R1 and R2 is it is \_R1\_ and \_R2\_ or _1 and _2 then Snaq will understand and differentate the R1 and R2 files. If any other identifiers are used, please prepare a manifest file and copy it to results/COHORT/ folder.\
Snaq will follow that manifest file if you provide it. Keep a copy of that manifest file somewhere outside the pipeline folder, because it could be overwritten by mistake.
* You need to send the snakemake command with basically two needed parameters ```--cores 10 --use-conda```, These two parameters are essential to run the analysis.
* After applied these two parameters you write the target. For example for an artifact of AB the target should be ```results/AB/AB.qza```. Snakemake will understand to import the data set saved in ```data/AB``` to ```AB.qza```; That will be done in two steps, first a manifest file is created, ```result/AB/AB_manifest.qza``` and then the real ```results/AB/AB.qza``` will be created.

* For a docker you need to send the command docker to run the pipeline. It is simple, you need to pull the image of snakemake, 

* Docker need few tricks, it runs inside a container, it is like running another computer inside your machine, this computer is isolated from the host if you don't tell the container to use your file system, moreover, if you use the file system inside the docker container, you will lose the data as soon as you stop your docker container. So, what we do is snakemake is already installed inside the image we are going to use. To run snakemake we need to give it access to our pipeline folder. so give it in the command line like this:
    - map your working directory to ```/work``` directory inside the contianer by using this command ```-v c:\snaq\:/work```
    - tell docker to use ```/work``` folder as your working directory. so when you go to ```/work``` directory, whatever you do will be saved in your working directory outside.
    - Start running commands after you set these parameters and other parameters like ```-it``` to make the container interactive and terminate it as soon as the code is finished. also we need to give the container name ```snakemake/snakemake``` This [image] (https://hub.docker.com/r/snakemake/snakemake) is created by snakemake developers. \ 


```
-v c:\snaq\:/work -w /work 
```

```
docker run -it -v d:\dustbox\snaq-test\:/work -w /work snakemake/snakemake snakemake --use-conda --cores 10 results/SRA/SRA+bb16t+fp-f17-r21crop+dd+cls-silva+rrf10000.zip
```

## Example dataset:
is available for testing [test data set](https://github.com/attayeb/snaq/releases/download/testing/AB.tar.gz); download it and extract it in data folder.