<p align="center">
<img src="logo_snaq.png">
</p>

# Snaq Pipeline

Snaq is a snakemake pipeline for Microbiome data analsysis using QIIME2.

This pipeline works in Linux, Mac and Windows (Ubuntu on Windows). It also can run using Docker container system.

## Citation:

Mohsen A, Chen Y-A, Allendes Osorio RS, Higuchi C and Mizuguchi K (2022) Snaq: A Dynamic Snakemake Pipeline for Microbiome Data Analysis With QIIME2. Front. Bioinform. 2:893933. doi: [10.3389/fbinf.2022.893933](https://doi.org/10.3389/fbinf.2022.893933)


## Installation

### Windows Subsystem for Linux.

* Install Ubuntu for windows 10 following the instructions in this [website](https://ubuntu.com/tutorials/ubuntu-on-windows#1-overview) 

* In "Ubuntu bash command line" install and test the pipeline following the same steps mentioned in Linux and Mac.

### Linux and Mac

All these steps should be executed in the terminal (linux and Mac) or (Ubuntu bash command line) in windows, 

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
download the latest release file from this repository and extract it in a new folder, or clone this repository:

### Docker (works for Windows, Mac and Linux):
* [Install docker](https://docs.docker.com/get-docker/)
* Using Terminal, command prompt, or windows [PowerShell](https://en.wikipedia.org/wiki/PowerShell) depending on your system, clone docker image for snakemake by sending this command:
```
docker pull snakemake/snakemake
```
* Download the latest release [Source code (zip)](https://github.com/attayeb/snaq/archive/refs/tags/v1.0.1.zip) from github repository, and extract it, or clone this repository.

* check the integrity of the pipeline by sending this command:
```
docker run --rm -it -v "$PWD":/snaq -w /snaq snakemake/snakemake snakemake -lt
```

## How to do the analysis:


* Create a new folder inside data folder use only letters in capital to name it, no spaces. Eg: AB, CONTROL, COHORTONE etc.


* Copy your paired-end fastq files to the folder, check the identifier of the R1 and R2 is it is \_R1\_ and \_R2\_ or _1 and _2 then Snaq will understand and differentate the R1 and R2 files. If any other identifiers are used, please prepare a manifest file and copy it to results/<DATASET>/ folder.\
Snaq will follow that manifest file if you provide it. Keep a copy of that manifest file somewhere outside the pipeline folder, because it could be overwritten by mistake.
* You need to send the snakemake command with basically two needed parameters ```--cores <number of cores> --use-conda```, These two parameters are essential to run the analysis.
* After these two parameters you type the analysis target. For example to import the data to QIIME2, an artifact will be created with .qza extension. To do that for a cohort names "AB" the target should be ```results/AB/AB.qza```. Snakemake will understand to import the data set saved in ```data/AB``` folder to ```results/AB/AB.qza``` artifact file; That will be done in two steps, first a manifest file is created, ```result/AB/AB_manifest.qza``` and then the files listed in that manifest files will be imported to ```results/AB/AB.qza```.

## Few important points about docker
* Docker creates a container depending on an image, The image can be created or downloaded. The command ```docker pull snakemake/snakemake``` will download the required image to run sanakemake.
* When ```docker run -it snakemake/snakemake``` is executed, a container will be built which basically means running a small virtual linux PC inside your host system, and whatever command you send after that will run inside that virtual PC.
* For example if you send command ```docker run -it snakemake/snakemake snakemake``` then the snakemake software will run inside that created container. 
* As soon as you stop running the docker container it resets back to original status, whatever modification you make are not permenant.
* If you want changes to remain, you can link a folder from host machine to the container, and the container will save, modify, read from that folder, for that we use ```-v``` parameter in docker command. To map the pipeline folder to a folder inside the container you can:
```docker run -it -v c:\snaq:/snaq -w /snaq snakemake/snakemake```

* This command will start a container using sankemake image, and linke c:\snaq folder in the windows host PC to /work folder inside the container, and it make the working directory /snaq. So what ever command you send will run inside /snaq folder, any change done will be permenant in that folder.

To run a basic task in docker then the the command should be like this

```
docker run --rm -it -v [snaq folder in host system]:/snaq -w /snaq snakemake/snakemake snakemake --use-conda --cores 10 results/AB/AB+bb-t16+fp-f17-r21+dd+cls-silva+rrf10000.zip
```

## Example dataset:
is available for testing: can be downloaded from [here](https://github.com/attayeb/snaq/releases/download/data/AB.tar.gz); download it and extract it in data folder.

Possible targets are:

```
AB_manifest.tsv
AB.qza
AB+bb-t16+fp-f17-r21+dd+cls-gg_asv.biom
AB+bb-t16+fp-f17-r21+dd+cls-gg+phyloseq.RDS
AB+bb-t16+fp-f17-r21+dd+cls-gg+rrf-d10000+beta_braycurtis.tsv
AB+bb-t16+fp-f17-r21+dd+cls-gg+rrf-d10000+beta_jaccard.tsv
AB+bb-t16+fp-f17-r21+dd+cls-gg+rrf-d10000+manta_tax.tsv
AB+bb-t16+fp-f17-r21+dd+cls-gg+rrf-d10000+manta.tsv
AB+bb-t16+fp-f17-r21+dd+cls-gg+rrf-d10000+otu_tax.biom
AB+bb-t16+fp-f17-r21+dd+cls-gg+rrf-d10000+otu_tax_biom.tsv
AB+bb-t16+fp-f17-r21+dd+cls-gg+rrf-d10000+otu_tax.qza
AB+bb-t16+fp-f17-r21+dd+cls-gg+rrf-d10000.zip
AB+bb-t16+fp-f17-r21+dd+cls-gg_taxonomy.csv
AB+bb-t16+fp-f17-r21+dd+cls-gg_taxonomy.qza
AB+bb-t16+fp-f17-r21+dd+fasttree.nwk
AB+bb-t16+fp-f17-r21+dd+fasttree_rooted.qza
AB+bb-t16+fp-f17-r21+dd+rrf-d10000+alphadiversity.tsv
AB+bb-t16+fp-f17-r21+dd+rrf-d10000+beta_unweightedunifrac.csv
AB+bb-t16+fp-f17-r21+dd+rrf-d10000+beta_unweightedunifrac.qza
AB+bb-t16+fp-f17-r21+dd+rrf-d10000+beta_weightedunifrac.csv
AB+bb-t16+fp-f17-r21+dd+rrf-d10000+beta_weightedunifrac.qza
AB+bb-t16+fp-f17-r21+dd_seq.csv
AB+bb-t16+fp-f17-r21+dd_seq.qza
AB+bb-t16+fp-f17-r21+dd_stats.qza
AB+bb-t16+fp-f17-r21+dd_table.qza
AB+bb-t16+fp-f17-r21+dd_table+rrf-d10000.csv
AB+bb-t16+fp-f17-r21+dd_table+rrf-d10000.qza
AB+bb-t16+fp-f17-r21.qza
```
For more details please check the paper.

## More detailed documentation is in preparaion
