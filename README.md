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
docker run -it -v [snaq folder in windows]:/snaq -w /snaq snakemake/snakemake snakemake --use-conda --cores 10 results/AB/AB+bb16t+fp-f17-r21crop+dd+cls-silva+rrf10000.zip
```

## Example dataset:
is available for testing [test data set](https://github.com/attayeb/snaq/releases/download/testing/AB.tar.gz); download it and extract it in data folder.

Possible targets are:

```

AB+bb12t.qza

AB+bb14t.qza
AB+bb16t+fp-f17-r21crop+dd+cls-gg_asv.biom
AB+bb16t+fp-f17-r21crop+dd+cls-gg+phyloseq.RDS
AB+bb16t+fp-f17-r21crop+dd+cls-gg+rrf10000+beta_braycurtis.tsv
AB+bb16t+fp-f17-r21crop+dd+cls-gg+rrf10000+beta_jaccard.tsv
AB+bb16t+fp-f17-r21crop+dd+cls-gg+rrf10000+manta_tax.tsv
AB+bb16t+fp-f17-r21crop+dd+cls-gg+rrf10000+manta.tsv
AB+bb16t+fp-f17-r21crop+dd+cls-gg+rrf10000+otu_tax.biom
AB+bb16t+fp-f17-r21crop+dd+cls-gg+rrf10000+otu_tax_biom.tsv
AB+bb16t+fp-f17-r21crop+dd+cls-gg+rrf10000+otu_tax.qza
AB+bb16t+fp-f17-r21crop+dd+cls-gg+rrf10000.zip
AB+bb16t+fp-f17-r21crop+dd+cls-gg_taxonomy.csv
AB+bb16t+fp-f17-r21crop+dd+cls-gg_taxonomy.qza
AB+bb16t+fp-f17-r21crop+dd+fasttree.nwk
AB+bb16t+fp-f17-r21crop+dd+fasttree_rooted.qza
AB+bb16t+fp-f17-r21crop+dd+rrf10000+alphadiversity.tsv
AB+bb16t+fp-f17-r21crop+dd+rrf10000+beta_unweightedunifrac.csv
AB+bb16t+fp-f17-r21crop+dd+rrf10000+beta_unweightedunifrac.qza
AB+bb16t+fp-f17-r21crop+dd+rrf10000+beta_weightedunifrac.csv
AB+bb16t+fp-f17-r21crop+dd+rrf10000+beta_weightedunifrac.qza
AB+bb16t+fp-f17-r21crop+dd_seq.csv
AB+bb16t+fp-f17-r21crop+dd_seq.qza
AB+bb16t+fp-f17-r21crop+dd_stats.qza
AB+bb16t+fp-f17-r21crop+dd_table.qza
AB+bb16t+fp-f17-r21crop+dd_table+rrf10000.csv
AB+bb16t+fp-f17-r21crop+dd_table+rrf10000.qza
AB+bb16t+fp-f17-r21crop.qza
AB+bb16t.qza
AB+bb18t.qza
AB+bb22t.qza
AB_manifest.tsv
AB.qza
```