"""
Snakemake pipeline for analyzing 16S RNA data using QIIME2

Usage
-----

The simplest command is:

snakemake --cores 10 --use-conda results/AB/AB+fp-f17-r21+bb-t18+cls-gg+rrf10000.zip
"""

from platform import system

_os = system()

if _os == "Linux":
     qiime_env = "envs/qiime2-2021.8-py38-linux-conda.yml"
elif _os == "Darwin":
     qiime_env = "envs/qiime2-2021.8-py38-osx-conda.yml"


rule export_artifact_2:
     """Export Artifact content to a folder"""
     message:
          "Extracting artifact"
     input:
          "results/{cohort, [A-Z]}/{cohort}.qza"
     output:
          directory("temp/{cohort}/{cohort}/")
     conda:
          qiime_env
     shell:
          "qiime tools export "
          "--input-path {input} "
          "--output-path {output}"



rule export_artifact:
     """Export Artifact content to a folder"""
     message:
          "Extracting artifact"
     input:
          "results/{cohort, [A-Z]}/{cohort}+{etc}.qza"
     output:
          directory("temp/{cohort}/{cohort}+{etc}")
     conda:
          qiime_env
     shell:
          "qiime tools export "
          "--input-path {input} "
          "--output-path {output}"


def get_allfile_names(wildcards):
     """Get all files from a folder"""
     input_folder = os.path.join("temp", wildcards.cohort, wildcards.id)
     return " ".join([os.path.join(input_folder, x) for x in os.listdir(input_folder) if "fastq" in x])


rule help:
     """Print out the information of rules on screen (docstrings)."""
     run:
          from termcolor import cprint
          for rule in workflow.rules:
               try:
                    if not rule.name.startswith('_'):
                         cprint(rule.name, "red", attrs=["bold"])
                         cprint(rule.docstring, "green")
                         cprint("")
               except:
                    pass
          
rule qza_fastqc:
     """Run fastqc quality control analysis"""
     message:
          "Applying FASTQC"
     input:
          "temp/{cohort}/{id}"          
     output:
          directory("results/{cohort}/quality/{id}/fastqc/")
     threads:
          20
     conda:
          "envs/quality.yml"
     params:
          get_allfile_names
     shell:
          "mkdir {output} && "
          "fastqc -o {output} -f fastq -t {threads} {params}"

rule qza_multiqc:
     """Combines multiple Fastqc reports using MultiQC. This rule works for one folder"""
     message:
          "MultiQC"
     input:
          "results/{cohort}/quality/{id}/fastqc/"
     output:
          directory("results/{cohort}/quality/{id}/multiqc/")
     conda:
          "envs/quality.yml"
     shell:
          "multiqc -o {output} {input}"

rule download_names_and_taxonpath:
     """Download database information from github"""
     output:
          taxonpath = "db/taxonpath.json",
          names = "db/names.json"
     shell:
          """
          wget https://github.com/attayeb/snaq/releases/download/testing/names.json
          wget https://github.com/attayeb/snaq/releases/download/testing/taxonpath.json
          """

rule dataset_multiqc:
     """Combines multiple Fastqc reports using Multiqc. This rule combines all the FastqC reports of one cohort"""
     message:
          "MultiQC"
     input:
          "results/{cohort}/quality/"
     output:
          directory("results/{cohort}/quality_summary/")
     conda:
          "envs/quality.yml"
     shell:
          "multiqc -dd 2 -d -o {output} {input}"

rule manifest:
     """Create manifest file: utilizing scripts/create_manifest_file.py script"""
     message:
          "Creating manifest file"
     input:
          "data/{cohort, [A-Z]}/"
     output:
          "results/{cohort}/{cohort}_manifest.tsv"
     conda:
          "envs/other.yml"
     shell:
          "python scripts/create_manifest_file.py -i {input} -o {output}"

rule import_data:
     """Import data: Import the raw fastq files to Qiime2 artifact with qza extension"""
     message:
          "Import data using manifest file"
     input:
          "results/{cohort}/{cohort}_manifest.tsv" 
     output:
          "results/{cohort}/{cohort}.qza"
     message:
          "Import data"
     conda: 
          qiime_env
     shell:
          "qiime tools import "
          "--type 'SampleData[PairedEndSequencesWithQuality]' "
          "--input-path {input} "
          "--input-format PairedEndFastqManifestPhred33V2 "
          "--output-path {output} "


rule trim_fastp:
     """Crop primers: Crop primers from both end reads."""
     input:
          qza="results/{cohort, [A-Z]}/{id}.qza"
     output:
          "results/{cohort}/{id}+fp-f{len1, \d+}-r{len2, \d+}.qza"
     message:
          "Trimming using fastp"
     conda: 
          qiime_env
     shell:
          "python scripts/fastp.py --inputf {input} "
          "--len1 {wildcards.len1} --len2 {wildcards.len2} "
          "--outputf {output}"


rule trim_bbduk:
     """Quality trimming using bbduk"""
     input:
          qza="results/{cohort, [A-Z]}/{id}.qza"
     output:
          "results/{cohort}/{id}+bb-t{threshold, \d+}.qza"
     message:
          "Trimming using bbduk"
     conda: 
          qiime_env
     shell:
          "python scripts/bbduk.py -i {input} "
          "-q {wildcards.threshold} -o {output}"


rule dada2:
     """Dada2 algorithm"""
     input:
          "results/{cohort, [A-Z]}/{id}.qza"
     output:
          table="results/{cohort}/{id}+dd_table.qza",
          stats="results/{cohort}/{id}+dd_stats.qza",
          repseq="results/{cohort}/{id}+dd_seq.qza"
     message:
          "Dada2 analysis"
     threads: 30
     conda: 
          qiime_env    
     shell:
          "qiime dada2 denoise-paired "
          "--p-trunc-len-f 0 --p-trunc-len-r 0 "
          "--i-demultiplexed-seqs {input} "
          "--o-table {output.table} "
          "--o-representative-sequences {output.repseq} "
          "--o-denoising-stats {output.stats} "
          "--verbose --p-n-threads {threads}"


rule rarefy:
     """rarefy feature table."""
     input:
          "results/{cohort}/{id}_table.qza"
     output:
          "results/{cohort}/{id}_table+rrf-d{r, \d+}.qza"
     conda:
          qiime_env
     shell:
          "qiime feature-table rarefy "
          "--i-table {input} "
          "--p-sampling-depth {wildcards.r} "
          "--o-rarefied-table {output}"

rule plot_dada_stats:
     """Plot dada2 stats results in form of distribution plot of the percentage of remained reads from the original data"""
     input:
          "results/{cohort}/{id}+dd_stats.qza"
     output:
          "results/{cohort}/{id}+dd_stats.jpg"
     conda:
          qiime_env
     shell:
          "python scripts/plot_dada.py --inp {input} --plot {output}"


def get_dada_jpgs(wildcards):
     """Get dada2 plot file names from a folder"""
     input_folder = os.path.join("results", wildcards.cohort)
     ret = [os.path.join(input_folder, x) for x in os.listdir(input_folder) if "+dd_stats.qza" in x]
     ret = [x.replace(".qza", ".jpg") for x in ret]
     return ret

def get_dada_jpgs_comma_separated(wildcards):
     """Get dada2 plot file names from a folder separated by comma"""
     input_folder = os.path.join("results", wildcards.cohort)
     ret = [os.path.join(input_folder, x) for x in os.listdir(input_folder) if "+dd_stats.qza" in x]
     ret = [x.replace(".qza", ".jpg") for x in ret]
     return ",".join(ret)

rule dada_stats_report:
     """Combines multiple dada2 stats plots"""
     input:
          get_dada_jpgs
     output:
          "results/{cohort}/{cohort}_dada_stats.pdf"
     conda:
          "envs/other.yml"
     params:
          get_dada_jpgs_comma_separated
     shell:
          "python scripts/report_stats.py --inp {params} --outp {output}"



rule download_silva_classifier:
     """Download pretrained SILVA taxonomy classifier"""
     output:
          "classifiers/silva-classifier.qza"
     conda:
          "envs/other.yml"
     shell:
          "cd classifiers && "
          "wget https://zenodo.org/record/5535616/files/silva-classifier.qza"

rule download_gg_classifier:
     """Download GreenGenes taxonomy classifier"""
     output:
          "classifiers/gg-classifier.qza"
     conda:
          "envs/other.yml"
     shell:
          "cd classifiers && "
          "wget https://zenodo.org/record/5535616/files/gg-classifier.qza"

rule download_silvav34_classifier:
     """Download V3-V4 region pretrained SILVA classifier."""
     output:
          "classifiers/silvaV34-classifier.qza"
     conda:
          "envs/other.yml"
     shell:
          "cd classifiers && "
          "wget https://zenodo.org/record/5535616/files/silvaV34-classifier.qza"


rule taxonomy:
     """Assign taxonomy to ASVs"""
     input:
          seq = "results/{cohort}/{id}_seq.qza",
          classifier = "classifiers/{cls}-classifier.qza"
     output:
          taxonomy= "results/{cohort}/{id}+cls-{cls}_taxonomy.qza",
     conda: 
          qiime_env
     threads: 30
     message:
          "Assign taxonomy using {wildcards.cls} database"
     shell:
          "qiime feature-classifier classify-sklearn "
          "--i-classifier {input.classifier} "
          "--i-reads {input.seq} --p-n-jobs {threads} "
          "--o-classification {output.taxonomy}"


rule mafft:
     """Creating phylogenetic tree step 1 """
     input:
          "results/{cohort}/{id}_seq.qza"
     output:
          "results/{cohort}/tree/{id}_seqaligned.qza"
     conda: 
          qiime_env
     message:
          "Mafft alignment"
     shell:
          "qiime alignment mafft "
          "--i-sequences {input} "
          "--o-alignment {output}"

rule alignment_mask:
     """Creating phylogenic tree step 2
     """
     input:
          "results/{cohort}/tree/{id}_seqaligned.qza"
     output:
          "results/{cohort}/tree/{id}_seqaligned_masked.qza"
     message:
          "Alignment mask"
     conda: 
          qiime_env     
     shell:
          "qiime alignment mask "
          "--i-alignment {input} "
          "--o-masked-alignment {output}"

rule fasttree:
     """Creating phylogenetic tree step 3 
     """
     input:
          "results/{cohort}/tree/{id}_seqaligned_masked.qza"
     output:
          "results/{cohort}/tree/{id}+fasttree.qza"
     conda: 
          qiime_env
     message:
          "creating tree"
     shell:
          "qiime phylogeny fasttree "
          "--i-alignment {input} "
          "--o-tree {output}"

rule midpoint_root:
     """Creating phylogenetic tree step 4
     """
     input:
          "results/{cohort}/tree/{id}+fasttree.qza"
     output:
          "results/{cohort}/{id}+fasttree_rooted.qza"
     message:
          "Midpoint rooting the tree"
     conda: 
          qiime_env
     shell:
          "qiime phylogeny midpoint-root "
          "--i-tree {input} "
          "--o-rooted-tree {output}"


rule export_tree:
     """Creating phylogenitic tree step 5: Exporting tree from Artifact file
     """
     input:
          "results/{cohort}/{id}+fasttree_rooted.qza"
     output:
          "results/{cohort}/{id}+fasttree.nwk"
     conda: 
          qiime_env
     message:
          "Save the tree.nwk file"
     params:
          "results/{cohort}/{id}+fasttree_rooted"
     shell:
          "qiime tools export "
          "--input-path {input} "
          "--output-path {params} && "
          "cp {params}/tree.nwk {output} && "
          "rm -r {params}"
    

rule make_biom:
     """Create Biom table"""
     input:
          table="results/{cohort}/{id}_table.qza",
          taxonomy="results/{cohort}/{id}+cls-{cls}_taxonomy.qza"
     output:
          "results/{cohort}/{id}+cls-{cls}_asv.biom"
     message:
          "Making biom table {output}"
     conda: 
          qiime_env
     shell:
          "python scripts/make_biom.py --tablef {input.table} "
          "--taxonomy {input.taxonomy} "
          "--output {output}"

rule extract_biom:
     """Auxillary rule to extract biom from Artifact (QZA)"""
     input:
          "results/{cohort}/{cohort}+{id}.qza"
     output:
          "results/{cohort}/{cohort}+{id}.biom"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype biom"

rule extract_sequence:
     """Auxillary rule to exract dada2 ASV sequences from (qza) Artifact to tsv file"""
     input:
          "results/{cohort}/{cohort}+{id}+dd_seq.qza"
     output:
          "results/{cohort}/{cohort}+{id}+dd_seq.tsv"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule extract_stats:
     """Auxillary rule to exract dada2 stats from (qza) Artifact to tsv file"""
     input:
          "results/{cohort}/{cohort}+{id}+dd_stats.qza"
     output:
          "results/{cohort}/{cohort}+{id}+dd_stats.tsv"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule export_phyloseq:
     """Create Phyloseq object RDS file"""
     input:
          biom="results/{cohort}/{id}+cls-{cls}_asv.biom",
          tree="results/{cohort}/{id}+fasttree.nwk"
     output:
          "results/{cohort}/{id}+cls-{cls}+phyloseq.RDS"
     conda:
          "envs/phyloseq.yml"
     shell:
          "Rscript scripts/export_phyloseq.R "
          "--biom {input.biom} --tree {input.tree} "
          "--outp {output}"
          

rule weighted_unifrac:
     """Computes beta diversity weighted unifrac"""
     input:
          table="results/{cohort}/{id}_table+rrf-d{r}.qza",
          tree="results/{cohort}/{id}+fasttree_rooted.qza"
     output:
          "results/{cohort}/{id}+rrf-d{r}+beta_weightedunifrac.qza"
     conda:
          qiime_env
     shell:
          "qiime diversity-lib weighted-unifrac "
          "--i-table {input.table} --i-phylogeny {input.tree} "
          "--o-distance-matrix {output}"

rule unweighted_unifrac:
     """Computes beta diversity weighted unifrac"""
     input:
          table="results/{cohort}/{id}_table+rrf-d{r}.qza",
          tree="results/{cohort}/{id}+fasttree_rooted.qza"
     output:
          "results/{cohort}/{id}+rrf-d{r}+beta_unweightedunifrac.qza"
     conda:
          qiime_env
     shell:
          "qiime diversity-lib unweighted-unifrac "
          "--i-table {input.table} --i-phylogeny {input.tree} "
          "--o-distance-matrix {output}"

rule extract_taxonomy_tsv:
     """converts taxonomy.qza to tsv file"""
     input:
          "results/{cohort}/{id}_taxonomy.qza"
     output:
          "results/{cohort}/{id}_taxonomy.tsv"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule extract_dadatable_tsv:
     """converts dada table.qza to tsv"""
     input:
          "results/{cohort}/{id}_table+rrf-d{r}.qza"           
     output:
          "results/{cohort}/{id}_table+rrf-d{r}.tsv"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule extract_unifrac_tsv:
     """converts unifrac qza artifact to tsv"""
     input:
          "results/{cohort}/{id}unifrac.qza"           
     output:
          "results/{cohort}/{id}unifrac.tsv"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype distance"




rule merge_dadatable:
     """Merge dada table from 2 datasets"""
     input:
          f1="results/{cohort1}/{cohort1}+{id}_table.qza",
          f2="results/{cohort2}/{cohort2}+{id}_table.qza"
     output:
          "results/{cohort1}-{cohort2}/{cohort1}-{cohort2}+{id}_table.qza"
     conda:
          qiime_env
     shell:
          "qiime feature-table merge "
          "--i-tables {input.f1} "
          "--i-tables {input.f2} "
          "--o-merged-table {output}"

rule merge_dadaseq:
     """Merge dada seq file from 2 datasets"""
     input:
          f1="results/{cohort1}/{cohort1}+{id}_seq.qza",
          f2="results/{cohort2}/{cohort2}+{id}_seq.qza"
     output:
          "results/{cohort1}-{cohort2}/{cohort1}-{cohort2}+{id}_seq.qza"
     conda:
          qiime_env
     shell:
          "qiime feature-table merge-seqs "
          "--i-data {input.f1} "
          "--i-data {input.f2} "
          "--o-merged-data {output}"

rule merge_taxonomy:
     priority:
          1
     input:
          f1="results/{cohort1}/{cohort1}+{id}_taxonomy.qza",
          f2="results/{cohort2}/{cohort2}+{id}_taxonomy.qza"
     output:
          "results/{cohort1}-{cohort2}/{cohort1}-{cohort2}+{id}_taxonomy.qza"
     conda:
          qiime_env
     shell:
          "qiime feature-table merge-taxa "
          "--i-data {input.f1} "
          "--i-data {input.f2} "
          "--o-merged-data {output}"

rule collapse_tax:
     """collapse taxonomy table to species level"""
     input:
          table="results/{cohort}/{cohort}+{id}+dd_table+rrf-d{r}.qza",
          tax="results/{cohort}/{cohort}+{id}+dd+{etc}_taxonomy.qza"
     output:
          "results/{cohort}/{cohort}+{id}+dd+{etc}+rrf-d{r}+otu_tax.qza"
     conda:
          qiime_env
     shell:
          "qiime taxa collapse --i-table {input.table} "
          "--p-level 7 "
          "--i-taxonomy {input.tax} "
          "--o-collapsed-table {output}"


rule create_metadata_file:
     input:
          "results/{cohort}/{cohort}_manifest.tsv"
     output:
          "results/{cohort}/{cohort}_metadata.tsv"
     conda:
          "envs/other.yml"
     shell:
          "python scripts/create_metadata_file.py -i {input} -o {output}"
     

rule core_metrics:
     input:
          table="results/{cohort}/{id}_table+rrf-d{r}.qza",
          metadata="results/{cohort}/{cohort}_metadata.tsv"
     output:
          directory("results/{cohort}/{id}+rrf-d{r}+coremetrics/")
     threads:
          30
     conda:
          qiime_env
     shell:
          "qiime diversity core-metrics "
          "--p-sampling-depth {wildcards.r} "
          "--i-table {input.table} "
          "--m-metadata-file {input.metadata} "
          "--output-dir {output} "
          "--p-n-jobs {threads}"

rule alpha_diversity:
     """compoutes alpha diversity"""
     input:
          "results/{cohort}/{id}_table+rrf-d{r}.qza"
     output:
          "results/{cohort}/{id}+rrf-d{r}+alphadiversity.tsv"
     conda:
          qiime_env          
     shell:
          "python scripts/alpha_diversity.py --inp {input} "
          "--outp {output}"

rule beta_diversity:
     """computes non-phylogenetic beta diversity"""
     input:
          "results/{cohort}/{id}+rrf-d{r}+otu_tax.qza"
     output:
          "results/{cohort}/{id}+rrf-d{r}+beta_braycurtis.tsv",
          "results/{cohort}/{id}+rrf-d{r}+beta_jaccard.tsv",
     params:
          "results/{cohort}/{id}+rrf-d{r}+beta.tsv"
     conda:
          qiime_env
     shell:
          "python scripts/beta_diversity.py --inp {input} --outp {params}"

rule biom_to_tsv:
     """converts biom table to tsv"""
     input:
          "results/{cohort}/{id}.biom"
     output:
          "results/{cohort}/{id}_biom.tsv"
     conda:
          qiime_env
     shell:
          "biom convert -i {input} -o {output} --to-tsv"

rule manta:
     """Produces manta output"""
     input:
          tsv="results/{cohort}/{id}+cls-{cls}+rrf-d{r}+otu_tax_biom.tsv",
          taxonpath="db/taxonpath.json",
          names="db/names.json"
     output:
          full="results/{cohort}/{id}+cls-{cls}+rrf-d{r}+manta.tsv",
          tax="results/{cohort}/{id}+cls-{cls}+rrf-d{r}+manta_tax.tsv"
     params:
          db=lambda wildcards: "1" if wildcards.cls=="gg" else "2"
     conda:
          "envs/other.yml"
     shell:
          "python scripts/manta.py "
          "-i {input.tsv} -o {output.full} -x {output.tax} "
          "-t {input.taxonpath} -n {input.names} "
          "-d {params.db} -r {wildcards.r}"

rule summary:
     """Produces summarized results in zipped file"""
     input:
          "results/{cohort}/{id}+cls-{cls}_taxonomy.tsv",
          "results/{cohort}/{id}_table+rrf-d{r}.tsv",
          "results/{cohort}/{id}+rrf-d{r}+beta_weightedunifrac.tsv",
          "results/{cohort}/{id}+rrf-d{r}+beta_unweightedunifrac.tsv",
          "results/{cohort}/{id}+cls-{cls}+rrf-d{r}+otu_tax.biom",
          "results/{cohort}/{id}+cls-{cls}+rrf-d{r}+otu_tax_biom.tsv",
          "results/{cohort}/{id}+cls-{cls}+phyloseq.RDS",
          "results/{cohort}/{id}+rrf-d{r}+alphadiversity.tsv",
          "results/{cohort}/{id}+cls-{cls}+rrf-d{r}+manta.tsv",
          "results/{cohort}/{id}+cls-{cls}+rrf-d{r}+manta_tax.tsv",
          "results/{cohort}/{id}_seq.tsv",
          "results/{cohort}/{id}+cls-{cls}+rrf-d{r}+beta_braycurtis.tsv",
          "results/{cohort}/{id}+cls-{cls}+rrf-d{r}+beta_jaccard.tsv"
     output:
          "results/{cohort}/{id}+cls-{cls}+rrf-d{r}.zip"
     conda:
          "envs/other.yml"
     shell:
          "zip -j {output} {input}"


ruleorder: merge_taxonomy > taxonomy > manifest
ruleorder: merge_dadatable > rarefy > manifest
ruleorder: export_phyloseq  > extract_biom
ruleorder: extract_biom > make_biom