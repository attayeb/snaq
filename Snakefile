"""
Snakemake pipeline for analyzing 16S RNA data using QIIME2

How to run it:

"""

rule explain:
     input:
          "explain.txt"
     log:
          "a.log"
     shell:
          "cat {input}"


def get_allfile_names(wildcards):
     """Get all files from a foler"""
     input_folder = os.path.join("data", wildcards.cohort)
     return [os.path.join(input_folder, x) for x in os.listdir(input_folder)]

rule print_help:
     """  Print help
     Input: None
     Output: None
     Action: print out the information of rules on screen
     """
     run:
          from termcolor import cprint
          for rule in workflow.rules:
               try:
                    if not rule.name.startswith('_'):
                         cprint(rule.name, "red", attrs=["bold"])
                         cprint(rule.docstring, "green")
               except:
                    pass
          
          
          


rule fastqc:
     """  Fastqc
     Input: Fastq files
     Output: fastq report html file.
     Action: Run fastqc quality control analysis
     """
     input:
          get_allfile_names
     output:
          directory("quality/{cohort}/fastqc/")
     threads:
          20
     conda:
          "envs/quality.yml"
     shell:
          "mkdir {output} && "
          "fastqc -o {output} -f fastq -t {threads} {input}"

rule multiqc:
     input:
          "quality/{cohort}/fastqc/"
     output:
          directory("quality/{cohort}/multiqc/")
     conda:
          "envs/quality.yml"
     shell:
          "multiqc -o {output} {input}"

rule import_data:
     """  Import data:
     Input: folder with fastq files (pair-ended)
     Output: QIIME2 data sequence artiact.
     Action: Import the raw fastq files to Qiime2 artifact with qza 
     """
     input:
          "data/{cohort}" 
     output:
          "qza/{cohort}/{cohort}_raw.qza"
     message:
          "Import data"
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime tools import "
	     "--type 'SampleData[PairedEndSequencesWithQuality]' "
	     "--input-path {input} "
	     "--input-format CasavaOneEightSingleLanePerSampleDirFmt "
	     "--output-path {output} "

rule cutadapt:
     input:
          "qza/{cohort}/{cohort}_raw.qza"
     output:
          "qza/{cohort}/{cohort}_cutadapt.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     threads:
          20
     shell:
          "qiime cutadapt trim-paired "
          "--i-demultiplexed-sequences {input} "
          "--p-cores {threads} "
          "--o-trimmed-sequences {output}"



rule trim_bbduk:
     input:
          qza="qza/{cohort}/{cohort}_raw.qza"
     output:
          "qza/{cohort}/{cohort}_bb{threshold}t.qza"
     message:
          "Trimming using bbduk"
     params:
          lambda wildcards, output: output[0].replace(".qza", ".log")
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "python scripts/bbduk.py -i {input} -q {wildcards.threshold} -o {output} > {params}"


rule dada2:
     input:
          "qza/{cohort}/{cohort}_{etc}.qza"
     output:
          table="qza/{cohort}/{cohort}_{etc}_f{f}_r{r}_dadatable_0r.qza",
	     stats="qza/{cohort}/{cohort}_{etc}_f{f}_r{r}_dadastats.qza",
	     ref_seq="qza/{cohort}/{cohort}_{etc}_f{f}_r{r}_dadaseq.qza"
     message:
          "Dada2 analysis"
     threads: 30
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"    
     shell:
          "qiime dada2 denoise-paired "
	     "--p-trim-left-f {wildcards.f} --p-trim-left-r {wildcards.r} "
	     "--p-trunc-len-f 0 --p-trunc-len-r 0 "
	     "--i-demultiplexed-seqs {input} "
	     "--o-table {output.table} "
	     "--o-representative-sequences {output.ref_seq} "
	     "--o-denoising-stats {output.stats} "
	     "--verbose --p-n-threads {threads}"

rule rarefy:
     input:
          "qza/{cohort}/{id}_dadatable_0r.qza"
     output:
          "qza/{cohort}/{id}_dadatable_{r}r.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime feature-table rarefy "
          "--i-table {input} "
          "--p-sampling-depth {wildcards.r} "
          "--o-rarefied-table {output}"

rule plot_dada_stats:
     input:
          "qza/{cohort}/{id}_dadastats.qza"
     output:
          "qza/{cohort}/plots/{id}_dadastats.pdf"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "python scripts/plot_dada.py --inp {input} --plot {output}"



rule download_silva_classifier:
     output:
          "classifiers/silva-classifier.qza"
     shell:
          "cd classifiers && wget https://zenodo.org/record/5535616/files/silva-classifier.qza"

rule download_gg_classifier:
     output:
          "classifiers/gg-classifier.qza"
     shell:
          "cd classifiers && "
          "wget https://zenodo.org/record/5535616/files/gg-classifier.qza"

rule download_silvav34_classifier:
     output:
          "classifiers/silvaV34-classifier.qza"
     shell:
          "cd classifiers && "
          "wget https://zenodo.org/record/5535616/files/silvaV34-classifier.qza"


rule taxonomy:
     input:
          refseq = "qza/{cohort}/{id}_dadaseq.qza",
	     classifier = "classifiers/{cls}-classifier.qza"
     output:
          taxonomy= "qza/{cohort}/{id}_{cls}_taxonomy.qza",
	#table = "qza/{resfseq}_{cls}_table.qzv"
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"
     threads: 30
     message:
          "Assign taxonomy using {wildcards.cls} database"
     shell:
          "qiime feature-classifier classify-sklearn "
	     "--i-classifier {input.classifier} "
	     "--i-reads {input.refseq} --p-n-jobs {threads} "
	     "--o-classification {output.taxonomy}"

#tree


rule mafft:
     input:
          "qza/{cohort}/{id}_dadaseq.qza"
     output:
          "qza/{cohort}/tree/{id}_dadaseqaligned.qza"
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"
     message:
          "Mafft alignment"
     shell:
          "qiime alignment mafft "
          "--i-sequences {input} "
          "--o-alignment {output}"

rule alignment_mask:
     input:
          "qza/{cohort}/tree/{id}_dadaseqaligned.qza"
     output:
          "qza/{cohort}/tree/{id}_dadaseqaligned_masked.qza"
     message:
          "Alignment mask"
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"     
     shell:
          "qiime alignment mask "
          "--i-alignment {input} "
          "--o-masked-alignment {output}"

rule fasttree:
     input:
          "qza/{cohort}/tree/{id}_dadaseqaligned_masked.qza"
     output:
          "qza/{cohort}/tree/{id}_fasttree.qza"
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"
     message:
          "creating tree"
     shell:
          "qiime phylogeny fasttree "
          "--i-alignment {input} "
          "--o-tree {output}"
          
          

rule midpoint_root:
     input:
          "qza/{cohort}/tree/{id}_fasttree.qza"
     output:
          "qza/{cohort}/{id}_fasttree_rooted.qza"
     message:
          "Midpoint rooting the tree"
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime phylogeny midpoint-root "
          "--i-tree {input} "
          "--o-rooted-tree {output}"


rule export_tree:
     input:
          "qza/{cohort}/{id}_fasttree_rooted.qza"
     output:
          "qza/{cohort}/{id}_fasttree.nwk"
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"
     message:
          "Save the tree.nwk file"
     params:
          "qza/{cohort}/{id}_fasttree_rooted"
     shell:
          """qiime tools export \
               --input-path {input} \
               --output-path {params} 
          cp {params}/tree.nwk {output}

          rm -r {params}"""


rule make_biom:
     input:
          table="qza/{cohort}/{id}_dadatable_{r}r.qza",
          taxonomy="qza/{cohort}/{id}_{cls}_taxonomy.qza"
     output:
          "qza/{cohort}/{id}_{cls}_{r}r.biom"
     message:
          "Making biom table {output}"
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "python scripts/make_biom.py --tablef {input.table} "
          "--taxonomy {input.taxonomy} "
          "--output {output}"

rule export_phyloseq:
     input:
          biom="qza/{cohort}/{id}_{cls}_{r}r.biom",
          tree="qza/{cohort}/{id}_fasttree.nwk"
     output:
          "qza/{cohort}/{id}_{cls}_{r}r_phyloseq.RDS"
     conda:
          "envs/phyloseq.yml"
     shell:
          "Rscript scripts/export_phyloseq.R --biom {input.biom} --tree {input.tree} --outp {output}"
          
rule weighted_unifrac:
     input:
          table="qza/{cohort}/{id}_dadatable_{r}r.qza",
          tree="qza/{cohort}/{id}_fasttree_rooted.qza"
     output:
          "qza/{cohort}/{id}_{r}r_weightedunifrac.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime diversity-lib weighted-unifrac "
          "--i-table {input.table} --i-phylogeny {input.tree} "
          "--o-distance-matrix {output}"

rule unweighted_unifrac:
     input:
          table="qza/{cohort}/{id}_dadatable_{r}r.qza",
          tree="qza/{cohort}/{id}_fasttree_rooted.qza"
     output:
          "qza/{cohort}/{id}_{r}r_unweightedunifrac.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime diversity-lib unweighted-unifrac "
          "--i-table {input.table} --i-phylogeny {input.tree} "
          "--o-distance-matrix {output}"

rule extract_taxonomy_csv:
     input:
          "qza/{cohort}/{id}_taxonomy.qza"
     output:
          "qza/{cohort}/{id}_taxonomy.csv"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule extract_dadatable_csv:
     input:
          "qza/{cohort}/{id}_dadatable_{r}r.qza"           
     output:
          "qza/{cohort}/{id}_dadatable_{r}r.csv"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule extract_unifrac_csv:
     input:
          "qza/{cohort}/{id}unifrac.qza"           
     output:
          "qza/{cohort}/{id}unifrac.csv"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype distance"


rule merge_dadatable:
     input:
          f1="qza/{cohort1}/{cohort1}_{id}_dadatable_{r}.qza",
          f2="qza/{cohort2}/{cohort2}_{id}_dadatable_{r}.qza"
     output:
          "qza/{cohort1}-{cohort2}/{cohort1}-{cohort2}_{id}_dadatable_{r}.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime feature-table merge "
          "--i-tables {input.f1} "
          "--i-tables {input.f2} "
          "--o-merged-table {output}"

rule merge_dadaseq:
     input:
          f1="qza/{cohort1}/{cohort1}_{id}_dadaseq.qza",
          f2="qza/{cohort2}/{cohort2}_{id}_dadaseq.qza"
     output:
          "qza/{cohort1}-{cohort2}/{cohort1}-{cohort2}_{id}_dadaseq.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime feature-table merge-seqs "
          "--i-data {input.f1} "
          "--i-data {input.f2} "
          "--o-merged-data {output}"

rule merge_taxonomy:
     priority:
          1
     input:
          f1="qza/{cohort1}/{cohort1}_{id}_taxonomy.qza",
          f2="qza/{cohort2}/{cohort2}_{id}_taxonomy.qza"
     output:
          "qza/{cohort1}-{cohort2}/{cohort1}-{cohort2}_{id}_taxonomy.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime feature-table merge-taxa "
          "--i-data {input.f1} "
          "--i-data {input.f2} "
          "--o-merged-data {output}"

rule alpha_diversity:
     input:
          "qza/{cohort}/{id}_dadatable_{r}r.qza"
     output:
          "qza/{cohort}/{id}_{r}r_alphadiversity.tsv"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"          
     shell:
          "python scripts/alpha_diversity.py --inp {input} "
          "--outp {output}"

rule manta:
     input:
          taxonomy="qza/{cohort}/{id}_{cls}_taxonomy.csv",
          abundancy="qza/{cohort}/{id}_dadatable_{r}r.csv",
          wunifrac="qza/{cohort}/{id}_{r}r_weightedunifrac.csv",
          uwunifrac="qza/{cohort}/{id}_{r}r_unweightedunifrac.csv"

     output:
          "qza/{cohort}/manta/{id}_{cls}_{r}r.zip"
     shell:
          "zip -j {output} {input}"


def get_kraken2_files(wildcards):
     files = os.listdir("data/{}".format(wildcards.cohort))
     files = [file.split("_")[0] for file in files if "_1.fastq.gz" in file]
     return files


KRAKEN = "/data/microbiome/kraken2/kraken2"


rule kraken2:
     input:
          r1 = "data/{cohort}/{id}_L001_R1_001.fastq.gz",
          r2 = "data/{cohort}/{id}_L001_R2_001.fastq.gz"
     output:
          report = "kmer/{cohort}/kraken/{id}_{db}_kraken.rep",
          output = "kmer/{cohort}/kraken/{id}_{db}_kraken.outp"
     shell:
          "{KRAKEN} --db /data/microbiome/kraken2/{wildcards.db} "
          "--threads 30 "
          "--output {output.output} "
          "--report {output.report} "
          "{input.r1} {input.r2}"

BRACKEN = "/data/microbiome/bracken/bracken"


rule bracken:
     input:
          "kmer/{cohort}/kraken/{id}_{db}_kraken.rep"
     output:
          report = "kmer/{cohort}/bracken/{id}_{db}_{l}_bracken.rep",
          output = "kmer/{cohort}/bracken/{id}_{db}_{l}_bracken.outp"
     shell:
          "{BRACKEN} -d /data/microbiome/kraken2/{wildcards.db} "
          "-r 300 -i {input} -l {wildcards.l} -o {output.output} -w {output.report}"

def get_kraken_files(wildcards):
     files = os.listdir("data/{}/".format(wildcards.cohort))
     files = [file.split("_L001_")[0] for file in files if "_R1_" in file]
     ret = [os.path.join("kmer/{}/kraken/".format(wildcards.cohort), "{}_{}_kraken.rep").format(file, wildcards.db) for file in files]
     
     return ret

rule merge_kraken:
     input:
          get_kraken_files
     output:
          "kmer/{cohort}/summarized/{cohort}_{db}_kraken.json"
     params:
          inputfiles=lambda wildcards, input: ",".join(input)
     shell:
          "python scripts/summarize.py --files {params.inputfiles} --outp {output}"



def get_bracken_files(wildcards):
     files = os.listdir("data/{}/".format(wildcards.cohort))
     files = [file.split("_L001_")[0] for file in files if "_R1_" in file]
     ret = [os.path.join("kmer/{}/bracken/".format(wildcards.cohort), "{}_{}_{}_bracken.rep").format(file, wildcards.db, wildcards.l) for file in files]
     
     return ret




rule merge_bracken:
     input:
          files=get_bracken_files,
          merged_kraken="kmer/{cohort}/summarized/{cohort}_{db}_kraken.json"
     output:
          "kmer/{cohort}/summarized/{cohort}_{db}_{l}_bracken.json"
     params:
          inputfiles=lambda wildcards, input: ",".join(input.files)
     shell:
          "python scripts/summarize.py --files {params.inputfiles} --outp {output}"


rule make_cami:
     input:
          "kmer/{cohort}/bracken/{id}_{db}_{l}_bracken.rep"
     output:
          "kmer/{cohort}/cami/{id}_{db}_{l}.cami"
     shell:
          "python scripts/kraken2cami.py -i {input} -o {output} -d db/ -s {wildcards.id}"


ruleorder: taxonomy > merge_taxonomy 
ruleorder: rarefy > merge_dadatable