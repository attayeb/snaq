"""
Snakemake pipeline for analyzing 16S RNA data using QIIME2

How to run it:

"""

rule explain:
     input:
          "explain.txt"
     shell:
          "cat {input}"


rule export_artifact:
     input:
          "results/{cohort}/{cohort}_{etc}.qza"
     output:
          directory("temp/{cohort}/{cohort}_{etc}/")
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime tools export "
          "--input-path {input} "
          "--output-path {output}"


def get_allfile_names(wildcards):
     """Get all files from a foler"""
     input_folder = os.path.join("temp", wildcards.cohort, wildcards.cohort+"_"+wildcards.etc)
     return " ".join([os.path.join(input_folder, x) for x in os.listdir(input_folder) if "fastq" in x])


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
          


rule qza_fastqc:
     """  Fastqc
     Input: Fastq files
     Output: fastq report html file.
     Action: Run fastqc quality control analysis
     """
     input:
          "temp/{cohort}/{cohort}_{etc}"          
     output:
          directory("quality/{cohort}/{cohort}_{etc}/fastqc/")
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
     input:
          "quality/{cohort}/{id}/fastqc/"
     output:
          directory("quality/{cohort}/{id}/multiqc/")
     conda:
          "envs/quality.yml"
     shell:
          "multiqc -o {output} {input}"

rule manifest:
     input:
          "data/{cohort}"
     output:
          "results/{cohort}/{cohort}_manifest.tsv"
     conda:
          "envs/other.yml"
     shell:
          "python scripts/create_manifest_file.py -i {input} -o {output}"

rule import_data:
     """  Import data:
     Input: folder with fastq files (pair-ended)
     Output: QIIME2 data sequence artiact.
     Action: Import the raw fastq files to Qiime2 artifact with qza 
     """
     input:
          "results/{cohort}/{cohort}_manifest.tsv" 
     output:
          "results/{cohort}/{cohort}_raw.qza"
     message:
          "Import data"
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime tools import "
	     "--type 'SampleData[PairedEndSequencesWithQuality]' "
	     "--input-path {input} "
	     "--input-format PairedEndFastqManifestPhred33V2 "
	     "--output-path {output} "

rule cutadapt:
     input:
          "results/{cohort}/{cohort}_raw.qza"
     output:
          "results/{cohort}/{cohort}_ca.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     threads:
          20
     shell:
          "qiime cutadapt trim-paired "
          "--i-demultiplexed-sequences {input} "
          "--p-cores {threads} "
          "--o-trimmed-sequences {output}"



rule trim_fastp:
     input:
          qza="results/{cohort}/{cohort}_raw.qza"
     output:
          "results/{cohort}/{cohort}_fp-f{len1}-r{len2}crop.qza"
     message:
          "Trimming using fastp"
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "python scripts/fastp.py --inputf {input} --len1 {wildcards.len1} --len2 {wildcards.len2} --outputf {output}"


rule trim_bbduk:
     input:
          qza="results/{cohort}/{cohort}_{etc}.qza"
     output:
          "results/{cohort}/{cohort}_{etc}_bb{threshold}t.qza"
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
          "results/{cohort}/{cohort}_{etc}.qza"
     output:
          table="results/{cohort}/{cohort}_{etc}_dd-table_rrf0.qza",
	     stats="results/{cohort}/{cohort}_{etc}_dd-stats.qza",
	     repseq="results/{cohort}/{cohort}_{etc}_dd-seq.qza"
     message:
          "Dada2 analysis"
     threads: 30
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"    
     shell:
          "qiime dada2 denoise-paired "
	     "--p-trunc-len-f 0 --p-trunc-len-r 0 "
	     "--i-demultiplexed-seqs {input} "
	     "--o-table {output.table} "
	     "--o-representative-sequences {output.repseq} "
	     "--o-denoising-stats {output.stats} "
	     "--verbose --p-n-threads {threads}"


rule rarefy:
     input:
          "results/{cohort}/{id}-table_rrf0.qza"
     output:
          "results/{cohort}/{id}-table_rrf{r}.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime feature-table rarefy "
          "--i-table {input} "
          "--p-sampling-depth {wildcards.r} "
          "--o-rarefied-table {output}"

rule plot_dada_stats:
     input:
          "results/{cohort}/{id}_dd-stats.qza"
     output:
          "results/{cohort}/plots/{id}_dd-stats.pdf"
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
          repseq = "results/{cohort}/{id}-seq.qza",
	     classifier = "classifiers/{cls}-classifier.qza"
     output:
          taxonomy= "results/{cohort}/{id}_cls-{cls}_taxonomy.qza",
	#table = "results/{resfseq}_{cls}_table.qzv"
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"
     threads: 30
     message:
          "Assign taxonomy using {wildcards.cls} database"
     shell:
          "qiime feature-classifier classify-sklearn "
	     "--i-classifier {input.classifier} "
	     "--i-reads {input.repseq} --p-n-jobs {threads} "
	     "--o-classification {output.taxonomy}"


rule mafft:
     input:
          "results/{cohort}/{id}-seq.qza"
     output:
          "results/{cohort}/tree/{id}-seqaligned.qza"
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
          "results/{cohort}/tree/{id}-seqaligned.qza"
     output:
          "results/{cohort}/tree/{id}-seqaligned-masked.qza"
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
          "results/{cohort}/tree/{id}-seqaligned-masked.qza"
     output:
          "results/{cohort}/tree/{id}_fasttree.qza"
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
          "results/{cohort}/tree/{id}_fasttree.qza"
     output:
          "results/{cohort}/{id}_fasttree_rooted.qza"
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
          "results/{cohort}/{id}_fasttree_rooted.qza"
     output:
          "results/{cohort}/{id}_fasttree.nwk"
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"
     message:
          "Save the tree.nwk file"
     params:
          "results/{cohort}/{id}_fasttree_rooted"
     shell:
          """qiime tools export \
               --input-path {input} \
               --output-path {params} 
          cp {params}/tree.nwk {output}

          rm -r {params}"""


rule make_biom:
     input:
          table="results/{cohort}/{id}-table_rrf{r}.qza",
          taxonomy="results/{cohort}/{id}_cls-{cls}_taxonomy.qza"
     output:
          "results/{cohort}/{id}_cls-{cls}_rrf{r}.biom"
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
          biom="results/{cohort}/{id}_cls-{cls}_rrf{r}.biom",
          tree="results/{cohort}/{id}_fasttree.nwk"
     output:
          "results/{cohort}/{id}_cls-{cls}_rrf{r}_phyloseq.RDS"
     conda:
          "envs/phyloseq.yml"
     shell:
          "Rscript scripts/export_phyloseq.R --biom {input.biom} --tree {input.tree} --outp {output}"
          
rule weighted_unifrac:
     input:
          table="results/{cohort}/{id}-table_rrf{r}.qza",
          tree="results/{cohort}/{id}_fasttree_rooted.qza"
     output:
          "results/{cohort}/{id}_rrf{r}_weightedunifrac.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime diversity-lib weighted-unifrac "
          "--i-table {input.table} --i-phylogeny {input.tree} "
          "--o-distance-matrix {output}"

rule unweighted_unifrac:
     input:
          table="results/{cohort}/{id}-table_rrf{r}.qza",
          tree="results/{cohort}/{id}_fasttree_rooted.qza"
     output:
          "results/{cohort}/{id}_rrf{r}_unweightedunifrac.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime diversity-lib unweighted-unifrac "
          "--i-table {input.table} --i-phylogeny {input.tree} "
          "--o-distance-matrix {output}"

rule extract_taxonomy_csv:
     input:
          "results/{cohort}/{id}_taxonomy.qza"
     output:
          "results/{cohort}/{id}_taxonomy.csv"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule extract_dadatable_csv:
     input:
          "results/{cohort}/{id}-table_rrf{r}.qza"           
     output:
          "results/{cohort}/{id}-table_rrf{r}.csv"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule extract_unifrac_csv:
     input:
          "results/{cohort}/{id}unifrac.qza"           
     output:
          "results/{cohort}/{id}unifrac.csv"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype distance"


rule merge_dadatable:
     input:
          f1="results/{cohort1}/{cohort1}_{id}-table_rrf{r}.qza",
          f2="results/{cohort2}/{cohort2}_{id}-table_rrf{r}.qza"
     output:
          "results/{cohort1}-{cohort2}/{cohort1}-{cohort2}_{id}-table_rrf{r}.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime feature-table merge "
          "--i-tables {input.f1} "
          "--i-tables {input.f2} "
          "--o-merged-table {output}"

rule merge_dadaseq:
     input:
          f1="results/{cohort1}/{cohort1}_{id}-seq.qza",
          f2="results/{cohort2}/{cohort2}_{id}-seq.qza"
     output:
          "results/{cohort1}-{cohort2}/{cohort1}-{cohort2}_{id}-seq.qza"
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
          f1="results/{cohort1}/{cohort1}_{id}_taxonomy.qza",
          f2="results/{cohort2}/{cohort2}_{id}_taxonomy.qza"
     output:
          "results/{cohort1}-{cohort2}/{cohort1}-{cohort2}_{id}_taxonomy.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime feature-table merge-taxa "
          "--i-data {input.f1} "
          "--i-data {input.f2} "
          "--o-merged-data {output}"

rule collapse_tax:
     input:
          table="results/{cohort}/{cohort}_{id}_dd-table_rrf{r}.qza",
          tax="results/{cohort}/{cohort}_{id}_dd_{etc}_taxonomy.qza"
     output:
          "results/{cohort}/{cohort}_{id}_dd_{etc}_rrf{r}_taxonomycollapsed.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
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
          table="results/{cohort}/{id}-table_rrf{r}.qza",
          metadata="results/{cohort}/{cohort}_metadata.tsv"
     output:
          directory("results/{cohort}/{id}_rrf{r}_coremetrics/")
     threads:
          30
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime diversity core-metrics "
          "--p-sampling-depth {wildcards.r} "
          "--i-table {input.table} "
          "--m-metadata-file {input.metadata} "
          "--output-dir {output} "
          "--p-n-jobs {threads}"

rule alpha_diversity:
     input:
          "results/{cohort}/{id}-table_rrf{r}.qza"
     output:
          "results/{cohort}/{id}_rrf{r}_alphadiversity.tsv"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"          
     shell:
          "python scripts/alpha_diversity.py --inp {input} "
          "--outp {output}"

rule summary:
     input:
          "results/{cohort}/{id}_cls-{cls}_taxonomy.csv",
          "results/{cohort}/{id}-table_rrf{r}.csv",
          "results/{cohort}/{id}_rrf{r}_weightedunifrac.csv",
          "results/{cohort}/{id}_rrf{r}_unweightedunifrac.csv",
          "results/{cohort}/{id}_cls-{cls}_rrf{r}.biom",
          "results/{cohort}/{id}_cls-{cls}_rrf{r}_phyloseq.RDS",
          "results/{cohort}/{id}_rrf{r}_alphadiversity.tsv"


     output:
          "results/{cohort}/{id}_cls-{cls}_rrf{r}.zip"
     shell:
          "zip -j {output} {input}"

ruleorder: trim_bbduk > trim_fastp
ruleorder: taxonomy > merge_taxonomy 
ruleorder: rarefy > merge_dadatable