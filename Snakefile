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
          "qza/{cohort}/{cohort}_ca.qza"
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
          qza="qza/{cohort}/{cohort}_{etc}.qza"
     output:
          "qza/{cohort}/{cohort}_{etc}_bb{threshold}t.qza"
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
          table="qza/{cohort}/{cohort}_{etc}_dd-f{f}-r{r}-table_rrf0.qza",
	     stats="qza/{cohort}/{cohort}_{etc}_dd-f{f}-r{r}-stats.qza",
	     ref_seq="qza/{cohort}/{cohort}_{etc}_dd-f{f}-r{r}-seq.qza"
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
          "qza/{cohort}/{id}-table_rrf0.qza"
     output:
          "qza/{cohort}/{id}-table_rrf{r}.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime feature-table rarefy "
          "--i-table {input} "
          "--p-sampling-depth {wildcards.r} "
          "--o-rarefied-table {output}"

rule plot_dada_stats:
     input:
          "qza/{cohort}/{id}-stats.qza"
     output:
          "qza/{cohort}/plots/{id}-stats.pdf"
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
          refseq = "qza/{cohort}/{id}-seq.qza",
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
          "qza/{cohort}/{id}_seq.qza"
     output:
          "qza/{cohort}/tree/{id}_seqaligned.qza"
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
          "qza/{cohort}/tree/{id}_seqaligned.qza"
     output:
          "qza/{cohort}/tree/{id}_seqaligned_masked.qza"
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
          "qza/{cohort}/tree/{id}_seqaligned_masked.qza"
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
          table="qza/{cohort}/{id}-table-rrf{r}.qza",
          taxonomy="qza/{cohort}/{id}_{cls}_taxonomy.qza"
     output:
          "qza/{cohort}/{id}_{cls}_rrf{r}.biom"
     message:
          "Making biom table {output}"
     conda: 
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "python scripts/make_biom.py --tablef {input.table} "
          "--taxonomy {input.taxonomy} "
          "--output {output}"

rule export_artifact:
     input:
          "qza/{cohort}/{id}.qza"
     output:
          "temp/{cohort}/{id}/"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime tools extract "
          "--input-path {input} "
          "--output-path {output}"


rule export_phyloseq:
     input:
          biom="qza/{cohort}/{id}_{cls}_rrf{r}.biom",
          tree="qza/{cohort}/{id}_fasttree.nwk"
     output:
          "qza/{cohort}/{id}_{cls}_{r}r_phyloseq.RDS"
     conda:
          "envs/phyloseq.yml"
     shell:
          "Rscript scripts/export_phyloseq.R --biom {input.biom} --tree {input.tree} --outp {output}"
          
rule weighted_unifrac:
     input:
          table="qza/{cohort}/{id}-table_rrf{r}.qza",
          tree="qza/{cohort}/{id}_fasttree_rooted.qza"
     output:
          "qza/{cohort}/{id}_rrf{r}_weightedunifrac.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime diversity-lib weighted-unifrac "
          "--i-table {input.table} --i-phylogeny {input.tree} "
          "--o-distance-matrix {output}"

rule unweighted_unifrac:
     input:
          table="qza/{cohort}/{id}-table_rrf{r}.qza",
          tree="qza/{cohort}/{id}_fasttree_rooted.qza"
     output:
          "qza/{cohort}/{id}_rrf{r}_unweightedunifrac.qza"
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
          "qza/{cohort}/{id}-table_rrf{r}.qza"           
     output:
          "qza/{cohort}/{id}-table_rrf{r}.csv"
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
          f1="qza/{cohort1}/{cohort1}_{id}-table_rrf{r}.qza",
          f2="qza/{cohort2}/{cohort2}_{id}-table_rrf{r}.qza"
     output:
          "qza/{cohort1}-{cohort2}/{cohort1}-{cohort2}_{id}_table-rrf{r}.qza"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"
     shell:
          "qiime feature-table merge "
          "--i-tables {input.f1} "
          "--i-tables {input.f2} "
          "--o-merged-table {output}"

rule merge_dadaseq:
     input:
          f1="qza/{cohort1}/{cohort1}_{id}-seq.qza",
          f2="qza/{cohort2}/{cohort2}_{id}-seq.qza"
     output:
          "qza/{cohort1}-{cohort2}/{cohort1}-{cohort2}_{id}-seq.qza"
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
          "qza/{cohort}/{id}-table_rrf{r}.qza"
     output:
          "qza/{cohort}/{id}_ttf{r}_alphadiversity.tsv"
     conda:
          "envs/qiime2-latest-py38-linux-conda.yml"          
     shell:
          "python scripts/alpha_diversity.py --inp {input} "
          "--outp {output}"

rule manta:
     input:
          taxonomy="qza/{cohort}/{id}_{cls}_taxonomy.csv",
          abundancy="qza/{cohort}/{id}_dadatable_rrf{r}.csv",
          wunifrac="qza/{cohort}/{id}_rrf{r}_weightedunifrac.csv",
          uwunifrac="qza/{cohort}/{id}_rrf{r}_unweightedunifrac.csv"

     output:
          "qza/{cohort}/manta/{id}_{cls}_rrf{r}.zip"
     shell:
          "zip -j {output} {input}"


ruleorder: taxonomy > merge_taxonomy 
ruleorder: rarefy > merge_dadatable