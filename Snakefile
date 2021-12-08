"""
Snakemake pipeline for analyzing 16S RNA data using QIIME2

Usage
-----

snakemake 

"""
from platform import system

_os = system()

if _os == "Linux":
     qiime_env = "envs/qiime2-2021.8-py38-linux-conda.yml"
elif _os == "Darwin":
     qiime_env = "envs/qiime2-2021.8-py38-osx-conda.yml"



rule explain:
     """
     
     """
     input:
          "explain.txt"
     shell:
          "cat {input}"

rule export_artifact_2:
     """Export Artifact content to a folder
       Input
       -----
          QIIME2 Artifact file (qza)
     
       Output
       ------
          Directory contains the content of the Artifact
     
       Tools
       -----
          export command from qiime tools
     
       Actions
       -------
          Export the artifact content to ouput folder"""
     message:
          "Extracting artifact"
     input:
          "results/{cohort}/{cohort}.qza"
     output:
          directory("temp/{cohort}/{cohort}/")
     conda:
          qiime_env
     shell:
          "qiime tools export "
          "--input-path {input} "
          "--output-path {output}"



rule export_artifact:
     """Export Artifact content to a folder
       Input
       -----
          QIIME2 Artifact file (qza)
     
       Output
       ------
          Directory contains the content of the Artifact
     
       Tools
       -----
          export command from qiime tools
     
       Actions
       -------
          Export the artifact content to ouput folder"""
     message:
          "Extracting artifact"
     input:
          "results/{cohort}/{cohort}+{etc}.qza"
     output:
          directory("temp/{cohort}/{cohort}+{etc}/")
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
     """  Print help
     
       Action 
       ------
          Print out the information of rules on screen (docstrings)."""
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
     Input
     -----
          Artifact of type  
          Fastq files
     Output: fastq report html file.
     Action: Run fastqc quality control analysis
     """
     message:
          "Applying FASTQC"
     input:
          "temp/{cohort}/{id}"          
     output:
          directory("quality/{cohort}/{id}/fastqc/")
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
     message:
          "MultiQC"
     input:
          "quality/{cohort}/{id}/fastqc/"
     output:
          directory("quality/{cohort}/{id}/multiqc/")
     conda:
          "envs/quality.yml"
     shell:
          "multiqc -o {output} {input}"

rule qza_cohort_multiqc:
     message:
          "MultiQC"
     input:
          "quality/{cohort}/fastqc/"
     output:
          directory("quality/{cohort}/multiqc/")
     conda:
          "envs/quality.yml"
     shell:
          "multiqc -o {output} {input}"

rule manifest:
     """Create manifest file
       Input
       -----
          folder name that contain the data (data/cohort/)
       
       Output
       ------
          manifest file (results/cohort/cohort_manifest.tsv)
       
       Tools
       -----
          utilizes scripts/create_manifest_file.py script
     """
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
     """Import data:
       input
       -----
          manifest file (pair-ended) [tsv].
       Output 
       ------
          QIIME2 data sequence artifact.
       tools
       -----
          utilizes "qiime tools import" command.
       action
       ------
          Import the raw fastq files to Qiime2 artifact with qza extension
     """
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
     """Crop primers:
       input
       -----
          ARTIFACT SampleData[PairedEndSequencesWithQuality] 
       Output 
       ------
          ARTIFACT SampleData[PairedEndSequencesWithQuality]
       tools
       -----
          utilizes "scripts/fastp.py" script which uses is a wrapper of fastp, takes input as a qiime2 artifact and produce qiime2 artifact output
       parameters
       ----------
          takes forward cropping length and backward cropping length
       action
       ------
          Crop primers from both end reads.
       usage
       -----
          fp-f17-r21crop : trims 17 bases from forward read and 21 bases form backward read.
     
          ex:
          snakemake --cores 10 --use-conda results/CH/CH_fp-f17-r21.qza
     """
     input:
          qza="results/{cohort, [A-Z]}/{id}.qza"
     output:
          "results/{cohort}/{id}+fp-f{len1}-r{len2}crop.qza"
     message:
          "Trimming using fastp"
     conda: 
          qiime_env
     shell:
          "python scripts/fastp.py --inputf {input} --len1 {wildcards.len1} --len2 {wildcards.len2} --outputf {output}"


rule trim_bbduk:
     """Quality trimming using bbduk:
       input
       -----
          original seq data [ARTIFACT SampleData[PairedEndSequencesWithQuality]]
       Output 
       ------
          quality trimmed seq data [ARTIFACT SampleData[PairedEndSequencesWithQuality]]
       tools
       -----
          utilizes "scripts/bbduk.py" script which uses is a wrapper of bbduk.sh, takes input as a qiime2 artifact and produce qiime2 artifact output
       parameters
       ----------
          takes quality trimming threshold as bb{threshold}t
       action
       ------
          quality trimming of fastq files (pair-ended).
       usage
       -----
          bb16t : trims the distal part of both reads (R1 and R2) when the average quality score of the distal segment is less than 16.
     
          ex:
          quality trimming of the raw files:
          - snakemake --cores 10 --use-conda results/CH/CH_bb16t.qza
          quality trimming of the primers croped data:
          - snakemake --cores 10 --use-conda results/CH/CH+fp-17f-r21+bb16t.qza
     """
     input:
          qza="results/{cohort, [A-Z]}/{id}.qza"
     output:
          "results/{cohort}/{id}+bb{threshold}t.qza"
     message:
          "Trimming using bbduk"
     params:
          lambda wildcards, output: output[0].replace(".qza", ".log")
     conda: 
          qiime_env
     shell:
          "python scripts/bbduk.py -i {input} -q {wildcards.threshold} -o {output} > {params}"


rule dada2:
     """Dada2 algorithm
       input
       -----
          ARTIFACT SampleData[PairedEndSequencesWithQuality] 
       Output 
       ------
          {dd-table_rrf0.qza} ARTIFACT FeatureTable[Frequency]
          {dd-seq.qza} ARTIFACT FeatureData[Sequence]
          {dd-stats.qza} ARTIFACT SampleData[DADA2Stats]
       tools
       -----
          utilizes "qiime dada2 denoise-paired" command
       action
       ------
          denoises paired-end sequences, dereplicates them, and filters chimeras.
       usage
       -----
          you can specify any output among the three output to get them all.
          ex: dd-seq.qza
           
          ex:
          dada2 of raw data without primer cropping or quality trimming:
          - sankemake --cores 10 --use-conda results/CH/CH_dd-table_rrf0.qza
          pair-end denoising of the data after cropping primers:
          - snakemake --cores 10 --use-conda results/CH/CH_bb16t_dd-seq.qza
          dada2 of data after quality trimming of the raw files after primers cropping:
          - snakemake --cores 10 --use-conda results/CH/CH_fp-17f-r21_bb16t_dd-stats.qza
        
        note
        ----
          rrf0 was added to distinguish the original input before rarefaction from rarefied table files.

     """ 
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
     """Dada2 table file rarefication
       input
       -----
          ARTIFACT FeatureTable[Frequency]
       Output 
       ------
          ARTIFACT FeatureTable[Frequency]
       tools
       -----
          utilizes "qiime feature-table rarefy" command
       action
       ------
          rarefy feature table.
       usage
       -----
           
          ex:
          rarefy the table to 10000
          - sankemake --cores 10 --use-conda results/CH/CH_dd-table_rrf10000.qza
          rarefy the table to 1000
          - snakemake --cores 10 --use-conda results/CH/CH_bb16t_dd-table_rrf1000.qza
     """
     input:
          "results/{cohort}/{id}_table.qza"
     output:
          "results/{cohort}/{id}_table+rrf{r}.qza"
     conda:
          qiime_env
     shell:
          "qiime feature-table rarefy "
          "--i-table {input} "
          "--p-sampling-depth {wildcards.r} "
          "--o-rarefied-table {output}"

rule plot_dada_stats:
     input:
          "results/{cohort}/{id}+dd_stats.qza"
     output:
          "results/{cohort}/plots/{id}+dd_stats.pdf"
     conda:
          qiime_env
     shell:
          "python scripts/plot_dada.py --inp {input} --plot {output}"


rule download_silva_classifier:
     """ Download pretrained SILVA taxonomy classifier
     output
     ------
          ARTIFACT TaxonomicClassifier
     action
     ------
          Download pretrained SILVA classifier.
     """
     output:
          "classifiers/silva-classifier.qza"
     conda:
          "envs/other.yml"
     shell:
          "cd classifiers && wget https://zenodo.org/record/5535616/files/silva-classifier.qza"

rule download_gg_classifier:
     """ Download GreenGenes taxonomy classifier
     output
     ------
          ARTIFACT TaxonomicClassifier
     action
     ------
          Download pretrained GreenGenes classifier.
     """
     output:
          "classifiers/gg-classifier.qza"
     conda:
          "envs/other.yml"
     shell:
          "cd classifiers && "
          "wget https://zenodo.org/record/5535616/files/gg-classifier.qza"

rule download_silvav34_classifier:
     """ Download V3-V4 region pretrained SILVA classifier.
     output
     ------
          ARTIFACT TaxonomicClassifier
     action
     ------
          Download pretrained SILVA V3-V4 region classifier.
     """
     output:
          "classifiers/silvaV34-classifier.qza"
     conda:
          "envs/other.yml"
     shell:
          "cd classifiers && "
          "wget https://zenodo.org/record/5535616/files/silvaV34-classifier.qza"


rule taxonomy:
     """Dada2 table file rarefication
       input
       -----
          reads = ARTIFACT FeatureTable[Frequency]
          classifier = ARTIFACT TaxonomicClassifier
       Output 
       ------
          taxonomy = ARTIFACT FeatureData[Taxonomy]
       tools
       -----
          utilizes "qiime feature-classifier classify-sklearn" command
       action
       ------
          Assign taxonomy to ASVs.
       usage
       -----
          classifier need ot be assigned after _cls-{classifier}
          available options:
               - GreenGenes: _cls-gg-taxonomy.qza
               - SILVA: _cls-silva-taxonomy.qza
               _ SILVA V3-V4: _cls-silvaV34-taxonomy.qza
          Because it is will understood which output form dada2 is going to be used, no need to add -seq after dd
          
          ex:
          - sankemake --cores 10 --use-conda results/CH/CH_dd_cls-gg-taxonomy.qza
     """
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
     """ Creating phylogeny tree step 1 
     """
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
     """ Creating phylogeny tree step 2
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
     """ Creating phylogeny tree step 3 
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
     """ Creating phylogeny tree step 4
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
     """ Creating phylogeny tree step 5:
          Exporting tree from Artifact file
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
          """qiime tools export \
               --input-path {input} \
               --output-path {params} 
          cp {params}/tree.nwk {output}

          rm -r {params}"""
    

rule make_biom:
     """Create Biom table
       Input
       -----
          - Abdundance table (ARTIFACT FeatureTable[Frequency])
          - Taxonomy table (ARTIFACT FeatureData[Taxonomy])
       Output
       ------
          biom Table
       Tools
       -----
          utilizes "scripts/make_biom.py" which merge the table frequency with the taxonomy
       Usage
       -----
          This rule produce biom table from FeatureTable[Frequency] without taxonomy collapsing. The rows of the presented biom table are the original ASVs

          ex:
          To produce the analysis of CH cohort with:
               * primer cropping at f:17 and r:21, 
               * quality trimming with 14 thereshold,
               * dada 2 denoising
               * using SILVA classifier for taxonomy assignment
               * rarefied at 1000
          snakemake --cores 10 --use-conda results/CH/CH_fp-17f-21rcrop_bb14t_dd_cls-silva_rrf10000.biom

          To produce the analysis of CH cohort with:
               * primer cropping at f:17 and r:21, 
               * quality trimming with 14 thereshold,
               * dada 2 denoising
               * using GreenGene classifier for taxonomy assignment
               * rarefied at 1000
          snakemake --cores 10 --use-conda results/CH/CH_fp-17f-21rcrop_bb14t_dd_cls-gg_rrf10000.biom
     """
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
     """Auxillary rule to exract dada2 ASV sequences from (qza) Artifact to csv file"""
     input:
          "results/{cohort}/{cohort}+{id}+dd_seq.qza"
     output:
          "results/{cohort}/{cohort}+{id}+dd_seq.csv"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule export_phyloseq:
     """Create Phyloseq object RDS file
       Input
       -----
          Biom table without rarefication [.biom file]
          Phylogeny tree [.nwk file]
       Output
       ------
          Phyloseq object saved in RDS file
       Tools
       -----
          utilizes "scripts/export_phyloseq.R" which is wriitn in R
       Usage
       -----
          ex:
          snakemake --core 10 --use-conda results/CH/CH_fp-f12-r22crop_bb12t_dd_cls-gg_rrf0_phyloseq.RDS
       Note
       ----
          It is not possible to create Phyloseq Object from data after taxonomy collapse because the tree is available only for ASV level,
          rarefication is possible in R after uploading the object to R
     """
     input:
          biom="results/{cohort}/{id}+cls-{cls}_asv.biom",
          tree="results/{cohort}/{id}+fasttree.nwk"
     output:
          "results/{cohort}/{id}+cls-{cls}+phyloseq.RDS"
     conda:
          "envs/phyloseq.yml"
     shell:
          "Rscript scripts/export_phyloseq.R --biom {input.biom} --tree {input.tree} --outp {output}"
          

rule weighted_unifrac:
     input:
          table="results/{cohort}/{id}_table+rrf{r}.qza",
          tree="results/{cohort}/{id}+fasttree_rooted.qza"
     output:
          "results/{cohort}/{id}+rrf{r}+beta_weightedunifrac.qza"
     conda:
          qiime_env
     shell:
          "qiime diversity-lib weighted-unifrac "
          "--i-table {input.table} --i-phylogeny {input.tree} "
          "--o-distance-matrix {output}"

rule unweighted_unifrac:
     input:
          table="results/{cohort}/{id}_table+rrf{r}.qza",
          tree="results/{cohort}/{id}+fasttree_rooted.qza"
     output:
          "results/{cohort}/{id}+rrf{r}+beta_unweightedunifrac.qza"
     conda:
          qiime_env
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
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule extract_dadatable_csv:
     input:
          "results/{cohort}/{id}_table+rrf{r}.qza"           
     output:
          "results/{cohort}/{id}_table+rrf{r}.csv"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype metadata"

rule extract_unifrac_csv:
     input:
          "results/{cohort}/{id}unifrac.qza"           
     output:
          "results/{cohort}/{id}unifrac.csv"
     conda:
          qiime_env
     shell:
          "python scripts/artifact_view.py --artifact {input} "
          "--filename {output} --filetype distance"




rule merge_dadatable:
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
     input:
          table="results/{cohort}/{cohort}+{id}+dd_table+rrf{r}.qza",
          tax="results/{cohort}/{cohort}+{id}+dd+{etc}_taxonomy.qza"
     output:
          "results/{cohort}/{cohort}+{id}+dd+{etc}+rrf{r}+otu_tax.qza"
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
          table="results/{cohort}/{id}_table+rrf{r}.qza",
          metadata="results/{cohort}/{cohort}_metadata.tsv"
     output:
          directory("results/{cohort}/{id}+rrf{r}+coremetrics/")
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
     input:
          "results/{cohort}/{id}_table+rrf{r}.qza"
     output:
          "results/{cohort}/{id}+rrf{r}+alphadiversity.tsv"
     conda:
          qiime_env          
     shell:
          "python scripts/alpha_diversity.py --inp {input} "
          "--outp {output}"

rule beta_diversity:
     input:
          "results/{cohort}/{id}+rrf{r}+otu_tax.qza"
     output:
          "results/{cohort}/{id}+rrf{r}+beta_braycurtis.tsv",
          "results/{cohort}/{id}+rrf{r}+beta_jaccard.tsv",
     params:
          "results/{cohort}/{id}+rrf{r}+beta.tsv"
     conda:
          qiime_env
     shell:
          "python scripts/beta_diversity.py --inp {input} --outp {params}"

rule biom_to_tsv:
     input:
          "results/{cohort}/{id}.biom"
     output:
          "results/{cohort}/{id}_biom.tsv"
     conda:
          qiime_env
     shell:
          "biom convert -i {input} -o {output} --to-tsv"

rule manta:
     input:
          tsv="results/{cohort}/{id}+cls-{cls}+rrf{r}+otu_tax_biom.tsv",
          taxonpath="db/taxonpath.json",
          names="db/names.json"
     output:
          full="results/{cohort}/{id}+cls-{cls}+rrf{r}+manta.tsv",
          tax="results/{cohort}/{id}+cls-{cls}+rrf{r}+manta_tax.tsv"
     params:
          db=lambda wildcards: "1" if wildcards.cls=="gg" else "2"
     conda:
          "envs/other.yml"
     shell:
          "python scripts/manta.py -i {input.tsv} -o {output.full} -x {output.tax} "
          "-t {input.taxonpath} -n {input.names} -d {params.db} -r {wildcards.r}"

rule summary:
     input:
          "results/{cohort}/{id}+cls-{cls}_taxonomy.csv",
          "results/{cohort}/{id}_table+rrf{r}.csv",
          "results/{cohort}/{id}+rrf{r}+beta_weightedunifrac.csv",
          "results/{cohort}/{id}+rrf{r}+beta_unweightedunifrac.csv",
          "results/{cohort}/{id}+cls-{cls}+rrf{r}+otu_tax.biom",
          "results/{cohort}/{id}+cls-{cls}+rrf{r}+otu_tax_biom.tsv",
          "results/{cohort}/{id}+cls-{cls}+phyloseq.RDS",
          "results/{cohort}/{id}+rrf{r}+alphadiversity.tsv",
          "results/{cohort}/{id}+cls-{cls}+rrf{r}+manta.tsv",
          "results/{cohort}/{id}+cls-{cls}+rrf{r}+manta_tax.tsv",
          "results/{cohort}/{id}_seq.csv",
          "results/{cohort}/{id}+cls-{cls}+rrf{r}+beta_braycurtis.tsv",
          "results/{cohort}/{id}+cls-{cls}+rrf{r}+beta_jaccard.tsv"
          #"results/{cohort}/plots/{id}_stats.pdf"


     output:
          "results/{cohort}/{id}+cls-{cls}+rrf{r}.zip"
     conda:
          "envs/other.yml"
     shell:
          "zip -j {output} {input}"

#ruleorder: trim_bbduk > trim_fastp
ruleorder: merge_taxonomy > taxonomy > manifest
ruleorder: merge_dadatable > rarefy > manifest
ruleorder: export_phyloseq  > extract_biom
ruleorder: extract_biom > make_biom