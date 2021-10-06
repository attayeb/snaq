import subprocess
import click
import qiime2
import biom
import os, zipfile
import pandas as pd
from qiime2 import Visualization, Artifact, Metadata
from os.path import join
from pandas import read_csv, unique
from os import listdir, mkdir
import shutil

def run_command(cmd, verbose=True):
    print("Running external command line application. This may print "
          "messages to stdout and/or stderr.")
    print("The command being run is below. This command cannot "
          "be manually re-run as it will depend on temporary files that "
          "no longer exist.")
    print("\nCommand:", end=' ')
    print(" ".join(cmd), end='\n\n')
    subprocess.run(cmd, check=True)

def fastp_trim(artifact, len1, len2):
    from q2_types.per_sample_sequences import (
        SingleLanePerSamplePairedEndFastqDirFmt,
        FastqManifestFormat,
        YamlFormat)
    art_ = artifact.view(SingleLanePerSamplePairedEndFastqDirFmt)
    manifest_o = pd.read_csv(os.path.join(
        str(art_), art_.manifest.pathspec), header=0, comment='#')
    manifest = manifest_o.copy()
    manifest.filename = manifest.filename.apply(
        lambda x: os.path.join(str(art_), x))
    id_to_fps = manifest.pivot(
        index='sample-id', columns='direction', values='filename')
    result = SingleLanePerSamplePairedEndFastqDirFmt()
    for _, (__, (fwd_fp, rev_fp)) in enumerate(id_to_fps.iterrows()):
        path1 = os.path.split(fwd_fp)[1]
        path2 = os.path.split(rev_fp)[1]

        p1 = str(os.path.join(result.path, path1))
        p2 = str(os.path.join(result.path, path2))

        cmd = ['fastp', '--in1 ' + fwd_fp, '--in2 ' + rev_fp, '--out1 ' + p1, 
        '-out2 ' + p2, '-Q ', "--trim_front1 "+ str(len1), "--trim_front2 "+str(len2)]
        run_command(cmd)
        os.remove(os.path.join(result.path, "fastp.json"))
        os.remove(os.path.join(result.path, "fastp.html"))


    result.manifest.write_data(art_.manifest.view(
        FastqManifestFormat), FastqManifestFormat)
    result.metadata.write_data(art_.metadata.view(YamlFormat), YamlFormat)
    return Artifact.import_data('SampleData[PairedEndSequencesWithQuality]', result)

@click.command()
@click.option("-i", "--inputf", required=True, type=str)
@click.option("--len1", required=True, type=int)
@click.option("--len2", required=True, type=int)
@click.option("-o", "--outputf", required=True, type=str)
def analyze(inputf, len1, len2, outputf):
    art = Artifact.load(inputf)
    trimmed = fastp_trim(art, len1, len2)
    trimmed.save(outputf)

if __name__ == "__main__":
    analyze()
