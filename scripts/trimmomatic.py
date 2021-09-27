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

def trimmomatic(artifact, trimming_threshold, sliding_window, headcrop=0):
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
        p1_u = str(os.path.join(result.path, "temp_"+path1))
        p2 = str(os.path.join(result.path, path2))
        p2_u = str(os.path.join(result.path, "temp_"+path2))

        cmd = ['trimmomatic', 'PE', fwd_fp, rev_fp, 
               p1, p1_u, p2, p2_u, "HEADCROP:"+str(headcrop),
               'SLIDINGWINDOW:'+str(sliding_window)+":"+str(trimming_threshold)]
        run_command(cmd)
        os.remove(p1_u)
        os.remove(p2_u)

    result.manifest.write_data(art_.manifest.view(
        FastqManifestFormat), FastqManifestFormat)
    result.metadata.write_data(art_.metadata.view(YamlFormat), YamlFormat)
    return Artifact.import_data('SampleData[PairedEndSequencesWithQuality]', result)

@click.command()
@click.option("-i", "file_name", type=str, required=True)
@click.option("-q", "quality_threshold", type=int, required=True)
@click.option("-w", "sliding_window", type=int, required=True)
@click.option("-h", "headcrop", type=int, required=True)
@click.option("-o", "output", type=str, required=True)
def analyze(file_name, quality_threshold, sliding_window, headcrop, output):
    art = Artifact.load(file_name)
    trimmed = trimmomatic(art, quality_threshold, sliding_window, headcrop)
    trimmed.save(output)


if __name__ == "__main__":
    analyze()
