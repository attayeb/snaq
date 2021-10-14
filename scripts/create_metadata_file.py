import pandas as pd
import click
import os
from os.path import join, abspath




identifiers = [{'R1':"_1", 'R2':"_2"},
               {'R1':"_R1_", "R2": "_R2_"}]


@click.command()
@click.option("-i", "manifest_file_name", required=True, type=str)
@click.option("-o", "metadata_file_name", required=True, type=str)
def create_metadata_file(manifest_file_name, metadata_file_name):
    df = pd.read_csv(manifest_file_name, sep="\t")
    df.drop(["forward-absolute-filepath", "reverse-absolute-filepath"], axis="columns")
    df.to_csv(metadata_file_name, sep="\t", index=False)

if __name__ == "__main__":
    create_metadata_file()