import click
from qiime2 import Artifact
from qiime2 import Metadata
from skbio import DistanceMatrix
import pandas as pd
import biom
        

@click.command()
@click.option("--artifact")
@click.option("--filename")
@click.option("--filetype")
def export(artifact, filename, filetype):
    
    if filetype=="metadata":        
        df = Artifact.load(artifact).view(Metadata).to_dataframe()
        df.to_csv(filename)
    
    if filetype=="distance":        
        df = Artifact.load(artifact).view(DistanceMatrix).to_data_frame()
        df.to_csv(filename)

    if filetype=="biom":
        with open(filename, "w") as f:
            Artifact.load(artifact).view(biom.Table).to_json(generated_by="QIIME2", direct_io=f)


if __name__ == "__main__":
    export()