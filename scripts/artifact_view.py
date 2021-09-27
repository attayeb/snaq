import click
from qiime2 import Artifact
from qiime2 import Metadata
from skbio import DistanceMatrix
import pandas as pd
        

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


if __name__ == "__main__":
    export()