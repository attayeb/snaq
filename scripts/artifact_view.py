import click
from qiime2 import Artifact
from qiime2 import Metadata
from skbio import DistanceMatrix
import pandas as pd
import biom
import json        

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
        taxonomy_levels = ['kingdum', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        art = Artifact.load(artifact).view(biom.Table)
        
        meta_ = {x['id']: x['id'].split(";") for x in json.loads(art.to_json("QIIME2"))['rows']}

        meta__ = {}
        
        for k, v in meta_.items():
            current = {}
            for i in range(len(v)):
                current[taxonomy_levels[i]] = v[i]
                meta__.update({k:{'taxonomy': current}})
        
        art.add_metadata(meta__, axis='observation')
        art.type = "OTU table"
        art = art.remove_empty()
        with open(filename, "w") as f:
            f.write(art.to_json(generated_by="QIIME2"))


if __name__ == "__main__":
    export()