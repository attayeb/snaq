import click
import pandas as pd
import biom
import qiime2
from qiime2 import Artifact

@click.command()
@click.option("--tablef")
@click.option("--taxonomy")
@click.option("--output")
def create_biom_table(tablef, taxonomy, output):

    message="QIIME2"
    biom_table = Artifact.load(tablef).view(biom.Table)
    biom_table.type = "OTU table"
    md = Artifact.load(taxonomy).view(pd.DataFrame)
    md.index.name = '#OTU ID'
    md.columns = ['taxonomy', 'confidence']
    biom_table.add_metadata(md.to_dict("index"), axis="observation")
    
    with open(output, "w") as f:
        f.write(biom_table.to_json(message))

if __name__ == "__main__":
    create_biom_table()