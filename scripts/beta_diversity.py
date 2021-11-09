import click
from qiime2 import Artifact
from qiime2 import Metadata
from skbio import DistanceMatrix
import pandas as pd
from qiime2.plugins.diversity.pipelines import beta

@click.command()
@click.option("--inp")
@click.option("--metrics", default='braycurtis,jaccard')
@click.option("--outp")
def beta_diversity(inp, metrics, outp):
    a = Artifact.load(inp)
    _metrics=metrics.split(",")

    for metric in _metrics:
        r = beta(a, metric)
        r.distance_matrix.view(DistanceMatrix).to_data_frame().to_csv(outp.split(".")[0]+"_" + 
            metric + "."+outp.split(".")[1])
if __name__ == "__main__":
    beta_diversity()
