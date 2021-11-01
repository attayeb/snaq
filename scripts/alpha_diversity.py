    """This script calculates Alpha diversity using QIIME2 functions.
    """



import click
from qiime2 import Artifact
from qiime2 import Metadata
import pandas as pd
from qiime2.plugins.diversity.pipelines import alpha

@click.command()
@click.option("--inp")
@click.option("--metrics", default='simpson,chao1,shannon,observed_features')
@click.option("--outp")
def alpha_diversity(inp, metrics, outp):
    """[summary]

    Args:
        inp (str): input file name
        metrics (str): metrics separated by comma, no spaces
        outp (str): output file names
    """
    a = Artifact.load(inp)
    _metrics=metrics.split(",")
    ret = []

    for metric in _metrics:
        r = alpha(a, metric)
        ret.append(r.alpha_diversity.view(Metadata).to_dataframe())

    x = pd.concat(ret, axis=1)
    Metadata(x).save(outp)
if __name__ == "__main__":
    alpha_diversity()
