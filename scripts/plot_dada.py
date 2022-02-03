import click
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns

from qiime2 import Artifact
from qiime2 import Metadata

@click.command()
@click.option("--inp")
@click.option("--plot")
def plot_dada(inp, plot):
    art = Artifact.load(inp)
    df = art.view(Metadata).to_dataframe()
    sns.displot(x='percentage of input non-chimeric', data=df)
    plt.savefig(plot, dpi=300)

if __name__ == "__main__":
    plot_dada()