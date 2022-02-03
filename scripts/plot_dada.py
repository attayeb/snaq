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
    plt.figure(figsize=(5, 5), dpi=100)
    sns.displot(x='percentage of input non-chimeric', data=df)
    plt.title(plot.split("/")[-1].replace("+dd_stats.jpg", ""))
    plt.xlim([0, 100])
    plt.tight_layout()
    plt.savefig(plot)

if __name__ == "__main__":
    plot_dada()