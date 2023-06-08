import pandas as pd
import json
import click
@click.command()
@click.option("-i", "input_file", required=True, type=str)
@click.option("-o", "output_file", required=True, type=str)
def prepare_manta_alpha_diversity(input_file, output_file):


    # alpha diversity for manta:
    alphadiversity_df = pd.read_csv(input_file,comment="#", sep="\t")
    alphadiversity_df.drop("observed_features", axis="columns", inplace=True)
    alphadiversity_df['method_id'] = 1

    alphadiversity_df.rename(columns={'Sample ID':'sample_id',
                                      "shannon_entropy":"shannon"}, inplace=True)
    alphadiversity_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    prepare_manta_alpha_diversity()

