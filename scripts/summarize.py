import click
import os
import json
import pandas as pd

@click.command()
@click.option("--files")
@click.option("--outp")
def summarize(files, outp):
    files = files.split(",")
    result = {}
    for file in files:
        _df = pd.read_csv(file, sep="\t")
        _df.columns=['percent', "count_accu", "count", "TL", "TID", "Taxonomy"]
        _df['index'] = _df.index 
        _df = _df.iloc[:,[6,3,4,5,2,1,0]]
        _df['Taxonomy'] = _df['Taxonomy'].apply(lambda x: x.strip())
        result[file] = _df.to_dict(orient="index")


    if os.path.isfile(outp):
        raise Exception(outp, "File exists")
    else:
        with open(outp, "w") as f:
            json.dump(result, f, indent=4)


if __name__ == "__main__":
    summarize()