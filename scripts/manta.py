import pandas as pd
import json
import click

def get_rank_from_ncbi(rank, taxonomy, taxonpath):
    ret = [[taxonomy.get(x[3:]) for x in y] for y in rank]
    ret1 = []
    for x in ret:
        txp = taxonpath[x[max([i for i in range(len(x)) if x[i] is not None])]]
        item = []
        
        for r in ['k', 'p', 'c', 'o', 'f', 'g', 's']:
            _item = txp[r]
            if _item == "":
                _item = "uc"
            item.append(_item)
        ret1.append(item)
    return ret1

@click.command()
@click.option("-i", "input_file", required=True, type=str)
@click.option("-o", "output_file", required=True, type=str)
@click.option("-t", "taxonpath", required=True, type=str)
@click.option("-n", "names", required=True, type=str)
@click.option("-d", "database", required=True, type=str)
def manta(input_file, output_file, taxonpath, names, database):
    df = pd.read_csv(input_file, sep="\t", skiprows=[0])
    with open(taxonpath) as f:
        taxonpath=json.load(f)
    with open(names) as f:
        names=json.load(f)
    taxonomy = {v:k for k, v in names.items()}
    taxs = [x.split(";") for x in df['#OTU ID']]
    d = pd.DataFrame(get_rank_from_ncbi(taxs, taxonomy, taxonpath))
    d.columns = ['0', '1', '2', '3', '4', '5', '6']
    df.drop("#OTU ID", axis="columns", inplace=True)
    df3 = d.join(df)
    
    df3m = df3.melt(id_vars=['0', '1', '2', '3', '4', '5', '6'])
    
    df3m = df3m.loc[:,['variable', '0', '1', '2', '3', '4', '5', '6', 'value']]
    
    df3m['pct'] = df3m['value']/10
    
    df3m['db'] = int(database)
    df3m['method'] = 1
    df3m = df3m[df3m['value']!= 0]
    df3m.to_csv(output_file, header=None, index=None, sep="\t")

if __name__ == "__main__":
    manta()
