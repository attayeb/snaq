import pandas as pd
import json
import click

empty_rank = {
    'k': '', 'p':'', 'c':'', 'o':'', 'f':'','g':'','s':''
}

def get_rank_from_ncbi(rank, taxonomy, taxonpath):
    ret = [[taxonomy.get(x[3:]) for x in y] for y in rank]
    ret1 = []    
    for x in ret:
        try:
            txp = taxonpath.get(x[max([i for i in range(len(x)) if x[i] is not None])], empty_rank)
            item = []
        except Exception as e:
            print e
            print(txp)
        
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
@click.option("-r", "rarefaction", required=True, type=int)
@click.option("-x", "output_taxonomy", required=True, type=str)
def manta(input_file, output_file, taxonpath, names, database, rarefaction, output_taxonomy):
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
    df4 = df3m.loc[:,['0', '1', '2', '3', '4', '5', '6']].copy().melt()

    df3m = df3m.loc[:,['variable', '0', '1', '2', '3', '4', '5', '6', 'value']]
    
    df3m['pct'] = (df3m['value']/rarefaction)*100
    
    df3m['db'] = int(database)
    df3m['method'] = 1
    df3m = df3m[df3m['value']!= 0]
    df3m.to_csv(output_file, header=None, index=None, sep="\t")

    df4['variable'] = df4['variable'].astype(int) + 1
    df4 = df4[df4['value']!= "uc"]
    df4['names'] = [names.get(x) for x in df4['value']]
    df4.iloc[:,[1,0,2]].drop_duplicates().to_csv(output_taxonomy, index=None, header=None, sep="\t")

if __name__ == "__main__":
    manta()
