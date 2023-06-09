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
        if x == [None, None, None, None, None, None, None]:
            next
        else:
            max_id = x[max([i for i in range(len(x)) if x[i] is not None])]
            txp = taxonpath.get(max_id, empty_rank)
            item = []


            for r in ['k', 'p', 'c', 'o', 'f', 'g', 's']:
                _item = txp[r]
                if _item == "":
                    _item = "uc"
                item.append(_item)
            ret1.append(item)
    return ret1

def top_taxons(df):
    ret = df.sort_values(9, ascending=False)
    ret['cumsum'] = ret[9].cumsum()
    ret = ret[ret['cumsum'] < 90]
    #print(ret)
    ret.drop(['cumsum', 8,9,10], axis="columns", inplace=True)
    ret = ret.melt(id_vars=[0, 11])
    ret.drop_duplicates(inplace=True)
    ret.columns = ['sample_id', 'method_id','rank_id', 'taxonomy_id']
    ret = ret[ret['taxonomy_id']!= 'uc']
    return ret[['sample_id', 'rank_id', 'taxonomy_id', 'method_id']]

@click.command()
@click.option("-i", "input_file", required=True, type=str)
#@click.option("-v", "alphadiversity", required=True, type=str)
@click.option("-o", "output_file", required=True, type=str)
@click.option("-t", "taxonpath", required=True, type=str)
@click.option("-s", "sample_file_name", required=True, type=str)
@click.option("-a", "abundant_taxonomy", required=True, type=str)
@click.option("-n", "names", required=True, type=str)
@click.option("-d", "database", required=True, type=str)
@click.option("-r", "rarefaction", required=True, type=int)
@click.option("-x", "output_taxonomy", required=True, type=str)
#@click.option("-p", "output_alphadiversity", required=True, type=str)
def manta(input_file, output_file, taxonpath, abundant_taxonomy, sample_file_name, names, database, rarefaction, output_taxonomy):
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
    df3m['method_id'] = 1
    df3m = df3m[df3m['value']!= 0]
    df3m.rename(columns={
        '0' : "kingdom_id", '1': 'phylum_id', '2':'class_id', '3':'order_id', '4':'family_id','5':'genus_id', '6':'species_id','value':'read_num', 'variable':'sample_id', 'pct':'read_pct', 'db':'reference_db_id'}, inplace=True)
    df3m['read_num'] = df3m['read_num'].astype(int)
    df3m.to_csv(output_file, index=None)

    df3m['sample_id'].rename("id").drop_duplicates().to_csv(sample_file_name, index=False)

    df4['variable'] = df4['variable'].astype(int) + 1
    df4 = df4[df4['value']!= "uc"]
    df4['names'] = [names.get(x) for x in df4['value']]
    df4=df4.iloc[:,[1,0,2]].drop_duplicates()
    df4.columns = ['id', "rank_id", "name"]
    df4.to_csv(output_taxonomy, index=None)

    # abundant_taxons
    df5 = df3m.copy()
    df5.columns = range(12)
    abundant_taxons = pd.DataFrame(df5.groupby(0).apply(lambda x: top_taxons(x))).reset_index(drop=True)
    abundant_taxons.to_csv(abundant_taxonomy, index=False)

    # alpha diversity for manta:
    #alphadiversity_df = pd.read_csv(alphadiversity,comment="#")
    #alphadiversity_df.to_csv(output_alphadiversity)



if __name__ == "__main__":
    manta()
