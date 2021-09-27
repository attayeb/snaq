import json
from os.path import isfile, join, basename
import argparse
import pandas as pd 

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Create CAMI file out of kaiju results")
    parser.add_argument("-i", dest="input", help="Input file name", required=True)
    parser.add_argument("-o", dest="output", help="Output CAMI file name", required=True)
    parser.add_argument("-d", dest="db", help="Database folder where [names.json, taxonpath.json]")
    parser.add_argument("-f", dest="filt", type=float, help="Percentage threshold below whicn will be ignored [default=0]", default=0)
    parser.add_argument("-s", dest="sid", help="Sample ID [default: file name]", default="")

    args = parser.parse_args()

    if args.sid == "":
        sid = basename(args.output).replace("-cami.txt", "")
    else:
        sid = args.sid

    with open(join(args.db, "names.json")) as f:
        names = json.load(f)

    names[''] = ''
    with open(join(args.db, "taxonpath.json")) as f:
        taxpath = json.load(f)
    permitted_ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum']
    ranklevel= {
        'phylum' : 1,
        'class'  : 2,
        'order'  : 3,
        'family' : 4,
        'genus'  : 5,
        'species': 6
        }
    rankorder = {
        'phylum' : 0,
        'class'  : 1,
        'order'  : 2,
        'family' : 3,
        'genus'  : 4,
        'species': 5
    }
    df = pd.read_csv(args.input, sep="\t", header=None)

    output=[]
    for rowi in df.iterrows():
        row = rowi[1]
        if str(row[4]) in taxpath:
            if row[0] < args.filt :
                continue
            outrow = {}
            tx = taxpath[str(row[4])]
            if tx['rank'] not in permitted_ranks:
                continue
            outrow['taxid'] = str(row[4])
            outrow['rank'] = tx['rank']
            number=ranklevel[tx['rank']]
            txids = list(tx.values())[1:number+2]
            outrow['taxpath'] = "|".join(txids)
            outrow['taxpathsn'] = "|".join([names[x] for x in txids])
            outrow['percentage'] = row[0]
            output.append(outrow)

    sorted_output = sorted(output, key=lambda d: (rankorder[d['rank']], d['percentage']), reverse=True)
    outputstring="# Taxonomic Profiling Output\n"
    outputstring+="@SampleID:{}\n".format(sid)
    outputstring+="@version:0.9.1\n"
    outputstring+="@Ranks:superkingdum|phylum|class|family|order|genus|species\n"
    outputstring+="@TaxonomyID:\n"
    outputstring+="@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n"

    for item in sorted_output:
        outputstring += "{}\t{}\t{}\t{}\t{}\n".format(item['taxid'], item['rank'], item['taxpath'], item['taxpathsn'], item['percentage'])


    if not isfile(args.output):
        with open(args.output, "w") as cami:
            cami.write(outputstring)
    else:
        print("File is exist... printing to screen...")
        print(outputstring)

