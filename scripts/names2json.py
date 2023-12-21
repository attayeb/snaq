import argparse
import pandas as pd
import json
from os.path import isfile
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Convert names.dmp to json dictionary keeping only scientific names")
    parser.add_argument("-i", dest="input", help="Input dmp file name", required=True)
    parser.add_argument("-o", dest="output", help="Output json file name", required=True)

    args = parser.parse_args()

    names = pd.read_csv(args.input, sep="|")
    names.columns = [x.strip() for x in names.columns]
    names = names.loc[names['synonym']=="\tscientific name\t"]
    names = names[['1','all']]
    names['all'] = [x.strip() for x in names['all']]
    names.set_index('1', inplace=True)
    names_dict = names['all'].to_dict()
    if not isfile(args.output):
        with open(args.output, "w") as j:
            json.dump(names_dict, j, indent=3)
        print("finish...")
    else:
        print("output file exist")
    

