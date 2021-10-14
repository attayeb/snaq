import pandas as pd
import click
import os
from os.path import join, abspath




identifiers = [{'R1':"_1", 'R2':"_2"},
               {'R1':"_R1_", "R2": "_R2_"}]


@click.command()
@click.option("-i", "input_folder", required=True, type=str)
@click.option("-o", "manifest_file_name", required=True, type=str)
def create_manifest(input_folder, manifest_file_name):
    files = [join(input_folder, x) for x in os.listdir(input_folder)]
    for iden in identifiers:
        # step 1
        R1_id = iden['R1']
        R2_id = iden['R2']
        
        R1 = [x for x in files if R1_id in x]
        R1.sort()
        
        R2 = [x for x in files if R2_id in x]
        R2.sort()
        sample_id1 = [x.split('_')[0].split("/")[-1] for x in R1]
        sample_id2 = [x.split('_')[0].split("/")[-1] for x in R2]
        condition2 = all([x==y for x,y in zip(sample_id1, sample_id2)])
        if len(R1) > 0:
            break
    
    condition3 = all([_r[0].replace(R1_id, "_")==_r[1].replace(R2_id, '_') for _r in zip(R1, R2)])
    if all([condition2, condition3]):
        R1 = [abspath(x) for x in R1]
        R2 = [abspath(x) for x in R2]
        res = {"sample-id": sample_id1,
               "forward-absolute-filepath": R1,
               "reverse-absolute-filepath": R2}
    pd.DataFrame(res).to_csv(manifest_file_name, sep="\t", index=False)

if __name__ == "__main__":
    create_manifest()