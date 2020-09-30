from Bio import Entrez
import argparse
import pandas as pd
from functools import lru_cache
from tqdm import tqdm

parser = argparse.ArgumentParser(description='Script takes a list of NCBI TaxIDs as an input and returns linages.')
parser.add_argument('-e', '--email', help='NCBI requires email to connect to the server.')
parser.add_argument('-l', '--list', help='File with a list o NCBI TaxIDs (one per line).')
parser.add_argument('-out', '--out', help='Output File')
args = parser.parse_args()

Entrez.email = args.email

with open(args.list, 'r') as list_of_taxids:
    taxids = [taxid.rstrip() for taxid in list_of_taxids.readlines()]

# taxid (input), 'ScientificName', 'Rank'
lineage_df = pd.DataFrame()

# full lineage
lineage_full = pd.DataFrame()


@lru_cache(maxsize=1000)
def get_taxonomy(input_taxid):
    handle = Entrez.efetch(db="Taxonomy", id=input_taxid, retmode="xml")
    records = Entrez.read(handle)
    lineage_len = len(records[0]['LineageEx'])
    taxid_to_df = [input_taxid] * lineage_len
    sn_to_df = [None] * lineage_len
    rank_to_df = [None] * lineage_len

    for i in range(lineage_len):
        sn_to_df[i] = records[0]['LineageEx'][i]['ScientificName']
        rank_to_df[i] = records[0]['LineageEx'][i]['Rank']

    return pd.DataFrame({
               'TaxID': taxid_to_df,
               'ScientificName': sn_to_df,
               'Rank': rank_to_df
           }), \
           pd.DataFrame({
               'TaxID': [input_taxid],
               'Full_Linage': [';'.join(sn_to_df)],
               'Full_Rank': [';'.join(rank_to_df)]
           })


print(f'Processing {len(taxids)} TaxIDs')
for taxid in tqdm(taxids):
    taxid_l_df, taxid_lf_df = get_taxonomy(taxid)
    lineage_df = lineage_df.append(taxid_l_df)
    lineage_full = lineage_full.append(taxid_lf_df)

lineage_full = lineage_full.set_index('TaxID')
lineage_df = lineage_df.pivot(index='TaxID', columns='Rank', values='ScientificName')
lineage_df = lineage_df.merge(lineage_full, how='left', on='TaxID')
lineage_df.to_csv('Test/out.tsv', sep='\t')
print('Finished')