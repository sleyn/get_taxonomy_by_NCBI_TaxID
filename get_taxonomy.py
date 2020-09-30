from Bio import Entrez
import argparse
import pandas as pd
from tqdm import tqdm

parser = argparse.ArgumentParser(description='Script takes a list of NCBI TaxIDs as an input and returns linages.')
parser.add_argument('-e', '--email', help='NCBI requires email to connect to the server.')
parser.add_argument('-l', '--list', help='File with a list o NCBI TaxIDs (one per line).')
parser.add_argument('-o', '--out', help='Output File')
parser.add_argument('-f', '--log', help='Log_file', default='log.log')
args = parser.parse_args()

Entrez.email = args.email

with open(args.list, 'r') as list_of_taxids:
    original_taxids = [taxid.rstrip() for taxid in list_of_taxids.readlines()]

taxids = original_taxids.copy()
taxids = list(set(taxids))

# Will be used to output data in the order as it was in the input file
original_order_df = pd.DataFrame({'TaxID': original_taxids})
original_order_df = original_order_df.set_index('TaxID')

# NCBI can process maximum 200 IDs at one request
# Splip initial set by lists size 200
if len(taxids) > 200:
    taxids_200 = [taxids[i:min(i+200, len(taxids) - 1)] for i in range(0, len(taxids), 200)]
else:
    taxids_200 = [taxids]

# Allowed ranks
allowed_ranks = [
    'superkingdom',
    'clade',
    'phylum',
    'subphylum',
    'class',
    'subclass',
    'order',
    'suborder',
    'family',
    'subfamily',
    'tribe',
    'genus',
    'subgenus',
    'species group',
    'species',
    'subspecies',
    'serotype',
    'forma specialis',
    'strain',
    'serogroup',
    'biotype'
]

# taxid (input), 'ScientificName', 'Rank'
lineage_df = pd.DataFrame({
    'TaxID': list(taxids),
    'superkingdom': ['-']*len(taxids),
    'clade': ['-']*len(taxids),
    'phylum': ['-']*len(taxids),
    'subphylum': ['-']*len(taxids),
    'class': ['-']*len(taxids),
    'subclass': ['-']*len(taxids),
    'order': ['-']*len(taxids),
    'suborder': ['-']*len(taxids),
    'family': ['-']*len(taxids),
    'subfamily': ['-']*len(taxids),
    'tribe': ['-']*len(taxids),
    'genus': ['-']*len(taxids),
    'subgenus': ['-']*len(taxids),
    'species group': ['-']*len(taxids),
    'species': ['-']*len(taxids),
    'subspecies': ['-']*len(taxids),
    'serotype': ['-']*len(taxids),
    'forma specialis': ['-']*len(taxids),
    'strain': ['-']*len(taxids),
    'serogroup': ['-']*len(taxids),
    'biotype': ['-']*len(taxids),
    'Full_Linage': ['-']*len(taxids),
    'Full_Rank': ['-']*len(taxids),
    'All_TaxIDs': ['-']*len(taxids),
    'Current_TaxID': ['-']*len(taxids)
})

lineage_df = lineage_df.set_index('TaxID')

print(f'Processing {len(taxids)} TaxIDs')

# Collect all taxids at once
for taxid_list in tqdm(taxids_200):
    handle = Entrez.efetch(db="Taxonomy", id=taxid_list, retmode="xml")
    records = Entrez.read(handle)

    for record in range(len(records)):
        if 'AkaTaxIds' in records[record]:
            all_taxids_for_record = records[record]['AkaTaxIds'].copy()
            all_taxids_for_record.append(records[record]['TaxId'])
            for aka_taxids in all_taxids_for_record:
                if aka_taxids in taxids:
                    taxid = aka_taxids
        else:
            all_taxids_for_record = [records[record]['TaxId']]
            taxid = records[record]['TaxId']
        # NCBI can output original query taxid in the AkaTaxid as the main TaxID could change.
        # Check where is the original taxid

        # if record was not fetched skip further operations
        if not records:
            with open(args.log, 'a') as log_file:
                log_file.write(f'Could not fetch {taxid}\n')
            continue

        for i in range(len(records[record]['LineageEx'])):
            # Check if rank exists in the dataframe
            if records[record]['LineageEx'][i]['Rank'] in allowed_ranks:
                lineage_df.loc[taxid, records[record]['LineageEx'][i]['Rank']] = records[record]['LineageEx'][i]['ScientificName']

        lineage_full = ';'.join([records[record]['LineageEx'][i]['ScientificName'] for i in range(len(records[record]['LineageEx']))])
        lineage_df.loc[taxid, 'Full_Linage'] = lineage_full
        rank_full = ';'.join([records[record]['LineageEx'][i]['Rank'] for i in range(len(records[record]['LineageEx']))])
        lineage_df.loc[taxid, 'Full_Rank'] = rank_full
        lineage_df.loc[taxid, 'All_TaxIDs'] = ';'.join([str(t) for t in all_taxids_for_record])
        lineage_df.loc[taxid, 'Current_TaxID'] = records[record]['TaxId']

not_fetched = lineage_df.query("Full_Linage == '-'").index
with open(args.log, 'w') as log_file:
    [log_file.write(taxid_nf + '\n') for taxid_nf in not_fetched]

original_order_df = original_order_df.merge(lineage_df, how='left', on='TaxID')
original_order_df.to_csv(args.out, sep='\t')
print('Finished')
