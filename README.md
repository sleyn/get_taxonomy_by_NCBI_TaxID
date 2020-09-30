# get_taxonomy_by_NCBI_TaxID
Get linage string by NCBI TaxID

The script takes a list of NCBI TaxIDs (one per line) and translate them to linage strings via BioPython Entrez module.

##Requires:
 
 - BioPython
 - Pandas
 - tqdm
 
Usage:
```python
$ get_taxonomy.py [-h] [-e EMAIL] [-l LIST] [-out OUT]

Arguments:
  -h, --help            show this help message and exit
  -e EMAIL, --email EMAIL  NCBI requires email to connect to the server.
  -l LIST, --list LIST     File with a list o NCBI TaxIDs (one per line).
  -o OUT, --out OUT        Output File
  -f LOG, --log LOG        Log file (filed TaxIDs). Default: log.log 
```

The script is processing 100 TaxIDs in 1-2 second.

**NOTE:** Sometimes NCBI randomly does not return records for query for unknow reason. Try again to run a list of not received TaxIDs from the log file. In the issue remains - likely you have wrong TaxID or TaxID that was deleted from the Taxonomy database.  

## Output

Output is a Tab Separated Values (TSV) table with the following columns:

- TaxID
- One column per each rank. Only selected ranks are in the output. The most notable is that 'no rank' ranks are omitted as it is hard to automatically evaluate if they could be in the same column.
  - The ranks that are in the output:
    - superkingdom
    - clade
    - phylum
    - subphylum
    - class
    - subclass
    - order
    - suborder
    - family
    - subfamily
    - tribe
    - genus
    - subgenus
    - species group
    - species
    - subspecies
    - serotype
    - forma specialis
    - strain
    - serogroup
    - biotype
- Full_Linage - the string with full lineage without any omissions.
- Full_Rank - the string with ranks names to the corresponding positions in the *Full_Lineage* string.
- All_TaxIDs - all TaxIDs associated with the organism.
- Current_TaxID - the TaxID that NCBI shows as the most current.
