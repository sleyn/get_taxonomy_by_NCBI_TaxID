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
  -out OUT, --out OUT      Output File
```