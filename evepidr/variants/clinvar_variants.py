import pandas as pd
import requests
import json

from evepidr.variants.utils import mutate_sequence

def get_canonical_sequence_from_uniprot(accessions: list) -> dict:
    gene_to_sequence = {}
    
    for accession in accessions:
        # Accessing Uniprot protein data for uniprot_id through Proteins REST API
        request_url = "https://www.ebi.ac.uk/proteins/api/proteins/" + uniprot_id
        response = requests.get(request_url, headers={"Accept": "application/json"})
        if response.ok:
            responseBody = json.loads(response.text)
            gene_to_sequence[responseBody["gene"][0]["name"]["value"]] = responseBody["sequence"]["sequence"]
        else:
            print(f"{accession} not found.")

    return gene_to_sequence

def get_missense_mutations_from_clinvar(gene: str) -> pd.DataFrame:
    pass

def set_up_variants_dataframe(protein_seq: str, clinvar_df: pd.DataFrame) -> pd.DataFrame:
    pass
