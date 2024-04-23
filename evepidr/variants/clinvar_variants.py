import pandas as pd
import requests
import json

from evepidr.variants.utils import mutate_sequence

def get_canonical_sequence_from_uniprot(accessions: list) -> dict:
    """
    """
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
    """
    """
    pass

def set_up_variants_dataframe(protein_seq: str, clinvar_df: pd.DataFrame) -> pd.DataFrame:
    """
    """
    pass

def clinvar_snp_missense_id_list(gene: str, retmax: int=10000) -> list:
    """
    """
    requestURL = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=(({gene}%5BGene%20Name%5D)%20AND%20%22single%20nucleotide%20variant%22%5BType%20of%20variation%5D)%20AND%20%22missense%20variant%22%5BMolecular%20consequence%5D&retmax={retmax}"
    response = requests.get(requestURL)
    if response.ok:
        root = ET.fromstring(response.text)
        ids = [id_element.text for id_element in root.findall(".//Id")]
        if root.find(".//Count").text == '0':
            print(f"Zero variants found for {gene}.")
        if len(ids) < int(root.find(".//Count").text):
            print(f'Only {len(ids)} out of {root.find(".//Count").text} variant IDs are listed. To list all variant IDs, adjust retmax.')
        return ids
    else:
        raise Exception('An error occurred', 'Request URL invalid', requestURL)

def clinvar_variant_info(variant_ids: list, file_path: str) -> str:
    """
    """
    id_string = ",".join(variant_ids)
    request_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={id_string}"
    response = requests.get(request_url)
    if response.ok:
        with open(file_path, "wb") as file:
            file.write(response.content)
        return response.text
    else:
        raise Exception('An error occurred', 'Request URL invalid', requestURL)
