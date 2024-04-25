import pandas as pd
import requests
import json
import xml.etree.ElementTree as ET

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

def clinvar_snp_missense_variants_id_list(gene: str, retmax: int=10000) -> list:
    """
    """
    request_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=(({gene}%5BGene%20Name%5D)%20AND%20%22single%20nucleotide%20variant%22%5BType%20of%20variation%5D)%20AND%20%22missense%20variant%22%5BMolecular%20consequence%5D&retmax={retmax}"
    response = requests.get(request_url)
    if response.ok:
        root = ET.fromstring(response.text)
        ids = [id_element.text for id_element in root.findall(".//Id")]
        if root.find(".//Count").text == '0':
            print(f"Zero variants found for {gene}.")
        if len(ids) < int(root.find(".//Count").text):
            print(f'Only {len(ids)} out of {root.find(".//Count").text} variant IDs are listed. To list all variant IDs, adjust retmax.')
        return ids
    else:
        response.raise_for_status()
        sys.exit()

def clinvar_variant_info(variant_ids: list) -> ET.Element:
    """
    """
    compiled_xml = ET.Element("EntrezESummaryResults")

    for i in range(0, len(variant_ids), 1):
        chunk = variant_ids[i:i+1]
        id_string = ",".join(chunk)
        request_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={id_string}"
        response = requests.get(request_url)
        if response.ok:
            chunk_xml = ET.fromstring(response.text)
            for doc_summary in chunk_xml.findall(".//DocumentSummary"):
                compiled_xml.append(doc_summary)
        else:
            response.raise_for_status()
            sys.exit()

    return compiled_xml

def clean_clinvar_xml_variants(protein_sequences: dict, clinvar_xml_root: ET.Element, csv_file_path: str) -> pd.DataFrame:
    """
    """
    # require the following information: sequence, gene, clinvar_id, aa_substitution, pathogenicity
    # clinvar_id from ????
    # gene and aa_substitution from title
    # pathogenicity from germline_classification/clinical_impact_classification/oncogenicity_classification
    clinvar_xml_root


