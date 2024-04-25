import pandas as pd
import requests
import json
import xml.etree.ElementTree as ET

from evepidr.variants.utils import mutate_sequence, three_to_one_aa_code

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

def clean_clinvar_xml_variants(protein_sequences: dict, clinvar_xml: ET.Element) -> pd.DataFrame:
    """
    """
    # require the following information: sequence, gene, clinvar_id, aa_substitution, pathogenicity
    # make remember to add canon
    sequences = []
    genes = []
    aa_substitutions = []
    pathogenicities = []
    clinvar_ids = []
    
    for variant_xml in x.findall(".//DocumentSummary"):
        germline_classification = variant_xml.find('.//germline_classification/description').text
        germline_classification = adjust_clinvar_classification(germline_classification)
        if germline_classification in ['Pathogenic', 'Benign']:
            clinvar_ids.append(variant_xml.get('uid'))
            title = variant_xml.find('title').text
            parts = title.split('(')
            genes.append(parts[1].split(')')[0])
            mutation = parts[-1].split(')')[0]
            aa_substitutions.append(three_to_one_aa_code(mutation[2:6]) + mutation[6:-3] + three_to_one_aa_code(mutation[-3:]))
            sequences.append(mutate_sequence(canon_sequence, mutation[1:-1], mutation[1:-1], mutation[-1]))
            pathogenicities.append(germline_classification)
    
    data = {
        'Sequence': sequences,
        'Gene': genes,
        'AA Substitution': aa_substitutions,
        'Pathogenicity': pathogenicities,
        'ClinVar ID': clinvar_ids
    }

    return pd.DataFrame()
        
