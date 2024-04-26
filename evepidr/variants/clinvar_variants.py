import pandas as pd
import requests, sys
import json
import time
import re
import xml.etree.ElementTree as ET

from evepidr.variants.utils import mutate_sequence, three_to_one_aa_code, adjust_clinvar_classification

def get_canonical_sequence_from_uniprot(accessions: list) -> dict:
    """
    """
    gene_to_sequence = {}
    gene_to_accession = {}
    
    for accession in accessions:
        # Accessing Uniprot protein data for uniprot_id through Proteins REST API
        request_url = "https://www.ebi.ac.uk/proteins/api/proteins/" + accession
        response = requests.get(request_url, headers={"Accept": "application/json"})
        if response.ok:
            responseBody = json.loads(response.text)
            gene_to_sequence[responseBody["gene"][0]["name"]["value"]] = responseBody["sequence"]["sequence"]
            gene_to_accession[responseBody["gene"][0]["name"]["value"]] = accession
        else:
            print(f"{accession} not found.")

    return gene_to_sequence, gene_to_accession

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

def clinvar_variant_info(variant_ids: list) -> ET.ElementTree:
    """
    """
    compiled_xml = ET.Element("EntrezESummaryResults")

    for i in range(0, len(variant_ids), 500):
        chunk = variant_ids[i:i+500]
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
        time.sleep(0.333333334)

    return ET.ElementTree(compiled_xml)

def clean_clinvar_xml_variants(protein_sequences: dict, gene_to_accession: dict, clinvar_xml: ET.Element) -> pd.DataFrame:
    """
    """
    sequences = []
    genes = []
    aa_substitutions = []
    pathogenicities = []
    clinvar_ids = []
    uniprot_ids = []
    
    for variant_xml in clinvar_xml.findall(".//DocumentSummary"):
        germline_classification = variant_xml.find('.//germline_classification/description').text
        germline_classification = adjust_clinvar_classification(germline_classification)
        
        if germline_classification in ['Pathogenic', 'Benign']:
            title = variant_xml.find('title').text
            parts = title.split('(')
            mutation = parts[-1].split(')')[0]
            if re.fullmatch(r"p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}", mutation):
                try:
                    mutation = three_to_one_aa_code(mutation[2:5]) + mutation[5:-3] + three_to_one_aa_code(mutation[-3:])
                except ValueError:
                    mutation = None
                if mutation:
                    gene = parts[1].split(')')[0]
                    canon_sequence = protein_sequences.get(gene)
                    if canon_sequence:
                        sequence = mutate_sequence(canon_sequence, int(mutation[1:-1]), int(mutation[1:-1]), mutation[-1])
                        id = variant_xml.get('uid')

                        sequences.append(sequence)
                        genes.append(gene)
                        aa_substitutions.append(mutation)
                        pathogenicities.append(germline_classification)
                        clinvar_ids.append(id)
                        uniprot_ids.append(gene_to_accession.get(gene))

    for gene in set(genes):
        sequences.append(protein_sequences[gene])
        genes.append(gene)
        aa_substitutions.append(None)
        pathogenicities.append("Canon")
        clinvar_ids.append(None)
        uniprot_ids.append(gene_to_accession.get(gene))
        
    data = {
        'Sequence': sequences,
        'Gene': genes,
        'AA Substitution': aa_substitutions,
        'Pathogenicity': pathogenicities,
        'ClinVar ID': clinvar_ids,
        'UniProt ID': uniprot_ids
    }
    
    return pd.DataFrame(data)
