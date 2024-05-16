import pandas as pd
import requests, sys
import time
import re
import xml.etree.ElementTree as ET


def clinvar_snp_missense_variants_id_list(gene: str, retmax: int=10000) -> list:
    """
    Retrieves a list of ClinVar variant IDs for missense single nucleotide polymorphisms (SNPs) associated with a specified gene.

    This function queries the ClinVar database via NCBI's E-utilities API to fetch variant IDs for a specified gene that are categorized as "missense variant" under the type of single nucleotide variant. It prints warnings if no variants are found or if the number of retrieved IDs is less than the total number available.
    
    Parameters:
    - gene (str): The gene name to query for variants.
    - retmax (int, optional): The maximum number of variant IDs to retrieve. Defaults to 10,000.
    
    Returns:
    - list: A list of ClinVar variant IDs as strings.
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
    Fetches detailed information for a list of ClinVar variant IDs and compiles them into an XML structure.

    This function retrieves detailed variant information from the ClinVar database for specified variant IDs. It processes the IDs in chunks to avoid query size limits and aggregates the results into a single XML element tree. The function includes a short delay between chunks to comply with API usage policies.
    
    Parameters:
    - variant_ids (list): A list of ClinVar variant IDs for which detailed information is requested.
    
    Returns:
    - ET.ElementTree: An XML element tree containing detailed information for the requested variants.
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

def clean_clinvar_xml_variants(gene_to_uniprot_id: dict, clinvar_xml: ET.Element) -> pd.DataFrame:
    """
    Extracts relevant data from a ClinVar XML document and formats it into a pandas DataFrame.

    This function parses an XML element tree of ClinVar variant data, filtering and transforming the data based on specific criteria (e.g., variants classified as either 'Pathogenic' or 'Benign'). It maps three-letter amino acid codes to their single-letter counterparts among other transformations and compiles the results into a structured DataFrame.
    
    Parameters:
    - gene_to_uniprot_id (dict): A dictionary mapping gene names to UniProt IDs.
    - clinvar_xml (ET.Element): An XML element representing the root of ClinVar data.
    
    Returns:
    - pd.DataFrame: A DataFrame containing columns for Gene, AA Substitution, Pathogenicity, ClinVar ID, and UniProt ID.
    """
    genes = []
    aa_substitutions = []
    pathogenicities = []
    clinvar_ids = []
    uniprot_ids = []

    for variant_xml in clinvar_xml.findall(".//DocumentSummary"):
        germline_classification = variant_xml.find('.//germline_classification/description').text
        germline_classification = _adjust_clinvar_classification(germline_classification)

        if germline_classification in ['Pathogenic', 'Benign']:
            title = variant_xml.find('title').text
            parts = title.split('(')
            mutation = parts[-1].split(')')[0]
            if re.fullmatch(r"p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}", mutation):
                try:
                    mutation = _three_to_one_aa_code(mutation[2:5]) + mutation[5:-3] + _three_to_one_aa_code(mutation[-3:])
                except ValueError:
                    mutation = None
                if mutation:
                    gene = parts[1].split(')')[0]
                    id = variant_xml.get('uid')
                    genes.append(gene)
                    aa_substitutions.append(mutation)
                    pathogenicities.append(germline_classification)
                    clinvar_ids.append(id)
                    uniprot_ids.append(gene_to_uniprot_id.get(gene))

    data = {
        'Gene': genes,
        'AA Substitution': aa_substitutions,
        'Pathogenicity': pathogenicities,
        'ClinVar ID': clinvar_ids,
        'UniProt ID': uniprot_ids
    }

    return pd.DataFrame(data)


## ============ Helper functions ===========================================

def _three_to_one_aa_code(code: str) -> str:
    """
    Converts a three-letter amino acid code to its corresponding single-letter code.

    Parameters:
    - code (str): The three-letter amino acid code.
    
    Returns:
    - str: The single-letter amino acid code.
    
    Raises:
    - ValueError: If the input code is not a recognized three-letter amino acid code.
    """
    three_to_single = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Glu": "E",
        "Gln": "Q",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V"
    }

    one_code = three_to_single.get(code)

    if not one_code:
        raise ValueError(f"{code} is not an amino acid code")
    else:
        return one_code

def _adjust_clinvar_classification(clinvar_classification: str) -> str:
    """
    Adjusts the classification of ClinVar entries to a simpler form ('Pathogenic', 'Benign', or 'Other') based on predefined groupings.

    Parameters:
    - clinvar_classification (str): The original classification from ClinVar.
    
    Returns:
    - str: A simplified classification as 'Pathogenic', 'Benign', or 'Other'.
    """
    groups = {
        ("Benign", "Likely benign", "protective"): 'Benign',
        ("Likely pathogenic", "Pathogenic", "Likely pathogenic, low penetrance", "Pathogenic, low penetrance", "Likely risk allele", "Established risk allele", "association"): 'Pathogenic'
    }
    for key, value in groups.items():
        if clinvar_classification in key:
            return value
    return "Other"
