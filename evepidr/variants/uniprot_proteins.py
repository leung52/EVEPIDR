import requests
import json

def get_canonical_sequence_from_uniprot(uniprot_ids: list) -> dict:
    """
    """
    gene_to_sequence = {}
    gene_to_uniprot_ids = {}

    for id in uniprot_ids:
        # Accessing Uniprot protein data for uniprot_id through Proteins REST API
        request_url = "https://www.ebi.ac.uk/proteins/api/proteins/" + id
        response = requests.get(request_url, headers={"Accept": "application/json"})
        if response.ok:
            responseBody = json.loads(response.text)
            gene_to_sequence[responseBody["gene"][0]["name"]["value"]] = responseBody["sequence"]["sequence"]
            gene_to_uniprot_ids[responseBody["gene"][0]["name"]["value"]] = id
        else:
            print(f"{id} not found.")

    return gene_to_sequence, gene_to_uniprot_ids

def save_sequences_as_fasta(gene_to_sequence: dict, file_path: str) -> None:
    """
    """
    output_file = open(file_path,'w')
    for gene, seq in gene_to_sequence.items():
        identifier_line = ">" + gene + "\n"
        output_file.write(identifier_line)
        sequence_line = seq + "\n"
        output_file.write(sequence_line)
    output_file.close()
