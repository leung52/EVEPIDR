import requests
import json

def get_canonical_sequence_from_uniprot(uniprot_ids: list) -> dict:
    """
    Retrieves canonical protein sequences from the UniProt database for given UniProt IDs.

    This function queries the UniProt REST API for each ID in the uniprot_ids list. It extracts the gene name and sequence if the sequence length is less than 1024 amino acids. It prints a message if a sequence is too long or if a UniProt ID is not found.

    Parameters:
    - uniprot_ids (list): A list of UniProt IDs as strings.

    Returns:
    - tuple of two dictionaries:
        gene_to_sequence (dict): A dictionary mapping protein gene names to their corresponding protein sequences.
        gene_to_uniprot_ids (dict): A dictionary mapping protein gene names to their UniProt IDs.
    """
    gene_to_sequence = {}
    gene_to_uniprot_ids = {}

    for id in uniprot_ids:
        # Accessing Uniprot protein data for uniprot_id through Proteins REST API
        request_url = "https://www.ebi.ac.uk/proteins/api/proteins/" + id
        response = requests.get(request_url, headers={"Accept": "application/json"})
        if response.ok:
            responseBody = json.loads(response.text)
            if len(responseBody["sequence"]["sequence"]) < 1024:
                gene_to_sequence[responseBody["gene"][0]["name"]["value"]] = responseBody["sequence"]["sequence"]
                gene_to_uniprot_ids[responseBody["gene"][0]["name"]["value"]] = id
            else:
                print(f'{id} too long. Length: {len(responseBody["sequence"]["sequence"])}')
        else:
            print(f"{id} not found.")

    return gene_to_sequence, gene_to_uniprot_ids

def save_sequences_as_fasta(gene_to_sequence: dict, fasta_file: str) -> None:
    """
    Saves a dictionary of gene names and their protein sequences into a file in FASTA format.

    This function writes each gene and sequence pair from the gene_to_sequence dictionary to a file specified by fasta_file. Each sequence is formatted according to the FASTA file standard, starting with a '>' followed by the gene name.

    Parameters:
    - gene_to_sequence (dict): A dictionary mapping gene names to their corresponding protein sequences.
    - fasta_file (str): The file path where the FASTA file will be saved.

    Returns:
    - None
    """
    output_file = open(fasta_file,'w')
    for gene, seq in gene_to_sequence.items():
        identifier_line = ">" + gene + "\n"
        output_file.write(identifier_line)
        sequence_line = seq + "\n"
        output_file.write(sequence_line)
    output_file.close()
