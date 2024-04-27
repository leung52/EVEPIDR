import requests
import json

def get_canonical_sequence_from_uniprot(uniprot_ids: list) -> dict:
    """
    """
    gene_to_sequence = {}
    gene_to_uniprot_ids = {}
    
    for id in uniprot_ids:
        # Accessing Uniprot protein data for uniprot_id through Proteins REST API
        request_url = "https://www.ebi.ac.uk/proteins/api/proteins/" + accession
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
    # Read existing contents of the file, if it exists
    existing_gene_to_sequence = {}
    if pathlib.Path(file_path).is_file():
        with open(file_path, "r") as f:
            gene_name = None
            sequence = ""
            for line in f:
                if line.startswith(">"):
                    # Store previous gene's sequence, if any
                    if gene_name is not None:
                        existing_gene_to_sequence[gene_name] = sequence
                    # Extract gene name from FASTA header
                    gene_name = line.strip()[1:]
                    sequence = ""
                else:
                    sequence += line.strip()
            # Store last gene's sequence
            if gene_name is not None:
                existing_gene_to_sequence[gene_name] = sequence

    # Update existing sequences with new ones
    existing_gene_to_sequence.update(gene_to_sequence)

    # Write the updated contents back to the file
    with open(file_path, "w") as f:
        for gene, sequence in existing_gene_to_sequence.items():
            f.write(f">{gene}\n")
            f.write(f"{sequence}\n")
