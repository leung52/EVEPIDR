import numpy as np
import pandas as pd
import h5py


def fasta_to_dict(file_path: str) -> dict:
    """
    """
    parsed_seqs = {}

    with open(file_path, 'r') as f:
        curr_seq_id = None
        curr_seq = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if curr_seq_id is not None:
                    parsed_seqs[curr_seq_id] = ''.join(curr_seq)
                curr_seq_id = line[1:]
                curr_seq = []
                continue
            curr_seq.append(line)

        # Add the last sequence to the dictionary
        if curr_seq_id is not None:
            parsed_seqs[curr_seq_id] = ''.join(curr_seq)

    return parsed_seqs

def add_sequences_to_variants_df(variants_df: pd.DataFrame, protein_sequences: dict) -> pd.DataFrame:
    """
    """
    df = variants_df.copy()
    df['Sequence'] = None

    for index, row in df.iterrows():
        gene = row['Gene']
        sequence = protein_sequences.get(gene)
        if sequence:
            aa_substitution = row['AA Substitution']
            if len(aa_substitution) > 2:
                sequence = _mutate(sequence, aa_substitution)
            df.at[index, 'Sequence'] = sequence
        else:
            df.drop(index, inplace=True)

    return df

def reduce(embeddings: np.ndarray) -> np.ndarray:
    """
    """
    reduced_embeddings = []
    for full_embedding in embeddings:
        reduced_embeddings.append(full_embedding.mean(axis=0))
    return np.array(reduced_embeddings)

def normalise(embeddings: np.ndarray) -> np.ndarray:
    """
    """
    return embeddings - np.mean(embeddings, axis=0)

def cosine_distance(reduced_embeddings: dict, variants_df_w_sequence: pd.DataFrame, gene_to_sequence: dict) -> pd.DataFrame:
    """
    """
    gene_to_embedding = {}
    for key, value in gene_to_sequence.items():
        if value in reduced_embeddings:
            gene_to_embedding[key] = reduced_embeddings[value]
    
    
    variants_df = variants_df_w_sequence.copy()

    for index, row in variants_df.iterrows():
        if row['Sequence'] in reduced_embeddings:
            canon_embedding = gene_to_embedding.get(row['Gene'])
            if canon_embedding:
                variant_embedding = reduced_embeddings[row['Sequence']]
                cosine_distance = distance.cosine(canon_embedding, variant_embedding)
                variants_df.at[index, 'ESM-1b Cosine Distance'] = cosine_distance
            else:
                variants_df.drop(index, inplace=True)
        else:
            variants_df.drop(index, inplace=True)

    variants_df.drop('Sequence', axis=1, inplace=True)

    return variants_df

def hdf5_to_dict(file_path: str) -> dict:
    """
    """
    embeddings = {}
    with h5py.File(file_path, 'r') as hf:
        for sequence in hf.keys():
          embeddings[sequence] = hf[sequence][:]
    return embeddings



## ==================== helper function ===============================
def _mutate(protein_sequence: str, mutation: str) -> str:
    """
    """
    # Adjusting for 1-based index
    n = int(mutation[1:-1])
    start_index = n - 1
    end_index = n

    # Ensure the start and end positions are within the original sequence
    if start_index < 0 or end_index > len(protein_sequence) or start_index > end_index:
        raise ValueError("Invalid start or end position for the sequence mutation.")

    assert protein_sequence[start_index] == mutation[0]
    # Replace the sequence between start and end with new_sequence
    mutated_sequence = protein_sequence[:start_index] + mutation[-1] + protein_sequence[end_index:]

    return mutated_sequence
