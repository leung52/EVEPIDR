import numpy as np
import pandas as pd
import h5py
from scipy.spatial.distance import cosine


def fasta_to_dict(file_path: str) -> dict:
    """
    Reads a FASTA file and returns a dictionary mapping sequence identifiers to sequences.

    Parameters:
    - file_path (str): Path to the FASTA file.
    
    Returns:
    - dict: A dictionary where keys are sequence identifiers (without the '>' character), and values are sequences.
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

def add_sequences_to_variants_df(variants_df: pd.DataFrame, gene_to_sequence: dict) -> pd.DataFrame:
    """
    Adds sequences to a DataFrame of variants, mutating them as specified in the DataFrame entries.

    This function takes a DataFrame containing variants and their associated genes, retrieves the canonical sequence for each gene, performs specified mutations, and adds these sequences to the DataFrame.
    
    Parameters:
    - variants_df (pd.DataFrame): DataFrame containing columns 'Gene' and 'AA Substitution'.
    - gene_to_sequence (dict): Dictionary mapping gene names to their canonical sequences.
    
    Returns:
    - pd.DataFrame: The updated DataFrame with a new 'Sequence' column containing the mutated sequences.
    """
    df = variants_df.copy()
    df['Sequence'] = None

    for index, row in df.iterrows():
        gene = row['Gene']
        sequence = gene_to_sequence.get(gene)
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
    Reduces embeddings by taking the mean across the first dimension.

    Parameters:
    - embeddings (np.ndarray): An array of embeddings, typically with multiple dimensions per embedding.
    
    Returns:
    - np.ndarray: An array of reduced embeddings.
    """
    reduced_embeddings = []
    for full_embedding in embeddings:
        reduced_embeddings.append(full_embedding.mean(axis=0))
    return np.array(reduced_embeddings)

def normalise(embeddings: np.ndarray) -> np.ndarray:
    """
    Normalises embeddings by subtracting the mean across the zeroth axis.

    Parameters:
    - embeddings (np.ndarray): An array of embeddings.
    
    Returns:
    - np.ndarray: The normalised array of embeddings.
    """
    return embeddings - np.mean(embeddings, axis=0)

def cosine_distance(reduced_embeddings: dict, variants_df_w_sequence: pd.DataFrame, gene_to_sequence: dict) -> pd.DataFrame:
    """
    Calculates the cosine distance between canonical and variant sequence embeddings for each entry in a DataFrame.

    This function uses precomputed embeddings for canonical and variant sequences to compute the cosine distance, then updates the DataFrame with these values.

    Parameters:
    - reduced_embeddings (dict): Dictionary mapping sequences to their reduced embeddings.
    - variants_df_w_sequence (pd.DataFrame): DataFrame containing variant sequences.
    - gene_to_sequence (dict): Dictionary mapping genes to canonical sequences.
    
    Returns:
    - pd.DataFrame: Updated DataFrame with a new 'ESM-1b Cosine Distance' column.
    """
    gene_to_embedding = {}
    for key, value in gene_to_sequence.items():
        if value in reduced_embeddings:
            gene_to_embedding[key] = reduced_embeddings[value]
    
    
    variants_df = variants_df_w_sequence.copy()

    for index, row in variants_df.iterrows():
        if row['Sequence'] in reduced_embeddings:
            canon_embedding = gene_to_embedding.get(row['Gene'])
            if canon_embedding is not None:
                variant_embedding = reduced_embeddings[row['Sequence']]
                cosine_distance = cosine(canon_embedding, variant_embedding)
                variants_df.at[index, 'ESM-1b Cosine Distance'] = cosine_distance
            else:
                variants_df.drop(index, inplace=True)
        else:
            variants_df.drop(index, inplace=True)

    variants_df.drop('Sequence', axis=1, inplace=True)

    return variants_df

def hdf5_to_dict(file_path: str) -> dict:
    """
    Reads embeddings from an HDF5 file into a dictionary.

    Parameters:
    - file_path (str): Path to the HDF5 file.
    
    Returns:
    - dict: A dictionary where keys are sequence identifiers and values are arrays of embeddings.
    """
    embeddings = {}
    with h5py.File(file_path, 'r') as hf:
        for sequence in hf.keys():
          embeddings[sequence] = hf[sequence][:]
    return embeddings



## ==================== helper function ===============================
def _mutate(protein_sequence: str, mutation: str) -> str:
    """
    Mutates a given protein sequence at a specified position to a new amino acid.

    Parameters:
    - protein_sequence (str): The original amino acid sequence of the protein.
    - mutation (str): Mutation description in the format 'W102M', where 'W' is the wild-type amino acid, '102' is the position, and 'M' is the mutant.
    
    Returns:
    - str: The mutated protein sequence.
    
    Raises:
    - ValueError: If the mutation position is invalid or the mutation specification does not match the sequence.
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
