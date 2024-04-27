import numpy as np
import pandas as pd


def fasta_to_dict(file_path: str) -> dict:
    """
    """
    sequences = {}
    with open(file_path, "r") as f:
        title = None
        sequence = ""
        for line in f:
            if line.startswith(">"):
                # If title is not None, store the previous sequence
                if title is not None:
                    sequences[title] = sequence
                # Extract title from FASTA header
                title = line.strip()[1:]
                sequence = ""
            else:
                # Append sequence line
                sequence += line.strip()
        # Store the last sequence
        if title is not None:
            sequences[title] = sequence
    return sequences

def prepare_sequences_for_lms(variants_df: pd.DataFrame, protein_sequnces: dict) -> list:
    """
    """
    ls = []
    for (gene, df_chunk) in variants_df.groupby('Gene'):
        variants = list(df_chunk['AA Substitution'])
        canon_sequence = protein_sequences.get(gene)
        variant_sequences = []
        for variant in variants:
            variant_sequences.append(_mutate(canon_sequence, variant))
        ls.append(variant_sequences)
    return ls

def reduce(embeddings: np.ndarray) -> np.ndarray:
    """
    
    """
    reduced_embeddings = []
    for full_embedding in all_embeddings:
        reduced_embeddings.append(full_embedding.mean(axis=0))
    return np.array(reduced_embeddings)

def normalise(embeddings: np.ndarray) -> np.ndarray:
    """
    """
    return embeddings - np.mean(embeddings, axis=0)

def cosine_distance(reduced_embeddings: dict, anno_df: pd.DataFrame, file_path: str) -> pd.DataFrame:
    """
    Calculates the cosine distance between sequence vectors in reduced_embeddings and the
    canonical sequence vector of the corresponding protein, excluding sequences without embeddings.
    
    Parameters:
    reduced_embeddings (dict): Dictionary with sequences as keys and 1D vectors as values.
    anno_df (pd.DataFrame): DataFrame with columns including 'Sequence', 'AA Substitution',
        'Clinical Significance', 'Gene', 'UniProt ID', 'Substitution Location'.
      
    Returns:
    pd.DataFrame: New DataFrame with columns 'Sequence' and 'Cosine Distance', excluding rows with missing values.
    """
    # Filter to get only canonical sequences
    canon_df = anno_df[anno_df['Clinical Significance'] == 'canon']
    
    # Initialize a list to store the cosine distances
    distances = []
    
    # Loop through the annotation DataFrame
    for index, row in anno_df.iterrows():
        if row['Clinical Significance'] != 'canon':
            # Find the canonical sequence for the current non-canonical sequence's protein
            canon_seq = canon_df[canon_df['UniProt ID'] == row['UniProt ID']]['Sequence'].values[0]
            
            # Retrieve the vector embeddings for the canonical and current sequence
            canon_vec = reduced_embeddings.get(canon_seq)
            current_vec = reduced_embeddings.get(row['Sequence'])
            
            # Calculate cosine distance if both vectors are found
            if canon_vec is not None and current_vec is not None:
                distance = cosine(canon_vec, current_vec)
                distances.append({'Sequence': row['Sequence'], 'Cosine Distance': distance})
    
    # Create a DataFrame from the distances list
    distance_df = pd.DataFrame(distances)
    
    # Remove rows with empty values in 'Cosine Distance'
    distance_df.dropna(subset=['Cosine Distance'], inplace=True)
    
    distance_df.to_csv(file_path, index=False)

def batchify(sequences: list, batch_size: int) -> list:
  """
  Split sequences into batches of a specified size.

  Parameters:
  sequences (list): List of sequences.
  batch_size (int): Size of each batch.

  Returns:
  List of batches, where each batch is a list of sequences.
  """
  batches = [sequences[i:i + batch_size] for i in range(0, len(sequences), batch_size)]
  return batches

def save_to_h5py(embeddings: np.ndarray, embedding_label: list, save_to_file_path: str) -> None:
    """
    """
    with h5py.File(file_path, 'a') as hf:
        for sequence, embedding in zip(sequences, reduced_embeddings):
          embedding_data = np.array(embedding)
          if sequence in hf:
              # If the dataset already exists, overwrite the data
              hf[sequence][...] = embedding_data
          else:
              # If the dataset doesn't exist, create a new one
              hf.create_dataset(sequence, data=embedding_data, maxshape=(None,))

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
