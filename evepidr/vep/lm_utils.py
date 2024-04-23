import numpy as np
import pandas as pd


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
