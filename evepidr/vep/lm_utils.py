import numpy as np

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
