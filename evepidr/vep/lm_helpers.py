import h5py
import numpy as np


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

def save_to_h5py(embeddings: np.ndarray, sequences: list, file_path: str) -> None:
    """
    """
    with h5py.File(file_path, 'a') as hf:
        for sequence, embedding in zip(sequences, embeddings):
          embedding_data = np.array(embedding)
          if sequence in hf:
              # If the dataset already exists, overwrite the data
              hf[sequence][...] = embedding_data
          else:
              # If the dataset doesn't exist, create a new one
              hf.create_dataset(sequence, data=embedding_data, maxshape=(None,))
