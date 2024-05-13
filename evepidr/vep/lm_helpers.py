import h5py
import numpy as np


def batchify(sequences: list, batch_size: int) -> list:
    """
    Splits a list of items into batches of a specified size.

    Parameters:
    - sequences (list): The list of items to be batched.
    - batch_size (int): The maximum number of items per batch.
    
    Returns:
    - list: A list of batches, where each batch is a list of items up to the specified batch size.
    """
    batches = [sequences[i:i + batch_size] for i in range(0, len(sequences), batch_size)]
    return batches

def save_to_h5py(embeddings: np.ndarray, sequences: list, output_file: str) -> None:
    """
    Saves embeddings to an HDF5 file, using the sequence as the dataset name.

    Parameters:
    - embeddings (np.ndarray): Array of embeddings to be saved.
    - sequences (list): List of sequence identifiers corresponding to the embeddings.
    - output_file (str): Path to the output HDF5 file where embeddings will be stored.
    
    Returns:
    - None
    """
    # Ensure embedding is on CPU before converting to Numpy array
    embedding_data = embeddings.cpu().numpy()
    
    # Continue with saving the data to h5py file
    with h5py.File(output_file, 'a') as h5f:
        for seq, embed in zip(sequences, embedding_data):
            if seq in h5f:
                del h5f[seq]  # Delete the existing dataset
            h5f.create_dataset(seq, data=embed, compression='gzip')
