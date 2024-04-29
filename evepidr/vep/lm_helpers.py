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

def save_to_h5py(embeddings: np.ndarray, sequences: list, output_file: str) -> None:
    """
    """
    # Ensure embedding is on CPU before converting to Numpy array
    embedding_data = embeddings.cpu().numpy()
    
    # Continue with saving the data to h5py file
    with h5py.File(output_file, 'a') as h5f:
        for seq, embed in zip(sequences, embedding_data):
            if seq in h5f:
                del h5f[seq]  # Delete the existing dataset
            h5f.create_dataset(seq, data=embed, compression='gzip')
