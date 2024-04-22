import torch
import esm
import h5py
import pandas as pd

from evepidr.vep.utils import batchify

# TO DO: write DocStrings
# TO DO: check sequences are <= 1024 in length

# https://github.com/facebookresearch/esm

def load_esm_1b() -> tuple:
  """
  """
  model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
  return (model, alphabet)

def load_esm_1v() -> tuple:
  """
  """
  # Load the first ESM-1v model of five models ensemble
  model, alphabet = esm.pretrained.esm1v_t33_650M_UR90S_1()
  return (model, alphabet)

def embed(model_alpha: tuple, sequences: list, file_path: str) -> None:
  """
  """
  model = model_alpha[0]
  alphabet = model_alpha[1]
  
  # Put the model in evaluation mode
  model.eval()  

  # Prepare batch for processing
  batch_converter = alphabet.get_batch_converter() # Function from esm
  sequences_as_batches = batchify(sequences, 10) # Helper function; Splits up sequences into manageable chunks

  for sequence_batch in sequences_as_batches:
    # Convert enumerate object to list
    enumerated_sequences = list(enumerate(sequence_batch)) 
    batch_labels, batch_strs, batch_tokens = batch_converter(enumerated_sequences)

    # Compute embeddings
    with torch.no_grad():
      results = model(batch_tokens, repr_layers=[33])  # Layer 33 most commonly used
      token_representations = results["representations"][33]

    # Save embeddings to h5py file
    with h5py.File(file_path, 'a') as hf:
      for sequence, embedding in zip(sequence_batch, token_representations):
        embedding_data = np.array(embedding)
        if sequence in hf:
          # If the dataset already exists, overwrite the data
          hf[sequence][...] = embedding_data
        else:
          # If the dataset doesn't exist, create a new one
          hf.create_dataset(sequence, data=embedding_data, maxshape=(None,))
