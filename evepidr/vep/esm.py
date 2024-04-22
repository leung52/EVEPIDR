import torch
import esm
import h5py
import pandas as pd

from evepidr.vep.lm_utils import batchify, save_to_h5py

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
    For ESM_1b only.
    """
    # TO DO: Check sequences are all same length and <= 1024

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
        save_to_h5py(token_representations, sequence_batch, file_path)

def zero_shot_variant_prediction(model_alpha: tuple, sequences: list, file_path: str) -> None:
    """
    Code for ESM_1v only. Not tested on ESM-1b.
    """
    pass
