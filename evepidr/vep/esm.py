import torch
import esm
import pandas as pd

from evepidr.vep.lm_helpers import batchify

def embed_with_esm_1b(variants_with_sequences_df: pd.DataFrame, gene_to_sequence: dict, output_file: str) -> dict:
    """
    """
    model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    model.eval()

    if torch.cuda.is_available():
        model = model.cuda()
        print("Transferred model to GPU")

    # Prepare batch for processing
    batch_converter = alphabet.get_batch_converter() # Function from esm

    all_sequences = []
    for gene, group in variants_with_sequences_df.groupby('Gene'):
        ls = list(group['Sequence'])
        ls.append(gene_to_sequence[gene])
        all_sequences.append(ls)

    for sequence_set in all_sequences:
        sequences_as_batches = batchify(sequence_set, 10) # Helper function; Splits up sequences into manageable chunks
        for sequence_batch in sequences_as_batches:
            # Convert enumerate object to list
            enumerated_sequences = list(enumerate(sequence_batch))
            batch_labels, batch_strs, batch_tokens = batch_converter(enumerated_sequences)

            # Compute embeddings
            with torch.no_grad():
                results = model(batch_tokens, repr_layers=[33])  # Layer 33 most commonly used
                token_representations = results["representations"][33]

            # Save embeddings to h5py file
            save_to_h5py(token_representations, sequence_batch, output_file)

    return hdf5_to_dict(output_file)

def wt_marginals_with_esm_1v(variants_df: pd.DataFrame, gene_to_sequence: dict) -> pd.DataFrame:
    """
    """
    model, alphabet = esm.pretrained.esm1v_t33_650M_UR90S_1()
    model.eval()

    if torch.cuda.is_available():
        model = model.cuda()
        print("Transferred model to GPU")

    batch_converter = alphabet.get_batch_converter()

    variants_groupby_genes = {group_name: group_data for group_name, group_data in variants_df.groupby('Gene')}

    for gene, sequence in gene_to_sequence.items():
        df = variants_groupby_genes[gene]
        data = [(1, sequence)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)

        with torch.no_grad():
          token_probs = torch.log_softmax(model(batch_tokens.cuda())["logits"], dim=-1)
          df['ESM-1v WT Marginals'] = df.apply(
              lambda row: _label_row(
                  row['AA Substitution'], 
                  sequence, 
                  token_probs, 
                  alphabet
                  ),
              axis=1
          )

    prediction_df = pd.concat(list(variants_groupby_genes.values()))
    prediction_df.reset_index(drop=True, inplace=True)

    return prediction_df

def _label_row(aa_substitution, sequence, token_probs, alphabet):
    wt, idx, mt = aa_substitution[0], int(aa_substitution[1:-1]) - 1, aa_substitution[-1]
    assert sequence[idx] == wt, "The listed wildtype does not match the provided sequence"

    wt_encoded, mt_encoded = alphabet.get_idx(wt), alphabet.get_idx(mt)

    # add 1 for BOS
    score = token_probs[0, 1 + idx, mt_encoded] - token_probs[0, 1 + idx, wt_encoded]
    return score.item()
