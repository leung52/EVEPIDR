import pandas as pd

from evepidr.vep.lm_utils import *
from evepidr.vep.esm import *
from evepidr.vep.alphamissense import *

variants_df = pd.read_csv("clinvar_data_patho_benign.csv", index=False)

# ESM-1b predictions
gene_to_sequence = fasta_to_dict("asd_linked_idps.fasta")
sequences_for_esm1b, sequence_to_gene_and_mutation = prepare_sequences_for_lms(variants_df, protein_sequnces)
model_alpha = load_esm_1b()
embeddings_file = 'esm_1b_per_residue_embeddings.py'
for sequences in sequences_for_esm1b:
  embed(model_alpha, sequences, embeddings_file)
embeddings_dict = hdf5_to_dict(embeddings_file)
reduced_embeddings = reduce(embeddings_dict.values())
