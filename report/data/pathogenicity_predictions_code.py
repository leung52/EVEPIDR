import pandas as pd

from evepidr.vep.lm_utils import *
from evepidr.vep.esm import *
from evepidr.vep.alphamissense import *


variants_df = pd.read_csv("report/data/clinvar_data_patho_benign.csv")
gene_to_sequence = fasta_to_dict("report/data/asd_linked_idps.fasta")

## AlphaMissense Pathogenicities
# need to download and add AlphaMissense_aa_substitutions.tsv to folder
# https://zenodo.org/records/10813168
am_tsv_file_path = 'report/data/AlphaMissense_aa_substitutions.tsv'
variants_df = alpha_missense_scores(variants_df, am_tsv_file_path)

## ESM-1b Cosine Distances
variants_df = variants_df = add_sequences_to_variants_df(variants_df, gene_to_sequence)
embeddings_dict = embed_with_esm_1b(variants_df, gene_to_sequence, 'report/data/esm_1b_per_residue_embeddings.hdf5')
reduced_embeddings = reduce(list(embeddings_dict.values()))
normalised_embeddings = normalise(reduced_embeddings)
# Save reduced and normalised embeddings in csv
rows = []
for sequence, reduced, normalised in zip(embeddings_dict, reduced_embeddings, normalised_embeddings):
    rows.append({'Sequence': sequence, 'Type': 'Reduced', **{f'Dim_{i+1}': value for i, value in enumerate(reduced)}})
    rows.append({'Sequence': sequence, 'Type': 'Normalised', **{f'Dim_{i+1}': value for i, value in enumerate(normalised)}})
embeddings_df = pd.DataFrame(rows)
embeddings_df.to_csv('report/data/esm_1b_reduced_normalised_embeddings.csv', index=False)
normalised_dict = dict(zip(embeddings_dict, normalised_embeddings))
# Calculate cosine distances
variants_df = cosine_distance(normalised_dict, variants_df, gene_to_sequence)

## ESM-1v WildType Marginals
variants_df = wt_marginals_with_esm_1v(variants_df, gene_to_sequence)

## Save pathogenicity predictions in csv
variants_df.to_csv('report/data/predicted_pathogenicities.csv', index=False)
