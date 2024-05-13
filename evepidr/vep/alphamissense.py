import pandas as pd

## Structure of AlphaMissense_aa_substitutions.tsv
# uniprot_id - UniProtKB accession number of the protein in which the variant induces a single amino-acid substitution (UniProt release 2021_02).
# protein_variant - Amino acid change induced by the alternative allele, in the format <Reference amino acid><POS_aa><Alternative amino acid> (e.g. V2L). POS_aa is the 1-based position of the residue within the protein amino acid sequence.
# am_pathogenicity - Calibrated AlphaMissense pathogenicity scores (ranging between 0 and 1), which can be interpreted as the predicted probability of a variant being clinically pathogenic.
# am_class - Classification of the protein_variant into one of three discrete categories: 'likely_benign', 'likely_pathogenic', or 'ambiguous'. These are derived using the following thresholds: 'likely_benign' if alphamissense_pathogenicity < 0.34; 'likely_pathogenic' if alphamissense_pathogenicity > 0.564; and 'ambiguous' otherwise.

# source: https://zenodo.org/records/10813168
# AlphaMissense file is too large for GitHub, ask user to download it locally to use

def alpha_missense_scores(variants_df: pd.DataFrame, am_tsv_file_path: str) -> pd.DataFrame:
    """
    Integrates AlphaFold Missense (AM) predictions from a TSV file with a DataFrame of variant data based on matching UniProt IDs and amino acid substitutions.

    This function reads a TSV file containing AlphaFold Missense predictions and filters it to include only the predictions that match the UniProt IDs and amino acid substitutions provided in the input DataFrame. It merges these predictions into the original DataFrame, adding a new column with the AlphaFold pathogenicity scores.
    
    Parameters:
    - variants_df (pd.DataFrame): A DataFrame containing at least the columns 'UniProt ID' and 'AA Substitution', representing variants to be scored.
    - am_tsv_file_path (str): The file path to the TSV file containing AlphaFold Missense predictions. The file is expected to have columns for 'uniprot_id', 'protein_variant', and 'am_pathogenicity', and may include a header and other irrelevant data.
    
    Returns:
    - pd.DataFrame: The original DataFrame enriched with a new column 'AM Pathogenicity', containing the pathogenicity scores from AlphaFold corresponding to each variant.
    """
    am_predictions_df = pd.DataFrame()

    # Create a set for faster lookup
    lookup_set = set(zip(variants_df['UniProt ID'], variants_df['AA Substitution']))

    # Process the TSV file in chunks
    for chunk in pd.read_csv(am_tsv_file_path, sep='\t', chunksize=1000000, skiprows=3):
        # Filter the chunk based on the lookup set
        chunk_filtered = chunk[chunk.apply(lambda x: (x['uniprot_id'], x['protein_variant']) in lookup_set, axis=1)]
        # Append the filtered chunk to the results DataFrame
        am_predictions_df = pd.concat([am_predictions_df, chunk_filtered], ignore_index=True)
    
    am_predictions_df.rename(columns={'uniprot_id': 'UniProt ID', 'protein_variant': 'AA Substitution', 'am_pathogenicity': 'AM Pathogenicity'}, inplace=True)
    am_predictions_df.drop('am_class', axis=1, inplace=True)

    merged_df = pd.merge(variants_df, am_predictions_df, on=['UniProt ID', 'AA Substitution'])

    return merged_df
