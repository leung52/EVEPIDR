# uniprot_id - UniProtKB accession number of the protein in which the variant induces a single amino-acid substitution (UniProt release 2021_02).
# protein_variant - Amino acid change induced by the alternative allele, in the format <Reference amino acid><POS_aa><Alternative amino acid> (e.g. V2L). POS_aa is the 1-based position of the residue within the protein amino acid sequence.
# am_pathogenicity - Calibrated AlphaMissense pathogenicity scores (ranging between 0 and 1), which can be interpreted as the predicted probability of a variant being clinically pathogenic.
# am_class - Classification of the protein_variant into one of three discrete categories: 'likely_benign', 'likely_pathogenic', or 'ambiguous'. These are derived using the following thresholds: 'likely_benign' if alphamissense_pathogenicity < 0.34; 'likely_pathogenic' if alphamissense_pathogenicity > 0.564; and 'ambiguous' otherwise.


# AlphaMissense file is too large for GitHub, ask user to download it locally to use

import pandas as pd



def alpha_missense_scores(annotations_df: pd.DataFrame, alpha_missense_tsv_file: str, file_path: str) -> pd.DataFrame:
    """
    Filter rows from a large TSV file where 'uniprot_id' and 'protein_variant' match those in the input DataFrame,
    processing the file in chunks to avoid memory issues.

    Parameters:
    alpha_missense_tsv_file (str): The file path to the TSV file.
    annotations_df (pd.DataFrame): A DataFrame containing 'uniprot_id' and 'protein_variant' columns.

    Returns:
    pd.DataFrame: A DataFrame with rows from the TSV file that match the input DataFrame's criteria.
    """
    # Initialize an empty DataFrame to store the filtered results
    filtered_alpha_miss_df = pd.DataFrame()

    # Create a set for faster lookup
    lookup_set = set(zip(annotations_df['UniProt ID'], annotations_df['AA Substitution']))

    # Define chunk size
    chunk_size = 10**5  # Adjust based on your system's memory

    # Process the TSV file in chunks
    for chunk in pd.read_csv(alpha_missense_tsv_file, sep='\t', chunksize=chunk_size, skiprows=3):
        # Filter the chunk based on the lookup set
        chunk_filtered = chunk[chunk.apply(lambda x: (x['uniprot_id'], x['protein_variant']) in lookup_set, axis=1)]
        # Append the filtered chunk to the results DataFrame
        filtered_alpha_miss_df = pd.concat([filtered_alpha_miss_df, chunk_filtered], ignore_index=True)

    # Save the filtered results to a new CSV file
    filtered_alpha_miss_df.to_csv(file_path, index=False)

    return filtered_alpha_miss_df
