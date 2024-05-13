import pandas as pd


def label_if_substitution_in_idr(variant_df: pd.DataFrame, idr_regions: dict) -> pd.DataFrame:
    """
    Labels each variant in a DataFrame as 'IDR', 'Folded', or 'Canon' based on amino acid substitution locations relative to predefined intrinsically disordered regions (IDRs).

    This function iterates over a DataFrame containing gene variants and their amino acid substitutions. It checks if the numeric position of each substitution falls within any IDR ranges specified for its gene. The function assigns 'IDR' if the substitution is within an IDR, 'Folded' if outside an IDR, and 'Canon' if there is no substitution data.
    
    Parameters:
    - variant_df (pd.DataFrame): A DataFrame with at least the columns 'Gene' and 'AA Substitution'. 'AA Substitution' should have substitutions in a format where the position is numeric (e.g., 'K103N').
    - idr_regions (dict): A dictionary mapping gene names to lists of tuples, each tuple representing the start and end positions of an IDR.
    
    Returns:
    - pd.DataFrame: The modified DataFrame with an additional 'Substitution Region' column indicating whether each substitution is in an 'IDR', 'Folded', or 'Canon' region.
    """
    def label_region(row):
        if pd.isna(row['AA Substitution']):
            return 'Canon'
        else:
            if any(start <= int(row['AA Substitution'][1:-1]) <= end for start, end in idr_regions.get(row['Gene'], [])):
                return 'IDR'
            else:
                return 'Folded'

    variant_df['Substitution Region'] = variant_df.apply(label_region, axis=1)
    return variant_df
