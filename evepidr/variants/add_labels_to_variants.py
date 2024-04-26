import pandas as pd


def label_if_substitution_in_idr(variant_df: pd.DataFrame, idr_regions: dict) -> pd.DataFrame:
    """
    """
    variant_df['Substitution Region'] = variant_df.apply(lambda row: 'IDR' if any(start <= int(row['AA Substitution'][1:-1]) <= end for start, end in idr_regions.get(row['Gene'], [])) else 'Folded', axis=1)
    return variant_df
