import pandas as pd


def label_if_substitution_in_idr(variant_df: pd.DataFrame, idr_regions: dict) -> pd.DataFrame:
    """
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
