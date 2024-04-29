import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_predictive_scores_on_histogram(pathogenicities_df: pd.DataFrame, models_columns: ls, file_name: str):
    """
    """
    # Loop through each substitution region
    for region in pathogenicities_df['Substitution Region'].unique():
        fig, axes = plt.subplots(1, len(models), figsize=(15, 5), sharey=True)
        fig.suptitle(f'Predictive Scores for {region} Region')

        for ax, model_col in zip(axes, models_columns):
            # Filter the data for the current region and model
            region_data = pathogenicities_df[pathogenicities_df['Substitution Region'] == region]
            path_scores = region_data[region_data['Pathogenicity'] == 'Pathogenic'][model_col]
            benign_scores = region_data[region_data['Pathogenicity'] == 'Benign'][model_col]
            
            # Define bins based on the combined range of pathogenic and benign scores
            combined_scores = np.concatenate((path_scores, benign_scores))
            bins = np.linspace(np.min(combined_scores), np.max(combined_scores), 60)

            ax.hist(path_scores, bins=bins, alpha=0.6, label='Pathogenic')
            ax.hist(benign_scores, bins=bins, alpha=0.6, label='Benign')
            
            ax.set_title(model_col)
            ax.set_xlabel('Pathogenicity Score')
            ax.set_ylabel('Count')
            ax.legend()
            ax.tick_params(axis='x', rotation=45)

        plt.rcParams.update({'font.size': 12})
        plt.tight_layout(rect=[0, 0, 1, 1])  # Adjust subplots to fit the figure title
        plt.savefig(file_name)
        plt.close()
