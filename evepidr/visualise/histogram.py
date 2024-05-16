import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_predictive_scores_on_histogram(pathogenicities_df: pd.DataFrame, models: list, save_file: str) -> None:
    """
    Plots histogram comparisons of pathogenic and benign predictive scores for different models, separated by substitution regions.

    This function takes a DataFrame containing pathogenicity scores from multiple models and plots histograms for pathogenic versus benign scores, grouped by substitution regions. Each region's data is plotted in separate figures with histograms for each model, illustrating the distribution of scores. The figures are saved as separate image files.
    
    Parameters:
    - pathogenicities_df (pd.DataFrame): A DataFrame containing pathogenicity scores and labels, with columns for 'Substitution Region' and 'Pathogenicity'.
    - models (list): A list of column names from the DataFrame representing the models whose scores are to be plotted.
    - save_file (str): Base path and filename prefix where the histogram plots will be saved. Each plot will be saved with an additional identifier for the region.
    
    Returns:
    - None
    """
    plt.rcParams.update({'font.size': 12})
    for region in pathogenicities_df['Substitution Region'].unique():
        fig, axes = plt.subplots(1, len(models), figsize=(15, 5), sharey=False)  # Set sharey to False
        fig.suptitle(f'Predictive Scores for {region} Region')

        for ax, model_col in zip(axes, models):
            # Filter the data for the current region and model
            region_data = pathogenicities_df[pathogenicities_df['Substitution Region'] == region]
            path_scores = region_data[region_data['Pathogenicity'] == 'Pathogenic'][model_col]
            benign_scores = region_data[region_data['Pathogenicity'] == 'Benign'][model_col]
            
            # Define bins based on the combined range of pathogenic and benign scores
            combined_scores = np.concatenate((path_scores, benign_scores))
            bins = np.linspace(np.min(combined_scores), np.max(combined_scores), 50)

            ax.hist(path_scores, weights=np.ones(len(path_scores)) / len(path_scores), bins=bins, alpha=0.6, label='Pathogenic')
            ax.hist(benign_scores, weights=np.ones(len(benign_scores)) / len(benign_scores),  bins=bins, alpha=0.6, label='Benign')
            
            ax.set_title(model_col)
            ax.set_xlabel('Pathogenicity Score')
            ax.set_ylabel('Percentage')
            ax.legend()
            ax.tick_params(axis='x', rotation=45)
            ax.autoscale(enable=True, axis='y', tight=True)  # Auto-scale the y-axis
        plt.tight_layout(rect=[0, 0, 1, 1])  # Adjust subplots to fit the figure title
        plt.savefig(save_file + '_' + region + '.png')  # Save each region to a separate file
        plt.close()
