import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from scipy.stats import norm


def plot_roc_curve_and_save(true_values: list, models_scores: list, model_labels: list, plot_title: str, save_file: str, n_bootstraps:int=10000, ci:int=95) -> None:
    """
    Plots the ROC curve for each model's scores against true binary classification outcomes and saves the plot.
    
    This function takes arrays of true classification results and predicted scores from multiple models, plots the ROC curve for each model, and displays the AUC (Area Under Curve) score with a confidence interval. It also prints the AUC and its confidence interval for each model to the console. The plot is saved to a file.
    
    Parameters:
    - true_values (list): A list of true binary labels for the classification task.
    - models_scores (list of lists): A list containing arrays of scores from different models. Each sublist corresponds to one model.
    - model_labels (list): A list of strings that are labels for each model, used in the legend.
    - plot_title (str): The title for the plot.
    - save_file (str): File path for saved plot.
    - n_bootstraps (int, optional): The number of bootstrap samples to use when estimating the confidence interval. Default is 10,000.
    - ci (int, optional): The confidence level for the AUC confidence interval as a percentage. Default is 95.
    
    Returns:
    - None
    """
    plt.rcParams.update({'font.size': 18})
    plt.figure(figsize=(10, 8))

    true_values = np.array(true_values)  # Convert true_values to a numpy array

    for model_scores, label in zip(models_scores, model_labels):
        model_scores = np.array(model_scores)  # Convert model_scores to a numpy array
        fpr, tpr, _ = roc_curve(true_values, model_scores)
        roc_auc = auc(fpr, tpr)

        # Calculate the Hanley and McNeil standard error for AUC
        n1 = sum(true_values)
        n2 = len(true_values) - n1
        q1 = roc_auc / (2 - roc_auc)
        q2 = 2 * roc_auc**2 / (1 + roc_auc)
        SE = np.sqrt((roc_auc * (1 - roc_auc) + (n1 - 1) * (q1 - roc_auc**2) + (n2 - 1) * (q2 - roc_auc**2)) / (n1 * n2))

        # Calculate the confidence interval
        z_score = norm.ppf(1 - (1 - ci/100) / 2)
        confidence_lower = roc_auc - z_score * SE
        confidence_upper = roc_auc + z_score * SE

        plt.plot(fpr, tpr, linewidth=2, label=f'{label} (AUC = {roc_auc:.3f}, {ci}% CI = [{confidence_lower:.3f}, {confidence_upper:.3f}])')

    plt.plot([0, 1], [0, 1], 'k--', linewidth=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(plot_title + ' - ROC')
    plt.legend(loc="lower right")
    plt.grid(False)
    plt.savefig(save_file)
    plt.show()
    plt.close()

def get_scores_and_true_values(pathogenicities_df: pd.DataFrame) -> dict:
    """
    Extracts true binary values and corresponding scores for different models from a DataFrame grouped by substitution region.
    
    This function processes a DataFrame that contains pathogenicity data along with various model scores for different regions of substitution. It groups the data by substitution region, converts the pathogenicity to binary format, and extracts scores for different models. It returns a dictionary with each region as a key and tuples of true values, model scores, and labels as values.
    
    Parameters:
    - pathogenicities_df (pd.DataFrame): A DataFrame containing the columns 'Substitution Region', 'Pathogenicity', and various columns for model scores.
    
    Returns:
    - dict: A dictionary where keys are region type and values are tuples containing lists of true values, lists of lists of model scores, and list of model labels.
    """
    scored_and_true_values = {}

    for region_str, grouped_df in pathogenicities_df.groupby('Substitution Region'):
        true_values = grouped_df['Pathogenicity'].to_list()
        true_values = [1 if item == 'Pathogenic' else 0 for item in true_values]
        alphamissense = grouped_df['AM Pathogenicity'].to_list()
        esm1b_cosine = grouped_df['ESM-1b Cosine Distance'].to_list()
        esm1v_marginal = grouped_df['ESM-1v WT Marginals'].to_list()

        model_scores = [alphamissense, esm1b_cosine, esm1v_marginal]
        labels = ['AM Pathogenicity', 'ESM-1b Cosine Distance', 'ESM-1v WT Marginals']

        if region_str == 'Folded':
            region_str = 'Folded region'

        scored_and_true_values[region_str] = (true_values, model_scores, labels)
    return scored_and_true_values
