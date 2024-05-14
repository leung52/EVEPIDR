import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from scipy.stats import norm


def plot_roc_curve_and_save(model_predictions_df: pd.DataFrame, true_value_column: str, model_columns: list, plot_title: str, save_file: str, ci: int=95) -> None:
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

    true_values = [1 if x == 'Pathogenic' else 0 for x in model_predictions_df[true_value_column]]
    true_values = np.array(true_values)

    for model in model_columns:
        model_scores = np.array(model_predictions_df[model])
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

        plt.plot(fpr, tpr, linewidth=2, label=f'{model} (AUC = {roc_auc:.3f}, {ci}% CI = [{confidence_lower:.3f}, {confidence_upper:.3f}])')

    plt.plot([0, 1], [0, 1], 'k--', linewidth=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(plot_title)
    plt.legend(loc="lower right")
    plt.grid(False)
    plt.savefig(save_file)
    plt.show()
    plt.close()
