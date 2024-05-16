import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc, confusion_matrix
from scipy.stats import norm
from math import sqrt


def confusion_matrix(model_predictions_df: pd.DataFrame, true_value_column: str, control_model_column: str, model_columns: list, upper_threshold: float, lower_threshold: float, save_txt_file: str, ci: int=95) -> None:
    """
    Generates confusion matrices for multiple models compared to a control model's thresholds and writes the results to a text file.

    This function calculates the confusion matrix for each model's predictions against true binary outcomes. It uses fixed thresholds for the control model and dynamically computes thresholds for other models to match the True Positive Rate (TPR) of the control model. The results, including metrics like TPR, FPR, Precision, and their confidence intervals, are written to a specified text file.
    
    Parameters:
    - model_predictions_df (pd.DataFrame): DataFrame containing the true values and model predictions.
    - true_value_column (str): The name of the column in `model_predictions_df` that contains the true binary labels.
    - control_model_column (str): The name of the column that contains the control model's predictions.
    - model_columns (list): A list of column names representing other models whose predictions are to be evaluated.
    - upper_threshold (float): The upper threshold for classifying predictions from the control model as positive.
    - lower_threshold (float): The lower threshold for classifying predictions from the control model as negative.
    - save_txt_file (str): Path to the text file where the output will be written.
    - ci (int, optional): The confidence interval for precision calculation in percentage (default is 95).
    
    Returns:
    - None
    """
    true_values = [1 if x == 'Pathogenic' else 0 for x in model_predictions_df[true_value_column]]
    true_values = np.array(true_values)

    control_model = np.array(model_predictions_df[control_model_column])
    mask_high = control_model > upper_threshold
    mask_low = control_model < lower_threshold
    control_model_scores_binary = np.where(mask_low, 0, np.where(mask_high, 1, np.nan)) # nan for ambiguous
    target_tpr = np.sum((control_model_scores_binary == 1) & (true_values == 1)) / np.sum(true_values == 1)

    for model in model_columns:
        model_scores = np.array(model_predictions_df[model])
        if model == control_model_columns:
            predictions_array = np.where(model_scores_array > upper_threshold, 1,
                                         np.where(model_scores_array < lower_threshold, 0, np.nan))
            str_threshold = str(upper_threshold) + ', ' str(lower_threshold)
        else:
            fpr, tpr, thresholds = roc_curve(true_values, model_scores)
            index = np.nanargmin(np.abs(tpr - target_tpr))
            threshold = thresholds[index]
            predictions_array = model_scores_array >= threshold
            str_threshold = str(threshold)
            
        valid_indices = ~np.isnan(predictions_array)
        
        conf_matrix = confusion_matrix(true_values[valid_indices], predictions_array[valid_indices], labels=[0, 1])
        tn, fp, fn, tp = conf_matrix.ravel()
        # Calculating rates
        tpr = tp / (tp + fn) if (tp + fn) > 0 else 0  # True Positive Rate
        fpr = fp / (fp + tn) if (fp + tn) > 0 else 0  # False Positive Rate
        tnr = tn / (tn + fp) if (tn + fp) > 0 else 0  # True Negative Rate
        fnr = fn / (fn + tp) if (fn + tp) > 0 else 0  # False Negative Rate
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0  # Precision
        se = sqrt((precision * (1 - precision)) / (tp + fp))
        z_score = norm.ppf(1 - (1 - ci/100) / 2)
        ci = z_score * se
        
        with open(save_txt_file, 'a') as file:
            file.write(f'Confusion Matrix for {model} (Threshold={str_threshold}):\n')
            file.write(f'True Positives (TP): {tp}\n')
            file.write(f'True Negatives (TN): {tn}\n')
            file.write(f'False Positives (FP): {fp}\n')
            file.write(f'False Negatives (FN): {fn}\n\n')
            file.write(f'True Positive Rate (Recall, Sensitivity): {tpr}\n')
            file.write(f'True Negative Rate (Specificity): {tnr}\n')
            file.write(f'False Positive Rate: {fpr}\n')
            file.write(f'False Negative Rate: {fnr}\n\n')
            file.write(f'Precision: {precision}\n\n')
            file.write('----------------------------------------------\n')
