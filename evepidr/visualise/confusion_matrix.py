def find_true_pos_rate(true_values: list, model_scores: list, lower_threshold: float, upper_threshold: float) -> float:
    true_values_array = np.array(true_values)
    model_scores_array = np.array(model_scores)
    mask_high = model_scores_array > upper_threshold
    mask_low = model_scores_array < lower_threshold
    model_scores_binary = np.where(mask_low, 0, np.where(mask_high, 1, np.nan)) # nan for ambiguous
    model_tpr = np.sum((model_scores_binary == 1) & (true_values_array == 1)) / np.sum(true_values_array == 1)
    return model_tpr

def find_threshold_for_tpr(true_values: list, model_scores: list, target_tpr: float) -> float:
    fpr, tpr, thresholds = roc_curve(true_values, model_scores)
    # Find the threshold nearest to the target TPR
    index = np.nanargmin(np.abs(tpr - target_tpr))
    return thresholds[index]

def print_confusion_matrix(true_values: list, model_scores: list, model_label: str, upper_threshold: float, lower_threshold: float=None) -> None:
    true_values_array = np.array(true_values)
    model_scores_array = np.array(model_scores)

    if lower_threshold is not None:
        predictions_array = np.where(model_scores_array > upper_threshold, 1,
                                     np.where(model_scores_array < lower_threshold, 0, np.nan))
    else:
        predictions_array = model_scores_array >= upper_threshold

    # Filter out NaN values from predictions for confusion matrix calculation
    valid_indices = ~np.isnan(predictions_array)
    if valid_indices.any():  # Check if there are any non-NaN predictions
        conf_matrix = confusion_matrix(true_values_array[valid_indices], predictions_array[valid_indices], labels=[0, 1])
        tn, fp, fn, tp = conf_matrix.ravel()

        # Calculating rates
        tpr = tp / (tp + fn) if (tp + fn) > 0 else 0  # True Positive Rate
        fpr = fp / (fp + tn) if (fp + tn) > 0 else 0  # False Positive Rate
        tnr = tn / (tn + fp) if (tn + fp) > 0 else 0  # True Negative Rate
        fnr = fn / (fn + tp) if (fn + tp) > 0 else 0  # False Negative Rate
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0  # Precision
        se = sqrt((precision * (1 - precision)) / (tp + fp))
        z_score = 1.96  # Z-score for 95% confidence
        ci = z_score * se

        print(f'Confusion Matrix for {model_label} (Threshold={upper_threshold}, {lower_threshold}):')
        print(f'True Positives (TP): {tp}               True Positive Rate (Recall, Sensitivity): {tpr}')
        print(f'True Negatives (TN): {tn}               True Negative Rate (Specificity): {tnr}')
        print(f'False Positives (FP): {fp}              False Positive Rate: {fpr}')
        print(f'False Negatives (FN): {fn}              False Negative Rate: {fnr}')
        print(f'Precision: {precision} ± {ci}')
        print()
    else:
        print(f'No valid predictions available for {model_label} due to all being NaN.')
