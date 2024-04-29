from evepidr.visualise.roc_and_confusion_matrix import *


pathogenicities_df = pd.read_csv('report/data/predicted_pathogenicities.csv')
pathogenicities_dict = get_scores_and_true_values(pathogenicities_df)

for region, (true_values, models_scores, labels) in pathogenicities_dict.items():
    title = region + ' located substitutions'

    # ROC plot
    plot_roc_curve_and_save(true_values, models_scores, labels, title)

    # TPR of AM and threshold values for ESM-1b and ESM-1v
    am_lower_threshold = 0.34
    am_upper_threshold = 0.564
    am_tpr = find_true_pos_rate(true_values, models_scores[0], am_lower_threshold, am_upper_threshold)
    esm_1b_threshold = find_threshold_for_tpr(true_values, models_scores[1], am_tpr)
    esm_1v_threshold = find_threshold_for_tpr(true_values, models_scores[2], am_tpr)

    print(title)
    print_confusion_matrix(true_values, models_scores[0], 'AlphaMissense', am_upper_threshold, am_lower_threshold)
    print_confusion_matrix(true_values, models_scores[1], 'ESM-1b', esm_1b_threshold)
    print_confusion_matrix(true_values, models_scores[2], 'ESM-1v', esm_1v_threshold)
