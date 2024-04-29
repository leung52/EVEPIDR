import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc


def plot_roc_curve(true_values: list, models_scores: list, model_labels: list, plot_title: str) -> plt.figure:
    plt.figure(figsize=(10, 8))
    for model_scores, label in zip(models_scores, model_labels):
        fpr, tpr, _ = roc_curve(true_values, model_scores)
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, linewidth=2, label=f'{label} (AUC = {roc_auc:.2f})')

    plt.plot([0, 1], [0, 1], 'k--', linewidth=2)
    plt.xlim([0.0, 1.05])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(plot_title + 'Receiver Operating Characteristic')
    plt.legend(loc="lower right")
    plt.grid(False)
    plt.show()
    return plt.figure
