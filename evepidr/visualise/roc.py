import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc


def plot_roc_curve(true_values: list, model_scores: list, labels: list, plot_title: str) -> plt.plot:
    plt.figure(figsize=(10, 8))
    for y_score, label in zip(model_scores, labels):
        fpr, tpr, _ = roc_curve(true_values, y_score)
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, linewidth=2, label=f'{label} (AUC = {roc_auc:.2f})')

    plt.plot([0, 1], [0, 1], 'k--', linewidth=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(plot_title + 'Receiver Operating Characteristic')
    plt.legend(loc="lower right")
    plt.grid(True)
    plt.show()
