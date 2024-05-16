import pandas as pd

from evepidr.visualise.histogram import *
from evepidr.visualise.roc_and_confusion_matrix import *


pathogenicities_df = pd.read_csv('report/data/predicted_pathogenicities.csv')
pathogenicities_df['ESM-1v WT Marginals'] = -pathogenicities_df['ESM-1v WT Marginals']

model_columns = ['AM Pathogenicity', 'ESM-1b Cosine Distance', 'ESM-1v WT Marginals']

idr_df = df[df['Substitution Region'] == 'IDR']
folded_df = df[df['Substitution Region'] == 'Folded']


# ROC
plot_roc_curve_and_save(idr_df, 'Pathogenicity', 'IDR variants ROC', 'roc_idr.png')
plot_roc_curve_and_save(folded_df, 'Pathogenicity', 'Folded variants ROC', 'roc_folded.png')

# Confusion matrix
confusion_matrix(idr_df, 'Pathogenicity', 'AM Pathogenicity', model_columns, 0.34, 0.564, 'confusion_matrix.txt')
confusion_matrix(folded_df, 'Pathogenicity', 'AM Pathogenicity', model_columns, 0.34, 0.564, 'confusion_matrix.txt')

# Histogram
plot_predictive_scores_on_histogram(pathogenicities_df, model_columns, 'histogram')
