import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import conorm
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import zscore

colors = ["Green", "k", "Red"]
nodes = [0, 0.5, 1.0]
cmap2 = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))

'''
TMM norm
'''
tissues = ['Kidney']
gene_list = ['S6K1', 'RPS6KB1', 'p70S6k', 'S6', '4EBP1', 'eIF4G', 'Akt', 'Ulk1', 'Lipin1', 'TFEB', 'Grb10', 'PKCa',
             'Sgk1',
             'Raptor', 'Rictor', 'mSin1', 'Pras40', 'mTOR',
             'TSC1', 'TSC2', 'NDRG1', 'IRS1']
gene_list = ['Raptor', 'Rictor', 'mSin1', 'Pras40', 'mTOR']
gene_list = ['IRS1']

gene_list = [x.capitalize() for x in gene_list]

for t in tissues:
    meta = pd.read_csv('./phospho_analysis/4_' + t + '_ref_normed.csv')
    gene = meta.columns.to_series().iloc[3:].str.split('_', expand=True)
    gene_in = gene[gene.iloc[:, 0].isin(gene_list)]
    gene_in = gene_in.sort_values(by=0)
    mtor_selected = meta.loc[:, ['Run', 'Channel', 'Condition'] + gene_in.index.to_list()]
    mtor_selected = mtor_selected.sort_values(by='Condition')
    mtor_selected = mtor_selected[
        mtor_selected.Condition.isin(['Normoxia', 'acute11p', 'acute8p', 'chronic11p', 'chronic8p'])]
    mtor_selected = mtor_selected.reset_index(drop=True)

    gene_filtered = np.log1p(mtor_selected.iloc[:, 3:])
    mtor_selected.iloc[:, 3:] = zscore(gene_filtered, axis=0)
    for i in range(3, mtor_selected.shape[1]):
        temp = mtor_selected.iloc[:, [0, 1, 2] + [i]]
        plt.figure(figsize=(5, 4))
        sns.boxplot(data=temp, x="Condition", y=temp.columns[-1], showfliers=False)
        plt.title(t + temp.columns[-1][:5])
        plt.show()
