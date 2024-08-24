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
tissues = ['Brain', 'Kidney']
gene_list = ['S6K1', 'RPS6KB1', 'p70S6k', 'S6', '4EBP1', 'eIF4G', 'Akt', 'Ulk1', 'Lipin1', 'TFEB', 'Grb10', 'PKCa',
             'Sgk1',
             'Raptor', 'Rictor', 'mSin1', 'Pras40', 'mTOR',
             'TSC1', 'TSC2', 'NDRG1', 'IRS1']
gene_list = [x.capitalize() for x in gene_list]

for t in tissues:
    meta = pd.read_csv('./phospho_analysis/4_' + t + '_ref_normed.csv')
    gene = meta.columns.to_series().iloc[3:].str.split('_', expand=True)
    gene_in = gene[gene.iloc[:, 0].isin(gene_list)]
    gene_in = gene_in.sort_values(by=0)
    mtor_selected = meta.loc[:, ['Run', 'Channel', 'Condition'] + gene_in.index.to_list()]
    mtor_selected = mtor_selected.sort_values(by='Condition')
    mtor_selected = mtor_selected.reset_index(drop=True)

    cmap = cm.get_cmap('Blues', 3)
    lut1 = dict(zip(['Plex1', 'Plex2', 'Plex3'], cmap(range(len(mtor_selected.Run.unique())))))
    row_colors1 = mtor_selected.Run.map(lut1)

    cmap = cm.get_cmap('rainbow', len(mtor_selected.Condition.unique()))
    lut2 = dict(zip(['Normoxia', 'acute11p', 'acute8p', 'chronic11p', 'chronic8p', 'precondition11p', 'precondition8p'],
                    cmap(range(len(mtor_selected.Condition.unique())))))
    row_colors2 = mtor_selected.Condition.map(lut2)

    gene_filtered = np.log1p(mtor_selected.iloc[:, 3:])
    gene_filtered = zscore(gene_filtered, axis=0)
    # gene_filtered = (gene_filtered - gene_filtered.min()) / (gene_filtered.max() - gene_filtered.min())
    row_colors = pd.concat([row_colors1, row_colors2], axis=1)

    gene = pd.DataFrame(data=gene_filtered.values, columns=gene_in.iloc[:, 0])
    plt.figure()
    sns.clustermap(gene,
                   row_cluster=False, col_cluster=False, row_colors=row_colors,
                   annot=False, figsize=(15, 10),
                   cmap=cmap2, vmin=-2, vmax=2)
    plt.title(t)
    plt.show()
