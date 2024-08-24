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

tissues = ['Plasma']

for t in tissues:
    meta = pd.read_csv('./IRS/meta_' + t + '_ref_normed.csv')

    meta = meta.sort_values(by='Condition')
    # meta = meta[
    #     meta.Condition.isin(['Normoxia', 'acute11p', 'acute8p', 'chronic11p', 'chronic8p'])]

    meta = meta.reset_index(drop=True)

    cmap = cm.get_cmap('Blues', 3)
    lut1 = dict(zip(['Plex1', 'Plex2', 'Plex3'], cmap(range(len(meta['R.BlockName'].unique())))))
    row_colors1 = meta['R.BlockName'].map(lut1)

    cmap = cm.get_cmap('rainbow', len(meta.Condition.unique()))
    lut2 = dict(zip(meta.Condition.unique(),
                    cmap(range(len(meta.Condition.unique())))))
    row_colors2 = meta.Condition.map(lut2)

    gene_filtered = np.log1p(meta.iloc[:, 7:])
    gene_filtered = zscore(gene_filtered, axis=0)
    # gene_filtered = (gene_filtered - gene_filtered.min()) / (gene_filtered.max() - gene_filtered.min())
    row_colors = pd.concat([row_colors1, row_colors2], axis=1)
    row_colors = row_colors.set_index(meta.BioReplicate)

    gene = pd.DataFrame(data=gene_filtered.values, index=meta.BioReplicate)
    plt.figure()
    sns.clustermap(gene,
                   row_cluster=False, col_cluster=True, row_colors=row_colors,
                   annot=False, figsize=(15, 15),
                   cmap=cmap2, vmin=-2, vmax=2)
    plt.title(t)
    plt.show()
    print(gene.shape)
