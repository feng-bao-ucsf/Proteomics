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

for t in tissues:
    mtor_selected = pd.read_csv('./phospho_analysis/3_' + t + '_ref_normed_no_impute.csv')

    mtor_selected = mtor_selected.sort_values(by='Condition')
    mtor_selected = mtor_selected[
        mtor_selected.Condition.isin(['Normoxia', 'acute11p', 'acute8p', 'chronic11p', 'chronic8p'])]

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

    gene = pd.DataFrame(data=gene_filtered.values, )
    plt.figure()
    sns.clustermap(gene,
                   row_cluster=False, col_cluster=True, row_colors=row_colors,
                   annot=False, figsize=(15, 10),
                   cmap=cmap2, vmin=-2, vmax=2)
    plt.title(t)
    plt.show()
    print(gene.shape)
