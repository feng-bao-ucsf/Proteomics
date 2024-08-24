import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
from sklearn import preprocessing
from matplotlib.colors import LinearSegmentedColormap

colors = ["Green", "k", "Red"]
nodes = [0, 0.5, 1.0]
cmap2 = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))

tissues = ['Plasma']

# removed_outlier = []
for tissue in tissues:

    meta = pd.read_csv('./IRS/meta_' + tissue + '_ref_normed.csv')
    meta = meta.sort_values(by='Condition')
    meta = meta.reset_index(drop=True)
    exp = meta.iloc[:, 7:]

    corr_thres = 0.95
    min_count_cutoff = 100

    gene_filtered = np.log1p(exp)

    corr = gene_filtered.T.corr()

    corr_same_cate = np.zeros(corr.shape[1])

    for iii in range(corr.shape[1]):
        corr_select = corr.iloc[:, iii]
        corr_select_same_condition = corr_select[meta.Condition == meta.Condition.iloc[iii]]
        # if one extreme sample affect the results
        # if np.sum(corr_select_same_condition < 0.9) == 1:
        #     corr_select_same_condition = corr_select_same_condition[corr_select_same_condition >= 0.9]

        corr_same_cate[iii] = (corr_select_same_condition.sum() - 1) / (corr_select_same_condition.size - 1)

    meta.insert(value=corr_same_cate, loc=0, column='Correlation')

    meta = meta.sort_values(by=['Condition', 'Correlation'])

    meta['Correlation'] = meta['Correlation'] >= corr_thres
    meta['Correlation'] = meta['Correlation'].astype(int)

    cmap = cm.get_cmap('rainbow', len(meta.Condition.unique()))
    lut2 = dict(zip(meta.Condition.unique(), cmap(range(len(meta.Condition.unique())))))
    row_colors2 = meta.Condition.map(lut2)

    cmap = cm.get_cmap('Greys', 2)
    lut4 = dict(zip(range(2), cmap(range(2))))
    row_colors4 = meta['Correlation'].map(lut4)

    row_colors = pd.concat([row_colors2, row_colors4], axis=1)

    gene = pd.DataFrame(data=gene_filtered.values, index=gene_filtered.index)
    plt.figure()
    sns.clustermap(gene.T.corr(), xticklabels=meta.BioReplicate, yticklabels=meta.BioReplicate,
                   row_cluster=False, col_cluster=False, row_colors=row_colors, annot=False, figsize=(10, 10),
                   vmin=0.8, vmax=1)
    plt.show()
    # sns.heatmap(gene.T.corr(), annot=True, square=True, vmin=0, vmax=1)
    # plt.savefig('./QC_correlation/' + tissue + '_all_genes_0.9.png', dpi=300, bbox_inches='tight')
