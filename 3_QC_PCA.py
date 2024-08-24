'''
Using clustering and PCA to visualize the batch bias
'''
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib import cm
from sklearn.decomposition import PCA
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import plotly.express as px

tissues = ['Plasma']

# removed_outlier = []
for tissue in tissues:

    meta = pd.read_csv('./IRS/meta_' + tissue + '_ref_normed.csv')
    print(meta.shape)

    gene = meta.iloc[:, 7:]

    gene = np.log1p(gene)

    '''
    using top 2 PCs, if replicate is great than 3 (or 4) std, remove it
    '''

    std_fold = 3
    bad_samples = []

    pca = PCA(n_components=2)
    gene_pc = pca.fit_transform(gene.fillna(0))

    condition = meta.Dose.astype(str) + '_' + meta.Duration.astype(str)

    pc_pd = pd.DataFrame(
        data={'condition': condition.values, 'PC-1': gene_pc[:, 0], 'PC-2': gene_pc[:, 1],
              'Plex': meta['R.BlockName'].str[4]},
        index=meta.index)
    pc_mean = pc_pd.groupby('condition').median()
    pc_std = pc_pd.groupby('condition').mad()

    for treat in condition.unique():
        pc_pd_select = pc_pd[pc_pd['condition'] == treat]

        pc1_max = pc_mean.loc[treat, 'PC-1'] + pc_std.loc[treat, 'PC-1'] * std_fold
        pc1_min = pc_mean.loc[treat, 'PC-1'] - pc_std.loc[treat, 'PC-1'] * std_fold

        pc2_max = pc_mean.loc[treat, 'PC-2'] + pc_std.loc[treat, 'PC-2'] * std_fold
        pc2_min = pc_mean.loc[treat, 'PC-2'] - pc_std.loc[treat, 'PC-2'] * std_fold

        pc_pd_select_bad = pc_pd_select[(pc_pd_select['PC-1'] > pc1_max) | (pc_pd_select['PC-1'] < pc1_min) |
                                        (pc_pd_select['PC-2'] > pc2_max) | (pc_pd_select['PC-2'] < pc2_min)]
        if pc_pd_select_bad.shape[0] > 0:
            bad_samples = bad_samples + pc_pd_select_bad.index.to_list()
    bad_samples = pd.Series(bad_samples)

    sample_quality = pd.Series(data=['Good'] * pc_pd.shape[0], index=pc_pd.index)
    sample_quality.loc[bad_samples] = 'Bad'
    pc_pd.insert(loc=0, value=sample_quality, column='QC')

    plt.figure(figsize=(7, 7))
    sns.scatterplot(data=pc_pd, x='PC-1', y='PC-2', hue='Plex')
    # fig.write_html("./fig/IRS_" + tissue + "_plex.html")
    plt.show()

    # plt.figure(figsize=(7, 7))
    # sns.scatterplot(data=pc_pd, x='PC-1', y='PC-2', color='condition')
    # # fig.write_html("./fig/IRS_" + tissue + "_condition.html")
    # plt.show()

    plt.figure(figsize=(7, 7))
    sns.scatterplot(data=pc_pd, x='PC-1', y='PC-2', hue='condition',
                    legend=True,
                    s=50)
    plt.title(tissue)
    plt.show()
    # plt.savefig('./QC_PCA/PCA_' + tissue + '.png', dpi=300, bbox_inches='tight')

    # plt.figure(figsize=(7, 7))
    # sns.scatterplot(data=pc_pd, x='PC-1', y='PC-2', hue='time',
    #                 hue_order=['variable', '24', '8d', '28d'], palette='Reds')
    # plt.title(tissue + ' Time')
    # plt.show()
    #
    # plt.figure(figsize=(7, 7))
    # sns.scatterplot(data=pc_pd, x='PC-1', y='PC-2', hue='dose',
    #                 hue_order=['Normoxia', '11%', '11%-8%', '8%'], palette='Blues')
    # plt.title(tissue + ' Concentration')
    # plt.show()
