import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import conorm
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import zscore
from sklearn.decomposition import PCA

colors = ["Green", "k", "Red"]
nodes = [0, 0.5, 1.0]
cmap2 = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))

'''
TMM norm
'''
tissues = ['Brain', 'Kidney']

for t in tissues:
    meta = pd.read_csv('./phospho_analysis/3_' + t + '_ref_normed_no_impute.csv')
    exp = meta.iloc[:, 3:]
    exp = np.log1p(exp)

    pca = PCA(n_components=2)
    gene_pc = pca.fit_transform(exp)

    pc_pd = pd.DataFrame(
        data={'Condition': meta.Condition, 'Run': meta.Run, 'PC-1': gene_pc[:, 0], 'PC-2': gene_pc[:, 1]},
        index=meta.index)

    plt.figure(figsize=(7, 7))
    sns.scatterplot(data=pc_pd, x='PC-1', y='PC-2', hue='Condition',
                    legend=True,
                    hue_order=['Normoxia', 'acute11p', 'acute8p', 'chronic11p', 'chronic8p', 'precondition11p',
                               'precondition8p'],
                    s=50)
    plt.title(t)
    plt.show()

    plt.figure(figsize=(7, 7))
    sns.scatterplot(data=pc_pd, x='PC-1', y='PC-2', hue='Run',
                    legend=True,
                    s=50)
    plt.title(t)
    plt.show()
