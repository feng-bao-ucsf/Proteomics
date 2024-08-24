import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import conorm
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import zscore

'''
TMM norm
'''
tissues = ['Liver', 'Heart', 'Muscle', 'Brain', 'Kidney']

size_all = pd.DataFrame(np.zeros([5, 2], dtype=int), index=tissues, columns=['Protein', 'Site'])
plt.figure()
for t in tissues:
    print(t)

    mtor_selected = pd.read_csv('./phospho_analysis/2_' + t + '_phosphosites.csv')
    mtor_selected = mtor_selected.replace({0: np.nan})
    mtor_selected = mtor_selected.dropna(axis=1, how='all')

    names = mtor_selected.columns.to_series().iloc[3:].str.split('_', expand=True)

    size_all.loc[t, 'Protein'] = names.iloc[:, 1].unique().size
    size_all.loc[t, 'Site'] = names.index.to_series().unique().size

    a = names.iloc[:, 1].value_counts()

    plt.plot(range(a.size), a)
plt.show()

size_all.T.to_csv('./phospho_analysis/numbers.csv')
size_all = size_all.stack().reset_index()

plt.figure(figsize=(4.5, 3.5))
sns.barplot(data=size_all, x='level_0', y=0, hue='level_1')
plt.show()

kinase_info = pd.read_csv('./phospho_analysis/Human_Kinase_Substrate_List_PhosphoSitePlus.csv')
k_name = kinase_info.KINASE.unique()
kinase_percent = pd.DataFrame(np.zeros([k_name.size, 5]), index=k_name,
                              columns=tissues)

for t in tissues:
    print(t)
    mtor_selected = pd.read_csv('./phospho_analysis/2_' + t + '_phosphosites.csv')
    mtor_selected = mtor_selected.replace({0: np.nan})
    mtor_selected = mtor_selected.dropna(axis=1, how='all')

    names = mtor_selected.columns.to_series().iloc[3:].str.split('_', expand=True)
    data_names = names.iloc[:, 1] + '_' + names.iloc[:, 2]

    for k in k_name:
        k_select = kinase_info[kinase_info.KINASE == k]
        k_select_name = k_select['SUB_ACC_ID'] + '_' + k_select['SUB_MOD_RSD']
        k_overlap = k_select_name[k_select_name.isin(data_names.values)]
        kinase_percent.loc[k, t] = (k_overlap.size + 0.0) / k_select_name.size

plt.figure()
sns.heatmap(kinase_percent)
plt.show()
