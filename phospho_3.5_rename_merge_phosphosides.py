import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import conorm

# '''
# TMM norm
# '''
# tissues = ['Liver', 'Heart', 'Muscle', 'Brain', 'Kidney']
#
# for t in tissues:
#     meta = pd.read_csv('./phospho_analysis/3_' + t + '_ref_normed.csv')
#     phosphosites = pd.read_csv('./phospho_analysis/' + t + '_imputed_Meta_Modification_site.csv')
#     exp = meta.iloc[:, 3:]
#     exp = exp.T
#
#     phosphosites = phosphosites[phosphosites.Probability > 0.75]
#     phosphosites = phosphosites[phosphosites.Original.isin(exp.index.to_list())]
#
#     new_name = phosphosites.Gene + '_' + phosphosites.Uniprot + '_' + phosphosites.Phos_Site
#     phosphosites.insert(loc=0, value=new_name, column='Merged_name')
#
#     new_phospho = pd.DataFrame()
#     for name in new_name.unique():
#         phosphosites_sub = phosphosites[phosphosites.Merged_name == name]
#         exp_sub = exp.loc[phosphosites_sub.Original, :]
#         if exp_sub.shape[0] == 1:
#             new_phospho.insert(loc=0, value=np.squeeze(exp_sub.values), column=name)
#         else:
#             prob = phosphosites_sub.Probability.values
#             prob = np.expand_dims(prob, axis=0)
#             exp_sub_merge = np.matmul(prob, exp_sub.values)
#             new_phospho.insert(loc=0, value=np.squeeze(exp_sub_merge), column=name)
#
#     meta_new = pd.concat([meta.iloc[:, :3], new_phospho], axis=1)
