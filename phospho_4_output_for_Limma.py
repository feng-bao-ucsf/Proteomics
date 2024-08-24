import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import conorm

'''
TMM norm
'''
tissues = ['Liver', 'Heart', 'Muscle', 'Brain', 'Kidney']

for t in tissues:
    meta = pd.read_csv('./phospho_analysis/3_' + t + '_ref_normed.csv')
    sample_info = pd.read_excel('./MS_data/TMT labeling scheme.xlsx', sheet_name=t)
    sample_info['%Oxygen'] = sample_info['%Oxygen'].astype(str)
    sample_info['treatment \nduration'] = sample_info['treatment \nduration'].astype(str)
    sample_info['treatment \nduration'] = sample_info['treatment \nduration'].replace({'8d': 'precondition',
                                                                                       '24': 'acute',
                                                                                       '28d': 'chronic',
                                                                                       'nan': ''})

    sample_info['%Oxygen'] = sample_info['%Oxygen'].replace({'11': '11p',
                                                             '11%-8%': '8p',
                                                             '8': '8p'})

    condi_mapping = pd.DataFrame(data={'Condition': 'Group_' + sample_info.Group.astype(str),
                                       'Treatment': sample_info['treatment \nduration'].astype(str) +
                                                    sample_info['%Oxygen'].astype(str)})
    condi_mapping = condi_mapping.drop_duplicates(subset=['Condition'])
    condi_mapping = dict(zip(condi_mapping.iloc[:, 0], condi_mapping.iloc[:, 1]))

    meta = meta.replace(condi_mapping)
    meta.to_csv('./phospho_analysis/4_' + t + '_ref_normed.csv', index=False)

    exp = meta.iloc[:, 3:].T
    exp = np.log(exp)

    cond_matrix = pd.DataFrame(np.zeros([meta.shape[0], meta.Condition.unique().size]),
                               index=exp.columns, columns=meta.Condition.unique(), dtype=int)

    for cod in meta.Condition.unique():
        id = meta[meta.Condition == cod].index
        cond_matrix.loc[id, cod] = 1

    exp.to_csv('./phospho_DE/Exp_' + t + '.csv', index=True)
    cond_matrix.to_csv('./phospho_DE/Cod_' + t + '.csv', index=True)
    print(t)
