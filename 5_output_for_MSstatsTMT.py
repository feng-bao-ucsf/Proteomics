import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import conorm

'''
TMM norm
'''
tissues = ['Brain', 'Kidney']

for t in tissues:
    meta = pd.read_csv('./IRS/5_meta_' + t + '_ref_normed_TMM.csv')
    meta = meta[meta['TMT_condition'] != 'Norm']
    ms_data = pd.DataFrame()

    meta['TMT_condition'] = meta['TMT_condition'].replace({'Group_4': 'Group_1'})

    protein_names = meta.columns.to_list()[5:]
    replicate_data = pd.DataFrame(np.zeros([meta.shape[0], 1], dtype=int), index=meta.index)

    for cond in meta['TMT_condition'].unique():
        cond_size = meta[meta['TMT_condition'] == cond].shape[0]
        replicate_data[meta['TMT_condition'] == cond] = np.expand_dims(range(cond_size), 1)

    for protein in protein_names:
        temp = pd.DataFrame(data={'Abundance': meta[protein].values,
                                  'Protein': [protein] * meta.shape[0],
                                  'Run': meta.Run.values,
                                  'Channel': meta.Channel.values,
                                  'Mixture': np.ones(meta.shape[0]),
                                  'TechRepMixture': np.ones(meta.shape[0]),
                                  'Condition': meta['TMT_condition'].values,
                                  'BioReplicate': replicate_data.values[:, 0]})
        ms_data = ms_data.append(temp, ignore_index=True)
    ms_data.to_csv('./IRS/6_MSTMT_' + t + '.csv', index=False)
    print(t)
    # [1]
    # "Mixture"        "TechRepMixture" "Run"            "Channel"
    # [5]
    # "Protein"        "Abundance"      "BioReplicate"   "Condition"
