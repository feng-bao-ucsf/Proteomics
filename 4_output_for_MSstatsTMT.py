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
    meta = pd.read_csv('./IRS/meta_' + t + '_ref_normed.csv')
    meta.Condition.replace({'Group_4': 'Group_1'})

    ms_data = pd.DataFrame()

    protein_names = meta.columns.to_list()[7:]
    conditions = meta.Dose.astype(str) + '_' + meta.Duration.astype(str)

    for protein in protein_names:
        temp = pd.DataFrame(data={'Abundance': meta[protein].values,
                                  'Protein': [protein] * meta.shape[0],
                                  'Run': meta['R.BlockName'].values,
                                  'Channel': meta.Channel.values,
                                  'Mixture': np.ones(meta.shape[0]),
                                  'TechRepMixture': np.ones(meta.shape[0]),
                                  'Condition': meta.Condition.values,
                                  'BioReplicate': meta.BioReplicate})
        ms_data = ms_data.append(temp, ignore_index=True)
    # ms_data.Condition=ms_data.Condition.replace({ 'Normoxia_nan':})
    ms_data.to_csv('./IRS/MSTMT_' + t + '.csv', index=False)
    print(t)
    # [1]
    # "Mixture"        "TechRepMixture" "Run"            "Channel"
    # [5]
    # "Protein"        "Abundance"      "BioReplicate"   "Condition"
