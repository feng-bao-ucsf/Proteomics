import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import conorm

'''
1. Within each run, calculate average pool channel
2. Across runs, calculate the geometric mean of average pool channels
'''
tissues = ['Plasma']

for t in tissues:
    meta = pd.read_csv('./MS_data/' + t + '1_meta.csv')
    exp = meta.iloc[:, 7:]
    exp = exp.replace({0: np.nan})
    exp = exp.dropna(axis=1)

    exp_1 = exp[meta['R.BlockName'] == 'Plex1']
    exp_2 = exp[meta['R.BlockName'] == 'Plex2']
    exp_3 = exp[meta['R.BlockName'] == 'Plex3']

    exp_1_ref = exp[(meta['R.BlockName'] == 'Plex1') & (meta['Condition'] == 'Bridge')]
    exp_2_ref = exp[(meta['R.BlockName'] == 'Plex2') & (meta['Condition'] == 'Bridge')]
    exp_3_ref = exp[(meta['R.BlockName'] == 'Plex3') & (meta['Condition'] == 'Bridge')]

    exp_1_ref = exp_1_ref.mean(axis=0)
    exp_2_ref = exp_2_ref.mean(axis=0)
    exp_3_ref = exp_3_ref.mean(axis=0)

    ref_ave = pd.concat([exp_1_ref, exp_2_ref, exp_3_ref], axis=1)
    ref_geo_mean = np.exp(np.log(ref_ave).mean(axis=1))

    scale_1 = ref_geo_mean / exp_1_ref
    scale_2 = ref_geo_mean / exp_2_ref
    scale_3 = ref_geo_mean / exp_3_ref

    exp_1_scaled = exp_1 * scale_1
    exp_2_scaled = exp_2 * scale_2
    exp_3_scaled = exp_3 * scale_3

    exp_scaled = pd.concat([exp_1_scaled, exp_2_scaled, exp_3_scaled], axis=0)

    exp_tmm = conorm.tmm(exp_scaled.T)
    exp_tmm = exp_tmm.T

    meta_scaled = pd.concat([meta.iloc[:, :7], exp_tmm], axis=1)
    meta_scaled = meta_scaled[meta_scaled.Condition != 'Bridge']
    meta_scaled.to_csv('./IRS/meta_' + t + '_ref_normed.csv', index=False)
