import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import conorm

'''
1. Within each run, calculate average pool channel
2. Across runs, calculate the geometric mean of average pool channels
'''
tissues = ['Liver', 'Heart', 'Muscle', 'Brain', 'Kidney']

for t in tissues:
    print(t)
    meta = pd.read_csv('./phospho_analysis/2_' + t + '_phosphosites.csv')
    exp = meta.iloc[:, 3:]
    # exp = exp.replace({0: np.nan})
    # for r in meta.Run.unique():
    #     for cond in meta.Condition.unique():
    #         a = exp[(meta.Run == r) & (meta.Condition == cond)]
    #         exp[(meta.Run == r) & (meta.Condition == cond)] = a.fillna(a.mean())

    exp = exp.fillna(exp.mean())
    exp = exp.dropna(axis=1)

    exp_1 = exp[meta['Run'] == 'Plex1']
    exp_2 = exp[meta['Run'] == 'Plex2']
    exp_3 = exp[meta['Run'] == 'Plex3']

    exp_1_ref = exp[(meta['Run'] == 'Plex1') & (meta['Condition'] == 'Bridge')]
    exp_2_ref = exp[(meta['Run'] == 'Plex2') & (meta['Condition'] == 'Bridge')]
    exp_3_ref = exp[(meta['Run'] == 'Plex3') & (meta['Condition'] == 'Bridge')]

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

    meta_scaled = pd.concat([meta.iloc[:, :3], exp_tmm], axis=1)
    meta_scaled = meta_scaled[meta_scaled.Condition != 'Bridge']
    meta_scaled.to_csv('./phospho_analysis/3_' + t + '_ref_normed_simple_impute.csv', index=False)
