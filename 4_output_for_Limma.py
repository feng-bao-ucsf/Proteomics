import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import conorm

'''
TMM norm
'''
tissues = ['Plasma']

for t in tissues:
    meta = pd.read_csv('./IRS/meta_' + t + '_ref_normed.csv')

    condition = meta.loc[:, ['Duration', 'Dose']].astype(str)
    condition['Duration'] = condition['Duration'].replace({'8d': 'precondition',
                                                           '24': 'acute',
                                                           '28d': 'chronic',
                                                           'nan': ''})

    condition['Dose'] = condition['Dose'].replace({'11': '11p',
                                                   '11%-8%': '8p',
                                                   '8': '8p'})
    condition.insert(value=condition['Duration'] + condition['Dose'], loc=0, column='Condition')
    exp = meta.iloc[:, 7:].T
    exp = np.log(exp)

    cond_matrix = pd.DataFrame(np.zeros([meta.shape[0], condition.Condition.unique().size]),
                               index=exp.columns, columns=condition.Condition.unique(), dtype=int)

    for cod in condition.Condition.unique():
        id = condition[condition.Condition == cod].index
        cond_matrix.loc[id, cod] = 1

    exp.to_csv('./limma/Exp_' + t + '.csv', index=True)
    cond_matrix.to_csv('./limma/Cod_' + t + '.csv', index=True)
    print(t)
