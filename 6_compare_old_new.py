import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

tissues = ['Brain', ]
p_thres = 0.05
lfc_thres = 0.5

name = []
for t in tissues:
    de_old = pd.read_csv(
        '../Exam-2022_oct/Proteomics/' + str.lower(t) + '_proteomics_log2FC_pval_padj.txt',
        delimiter='\t', index_col=1)

    de_new = pd.read_csv('./DE/' + t + '/limma_acute11p.csv', index_col=0)

    de_old_select = de_old['11per_24h_log2FC']
    de_new_select = de_new['logFC']

    shared_name = de_old_select.index.to_list() + de_new_select.index.to_list()
    shared_name = pd.Series(shared_name)
    shared_name = shared_name.value_counts()

    shared_name = shared_name[shared_name == 2]

    plt.figure()
    plt.scatter(de_old_select[shared_name.index], de_new_select[shared_name.index], s=5)
    plt.axis('square')
    plt.xlim([-1, 1])
    plt.ylim([-1, 1])
    plt.show()
