import pandas as pd
import numpy as np

tissues = ['Plasma']

for tissue in tissues:
    data = pd.read_excel('./MS_data/' + tissue + '/Protein_TMT_Quant.xlsx')
    annot = pd.read_excel('./MS_data/' + tissue + '/Sample-Info.xlsx')

    sample_name = annot['BioReplicate']
    for sample in sample_name:
        meta = annot[annot.BioReplicate == sample]
