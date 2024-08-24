import pandas as pd
import numpy as np

tissues = ['Plasma']

for tissue in tissues:
    data = pd.read_excel('./MS_data/' + tissue + '/Protein_TMT_Quant.xlsx')
    annot = pd.read_excel('./MS_data/' + tissue + '/Sample-Info.xlsx')
    sample_info = pd.read_excel('./MS_data/TMT labeling scheme.xlsx', sheet_name=tissue)
    sample_info = sample_info.iloc[:48, :]

    data = data.set_index('PG.Genes')
    annot = annot.set_index('BioReplicate')
    data = data.T

    meta = pd.concat([annot, data.loc[annot.index, :]], axis=1)
    run_chan = meta['R.BlockName'] + meta.Channel.str[5:]
    condition = pd.DataFrame(data=np.zeros([run_chan.size, 2]), columns=['dose', 'duration'], index=run_chan.values)

    sample_info_condition = 'Plex' + sample_info.Plex.astype(int).astype(str) + \
                            '_' + sample_info['TMT channel'].astype(str)
    condition.loc[sample_info_condition, 'dose'] = sample_info['%Oxygen'].values
    try:
        condition.loc[sample_info_condition, 'duration'] = sample_info['treatment \nduration'].values
    except:
        condition.loc[sample_info_condition, 'duration'] = sample_info['treatment\nduration'].values
    condition = condition.replace({0: 'Bridge'})

    meta.insert(value=condition['dose'].values, loc=0, column='Dose')
    meta.insert(value=condition['duration'].values, loc=0, column='Duration')

    meta.to_csv('./MS_data/' + tissue + '1_meta.csv')

    print(tissue)
