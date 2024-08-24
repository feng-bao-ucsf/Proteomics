import pandas as pd
from tqdm import tqdm

tissues = ['Brain', 'Kidney']
channels = ['126', '127N', '127C', '128N', '128C',
            '129N', '129C', '130N', '130C', '131N', '131C',
            '132N', '132C', '133N', '133C', '134N', '134C',
            '135N']
channels = pd.Series(channels)

select_col = 'PSM.TMT18_' + channels
for t in tissues:
    psm = pd.read_csv('./phospho/' + t + '_Report_PhosphoTMT PSM Quant (Normal).csv')
    reorder_data = pd.DataFrame()
    for run in psm['R.BlockName'].unique():
        print(run)
        psm_select = psm[psm['R.BlockName'] == run]

        pep_exp = psm_select[select_col.to_list() + ['PSM.Phospho (STY)']]
        pep_exp = pep_exp.dropna(axis=0, subset=['PSM.Phospho (STY)'])
        pep_merge = pep_exp.groupby(by='PSM.Phospho (STY)').mean()

        for pp in tqdm(range(pep_merge.shape[0])):
            pp_select = pep_merge.iloc[pp, :]
            pp_name = pp_select.name
            meta_select = psm_select[psm_select['PSM.Phospho (STY)'] == pp_name].iloc[0, :]
            cond_select = meta_select.loc['Condition' + pp_select.index.to_series().str[3:]]
            temp = pd.DataFrame(data={'Channel': pp_select.index,
                                      'Abundance': pp_select.values,
                                      'Run': [run] * pp_select.size,
                                      'Condition': cond_select.values,
                                      'Gene': [meta_select['PG.Genes']] * pp_select.size,
                                      'Peptide': [pp_name] * pp_select.size,
                                      'Uniprot': [meta_select['PG.UniprotIds']] * pp_select.size})
            reorder_data = pd.concat([reorder_data, temp], axis=0, ignore_index=True)

    reorder_data = reorder_data.dropna()
    reorder_data.to_csv('./phospho_analysis/1_' + t + '_raw_data.csv')
