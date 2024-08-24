import pandas as pd
from tqdm import tqdm

tissues = ['Liver', 'Heart', 'Muscle', 'Brain', 'Kidney']

for t in tissues:
    pep = pd.read_csv('./phospho_analysis/1_' + t + '_raw_data.csv', index_col=0)
    pep = pep.dropna(subset=['Peptide'], axis=0)

    pep_unique = pep.Gene + '_' + pep.Uniprot + '_' + pep.Peptide
    pep_unique = pep_unique.unique()

    merged_data = pd.DataFrame(columns=['Run', 'Channel', 'Condition'] + list(pep_unique))

    for run in pep.Run.unique():
        print(run)
        for chan in pep.Channel.unique():
            pep_select = pep[(pep.Run == run) & (pep.Channel == chan)]
            merged_temp = pd.DataFrame(columns=['Run', 'Channel', 'Condition'] + list(pep_unique), index=[0])

            pep_name_select = pep_select.Gene + '_' + pep_select.Uniprot + '_' + pep_select.Peptide
            pep_select = pep_select.rename(index=dict(zip(pep_select.index, pep_name_select.values)))

            merged_temp.loc[0, pep_name_select] = pep_select.Abundance
            merged_temp.loc[0, 'Run'] = run
            merged_temp.loc[0, 'Channel'] = chan
            merged_temp.loc[0, 'Condition'] = pep_select.Condition[0]

            merged_data = pd.concat([merged_data, merged_temp], axis=0, ignore_index=True)

    merged_data.to_csv('./phospho_analysis/2_' + t + '_raw_data_widefmt.csv', index=False)
