import pandas as pd

#load up subjects tables
datasets = snakemake.params.datasets
tabular_tsvs = snakemake.input.tabular_tsvs
subjects_tsvs = snakemake.input.subjects_tsvs

df_list=[]
for tabular_tsv,subjects_tsv,dataset in zip(tabular_tsvs,subjects_tsvs,datasets):
    df_ = pd.read_csv(subjects_tsv,
                                       dtype={"participant_label": str},
                                       sep='\t')
    df_['participant_id'] = 'sub-' + df_['participant_label']
    
    #add a new dataset column
    df_['dataset'] = dataset

    df_tabular_ = pd.read_csv(tabular_tsv,
                                       dtype={"participant_label": str},
                                       sep='\t')

    if 'participant_label' in df_tabular_.columns:
        df_tabular_['participant_id'] = 'sub-' + df_tabular_['participant_label']
        df_tabular_ = df_tabular_.drop('participant_label',axis=1)

    if 'participant_label' in df_:
        df_ = df_.drop('participant_label',axis=1)
        
    df_merged_ = pd.merge(df_,df_tabular_,how='left',on='participant_id')

    df_list.append(df_merged_)
    

df_subjects = pd.concat(df_list)


df_subjects.to_csv(snakemake.output.tsv,sep='\t',index=False)
