import pandas as pd

df_subjects = pd.DataFrame({'participant_label': snakemake.params.subjects})

df_subjects['participant_id'] = 'sub-' + df_subjects['participant_label']
df_clinical = pd.read_csv(snakemake.input.tsv,sep='\t',dtype={"participant_label": str})

df = pd.merge(left=df_subjects,right=df_clinical,left_on='participant_id',right_on='participant_id')

columns = list(df.columns)
columns.remove('participant_label')

df[columns].to_csv(snakemake.output.tsv,sep='\t',index=False)
