import pandas as pd
import nibabel as nib
from glob import glob
import numpy as np
import statsmodels.api as sp
#imported the full paths since it wasn't working otherwise - not sure why 
from statsmodels.stats.weightstats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection

df_tabular = pd.read_csv(snakemake.input.tsv,sep='\t')

group1 = df.query(snakemake.params.group1_query)
group2 = df.query(snakemake.params.group2_query)

column_name = snakemake.params.column_name


