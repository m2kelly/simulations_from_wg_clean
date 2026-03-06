import pandas as pd
qval1='/home/maria/run_simulations_cactus_clean/scripts/2_run_cbase/hg38_Anc4/hg38_Anc4/Anc4->hg38/q_values_hg38_Anc4.txt'
qval2='/home/maria/run_simulations_cactus_clean/playground/q_values_cactus.txt'


df1=pd.read_csv(qval1,sep='\t', index_col=0, skiprows=1)
print(len(df1[df1['q_phi_neg'] < 0.05]))

df2=pd.read_csv(qval2,sep='\t', index_col=0, skiprows=1)
print(len(df2[df2['q_phi_neg'] < 0.05]))

