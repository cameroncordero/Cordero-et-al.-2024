import pandas as pd
from scipy.spatial.distance import cosine

file1 = '/media/cam/Working/8-oxodG/lesion_files/comparison_to_sbs/SRR_treated_cellular_69-70_fixed_context_percentages.txt'
file2 = '/media/cam/Working/8-oxodG/lesion_files/comparison_to_sbs/sbs96.txt'

df1 = pd.read_csv(file1, sep='\t', index_col=0, header=0)
df2 = pd.read_csv(file2, sep='\t', index_col=0, header=0)

# print(df1.head())
# print(df2.head())

# Select specific signature
array1 = df1.loc[:, 'SBS96A'].to_numpy().flatten()
array2 = df2.loc[:, 'SBS96A'].to_numpy().flatten()

print(sum(array1))
print(sum(array2))

print(1-cosine(array1, array2))