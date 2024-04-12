from pathlib import Path
import pandas as pd
from scipy.stats import chi2_contingency

def calculate_chi2_from_matrix(matrix_1, matrix_2):
    matrix_1 = Path(matrix_1)
    data_1 = pd.read_csv(matrix_1, sep='\t', index_col='MutationType')
    matrix_2 = Path(matrix_2)
    data_2 = pd.read_csv(matrix_2, sep='\t', index_col='MutationType')
    
    # Calculate the sum across rows for each data
    data_1['Sum'] = data_1.sum(axis=1)
    data_2['Sum'] = data_2.sum(axis=1)
    
    # Combine the two dataframes into one
    combined_data = pd.concat([data_1['Sum'], data_2['Sum']], axis=1, keys=['Group1_Sum', 'Group2_Sum'])
    
    # Filter out rows where the sum across both groups is zero
    filtered_data = combined_data[(combined_data['Group1_Sum'] > 0) | (combined_data['Group2_Sum'] > 0)]
    
    # Prepare the contingency table
    contingency_table = filtered_data[['Group1_Sum', 'Group2_Sum']]
    
    # Perform the Chi-squared test
    chi2, p, dof, expected = chi2_contingency(contingency_table)
    
    print(f"Chi-squared: {chi2}")
    print(f"P-value: {p}")
    print(f"Degrees of freedom: {dof}")
    print("Expected frequencies:")
    print(expected)

calculate_chi2_from_matrix('/media/cam/Working/8-oxodG/hmces/SigProfiler/wt_input_vcfs/output/ID/wt.ID83.all', '/media/cam/Working/8-oxodG/hmces/SigProfiler/hmces_input_vcfs/output/ID/hmces.ID83.all')