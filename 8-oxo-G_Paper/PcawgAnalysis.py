import seaborn as sns
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import math
import numpy as np

def convert_tsv_to_dict(tsv_file):
    """Convert a tsv file to a dictionary"""
    tsv_file = Path(tsv_file)
    # Read the TSV file into a DataFrame
    df = pd.read_csv(tsv_file, sep='\t', index_col=0)
    # Convert the DataFrame to a dictionary
    mutation_dict = df.to_dict()
    return mutation_dict

def create_sig_dict(indel_dict, signature = "ID83", percentage = 0.15):
    """Create a dict of specific signature"""
    sig_dict = {}
    for key, value in indel_dict.items():
        if value[signature]/sum(value.values()) > percentage:
            sig_dict[key] = value[signature]
    return sig_dict

def create_filter_sig_dict(dict, signature = "SBS18"):
    """Filter a dict to a specific signature"""
    sig_dict = {}
    for key, value in dict.items():
        if value[signature]/sum(value.values()) > 1:
            sig_dict[key] = value[signature]
    return sig_dict

def plot_pearson_scatter_sigs(sbs18_dict, id83_dict):
    signature_correlation_dict = {}
    # Get the keys of one of the nested dictionaries
    signatures = list(next(iter(sbs18_dict.values())).keys())
    for signature in signatures:
        sig_dict = {}
        for key, value in sbs18_dict.items():
            sig_dict[key] = value[signature]
        x_values = []
        y_values = []
        for key in sig_dict:
            try:
                id83_dict[key]
            except KeyError:
                continue
            if sig_dict[key] < 20:
                continue
            if id83_dict[key] < 20:
                continue
            x_values.append(sig_dict[key])
            y_values.append(id83_dict[key])
        
        # Calculate the Pearson correlation coefficient and the p-value
        corr, p_value = pearsonr(x_values, y_values)
        signature_correlation_dict[signature] = [corr, p_value*len(signatures), len(x_values)]
    
    return signature_correlation_dict

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib import colors as mcolors

def plot_correlation_coefficients(signature_correlation_dict, filter_p_value=False):
    """Plot the correlation coefficients of different signatures"""
    # Create a DataFrame from the signature_correlation_dict
    data = pd.DataFrame.from_dict(signature_correlation_dict, orient='index', columns=['Correlation Coefficient', 'p-value', 'number of samples'])
    data.reset_index(inplace=True)
    data.columns = ['Signature', 'Correlation Coefficient', 'p-value', 'number of samples']
    
    # Filter out rows with p-value greater than 0.01 if filter_p_value is True
    if filter_p_value:
        data = data[data['p-value'] <= 0.01]
    
    # Set the figure size
    fig, ax = plt.subplots(figsize=(6, 8))
    
    # Sort the data by the 'Correlation Coefficient' column
    data.sort_values('Correlation Coefficient', ascending=False, inplace=True)
    
    # Create the strip plot with smaller dots and some jitter
    top_10_signatures = data.nlargest(20, 'Correlation Coefficient')
    data['color'] = np.where(data['Signature'].isin(top_10_signatures['Signature']), data['Signature'], 'Other')
    
    # Create a color cycle with darker colors
    color_cycle = ['#FF0000', '#FFA500', '#008000', '#0000FF', '#800080']
    
    # Create a color dictionary based on the order of the signatures in the 'Signature' column
    ordered_signatures = data['Signature'].unique()
    num_colors = len(ordered_signatures)
    my_palette = dict(zip(ordered_signatures, [color_cycle[i % len(color_cycle)] for i in range(num_colors)]))
    
    # Check if 'Other' is in my_palette, if not, add it with a default color
    if 'Other' not in my_palette:
        my_palette['Other'] = 'gray'  # This is a shade of gray
    
    sns.stripplot(x=np.zeros(data.shape[0]), y='Correlation Coefficient', hue='color', data=data, size=10, jitter=0.2, ax=ax, palette=my_palette, edgecolor='black', linewidth=1)
    
    # Create the legend
    legend_patches = [patches.Patch(color=my_palette[sig], label=f'{sig} (Adjusted p={data.loc[data["Signature"]==sig, "p-value"].values[0]:.1e}, n={data.loc[data["Signature"]==sig, "number of samples"].values[0]})') for sig in top_10_signatures['Signature']]
    ax.legend(handles=legend_patches, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Set the title and labels
    ax.set_title('Correlation Coefficients of Different Signatures')
    ax.set_xlabel('')
    ax.set_ylabel('Correlation Coefficient')
    
    # Remove the x-axis ticks and set x-limits
    ax.set_xticks([])
    ax.set_xlim(-0.5, 0.5)
    
    # Tighten the layout
    plt.tight_layout()
    
    plt.show()

def plot_dictionaries(sbs18_dict, id83_dict, min_sbs18_mutations = 20, min_id83_mutations = 20):
    """Plot the values of two dictionaries"""
    # Create lists for the x and y values
    x_values = []
    y_values = []
    x_values_blue = []
    y_values_blue = []
    unaligned = 0
    aligned = 0
    for key in sbs18_dict:
        try:
            id83_dict[key]
        except KeyError:
            unaligned += 1
            continue
        if sbs18_dict[key] < min_sbs18_mutations or id83_dict[key] < min_id83_mutations:
            x_values.append(sbs18_dict[key])
            y_values.append(id83_dict[key])
        else:
            x_values_blue.append(sbs18_dict[key])
            y_values_blue.append(id83_dict[key])
            aligned += 1
    
    # Calculate the Pearson correlation coefficient and the p-value for the blue points
    corr, p_value = pearsonr(x_values_blue, y_values_blue)
    print(f'Pearson correlation: {corr}')
    print(f'p-value: {p_value}')
    print(f'Number of aligned samples: {aligned}')
    
    # log2 transform x and y
    x_values = [math.log2(x) for x in x_values]
    y_values = [math.log2(y) for y in y_values]
    x_values_blue = [math.log2(x) for x in x_values_blue]
    y_values_blue = [math.log2(y) for y in y_values_blue]
    
    # Fit a line to the data
    m, b = np.polyfit(x_values_blue, y_values_blue, 1)
    
    # Create the plot
    medium_blue = (0.0, 0.5, 1.0)  # RGB values
    sns.set(style="whitegrid")
    plt.scatter(x_values, y_values, color='gray', alpha=0.6, s=20)  # Plot the gray points
    plt.scatter(x_values_blue, y_values_blue, color=medium_blue, alpha=0.6, s=20)  # Plot the blue points
    plt.plot(x_values_blue, [m*x + b for x in x_values_blue], color='red')  # Plot the line
    plt.xlabel('SBS18', fontsize=14)
    plt.ylabel('ID83', fontsize=14)
    
    # Display the Pearson correlation, p-value, and number of aligned samples on the plot
    textstr = f'Pearson r: {corr:.2f}\np-value: {p_value:.2e}\nn={aligned}'
    props = dict(boxstyle='round', facecolor='lightgray', alpha=0.4)
    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=12,
        verticalalignment='top', bbox=props)
    
    # Make the plot square and keep the axes the same
    plt.gca().set_aspect('equal', adjustable='box')
    
    # Set the x and y ticks to be the same
    max_limit = max(max(x_values_blue + x_values), max(y_values_blue + y_values)) * 1.1
    min_limit = min(min(x_values_blue + x_values), min(y_values_blue + y_values)) * 0.9
    plt.xlim([min_limit, max_limit])
    plt.ylim([min_limit, max_limit])

    # Set the tick marks
    ticks = np.linspace(min_limit, max_limit, num=6)
    plt.xticks(ticks)
    plt.yticks(ticks)

    plt.show()

snp = convert_tsv_to_dict("/media/cam/Working/PCAWG/PCAWG_SBS18_Analysis/SNPs.tsv")
indel = convert_tsv_to_dict("/media/cam/Working/PCAWG/PCAWG_SBS18_Analysis/INDELs.tsv")

sbs18_dict = create_sig_dict(snp, signature="SBS18", percentage=0.00)
id83_dict = create_sig_dict(indel, signature="ID83", percentage=0.00)

# Call the function with your dictionaries
plot_dictionaries(sbs18_dict, id83_dict, min_sbs18_mutations=20, min_id83_mutations=20)
signature_correlation_dict = plot_pearson_scatter_sigs(snp, id83_dict)
# plot_correlation_coefficients(signature_correlation_dict, filter_p_value=True)