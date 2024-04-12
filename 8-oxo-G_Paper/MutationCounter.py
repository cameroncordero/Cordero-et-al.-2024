import pandas as pd
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def parse_files(directory_path):
    # Define your groups
    groups = {
        "WT Treated": ['_4_', '_5_', '_9_', '_10_'],
        "WT Untreated": ['_2_', '_3_', '_7_', '_8_'],
        "hmces-/- Treated": ['_14_', '_15_', '_19_', '_20_', '_24_', '_25_'],
        "hmces-/- Untreated":['_12_', '_13_', '_17_', '_18_', '_22_', '_23_']
    }

    # Resulting data frame
    all_data = []

    # Iterate through all files in the specified directory
    for file_path in Path(directory_path).glob('*'):
        file_str = str(file_path)  # Convert file path to string to check conditions
        
        # Determine the group of the current file
        group = None
        for key, values in groups.items():
            if any(num in file_str for num in values):
                group = key
                break
        
        if group:
            snp_count = 0
            indel_count = 0
            with open(file_path, 'r') as file:
                for line in file:
                    if line.startswith('#'):
                        continue  # Skip comments
                    else:
                        # Assuming columns are separated by tabs
                        columns = line.strip().split('\t')
                        ref, alt = columns[2], columns[3]  # Extract relevant data
                        
                        # Check if SNP or INDEL
                        if len(ref) == 1 and len(alt) == 1:
                            snp_count += 1
                        else:
                            indel_count += 1
            
            # Add file's data to the collection
            all_data.append([file_str, group, snp_count, indel_count])

    # Create DataFrame
    df = pd.DataFrame(all_data, columns=['Filename', 'Group', 'SNV', 'INDEL'])
    return df

# Using the function
# data = parse_files('/media/cam/Working/8-oxodG/mutation_counter')  # replace with your actual directory path
# data.to_csv('/media/cam/Working/8-oxodG/mutation_counter/parsed_data.tsv', index=False, sep='\t')  # Saving the data as a CSV file

def plot_data(dataframe):
    # Setting the overarching theme of the plots
    sns.set(style="whitegrid")

    # Define the order and custom color palette
    category_order = ["WT Untreated", "WT Treated", "hmces-/- Untreated", "hmces-/- Treated"]
    custom_palette = {
        "WT Untreated": "skyblue",  # or any color you prefer
        "WT Treated": "skyblue",
        "hmces-/- Untreated": "salmon",  # or any color you prefer
        "hmces-/- Treated": "salmon"
    }

    # Create a figure to specify the layout
    fig, axs = plt.subplots(2, 1, figsize=(6, 6))  # Adjusting the figure size for better clarity

    # Violin plot for SNP
    violin1 = sns.violinplot(x='Group', y='SNP', data=dataframe, order=category_order, palette=custom_palette, inner="point", ax=axs[0])
    axs[0].set_title('SNP counts')

    # Violin plot for INDEL
    violin2 = sns.violinplot(x='Group', y='INDEL', data=dataframe, order=category_order, palette=custom_palette, inner="point", ax=axs[1])
    axs[1].set_title('INDEL counts')

    # Set labels and other aesthetics
    for ax in axs.flat:
        ax.set_ylabel('Count')  # y-axis label
        ax.set_xlabel('Group')  # x-axis label
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')  # Rotate x-axis labels if they're long

    # Enhance legibility
    plt.tight_layout()
    plt.show()

def plot_data_with_horizontal_mean(dataframe):
    # Setting the overarching theme of the plots
    sns.set(style="whitegrid")

    # Define the order and custom color palette
    category_order = ["WT Untreated", "WT Treated", "hmces-/- Untreated", "hmces-/- Treated"]
    custom_palette = {
        "WT Untreated": "skyblue",
        "WT Treated": "skyblue",  # Same color for WT
        "hmces-/- Untreated": "salmon",
        "hmces-/- Treated": "salmon"  # Same color for hmces-/-
    }

    # Create a figure to specify the layout
    fig, axs = plt.subplots(2, 1, figsize=(4, 8))  # Adjusting the figure size for better clarity

    # Create a dot plot for SNPs, with consistent colors for each condition
    sns.stripplot(x='Group', y='SNP', data=dataframe, order=category_order, palette=custom_palette, jitter=0.25, size=8, ax=axs[0], edgecolor='black')

    # Calculate means and plot horizontal lines for SNPs
    means = dataframe.groupby('Group')['SNP'].mean().reindex(category_order)
    for i, mean in enumerate(means):
        # Draw each mean line with increased thickness (e.g., 2.5) and on top of the other elements (e.g., zorder=3)
        axs[0].hlines(mean, i - 0.25, i + 0.25, color='black', linewidth=2.5, zorder=3)
    axs[0].set_title('SNP counts')

    # Set y-axis to log scale for SNP plot
    axs[0].set_yscale('log')

    # Set major ticks at 1,000, 10,000, and 100,000
    axs[0].set_yticks([1000, 10000, 100000], minor=False)

    # Set minor ticks at each order of magnitude between 1,000 and 100,000
    axs[0].set_yticks([2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
                       20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000,
                       200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000], minor=True)

    # Label major ticks with the respective values and minor ticks with empty strings
    axs[0].set_yticklabels(['1,000', '10,000', '100,000'], minor=False)
    axs[0].set_yticklabels([''] * len(axs[0].get_yticks(minor=True)), minor=True)
    # Create a dot plot for INDELs, with consistent colors for each condition
    sns.stripplot(x='Group', y='INDEL', data=dataframe, order=category_order, palette=custom_palette, jitter=0.25, size=8, ax=axs[1], edgecolor='black')

    # Calculate means and plot horizontal lines for INDELs
    means = dataframe.groupby('Group')['INDEL'].mean().reindex(category_order)
    for i, mean in enumerate(means):
        # Draw each mean line with increased thickness (e.g., 2.5) and on top of the other elements (e.g., zorder=3)
        axs[1].hlines(mean, i - 0.25, i + 0.25, color='black', linewidth=2.5, zorder=3)
    axs[1].set_title('INDEL counts')

    # Create legend handles manually
    handles = [
        mpatches.Patch(color='skyblue', label='WT'),
        mpatches.Patch(color='salmon', label='hmces-/-')
    ]

    # Set the legend on the figure with custom handles
    fig.legend(handles=handles, loc='upper left', title="Genotype")

    # Set labels and other aesthetics
    for ax in axs.flat:
        ax.set_ylabel('Count')  # y-axis label
        ax.set_xlabel('')  # Remove x-axis label
        # Explicitly set the x-tick labels based on your category_order
        ax.set_xticklabels(['NT', 'KBrO3', 'NT', 'KBrO3'])
    
    # Enhance legibility
    plt.tight_layout()
    plt.show()

def count_mutations_in_vcfs(directory_path):

    # Iterate through all files in the specified directory
    data = {}
    for file_path in Path(directory_path).glob('*.vcf'):
        snv_count = 0
        indel_count = 0
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue  # Skip comments
                else:
                    # Assuming columns are separated by tabs
                    columns = line.strip().split('\t')
                    ref, alt = columns[3], columns[4]  # Extract relevant data
                    
                    # Check if SNP or INDEL
                    if len(ref) == 1 and len(alt) == 1:
                        snv_count += 1
                    else:
                        indel_count += 1
        # Add file's data to the collection
        data[file_path.name] = {'SNV': snv_count, 'INDEL': indel_count}
    
    data = pd.DataFrame(data).T

    # Add a total row
    data.loc['Total'] = data.sum()

    data.to_csv(Path(directory_path)/'mutation_counts.txt', index=True, sep='\t')  # Saving the data as a CSV file
    return data

count_mutations_in_vcfs('/media/cam/Working/8-oxodG/polh/XP_polh')  # replace with your actual directory path