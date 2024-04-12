import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

def sig_to_df(file, mut_number):
    df = pd.read_csv(file, sep='\t', index_col=0)
    
    # Drop the last column
    df = df.drop(df.columns[-1], axis=1)
    
    # Multiply every row by mut_number and round
    df = (df * mut_number).round(0)
    
    return df

def df_difference(df1, df2):
    return (df1 - df2)

def plot_indel(df):
    # Extract the mutation details from the MutationType
    df['MutationDetail'] = df.index.str.split(':').str[1:3].str.join(':')

    df['ID83A'] = (df['ID83A'] / df['ID83A'].sum()) * 100
    
    # Define a color map for different mutation details
    color_map = {
        'Del:C': 'blue',
        'Del:T': 'cyan',
        'Ins:C': 'red',
        'Ins:T': 'pink',
        'Del:R': 'green',
        'Ins:R': 'yellow',
        'Del:M': 'purple'
    }
    
    # Set the size of the figure
    plt.figure(figsize=(15, 10))
    
    # Plot each group with its corresponding color
    for mutation_type in df.index:
        mutation_detail = df.loc[mutation_type, 'MutationDetail']
        plt.bar(mutation_type, df.loc[mutation_type, 'ID83A'], color=color_map[mutation_detail], edgecolor='black', label=mutation_detail)
    
    # Rotate x labels for better visibility
    plt.xticks(rotation=90)
    
    # Set title, labels, and legend with larger fonts
    plt.title('Mutation Type Differences', fontsize=18)
    plt.xlabel('Mutation Type', fontsize=16)
    plt.ylabel('Percentage (%)', fontsize=16)
    
    # Create a custom legend without duplicate labels
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), title='Mutation Detail', fontsize=12)
    
    # Display the plot
    plt.tight_layout()
    plt.show()

def extract_mutation_types(df):
    # Extract the mutation type directly from the index
    mutation_types = df.index.str.extract(r'\[(.*?)\]')[0]
    print(mutation_types)

def plot_sbs(df):
    fig, ax = plt.subplots()  # Define the figure and axis here

    # Extract the mutation type directly from the index
    df['MutationType'] = df.index

    df['SBS96A'] = (df['SBS96A'] / df['SBS96A'].sum()) * 100

    mutation_types = df['MutationType'].str.extract(r'\[(.*?)\]')[0]
    
    # Create a new column for sorting purposes
    df['SortOrder'] = df['MutationType'].str.replace(r'\[.*?\]', '', regex=True)
    
    # Prepend a prefix based on the major mutation type
    prefix_map = {
        'C>A': '1',
        'C>G': '2',
        'C>T': '3',
        'T>A': '4',
        'T>C': '5',
        'T>G': '6'
    }
    df['MajorSort'] = mutation_types.map(prefix_map)
    df['SortOrder'] = df['MajorSort'] + df['SortOrder']
    
    # Sort the dataframe
    df = df.sort_values(by=['SortOrder'])
    
    # Define a color map for different mutation types
    color_map = {
        'C>A': 'cornflowerblue',
        'C>G': 'black',
        'C>T': 'red',
        'T>A': 'gray',
        'T>C': 'green',
        'T>G': 'pink'
    }
    
    # Create a set to keep track of labels that have been added
    added_labels = set()

    # Modify x-axis tick labels to include colors
    tick_labels = []
    for mutation in df['MutationType']:
        context = mutation.split('[')[0] + mutation.split('>')[0][-1] + mutation.split(']')[1]
        tick_labels.append(context)

        mutation_type = mutation.split('[')[1].split(']')[0]
        label = mutation_type if mutation_type not in added_labels else ""
        ax.bar(mutation, df.loc[mutation, 'SBS96A'], color=color_map[mutation_type], label=label)
        added_labels.add(mutation_type)

    ax.set_xticks(range(len(df.index)))
    ax.set_xticklabels(tick_labels, rotation=90, ha='center', fontname="Courier New", fontweight='bold', fontsize=8)
        
    # Add legend, title, and labels
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys())
    ax.set_title('SBS96A Mutations')
    ax.set_xlabel('Mutation Type')
    ax.set_ylabel('Percentage (%)')
    plt.tight_layout()
    plt.show()

polh = '/home/cam/Documents/repos/SigProfiler/polh/ID83/Suggested_Solution/ID83_De-Novo_Solution/Signatures/ID83_De-Novo_Signatures.txt'
wt = '/home/cam/Documents/repos/SigProfiler/wt/ID83/Suggested_Solution/ID83_De-Novo_Solution/Signatures/ID83_De-Novo_Signatures.txt'
plot_indel(df_difference(sig_to_df(polh, 1742), sig_to_df(wt, 1074)))

# polh = '/home/cam/Documents/repos/SigProfiler/polh/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt'
# wt = '/home/cam/Documents/repos/SigProfiler/wt/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt'
# plot_sbs(df_difference(sig_to_df(polh, 56496), sig_to_df(wt, 45881)))

hmces = '/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/SigProfiler/hmces_output/ID83/Suggested_Solution/ID83_De-Novo_Solution/Signatures/ID83_De-Novo_Signatures.txt'
wt = '/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/SigProfiler/wt_output/ID83/Suggested_Solution/ID83_De-Novo_Solution/Signatures/ID83_De-Novo_Signatures.txt'
plot_indel(df_difference(sig_to_df(hmces, 1048), sig_to_df(wt, 837)))

# hmces = '/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/SigProfiler/hmces_output/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt'
# wt = '/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/SigProfiler/wt_output/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt'
# plot_sbs(df_difference(sig_to_df(hmces, 41115), sig_to_df(wt, 31254)))
