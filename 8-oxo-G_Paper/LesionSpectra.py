from pathlib import Path
import pandas as pd

def reverse_complement(sequence):
    """Return the reverse complement of a sequence."""
    comp_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(comp_map[b] for b in reversed(sequence))

def count_lesion_context(lesion_bed, use_percentage=True):
    """
    Count the number of lesions in each context (5'-N-3').

    :param lesion_bed: Path to the BED file of lesions
    :param use_percentage: Boolean to decide between percentages (True) or raw counts (False)
    :return: None
    """
    lesion_bed = Path(lesion_bed)
    context_dict = {}
    
    with open(lesion_bed) as f:
        for line in f:
            tsv = line.strip().split('\t')
            context = tsv[6]
            if 'N' in context:
                continue
            strand = tsv[5]
            context = reverse_complement(context)
            mutation = 'C>A'
            key = context + '_' + mutation
            context_dict[key] = context_dict.setdefault(key, 0) + 1

    suffix = '_context_percentages.txt' if use_percentage else '_context_counts.txt'
    dictionary_to_file(context_dict, lesion_bed.with_name(lesion_bed.stem + suffix), use_percentage)

def dictionary_to_file(context_dict, path, use_percentage):
    """
    Write the lesion context data (percentages or raw counts) to a file.

    :param context_dict: A dictionary of lesion counts for each context
    :param use_percentage: Boolean to decide between percentages (True) or raw counts (False)
    :return: None
    """
    headers = ["C>A"]
    contexts = ["ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", 
                "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT"]

    df = pd.DataFrame(columns=headers, index=contexts)

    if use_percentage:
        total_counts = sum(context_dict.values())

        for key, count in context_dict.items():
            trinucleotide, mutation = key.split('_')
            percentage = count / total_counts * 100  # Calculate percentage
            df.at[trinucleotide, mutation] = percentage
    else:
        for key, count in context_dict.items():
            trinucleotide, mutation = key.split('_')
            df.at[trinucleotide, mutation] = count

    with open(path, 'w') as f:
        f.write(df.to_csv(sep='\t', index_label=False))

def compute_difference(df1_path, df2_path, output_path):
    """
    Compute the difference between the entries of two dataframes with the same keys and save the result.

    :param df1_path: Path to the first dataframe file.
    :param df2_path: Path to the second dataframe file.
    :param output_path: Path to save the resulting difference dataframe.
    :return: None
    """
    
    # Read the dataframes
    df1 = pd.read_csv(df1_path, sep='\t', index_col=0)
    df2 = pd.read_csv(df2_path, sep='\t', index_col=0)

    # Ensure they have the same structure
    if df1.columns.equals(df2.columns) and df1.index.equals(df2.index):
        # Compute the difference
        diff_df = df1 - df2
        # Save the difference dataframe
        diff_df.to_csv(output_path, sep='\t')
    else:
        print("The dataframes have different structures. Cannot compute difference.")

def filter_sbs96(input_filepath, save_path):
    """
    Filters the SBS96 dataset to retain only [C>A] mutations from the SBS96A column.
    The percentages are then recalculated to sum up to 100 and saved to a specified file.

    Args:
    - input_filepath (str): Path to the tab-separated SBS96 dataset.
    - save_path (str): Path where the resulting data will be saved.

    Returns:
    - None
    """
    # Read the data into a dataframe
    df = pd.read_csv(input_filepath, sep='\t', index_col=0)

    # Filter for [C>A] mutations and grab the SBS96A column
    ca_mutations = df[df.index.str.contains('\[C>A\]')].iloc[:, 0]

    # Rename the index to the format XCY
    ca_mutations.index = ca_mutations.index.str[0] + 'C' + ca_mutations.index.str[-1]

    # Normalize the percentages
    ca_mutations_percentage = ca_mutations / ca_mutations.sum() * 100

    # Save the resulting data
    ca_mutations_percentage.to_csv(save_path, sep='\t', header=True)

# Example usage
# count_lesion_context('/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/64-66_untreated_cellular/SRR_64-65-66.bed', use_percentage=True)
# count_lesion_context('/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/67-68/SRR_67-68.bed', use_percentage=True)
# count_lesion_context('/media/cam/Working/8-oxodG/lesion_files/bed/SRR_untreated_cellular_64-65-66_fixed.bed', use_percentage=True)

compute_difference( df1_path= '/media/cam/Working/8-oxodG/lesion_files/comparison_to_sbs/sbs96.txt',
                    df2_path= '/media/cam/Working/8-oxodG/lesion_files/comparison_to_sbs/SRR_treated_cellular_69-70_fixed_context_percentages.txt',
                    output_path= '/media/cam/Working/8-oxodG/lesion_files/comparison_to_sbs/kbro3-lesions.txt')

# filter_sbs96(input_filepath='/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/comparison_to_sbs/cog.sanger.ac.uk_cosmic-signatures-production_documents_v3.3_SBS18_PROFILE.txt',
#              save_path='/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/comparison_to_sbs/sbs18.txt')
# filter_sbs96(input_filepath='/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/comparison_to_sbs/SBS96_De-Novo_Signatures.txt',
#              save_path='/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/comparison_to_sbs/sbs96.txt')