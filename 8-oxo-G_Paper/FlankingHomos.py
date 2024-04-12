import re
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pathlib import Path
import glob
from scipy.stats import fisher_exact
import numpy as np
import os
from itertools import groupby
import pandas as pd
from itertools import product

def process_chromosome(context):
    counts = {'A': {'5_prime': {}, '3_prime': {}}, 'T': {'5_prime': {}, '3_prime': {}}, 'C': {'5_prime': {}, '3_prime': {}}, 'G': {'5_prime': {}, '3_prime': {}}}
    context = context.upper()
    for i in range(5, len(context) - 5):
        base = context[i]
        if base in ['A', 'T', 'C', 'G'] and 'N' not in context[i-5:i+6]:
            five_prime = context[i-5:i]
            three_prime = context[i+1:i+6]
            counts[base]['5_prime'].setdefault(five_prime, 0)
            counts[base]['5_prime'][five_prime] += 1
            counts[base]['3_prime'].setdefault(three_prime, 0)
            counts[base]['3_prime'][three_prime] += 1
    return counts

def calculate_background_context(fasta_file):
    fasta_file = Path(fasta_file)
    total_counts = {'A': {'5_prime': {}, '3_prime': {}}, 'T': {'5_prime': {}, '3_prime': {}}, 'C': {'5_prime': {}, '3_prime': {}}, 'G': {'5_prime': {}, '3_prime': {}}}
    with open(fasta_file, 'r') as f:
        context = ''
        for line in f:
            if line.startswith('>'):
                if context:  # if context is not empty, process the previous chromosome
                    counts = process_chromosome(context)
                    for base in ['A', 'T', 'C', 'G']:
                        for flank in ['5_prime', '3_prime']:
                            for sequence, count in counts[base][flank].items():
                                total_counts[base][flank].setdefault(sequence, 0)
                                total_counts[base][flank][sequence] += count
                    context = ''  # reset context for the next chromosome
            else:
                context += line.strip()
        if context:  # process the last chromosome
            counts = process_chromosome(context)
            for base in ['A', 'T', 'C', 'G']:
                for flank in ['5_prime', '3_prime']:
                    for sequence, count in counts[base][flank].items():
                        total_counts[base][flank].setdefault(sequence, 0)
                        total_counts[base][flank][sequence] += count

    # write the totals to a file
    for base in ['A', 'T', 'C', 'G']:
        with open(fasta_file.parent / f'{base}_flanking_context_counts.txt', 'w') as f:
            for flank in ['5_prime', '3_prime']:
                for sequence, count in total_counts[base][flank].items():
                    f.write(f"{flank}\t{sequence}\t{count}\n")

def reformat_apobec_data(input_file):
    input_file = Path(input_file)
    output_file = input_file.parent / (input_file.stem + '_reformatted.bed')

    with open(input_file, 'r') as f, open(output_file, 'w') as out:
        header = f.readline().strip().split('\t')  # Read the header
        indices = {name: index for index, name in enumerate(header)}  # Create a dictionary of field indices

        for line in f:
            fields = line.strip().split('\t')
            chrom = fields[indices['Chromosome']]
            start = int(fields[indices['Start_position']]) - 1  # 0-based start
            end = fields[indices['End_position']]
            ref = fields[indices['Reference_Allele']]
            alt = fields[indices['Tumor_Seq_Allele2']]
            context = fields[indices['"CONTEXT(+/-100)"']]
            context = context[:99] + context[99].upper() + context[100:]  # Capitalize the base before the deleted one
            var_type = fields[indices['Variant_Type']]
            if alt == '-':  # Deletion
                alt = '<DEL>'
            elif len(ref) > 1:  # Insertion
                alt = '<INS>'
            else:  # SNP
                alt = f'{ref}>{alt}'
            out.write(f'{chrom}\t{start}\t{end}\t.\t0\t+\t{context}\t{alt}\t{var_type}\t\n')

def reverse_complement(seq: str):
    '''
    Takes a nucleotide string in IUPAC and regular format and returns the reverse complement
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                        "R":"Y", "Y":"R", "S":"S", "W":"W", "K":"M",
                        "M":"K","B":"V", "D":"H", "H":"D", "V":"B",
                        "N":"N"}
    return "".join(complement.get(base, base) for base in reversed(seq))

def is_repeat(sequence):
    """Return the length of the homopolymer run at the front end of the sequence."""
    return max(len(list(g)) for _, g in groupby(sequence))

def calculate_homopolymer_run(sequence, flank):
    """Calculate the length of the homopolymer run at the start of the sequence for the 3' side and at the end for the 5' side."""
    if not sequence:
        return 0
    if flank == '5_prime':
        sequence = sequence[::-1]  # Reverse the sequence for 5' side
    run_length = 1
    for i in range(1, len(sequence)):
        if sequence[i] == sequence[0]:
            run_length += 1
        else:
            break
    return run_length

# This function loads the context counts from the files in the specified directory.
def load_context_counts(base_of_interest):
    # Initialize an empty dictionary to store the total counts for each context.
    total_counts = {base_of_interest: {'5_prime': {}, '3_prime': {}}}
    # Initialize an empty dictionary to store the repeat contexts for each flank and homopolymer run length.
    repeat_contexts = {'5_prime': {i: set() for i in range(1, 6)}, '3_prime': {i: set() for i in range(1, 6)}}

    # Loop over each file in the specified directory.
    for filename in glob.glob('/home/cam/Documents/genome_files/hg19/*_flanking_context_counts.txt'):
        # Extract the base from the filename.
        base = os.path.basename(filename).split('_')[0].upper()

        # If the base is 'G' or 'C', open the file and read the counts.
        base_list = [base_of_interest, reverse_complement(base_of_interest)]
        if base in base_list:
            with open(filename, 'r') as f:
                for line in f:
                    # Split the line into flank, sequence, and count.
                    flank, sequence, count = line.strip().split('\t')
                    sequence = sequence.upper()

                    # If the base is 'C', reverse complement the sequence, swap the flank, and change the base to 'G'.
                    if base == base_list[1]:
                        sequence = reverse_complement(sequence)
                        flank = '3_prime' if flank == '5_prime' else '5_prime'  # Swap the flank
                        base = base_of_interest

                    # Update the total counts dictionary with the count for the current base, flank, and sequence.
                    total_counts[base][flank][sequence] = int(count)

                    # Calculate the homopolymer run length for the current sequence and flank.
                    repeat_length = calculate_homopolymer_run(sequence, flank)

                    # If the homopolymer run length is between 1 and 5, add the sequence to the repeat contexts dictionary.
                    if 1 <= repeat_length <= 5:
                        repeat_contexts[flank][repeat_length].add(sequence)

    # Return the total counts and repeat contexts.
    return total_counts, repeat_contexts

####################
# left off below here
# need to change it do that if the base is in a homopolymer it willstill count it, however it needs to shift the direction of it a different way so it left justifies it the same


def handle_mutation(line, observed_counts, total_deletions, base_of_interest):
    tsv = line.strip().split('\t')
    var_type = tsv[8].upper()
    context = tsv[6]
    base_list = [base_of_interest, reverse_complement(base_of_interest)]
    if 'DEL' in var_type:
        del_refs = re.findall(r'[A-Z]+', context)
        deleted_base = del_refs[0][1]
        # if deletion length is greater than 1, skip
        if len(del_refs[0][1:]) > 1: # skip if there are more than 1 deletions
            return observed_counts, total_deletions
        # if reference base for deletion is same as deleted base, skip (shouldnt ever happen because of 5' justification)
        if del_refs[0][0] == deleted_base:
            return observed_counts, total_deletions
        # make sure deletion isn't part of a homopolymer by checking the right half of the string
        if re.split(r'[A-Z]+', context)[1][0].upper() == deleted_base:
            return observed_counts, total_deletions
        del_start_index = context.index(deleted_base)
        del_end_index = del_start_index + 1
        context = context.upper() # convert to uppercase so that the reference base is also uppercase and everything is consistent
        
        # Check if the deletion is at the end of the context string. If it is not, and the base following the deletion is different from the deleted base, then increment the total deletions count and set del_base to the deleted base.
        if del_end_index < len(context) and context[del_end_index] != deleted_base:
            total_deletions += 1
            del_base = deleted_base

            # If the deleted base is not 'G' or 'C', return the current counts and total deletions without making any changes.
            if del_base not in base_list:
                return observed_counts, total_deletions

            # If the deleted base is 'C', reverse and complement the context string, and set del_base to 'G'. Also, recalculate the start and end indices of the deletion in the reversed context.
            if del_base == base_list[1]:
                context = reverse_complement(context)
                del_base = base_of_interest
                del_start_index = len(context) - del_end_index
                del_end_index = del_start_index + 1

            # Get the 5 bases before and after the deletion in the context string.
            pre_del_context = context[:del_start_index][-5:]
            post_del_context = context[del_end_index:][:5]

            # For each flank (5' and 3'), calculate the length of the homopolymer run and set context_part to the sequence of bases in the flank.
            for flank in ['5_prime', '3_prime']:
                if flank == '5_prime':
                    homopolymer_run = calculate_homopolymer_run(pre_del_context[::-1], flank)
                    context_part = pre_del_context[::-1]
                else:
                    homopolymer_run = calculate_homopolymer_run(post_del_context, flank)
                    context_part = post_del_context

                # Update the observed_counts dictionary with the count of the current deletion. The dictionary is nested with the structure: base -> flank -> homopolymer run -> context -> count.
                observed_counts[del_base][flank].setdefault(homopolymer_run, {})
                observed_counts[del_base][flank][homopolymer_run].setdefault(context_part, 0)
                observed_counts[del_base][flank][homopolymer_run][context_part] += 1

        # Return the updated counts and total deletions.
        return observed_counts, total_deletions
    return observed_counts, total_deletions

def calculate_observed_counts(mutation_file, base_of_interest):
    observed_counts = {base_of_interest: {'5_prime': {}, '3_prime': {}}, 'C': {'5_prime': {}, '3_prime': {}}}
    total_deletions = 0
    with open(mutation_file, 'r') as f:
        for line in f:
            observed_counts, total_deletions = handle_mutation(line, observed_counts, total_deletions, base_of_interest)
    # print(observed_counts)
    return observed_counts, total_deletions

# This function calculates the enrichment scores for each context.
def calculate_enrichment_scores(observed_counts, total_counts, total_deletions, repeat_contexts, base_of_interest):
    # Initialize an empty dictionary to store the enrichment scores.
    enrichment_scores = {}

    # Loop over each base. In this case, only 'G' is considered.
    for base in base_of_interest:
        # Loop over each flank (5' and 3').
        for flank in ['5_prime', '3_prime']:
            # Loop over each possible homopolymer run length (1 to 5).
            for homopolymer_run in range(1, 6):
                # Calculate the observed count (O) for the current base, flank, and homopolymer run length.
                # This is the sum of the observed counts for all runs of length homopolymer_run or more.
                O = sum(sum(observed_counts[base][flank].get(run, {}).values()) for run in range(homopolymer_run, 6))

                # Calculate the total count (B) for the current base, flank, and homopolymer run length.
                # This is the sum of the total counts for all contexts in repeat_contexts[flank][homopolymer_run].
                B = sum(total_counts[base][flank].get(context, 0) for context in repeat_contexts[flank][homopolymer_run])

                # Calculate the observed count (O_non) for the current base, flank, and runs of length less than homopolymer_run.
                O_non = sum(sum(observed_counts[base][flank].get(run, {}).values()) for run in range(1, homopolymer_run))

                # Calculate the total count (B_non) for the current base, flank, and all contexts not in repeat_contexts[flank][homopolymer_run].
                B_non = sum(total_counts[base][flank].get(context, 0) for context in set(total_counts[base][flank].keys()) - set(repeat_contexts[flank][homopolymer_run]))

                # Calculate the total number of repeat and non-repeat sites.
                total_repeat_sites = B
                total_non_repeat_sites = B_non

                # Calculate the expected counts (E and E_non) for the current base, flank, and homopolymer run length.
                E = (B / (total_repeat_sites + total_non_repeat_sites)) * total_deletions
                E_non = (B_non / (total_repeat_sites + total_non_repeat_sites)) * total_deletions

                # Calculate the enrichment scores (ES and ES_non) for the current base, flank, and homopolymer run length.
                ES = O / E if E else 0
                ES_non = O_non / E_non if E_non else 0

                # Perform a Fisher's exact test to calculate the p-value for the observed and expected counts.
                _, p_value = fisher_exact([[O, B], [O_non, B_non]])

                # Store the enrichment scores and p-value in the enrichment_scores dictionary.
                enrichment_scores[(base, flank, homopolymer_run)] = {'Repeat_Enrichment_Score': ES, 'Non_Repeat_Enrichment_Score': ES_non, 'P_Value': p_value}

    # Return the enrichment scores.
    return enrichment_scores

def count_nucleotide_contributions(observed_counts, base_of_interest):
    # Initialize a dictionary to store the counts for each base, flank, and homopolymer run length.
    counts = {base_of_interest: {'5_prime': {i: {} for i in range(1, 6)}, '3_prime': {i: {} for i in range(1, 6)}}}

    # Loop over the observed counts dictionary.
    for base, flanks in observed_counts.items():
        for flank, runs in flanks.items():
            for run, contexts in runs.items():
                # For each context, get the first base (which is the base of the homopolymer run) and add the count to the corresponding entry in the counts dictionary.
                for context, count in contexts.items():
                    homopolymer_base = context[0]
                    counts[base][flank][run].setdefault(homopolymer_base, 0)
                    counts[base][flank][run][homopolymer_base] += count

    # Return the counts dictionary.
    return counts

def plot_enrichment_scores(enrichment_scores, nucleotide_counts, base_of_interest):
    # Create a figure with 4 subplots
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    axs = axs.ravel()

    # Flatten the enrichment scores for easier handling
    flat_scores = []
    for (base, flank, homopolymer_run), data in enrichment_scores.items():
        flat_scores.append({
            'base': base,
            'flank': flank,
            'homopolymer_run': homopolymer_run,
            'Repeat_Enrichment_Score': data['Repeat_Enrichment_Score'],
            'Non_Repeat_Enrichment_Score': data['Non_Repeat_Enrichment_Score'],
            'P_Value': data['P_Value']
        })

    # Convert to DataFrame for easier handling
    df = pd.DataFrame(flat_scores)

    # Define the colors for each nucleotide
    nucleotide_colors = {'a': 'tab:blue', 't': 'tab:orange', 'c': 'tab:green', 'g': 'tab:red'}

    # Plot each flank for 'G'
    for i, flank in enumerate(['5_prime', '3_prime']):
        ax = axs[i]
        df_subset = df[(df['base'] == base_of_interest) & (df['flank'] == flank)]
        x = np.arange(len(df_subset))
        width = 0.5  # Increase the width of the bars

        # Plot the enrichment scores
        ax.bar(x - width/2, df_subset['Repeat_Enrichment_Score'], width, label='Repeat', edgecolor='black', linewidth=1.2)
        ax.bar(x + width/2, df_subset['Non_Repeat_Enrichment_Score'], width, label='Non-Repeat', edgecolor='black', linewidth=1.2)

        ax.set_xlabel('Homopolymer Run Length')
        ax.set_ylabel('Enrichment Score')
        ax.set_title(f'Enrichment Scores for {base_of_interest} at {flank}')
        ax.set_xticks(x)

        # Change the labels of the x-axis ticks
        labels = [str(run) if run < 5 else '5+' for run in df_subset['homopolymer_run']]
        ax.set_xticklabels(labels)

        ax.legend()

        # Increase the thickness of the axis lines
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.5)

    # Plot the nucleotide counts
    for i, flank in enumerate(['5_prime', '3_prime']):
        ax = axs[i+2]
        # No filtering for 'G' base here
        df_subset = df[df['flank'] == flank]
        x = np.arange(len(df_subset))
        width = 0.5  # Increase the width of the bars

        # Plot the nucleotide counts
        bottom = np.zeros(len(df_subset))
        for nucleotide, color in nucleotide_colors.items():
            counts = [nucleotide_counts[base_of_interest][flank][homopolymer_run].get(nucleotide.upper(), 0) for homopolymer_run in df_subset['homopolymer_run'].values]
            ax.bar(x, counts, width, bottom=bottom, label=nucleotide.upper(), color=color, edgecolor='black', linewidth=1.2)
            bottom += counts

        ax.set_xlabel('Homopolymer Run Length')
        ax.set_ylabel('Count')
        ax.set_title(f'Nucleotide Counts for {base_of_interest} at {flank}')
        ax.set_xticks(x)

        # Change the labels of the x-axis ticks
        labels = [str(run) if run < 5 else '5+' for run in df_subset['homopolymer_run']]
        ax.set_xticklabels(labels)

        ax.legend()

        # Increase the thickness of the axis lines
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(1.5)

        # Increase the y-axis limit for more space at the top
        ax.set_ylim(0, max(bottom) * 1.1)

        # Add the total number of events at the top of each bar
        for i, total in enumerate(bottom):
            ax.text(i, total, str(int(total)), ha='center', va='bottom')

    fig.tight_layout()
    plt.show()

def calculate_context_enrichment(mutation_file, base_of_interest):
    # Convert the mutation_file to a Path object for easier handling.
    mutation_file = Path(mutation_file)

    # Load the total counts and repeat contexts from the context counts files.
    total_counts, repeat_contexts = load_context_counts(base_of_interest)

    # Calculate the observed counts and total number of deletions from the mutation file.
    observed_counts, total_deletions = calculate_observed_counts(mutation_file, base_of_interest)

    # Calculate the enrichment scores for each base, flank, and homopolymer run length.
    enrichment_scores = calculate_enrichment_scores(observed_counts, total_counts, total_deletions, repeat_contexts, base_of_interest)

    nucleotide_counts = count_nucleotide_contributions(observed_counts, base_of_interest)

    # Convert the enrichment scores to a DataFrame for easier handling.
    df = pd.DataFrame(enrichment_scores).T.reset_index()

    # Rename the columns of the DataFrame.
    df.columns = ['Base', 'Flank', 'Homopolymer_Run_Length', 'Repeat_Enrichment_Score', 'Non_Repeat_Enrichment_Score', 'P_Value']

    # Write the DataFrame to a TSV file.
    df.to_csv(mutation_file.with_name(mutation_file.stem + '_enrichment_scores.tsv'), sep='\t', index=False)

    # Return the enrichment scores.
    return enrichment_scores, nucleotide_counts

def main(mutation_file, base_of_interest):
    # Calculate enrichment scores
    enrichment_scores, nucleotide_counts = calculate_context_enrichment(mutation_file, base_of_interest)
    plot_enrichment_scores(enrichment_scores, nucleotide_counts, base_of_interest)

if __name__ == "__main__":
    # Assume the mutation file is stored in the variable mutation_file
    mutation_file = "/media/cam/Working/8-oxodG/hmces/flanks/HMCES_KBr.bed"
    main(mutation_file, 'C')