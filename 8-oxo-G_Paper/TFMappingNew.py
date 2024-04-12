import math
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.ndimage import gaussian_filter1d
import seaborn as sns
from scipy.interpolate import interp1d

def intersect_tf_sites_with_dnase(tf_binding_bed, dnase_bed):
    tf_binding_bed, dnase_bed = Path(tf_binding_bed), Path(dnase_bed)
    with open(tf_binding_bed) as tf, open(dnase_bed) as dnase:
        pass

def filter_intersected_tf_maps(intersect_1, intersect_2):
    tsv_sites = set()
    intersect_path = Path(intersect_1)
    with open(intersect_1) as f1, open(intersect_2) as f2:
        for line in f1:
            tsv = line.strip().split('\t')
            if int(tsv[2])-int(tsv[1]) == int(tsv[10]):
                tsv_sites.add('\t'.join(tsv[:7]))
        for line in f2:
            tsv = line.strip().split('\t')
            if int(tsv[2])-int(tsv[1]) == int(tsv[10]):
                tsv_sites.add('\t'.join(tsv[:7]))
    with open(intersect_path.with_name('DNAse_filtered.txt'), 'w') as o:
        for item in tsv_sites:
            o.write(item+'\n')
    return(intersect_path.with_name('DNAse_filtered.txt'))

def expand_motif_positions(tf_map):
    path = Path(tf_map)
    with open(tf_map, 'r') as f:
        with open(path.with_name(path.stem+'_expanded.bed'), 'w') as o:
            for line in f:
                tsv = line.strip().split('\t')
                pos0 = int(tsv[1])
                pos1 = int(tsv[2])
                rest = '\t'.join(tsv[3:])
                o.write(f'{tsv[0]}\t{str(pos0-1000)}\t{str(pos1+1000)}\t{rest}\n')
    return path.with_name(path.stem+'_expanded.bed')

def moving_average(data, window_size):
    kernel = np.ones(window_size) / window_size
    return np.convolve(data, kernel, mode='same')

def midpoint_tf_analysis(input_file, tf_map, window_size=5):
    sns.set_style("whitegrid")
    plt.rcParams['font.family'] = 'Arial'
    colors = sns.color_palette("Paired")

    input_path, tf_path = Path(input_file), Path(tf_map)
    output = tf_path.with_name(input_path.stem + '_' + tf_path.stem + 'intersect.txt')
    with subprocess.Popen(args=[f'bedtools intersect -wa -wb -a {input_file} -b {tf_map} > {output}'],
                          stdout=subprocess.PIPE, shell=True) as p:
        for text in p.stdout:
            print(text)
    print(f'bedtools intersect -wa -wb -a {input_file} -b {tf_map} > {output}')
    with open(output) as f:
        position_dict = {}
        for line in f:
            tsv = line.strip().split('\t')
            mut_start = int(tsv[1])
            tf_start = int(tsv[9])
            tf_end = int(tsv[10])
            mid_pt = math.floor((tf_end + tf_start) / 2)
            mut_position = mut_start - mid_pt
            position_dict[mut_position] = position_dict.setdefault(mut_position, 0) + 1

        sorted_data = sorted(position_dict.items())
        positions, counts = zip(*sorted_data)
        positions = np.array(positions)
        counts = np.array(counts)

        mask = (-1000 <= positions) & (positions <= 1000)
        filtered_positions = positions[mask]
        filtered_counts = counts[mask]

        counts_smooth = moving_average(filtered_counts, window_size)

        fig, ax = plt.subplots()

        ax.plot(filtered_positions, counts_smooth, color=colors[1], linewidth=3, label='Smoothed data')
        ax.scatter(filtered_positions, filtered_counts, color=colors[0], label='Original data', alpha=0.5)

        ax.set_xlabel('Position', fontsize=36)
        ax.set_ylabel('Counts', fontsize=36)
        ax.legend(loc='upper left', fontsize=24)
        ax.set_title('Comparison of Original and Smoothed Data', fontsize=36)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both', which='major', labelsize=24)
        ax.set_facecolor((0.95, 0.95, 0.95))
        ax.grid(True)

        fig.set_size_inches(16, 9)
        plt.savefig(tf_path.with_name(input_path.stem + '_' + tf_path.stem + '_intersect_smoothed.svg'), dpi=1000)
        plt.show()

def midpoint_tf_analysis_polyfit(input_file, tf_map, degree = 20):
    sns.set_style("whitegrid")
    plt.rcParams['font.family'] = 'Arial'
    colors = sns.color_palette("Paired")

    input_path, tf_path = Path(input_file), Path(tf_map)
    output = tf_path.with_name(input_path.stem + '_' + tf_path.stem + 'intersect.txt')
    with subprocess.Popen(args=[f'bedtools intersect -wa -wb -a {input_file} -b {tf_map} > {output}'],
                          stdout=subprocess.PIPE, shell=True) as p:
        for text in p.stdout:
            print(text)
    print(f'bedtools intersect -wa -wb -a {input_file} -b {tf_map} > {output}')
    with open(output) as f:
        position_dict = {}
        for line in f:
            tsv = line.strip().split('\t')
            mut_start = int(tsv[1])
            tf_start = int(tsv[9])
            tf_end = int(tsv[10])
            mid_pt = math.floor((tf_end + tf_start) / 2)
            mut_position = mut_start - mid_pt
            position_dict[mut_position] = position_dict.setdefault(mut_position, 0) + 1

        sorted_data = sorted(position_dict.items())
        positions, counts = zip(*sorted_data)
        positions = np.array(positions)
        counts = np.array(counts)

        mask = (-1000 <= positions) & (positions <= 1000)
        filtered_positions = positions[mask]
        filtered_counts = counts[mask]

        coeffs = np.polyfit(filtered_positions, filtered_counts, degree)
        poly = np.poly1d(coeffs)
        fitted_counts = poly(filtered_positions)

        fig, ax = plt.subplots()

        ax.plot(filtered_positions, fitted_counts, color=colors[1], linewidth=3, label=f'Polynomial fit (n={degree})')
        ax.scatter(filtered_positions, filtered_counts, color=colors[0], label='Original data', alpha=0.5)

        ax.set_xlabel('Position Relative to TF Binding Site Midpoint (bp)', fontsize=36)
        ax.set_ylabel('Mutation Counts', fontsize=36)
        ax.legend(loc='upper left', fontsize=24)
        # ax.set_title('IDK', fontsize=36)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both', which='major', labelsize=24)
        ax.set_facecolor((0.95, 0.95, 0.95))
        ax.grid(True)

        fig.set_size_inches(16, 9)
        plt.savefig(tf_path.with_name(input_path.stem + '_' + tf_path.stem + 'intersect_polyfit.svg'), dpi=1000)
        plt.show()

def melanoma_subset(input_file):
    path = Path(input_file)
    output = path.with_name(path.stem+'_subset.bed')
    with open(path) as f, open(output, 'w') as o:
        for line in f:
            pass # left off here

dnase1 = '/media/cam/Working/8-oxodG/Cam_calls_from_Data9/Analysis/TF_Binding/new_tf_map/TF_ENCFF278DHP_intersect.txt'
dnase2 = '/media/cam/Working/8-oxodG/Cam_calls_from_Data9/Analysis/TF_Binding/new_tf_map/TF_ENCFF742YZA_intersect.txt'
mutations = '/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/vcf_files/genotype_split/UV_nucleomutics/UV.mut'
old_map = '/media/cam/Data9/CortezAnalysis/Cam_calls/Analysis/TF_Binding/new_tf_map/non_filtered_TF_map_1000.bed'
melanoma = '/home/cam/Documents/UV_Data/MELA-AU_trinuc_context_mutations_filtered.bed6'

# with filtered data
filtered_file = filter_intersected_tf_maps(dnase1, dnase2)
expanded_file = expand_motif_positions(filtered_file)
# midpoint_tf_analysis(mutations, expanded_file, window_size=50)
# midpoint_tf_analysis_polyfit(mutations, expanded_file)


# # # with old map
# midpoint_tf_analysis(mutations, old_map, window_size= 50)
# midpoint_tf_analysis_polyfit(mutations, old_map, degree = 30)

# with melanoma data
midpoint_tf_analysis(mutations, expanded_file, window_size=50)
midpoint_tf_analysis_polyfit(mutations, expanded_file, degree=30)