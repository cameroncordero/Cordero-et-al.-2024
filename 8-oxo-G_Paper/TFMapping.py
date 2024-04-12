import re
import subprocess
import os
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.optimize import curve_fit
import random

def shorten_file(input_file, output_file):
    with open(input_file, 'r') as f:
        with open(output_file, 'w') as o:
            for line in f:
                tsv = line.strip().split('\t')
                new_line = '\t'.join(tsv[:8])
                o.write(new_line+'\n')

def mutations_in_sites(vcf_file, tf_map: dict):
    with open(vcf_file, 'r') as f:
        tf_pos_dict = {}
        for key in tf_map.keys():
            tf_pos_dict.setdefault(key)
        for line in f:
            if '#' in line: continue
            tsv = line.strip().split('\t')
            chrom = tsv[0]
            pos1 = int(tsv[1])-1
            pos2 = int(tsv[1])
            for key in tf_map.keys():
                if chrom in tf_map[key].keys():
                    for values in tf_map[key][chrom]:
                        for positions in values:
                            if pos1 - int(positions[0]) >= 0 and pos2 - int(positions[2]) <= 0:
                                mut_pos = pos1 - int(positions[0])
                                tf_pos_dict[key] = tf_pos_dict[key].setdefault(mut_pos, 0) +1
        return tf_pos_dict

def tf_mapping_proper(vcf_file, tf_map, tf_map_dict):
    tf_pos_dict = {}
    processed_lines = 0
    skipped_lines = 0
    map_lines = 0
    mutations_in_tf = 0
    for key in tf_map_dict.keys():
        tf_pos_dict[key] = {}
    with open(vcf_file, 'r') as vcf:
        with open(tf_map, 'r') as map:
            map_line = map.readline()
            map_lines += 1
            map_tsv = map_line.strip().split('\t')
            map_chrom = map_tsv[0][3:]
            map_pos0 = int(map_tsv[1])
            map_pos1 = int(map_tsv[2])
            map_tf = map_tsv[6]
            for vcf_line in vcf:
                if '#' in vcf_line[0]:
                    processed_lines += 1
                    continue
                vcf_tsv = vcf_line.strip().split('\t')
                vcf_chrom = vcf_tsv[0][3:]
                vcf_pos0 = int(vcf_tsv[1])-1
                vcf_pos1 = int(vcf_tsv[1])
                if '_' in vcf_chrom:
                    processed_lines += 1
                    continue
                
                ##### File Processing #####
                if vcf_chrom < map_chrom:
                    processed_lines += 1
                    continue
                elif vcf_chrom > map_chrom:
                    for line in map:
                        map_lines += 1
                        map_tsv = line.strip().split('\t')
                        map_chrom = map_tsv[0][3:]
                        map_pos0 = int(map_tsv[1])
                        map_pos1 = int(map_tsv[2])
                        map_tf = map_tsv[6]
                        if not vcf_chrom > map_chrom:
                            break
                if vcf_chrom == map_chrom:
                    if vcf_chrom == map_chrom and vcf_pos1 > map_pos1:
                        for line in map:
                            map_lines += 1
                            map_tsv = line.strip().split('\t')
                            map_chrom = map_tsv[0][3:]
                            map_pos0 = int(map_tsv[1])
                            map_pos1 = int(map_tsv[2])
                            map_tf = map_tsv[6]
                            if vcf_chrom != map_chrom:
                                break
                            if not vcf_pos1 > map_pos1:
                                break
                    if vcf_chrom == map_chrom and vcf_pos0 < map_pos0:
                        processed_lines += 1
                        continue
                    if vcf_chrom == map_chrom and vcf_pos0 >= map_pos0 and vcf_pos1 <= map_pos1:
                        mutations_in_tf += 1
                        tf_mut_pos = vcf_pos0 - map_pos0 + 1
                        tf_pos_dict[map_tf][tf_mut_pos] = tf_pos_dict[map_tf].setdefault(tf_mut_pos, 0) + 1
                        processed_lines += 1
                        continue
                    if vcf_chrom < map_chrom:
                        processed_lines += 1
                        continue
                skipped_lines += 1
                print(vcf_line)
                print(line)

    print(f'vcf processed lines = {processed_lines}, skipped lines = {skipped_lines}, total = {processed_lines+skipped_lines}')
    print(f'map lines = {map_lines}')
    print(mutations_in_tf)
    return tf_pos_dict

def count_mutations_in_each_site(tf_map_overlap: str, output_file: str):
    with open(tf_map_overlap, 'r') as f:
        with open(output_file, 'w') as o:
            motif_dict = {}
            for line in f:
                tsv = line.strip().split('\t')
                motif_dict[tsv[10]] = motif_dict.setdefault(tsv[10], 0) + 1
            o.write(f'Motif\tMutation Count\n')
            for key, value in motif_dict.items():
                o.write(f'{key}\t{value}\n')

def expand_motif_positions(tf_map, output_file):
    with open(tf_map, 'r') as f:
        with open(output_file, 'w') as o:
            for line in f:
                tsv = line.strip().split('\t')
                pos0 = int(tsv[1])
                pos1 = int(tsv[2])
                rest = '\t'.join(tsv[3:])
                o.write(f'{tsv[0]}\t{str(pos0-1000)}\t{str(pos1+1000)}\t{rest}\n')

def concat_pos_fasta_to_bed(tf_bed, tf_fasta, combined_bed):
    with open(tf_bed, 'r') as bed:
        with open(tf_fasta, 'r') as fa:
            with open(combined_bed, 'w') as o:
                for line in bed:
                    tsv = line.strip().split('\t')
                    fa_stuff = fa.readline()
                    faaaaaa = re.split(':|-', fa_stuff.strip().split('\t')[0])
                    fa_chrom_pos = '_'.join(faaaaaa)
                    bed_chrom_pos = '_'.join(tsv[:3])
                    if fa_chrom_pos == bed_chrom_pos:
                        o.write(line.strip()+'\t'+fa_stuff.strip().split('\t')[1].upper()+'\n')
                    else:
                        print('BROKE AF')

def calc_different_rates(file1, file2, output_file):
    with open(file1, 'r') as f1:
        with open(file2, 'r') as f2:
            with open(output_file, 'w') as o:
                mutation_counts_internal = {}
                mutation_counts_external = {}
                f1.readline()
                f2.readline()
                for line in f1:
                    tsv = line.strip().split('\t')
                    mutation_counts_internal[tsv[0]] = (tsv[1])
                for line in f2:
                    tsv = line.strip().split('\t')
                    mutation_counts_external[tsv[0]] = (tsv[1])
                o.write(f'Motif\tMutation rate internal\tMutation rate external\n')
                for key, internal in mutation_counts_internal.items():
                    external = (float(mutation_counts_external[key]) - float(internal))/1000
                    internal = float(internal)/20
                    o.write(f'{key}\t{str(internal)}\t{str(external)}\n')

def melanoma_subsetting(input_file, tf_map, output_file):
    with subprocess.Popen(args=[f'bedtools intersect -wa -wb -a {input_file} -b {tf_map} > {os.path.dirname(input_file)}/{output_file}_overlap.tsv'],
                        stdout=subprocess.PIPE, shell=True) as p:
        for text in p.stdout:
            print(text)
    with open(f'WorkingDir/{output_file}_overlap.tsv', 'r') as f:
        position_dict = {}
        for line in f:
            rando = random.randint(1,22961206)
            if rando > 413337: continue
            tsv = line.strip().split('\t')
            mut_start = int(tsv[1])
            tf_start = int(tsv[7])
            tf_end = int(tsv[8])
            mid_pt = math.floor((tf_end+tf_start)/2)
            mut_position = mut_start-mid_pt
            position_dict[mut_position] = position_dict.setdefault(mut_position, 0) + 1
        positions = list(position_dict.keys())
        counts = list(position_dict.values())
        positions = np.array(positions)
        counts = np.array(counts)
        model1 = np.poly1d(np.polyfit(positions, counts, 100))
        return model1(0)-model1(250)
        # polyline = np.linspace(positions.min(), positions.max())
        # plt.plot(polyline, model1(polyline), color='red', linewidth=3)
        # plt.scatter(positions, counts)
        # plt.show()

def midpoint_tf_analysis(input_file, tf_map, output_file):
    with subprocess.Popen(args=[f'bedtools intersect -wa -wb -a {input_file} -b {tf_map} > {os.path.dirname(input_file)}/{output_file}_overlap.tsv'],
                        stdout=subprocess.PIPE, shell=True) as p:
        for text in p.stdout:
            print(text)
    with open(f'WorkingDir/{output_file}_overlap.tsv', 'r') as f:
        position_dict = {}
        for line in f:
            tsv = line.strip().split('\t')
            mut_start = int(tsv[1])
            tf_start = int(tsv[9])
            tf_end = int(tsv[10])
            mid_pt = math.floor((tf_end+tf_start)/2)
            mut_position = mut_start-mid_pt
            position_dict[mut_position] = position_dict.setdefault(mut_position, 0) + 1
        positions = list(position_dict.keys())
        counts = list(position_dict.values())
        positions = np.array(positions)
        counts = np.array(counts)
        model1 = np.poly1d(np.polyfit(positions, counts, 100))
        polyline = np.linspace(positions.min(), positions.max())
        plt.plot(polyline, model1(polyline), color='red', linewidth=3)
        plt.scatter(positions, counts)
        plt.show()

def find_range(list_of_numbers):
    n_min = min(list_of_numbers)
    n_max = max(list_of_numbers)
    n_range = n_max - n_min
    return n_min, n_max, n_range

def find_median(list_of_numbers):
    list_of_numbers.sort()
    length = len(list_of_numbers)
    length_is_odd = True if length % 2 == 0 else False
    if length_is_odd:
        index = length//2
        median = list_of_numbers[index]
    else:
        index_1 = length//2
        index_2 = index_1 + 1
        median = (list_of_numbers[index_1] + list_of_numbers[index_2]) / 2
    return median

def find_mean(list_of_numbers):
    sum_n = sum(list_of_numbers)
    len_n = len(list_of_numbers)
    mean = sum_n/len_n
    return mean

# peak_height_list = []
# for _ in range(100):
#     peak_height_list.append(melanoma_subsetting(input_file='WorkingDir/MELA-AU_trinuc_context_mutations.bed',
#                                                 tf_map='WorkingDir/TF_map_1000.bed',
#                                                 output_file='MELA-AU_SNP_TF_map'
#                                                 )
#                             )

# print(find_range(peak_height_list))
# print(find_mean(peak_height_list))
# print(find_median(peak_height_list))


# midpoint_tf_analysis(   input_file='WorkingDir/KM_SNPs.bed',
#                         tf_map='WorkingDir/TF_map_1000.bed',
#                         output_file='KM_SNP_TF_map'
#                     )

# melanoma_subsetting(    input_file='WorkingDir/MELA-AU_trinuc_context_mutations.bed',
#                         tf_map='WorkingDir/TF_map_1000.bed',
#                         output_file='MELA-AU_SNP_TF_map'
#                     )