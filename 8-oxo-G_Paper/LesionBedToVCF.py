from pathlib import Path
import subprocess

def bed_to_vcf(bed_file):
    sample_vcf_data = 'PASS	SOMATIC;QSS=56;TQSS=1;NT=ref;QSS_NT=56;TQSS_NT=1;SGT=GG->GT;DP=76;MQ=60.00;MQ0=0;ReadPosRankSum=-0.82;SNVSB=0.00;SomaticEVS=19.83	DP:FDP:SDP:SUBDP:AU:CU:GU:TU	26:1:0:0:0,0:0,0:25,26:0,0	49:0:0:0:0,0:0,0:25,26:24,24'
    bed_file = Path(bed_file)
    vcf_file = bed_file.with_suffix('.vcf')
    with open(bed_file, 'r') as bed, open(vcf_file, 'w') as vcf:
        for line in bed:
            line = line.strip().split('\t')
            chrom = line[0]
            start_pos = line[1]
            id = line[3]
            end_pos = line[2]
            strandedness = line[5]
            ref = line[6][1]
            if strandedness == '+':
                alt = 'T'
            else:
                ref = 'C'
                alt = 'A'
            vcf_line = '\t'.join([chrom,str(int(end_pos)-1), '.', ref, alt, '.', sample_vcf_data])
            vcf.write(vcf_line + '\n')


bed_to_vcf('/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/proofed_files/SRR_untreated_cellular_64-65-66.bed')
bed_to_vcf('/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/proofed_files/SRR_treated_naked_67-68.bed')
bed_to_vcf('/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/proofed_files/SRR_treated_cellular_69-70.bed')

subprocess.run(['sort', '-k1,1', '-k2,2n', '/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/proofed_files/SRR_untreated_cellular_64-65-66.vcf', '-o', '/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/proofed_files/SRR_untreated_cellular_64-65-66.vcf'])
subprocess.run(['sort', '-k1,1', '-k2,2n', '/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/proofed_files/SRR_treated_naked_67-68.vcf', '-o', '/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/proofed_files/SRR_treated_naked_67-68.vcf'])
subprocess.run(['sort', '-k1,1', '-k2,2n', '/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/proofed_files/SRR_treated_cellular_69-70.vcf', '-o', '/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/proofed_files/SRR_treated_cellular_69-70.vcf'])

def one_line_only(file):
    file = Path(file)
    previous_line = ''
    with open(file, 'r') as f, open(file.with_suffix('.filtered.vcf'), 'w') as o:
        for line in f:
            if line == previous_line:
                continue
            else:
                previous_line = line
                o.write(line)

one_line_only('/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/proofed_files/SRR_untreated_cellular_64-65-66.vcf')
one_line_only('/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/proofed_files/SRR_treated_naked_67-68.vcf')
one_line_only('/media/cam/Working/8-oxodG/8-oxodG_Final_Analysis/lesion_files/proofed_files/SRR_treated_cellular_69-70.vcf')