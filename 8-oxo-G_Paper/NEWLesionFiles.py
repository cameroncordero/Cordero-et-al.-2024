from pathlib import Path
import subprocess

def bed_to_vcf(bed_file):
    with open(bed_file, 'r') as bed, open(bed_file.with_suffix('.vcf'), 'w') as vcf:
        sample_vcf_data = 'PASS	SOMATIC;QSS=56;TQSS=1;NT=ref;QSS_NT=56;TQSS_NT=1;SGT=GG->GT;DP=76;MQ=60.00;MQ0=0;ReadPosRankSum=-0.82;SNVSB=0.00;SomaticEVS=19.83	DP:FDP:SDP:SUBDP:AU:CU:GU:TU	26:1:0:0:0,0:0,0:25,26:0,0	49:0:0:0:0,0:0,0:25,26:24,24'
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
            vcf_line = '\t'.join([chrom, str(int(end_pos)-1), '.', ref, alt, '.', 'PASS'])
            vcf.write(vcf_line + '\n')