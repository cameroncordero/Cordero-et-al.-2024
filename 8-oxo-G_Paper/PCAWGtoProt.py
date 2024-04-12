from pathlib import Path
import subprocess

def convert_files_to_vcf(directory, signature, cutoff=0.8):
    """Convert all decomposed files in a directory to a single VCF file"""
    directory = Path(directory)
    decomposed_files = directory.glob('*.txt')  # Adjust the file pattern as needed

    with open(directory / 'combined.vcf', 'w') as output:
        output.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'SAMPLE']) + '\n')

        for decomposed_file in decomposed_files:
            with open(decomposed_file, 'r') as file:
                header = file.readline().strip().split('\t')
                chrom_index = header.index('Chr')
                pos_index = header.index('Pos')
                mutation_type_index = header.index('MutationType')
                signature_index = header.index(signature)

                next(file)  # Skip header
                for line in file:
                    tsv = line.strip().split('\t')
                    if float(tsv[signature_index]) > cutoff:
                        chrom = tsv[chrom_index]
                        pos = tsv[pos_index]
                        mutation_type = tsv[mutation_type_index]
                        ref = mutation_type.split('>')[0][-1]
                        alt = mutation_type.split('>')[1][0]
                        output.write(f"chr{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\t{decomposed_file.stem}\n")

    # After the convert_files_to_vcf function
    subprocess.run(['sort', '-k1,1V', '-k2,2n', str(directory / 'combined.vcf')], stdout=open(directory / 'combined_sorted.vcf', 'w'))


convert_files_to_vcf('/media/cam/Working/8-oxodG/pcawg/cosmic/Assignment_Solution/Activities/Decomposed_Mutation_Probabilities', 'SBS18', 0.5)