from pathlib import Path
import tempfile
import subprocess
import os
import itertools

class VCFtoCustom():

    def __init__(self, vcfpath, mutation_type=None, context = None, expansion_length=1, assembly_path='/home/cam/Documents/genome_files/hg19/hg19.fa', split=True) -> None:
        self.vcfpath = Path(vcfpath)
        self.bedpath = self.vcfpath.with_suffix('.bed')
        if mutation_type:
            self.bedpath = self.bedpath.with_name(self.bedpath.stem + f'_{mutation_type}.bed')
        if context:
            self.bedpath = self.bedpath.with_name(self.bedpath.stem + f'_{context}.bed')
        if mutation_type and context:
            self.bedpath = self.bedpath.with_name(self.bedpath.stem + f'_{mutation_type}_{context}.bed')
        self.context_length = expansion_length
        self.assembly_path = assembly_path
        self.split = split

        # If a mutation type is specified, expand it to all possible specific mutations
        if mutation_type:
            self.mutation_types = self.mutation_combinations(mutation_type)
        else:
            self.mutation_types = self.mutation_combinations('N>N')

        # If a context is specified, expand it to all possible contexts
        if context:
            self.contexts = self.contexts_in_IUPAC(context)
        else:
            self.contexts = self.contexts_in_IUPAC('NNN')


    def convert(self) -> None:
        temp_files = []
        try:
            with open(self.vcfpath, 'r') as vcf:
                with tempfile.NamedTemporaryFile('w+t', delete=False) as temp_bed:
                    temp_files.append(temp_bed.name)
                    for line in vcf:
                        if line.startswith('#'):
                            continue
                        else:
                            line = line.strip().split('\t')
                            chrom = line[0]
                            if 'chr' not in chrom:
                                # chrom = 'chr' + chrom
                                pass
                            start = str(int(line[1]) - 1 - self.context_length)
                            end = str(int(line[1]) + len(line[3]) - 1 + self.context_length)
                            name = '.'
                            value = '0'

                            ref_base = line[3]  # reference base
                            alt_base = line[4]  # alternative base

                            sample_name = line[-1]  # sample name

                            strand = '+'  # Default strand

                            # Determine the type of mutation
                            if len(ref_base) == 1 and len(alt_base) == 1:
                                mutation_type = 'SNP'
                            elif len(ref_base) < len(alt_base):
                                mutation_type = 'INS'
                            elif len(ref_base) > len(alt_base):
                                mutation_type = 'DEL'
                            else:
                                mutation_type = 'MNP'  # Or other suitable type for your case

                            mutation = ref_base + '>' + alt_base  # Default mutation representation

                            # Apply special consideration only for SNPs that match the desired mutation type.
                            if mutation_type == 'SNP':
                                mutation = ref_base + '>' + alt_base

                            if int(start) < 0:
                                continue
                            # Write the current entry into the temporary BED file.
                            temp_bed.write(f'{chrom}\t{start}\t{end}\t{name}\t{value}\t{strand}\t{mutation}\t{mutation_type}\t{sample_name}\n')

                    temp_bed.flush()
                    temp_bed.seek(0)
                    # command = [
                    #     'sort', '-k1,1', 
                    #     '-k2,2n', '-k3,3n',
                    #     temp_bed.name, '-o', temp_bed.name
                    #     ]
                    # subprocess.run(command)
                    temp_fasta = self.get_fasta(temp_bed.name)
                    temp_files.append(temp_fasta)
                    contexts = self.contexts[0]
                    reverse_complement_contexts = self.contexts[1]
                    with open(self.bedpath, 'w') as bed, open(temp_fasta, 'r') as temp_fasta:
                        for line in temp_bed:
                            fasta_coordinates = temp_fasta.readline()
                            context = temp_fasta.readline().strip()
                            line = line.strip().split('\t')
                            strand = line[5]
                            mutation_type = line[7]
                            mutation = line[6]
                            if mutation_type == 'SNP':    
                                # Check if the mutation is in the list of specific mutations we're interested in
                                if mutation in self.mutation_types:
                                    # Check if the context is in the list of contexts we're interested in
                                    if context[self.context_length-1:self.context_length+2].upper() in contexts:
                                        strand = '+'
                                    else:
                                        continue  # Skip this mutation if the context is not in the list
                                elif (self.reverse_complement(ref_base) + '>' + self.reverse_complement(alt_base)) in self.mutation_types:
                                    # Check if the reverse complement of the context is in the list of contexts we're interested in
                                    if self.reverse_complement(context) in reverse_complement_contexts:
                                        mutation = self.reverse_complement(ref_base) + '>' + self.reverse_complement(alt_base)
                                        strand = '-'
                                    else:
                                        continue  # Skip this mutation if the reverse complement of the context is not in the list
                                else:
                                    continue  # Skip this mutation if it's not in the list of mutation types
                            if line[1] not in fasta_coordinates and line[2] not in fasta_coordinates:
                                print('Error: BED and FASTA coordinates do not match')
                                print(fasta_coordinates, context, line)
                                temp_bed.readline()
                                if fasta_coordinates == '':
                                    print('Error: FASTA coordinates are empty')
                                    continue
                                if context == '':
                                    print('Error: Context is empty')
                                    continue
                                # exit(1)
                            line[1] = str(int(line[1]) + self.context_length)
                            line[2] = str(int(line[2]) - self.context_length)
                            new_line = '\t'.join(line[:6]) + '\t' + context.strip()[:self.context_length].lower()+ context.strip()[self.context_length:-self.context_length].upper() + context.strip()[-self.context_length:].lower() + '\t' + '\t'.join(line[6:]) + '\n'
                            bed.write(new_line)
                    self.sort_bed(self.bedpath)
                    if self.split:
                        self.split_by_type()
        finally:
            # pass
            for file in temp_files:
                try:
                    os.remove(file)  # Delete each temporary file
                except FileNotFoundError:
                    pass  # File already deleted, do nothing.

    def mutation_combinations(self, mut_type: str):
        """Takes a mutation type in the format N>N and returns a list of all possible nucleotide mutations.

        Args:
            mut_type (str): the mutation type in the format N>N

        Returns:
            possible_mutations (list): a list of all possible nucleotide mutations
        """
        iupac_trans = {
            "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT",
            "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
            "N": "ACGT", 'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G'
        }
        
        # Split the mutation type at the '>' character
        from_nuc, to_nuc = mut_type.split('>')
        
        # Get the possible nucleotides for each side of the '>'
        from_nucs = iupac_trans[from_nuc]
        to_nucs = iupac_trans[to_nuc]
        
        # Generate all possible combinations
        possible_mutations = [f"{f}>{t}" for f in from_nucs for t in to_nucs]
        
        return possible_mutations

    def reverse_complement(self, seq: str):
        """returns the reverse complement of nucleotide sequences in standard or IUPAC notation

        Args:
            seq (str): sequence of DNA in standard or IUPAC form that

        Returns:
            reverse_complement (str): the reverse complement of the input sequence
        """
        # make a lookup table
        complement_table = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "R": "Y",
            "Y": "R",
            "S": "S",
            "W": "W",
            "K": "M",
            "M": "K",
            "B": "V",
            "D": "H",
            "H": "D",
            "V": "B",
            "N": "N"
        }

        seq_rev = seq[::-1]
        complement_seq = "".join(complement_table.get(base, base) for base in seq_rev)
        return complement_seq
    
    def contexts_in_IUPAC(self, context: str):
        # takes an IUPAC context and returns a list of all possible contexts
        iupac_trans = {
            "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT",
            "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
            "N": "ACGT", 'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G'
        }
        context_list = []
        rev_comp_list = []
        for bases in itertools.product(*(iupac_trans[base] for base in context)):
            new_context = ''.join(bases)
            context_list.append(new_context)
            rev_comp_list.append(self.reverse_complement(new_context))
        return context_list, rev_comp_list

    def get_fasta(self, bed_file) -> str:
        with tempfile.NamedTemporaryFile('w', delete=False) as temp_fasta:
            command = [
                'bedtools', 'getfasta', '-s',
                '-fi', self.assembly_path,  # Input FASTA file
                '-bed', bed_file,  # Input BED file
                '-fo', temp_fasta.name  # Output file to write the extracted sequences
                ]
            subprocess.run(command)
            return temp_fasta.name


    def sort_bed(self, bed_file) -> str:
        command = [
            'sort', '-k1,1', 
            '-k2,2n', '-k3,3n',
            bed_file, '-o', bed_file
            ]
        subprocess.run(command)
        return bed_file

    def split_by_type(self):
        with open(self.bedpath) as bed:
            lines = bed.readlines()

        # Initialize mutation lists
        mutations = {'SNP': [], 'DNP': [], 'ONP': [], 'INS': [], 'DEL': []}

        prev_line = None
        combined_mutation = None  # To handle ongoing combination of SNPs

        for line in lines:
            line = line.strip().split('\t')
            chromosome, start, end, *rest, context, mutation, mutation_type, sample = line

            if mutation_type in ['INS', 'DEL']:
                mutations[mutation_type].append(line)
                continue

            if prev_line:
                prev_chromosome, prev_start, prev_end, *prev_rest, prev_context, prev_mutation, prev_mutation_type, prev_sample = prev_line

                same_chromosome = chromosome == prev_chromosome
                adjacent_position = int(prev_end) == int(start)
                same_sample = sample == prev_sample

                if same_chromosome and adjacent_position and same_sample:
                    # This is a continuation of a tandem mutation
                    combined_start = prev_start
                    combined_end = end

                    # Extract uppercase characters from both contexts
                    prev_upper = ''.join(filter(str.isupper, prev_context))
                    next_upper = ''.join(filter(str.isupper, context))

                    # Combine the contexts
                    combined_context = prev_context.split(prev_upper)[0] + prev_upper + next_upper + context.split(next_upper)[-1]

                    combined_mutation = prev_mutation.split('>')[0] + mutation[0] + '>' + prev_mutation.split('>')[-1] + mutation[-1]

                    # Determine new mutation type based on combined length
                    mutation_length = len(combined_mutation.split('>')[0])
                    if mutation_length == 2:
                        new_mutation_type = 'DNP'
                    else:
                        new_mutation_type = 'ONP'

                    # Update prev_line to represent combined mutation
                    prev_line = [chromosome, combined_start, combined_end, *rest, combined_context, combined_mutation, new_mutation_type, sample]
                else:
                    # Not combinable, handle previous combined mutation or standalone SNP
                    if combined_mutation:
                        mutations[new_mutation_type].append(prev_line)
                        combined_mutation = None  # Reset combined mutation tracker
                    else:
                        # Previous line was a standalone SNP, not part of a combination
                        mutations[prev_mutation_type].append(prev_line)

                    prev_line = line  # Move to the next line as the current line
            else:
                prev_line = line

        # Handle the last line after loop completion
        if combined_mutation:
            mutations[new_mutation_type].append(prev_line)
        elif prev_line:
            *_, mutation_type, _ = prev_line
            mutations[mutation_type].append(prev_line)

        # Write collected mutations to their respective files
        for mutation_type, lines in mutations.items():
            with open(self.bedpath.with_name(self.bedpath.stem + f'_{mutation_type}.bed'), 'w') as file:
                for line in lines:
                    file.write('\t'.join(line) + '\n')

def reverse_complement(seq: str):
    """returns the reverse complement of nucleotide sequences in standard or IUPAC notation

    Args:
        seq (str): sequence of DNA in standard or IUPAC form that

    Returns:
        reverse_complement (str): the reverse complement of the input sequence
    """
    # make a lookup table
    complement_table = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N"
    }

    seq_rev = seq[::-1]
    complement_seq = "".join(complement_table.get(base, base) for base in seq_rev)
    return complement_seq

def bed_to_vcf(bed_path):
    # artificially create a vcf file from the bed file
    bed_path = Path(bed_path)
    vcf_output_path = bed_path.with_suffix('.vcf')
    vcf_header = f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n"
    with open(vcf_output_path, 'w') as vcf:
        vcf.write(vcf_header)
        with open(bed_path) as bed:
            for line in bed:
                chrom, start, end, name, value, strand, context, mutation_type = line.strip().split('\t')
                ref = mutation_type[0]
                alt = mutation_type[-1]
                if strand == '-':
                    ref = reverse_complement(ref)
                    alt = reverse_complement(alt)
                pos = int(start) + 1
                vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\t.\t.\t.\n"
                vcf.write(vcf_line)

# bed_to_vcf('/media/cam/Working/UV_test_data/UV_nucleomutics/UV.mut')

conversion = VCFtoCustom('/media/cam/Working/APOBEC_Deletions/Strelka/Strelka_results/unique_mutations_pass_filtered.vcf', assembly_path='/home/cam/Documents/genome_files/SacCer3/SacCer3.fa', split=True, expansion_length=30)
conversion.convert()