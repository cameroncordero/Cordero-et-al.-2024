import csv

def process_tsv(input_filename, output_filename):
    # Open the input file in read mode
    with open(input_filename, 'r') as input_file:
        # Create a CSV reader object for the input file
        reader = csv.reader(input_file, delimiter='\t')

        # Open the output file in write mode
        with open(output_filename, 'w', newline='') as output_file:
            # Create a CSV writer object for the output file
            writer = csv.writer(output_file, delimiter='\t')

            # Iterate over the rows in the input file
            for row in reader:
                # Write the first three columns of the row to the output file
                writer.writerow(row[:3])


process_tsv( input_filename='WorkingDir/ENCFF742YZA.bed',
                output_filename='WorkingDir/ENCFF742YZA_truncated.bed')