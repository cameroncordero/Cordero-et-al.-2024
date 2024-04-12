from pathlib import Path

def remove_double_quotes(file_path, output_dir):
    # Read the lines from the file
    lines = file_path.read_text().splitlines()
    
    # Remove double quotes from each line
    new_lines = [line.replace('"', '') for line in lines]
    
    # Write the processed lines to a new file in the output directory
    new_file_path = output_dir / (file_path.stem + "_processed" + file_path.suffix)
    new_file_path.write_text('\n'.join(new_lines))

def main():
    # Automatically use the directory where the script is running
    input_dir = Path(input())

    # Construct the output directory path
    output_dir = input_dir / "processed_files"
    output_dir.mkdir(exist_ok=True)

    # Process all text files in the input directory
    for file_path in input_dir.glob('*.vcf'):  # Assuming text files, change the extension if required
        remove_double_quotes(file_path, output_dir)

if __name__ == "__main__":
    main()