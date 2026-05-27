import csv
import argparse
import textwrap

# to call this script: python csv_to_fasta.py <input_filepath> <output_filepath>
#for the file path include the .csv and .fasta, the script will creat a fasta file for you
# Set up argument parsing
parser = argparse.ArgumentParser(description="Convert a CSV file to a FASTA file with custom headers.")
parser.add_argument('input_file', type=str, help="Path to the input CSV file.")
parser.add_argument('output_file', type=str, help="Path to the output FASTA file.")
args = parser.parse_args()

# Open the CSV file for reading and the FASTA file for writing
with open(args.input_file, mode='r') as csv_file, open(args.output_file, mode='w') as fasta_file:
    csv_reader = csv.DictReader(csv_file)
    
    for row in csv_reader:
        # Extract the necessary fields from the CSV row
        ID = row['ID']
        contig = row['Contig']
        e_value = row['E-value']
        align_len = row['Alignment Length']
        seq_identity = row['Sequence Identity']
        query_start = row['Query Start']
        subject_start = row['Subject Start']
        subject_end = row['Subject End']
        sequence = row['Sequence']
        
        # Create the FASTA header
        header = f">{ID}|contig={contig}|E-value={e_value}|AlignLen={align_len}|SeqIdentity={seq_identity}|SubjectStart={subject_start}|SubjectEnd={subject_end}"

        # Format sequence into 60-character lines for FASTA format
        wrapped_sequence = "\n".join(textwrap.wrap(sequence, 60))

        # Write to FASTA file
        fasta_file.write(header + '\n' + wrapped_sequence + '\n')


print(f"FASTA file '{args.output_file}' has been created.")