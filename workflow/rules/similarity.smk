import pandas as pd

import warnings
warnings.filterwarnings("ignore")

rule extract_file:
    input:
        fasta = f"{sys.path[0]}/../../results/fasta_filter.fasta",


# Load the DataFrame
df = pd.read_csv('dataframe_args.csv')  # Assuming your DataFrame is stored in a CSV file

# Extract the values from the desired column
region_names = df['Folder'].tolist()

# Now you have the region names from the DataFrame column

# Define the function to filter the FASTA file
def filter_fasta(input_file, output_file, region_names):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        current_header = None
        sequence = ''
        for line in f_in:
            if line.startswith('>'):
                if current_header and current_header in region_names:
                    f_out.write(f'>{current_header}\n')
                    f_out.write(f'{sequence}\n')
                current_header = line.strip().split()[0][1:]
                sequence = ''
            else:
                sequence += line.strip()
        # Write the last sequence
        if current_header and current_header in region_names:
            f_out.write(f'>{current_header}\n')
            f_out.write(f'{sequence}\n')


# Call the function to filter the FASTA file
input_fasta_file = '/home/fernanda/abaumannii_prophages/results/fasta_filter.fasta'
output_fasta_file = 'filtered_args.fasta'
filter_fasta(input_fasta_file, output_fasta_file, region_names)

print("Filtered FASTA file created successfully.")

#%%
