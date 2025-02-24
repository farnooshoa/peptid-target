import pandas as pd
from Bio import SeqIO

# Load mutation data from CSV
mutation_data = pd.read_csv("her2_mutations.csv")

# Inspect the first few rows to understand the data structure
print("Initial Mutation Data:")
print(mutation_data.head())

# Filter mutations based on the 'start' and 'end' positions within the extracellular domain (ECD)
ecd_mutations = mutation_data[(mutation_data['start'] >= 1) & (mutation_data['end'] <= 645)]

# Print filtered mutations to verify
print("Filtered ECD Mutations:")
print(ecd_mutations)

# Extract mutation details from the 'attributes' column
ecd_mutations['mutation_detail'] = ecd_mutations['attributes'].str.extract(r'([a-zA-Z]+\d+[a-zA-Z]+)')

# Print mutation details to verify extraction
print("Mutation Details:")
print(ecd_mutations[['seqid', 'start', 'end', 'mutation_detail']])

# Clean the data by dropping rows with missing mutation details
ecd_mutations = ecd_mutations.dropna(subset=['mutation_detail'])

# Print out the number of missing entries and remove them
print("Missing Mutation Details (Before Cleaning):", ecd_mutations['mutation_detail'].isna().sum())

# Save filtered mutations to a new CSV file
ecd_mutations.to_csv("filtered_ecd_mutations.csv", index=False)

# Define HER2 protein sequence
fasta_file = "P04626.fasta"

# Read the sequence from the FASTA file
def read_fasta(file_path):
    try:
        for record in SeqIO.parse(file_path, "fasta"):
            return str(record.seq)  # Convert to a string
    except FileNotFoundError:
        print(f"Error: The file {file_path} was not found.")
        return None
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return None

# Load the HER2 sequence
her2_sequence = read_fasta(fasta_file)

# Check if the sequence was loaded successfully
if her2_sequence is None:
    raise ValueError("Failed to load HER2 sequence. Check the FASTA file path and format.")

# Print the sequence length to verify
print(f"Length of HER2 Sequence: {len(her2_sequence)}")

# Map mutation positions to HER2 sequence
def map_mutations_to_sequence(mutations, sequence):
    for _, row in mutations.iterrows():
        pos = row['start']  # We are using 'start' position for mutation
        ref_allele = row['mutation_detail'][0]  # Assuming the first letter is the reference allele
        alt_allele = row['mutation_detail'][-1]  # Assuming the last letter is the mutated allele
        
        # Check if the position is within the sequence length
        if pos > len(sequence):
            print(f"Warning: Position {pos} is outside the sequence length ({len(sequence)}).")
            continue
        
        # Check the mutation against the sequence at the given position
        if sequence[pos-1] != ref_allele:  # Account for Python's 0-based indexing
            print(f"Warning: Sequence mismatch at position {pos}. Expected: {ref_allele}, Found: {sequence[pos-1]}")
        else:
            print(f"Position: {pos}, Wild-Type: {ref_allele}, Mutant: {alt_allele}, Sequence at position: {sequence[pos-1]}")

# Run mutation mapping
map_mutations_to_sequence(ecd_mutations, her2_sequence)

# Save final cleaned mutation data
ecd_mutations.to_csv("final_filtered_ecd_mutations.csv", index=False)

# Define functional regions of HER2
functional_regions = [
    (1, 195),    # Subdomain I (Extracellular)
    (196, 320),  # Subdomain II (Dimerization Domain)
    (321, 488),  # Subdomain III (Extracellular)
    (489, 645),  # Subdomain IV (Extracellular)
    (712, 987)   # Kinase Domain (Intracellular)
]

# Filter mutations in functional regions
key_mutations = ecd_mutations[
    ecd_mutations.apply(lambda row: any(start <= row['start'] <= end for start, end in functional_regions), axis=1)
]

print("Key Mutations in Functional Regions:")
print(key_mutations)

# Extract peptide sequences around mutation sites
def extract_peptides(sequence, mutations, window_size=10):
    peptides = []
    for _, row in mutations.iterrows():
        pos = row['start'] - 1  # Convert to 0-based indexing
        start = max(0, pos - window_size)
        end = min(len(sequence), pos + window_size + 1)
        peptide = sequence[start:end]
        peptides.append((row['mutation_detail'], peptide))
    return peptides

# Generate peptide candidates
peptide_candidates = extract_peptides(her2_sequence, key_mutations, window_size=10)

print("Peptide Candidates:")
for mutation, peptide in peptide_candidates:
    print(f"Mutation: {mutation}, Peptide: {peptide}")
    # Save peptide candidates to a CSV file
peptide_df = pd.DataFrame(peptide_candidates, columns=["Mutation", "Peptide"])
peptide_df.to_csv("peptide_candidates.csv", index=False)