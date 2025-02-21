import pandas as pd

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



# Define HER2 protein sequence (UniProt: P04626)
her2_sequence = "MRLTVLQLLLLWVSGGLRAEAWTNLGR..." 

# Map mutation positions to HER2 sequence
def map_mutations_to_sequence(mutations, sequence):
    for _, row in mutations.iterrows():
        pos = row['start']  # We are using 'start' position for mutation
        ref_allele = row['mutation_detail'][0]  # Assuming the first letter is the reference allele
        alt_allele = row['mutation_detail'][-1]  # Assuming the last letter is the mutated allele
        
        # Check the mutation against the sequence at the given position
        if sequence[pos-1] != ref_allele:  # Account for Python's 0-based indexing
            print(f"Warning: Sequence mismatch at position {pos}. Expected: {ref_allele}, Found: {sequence[pos-1]}")
        else:
            print(f"Position: {pos}, Wild-Type: {ref_allele}, Mutant: {alt_allele}, Sequence at position: {sequence[pos-1]}")

# Run mutation mapping
map_mutations_to_sequence(ecd_mutations, her2_sequence)

# Save final cleaned mutation data
ecd_mutations.to_csv("final_filtered_ecd_mutations.csv", index=False)

