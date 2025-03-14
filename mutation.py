import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re
import argparse
from multiprocessing import Pool
import matplotlib.pyplot as plt
import os

# --- Utility Functions ---

def load_mutation_data(file_path):
    """Load and validate mutation data from CSV."""
    try:
        df = pd.read_csv(file_path)
        required_cols = ['start', 'end', 'attributes']  
        if not all(col in df.columns for col in required_cols):
            raise ValueError(f"Missing required columns in {file_path}: {required_cols}")
        return df
    except FileNotFoundError:
        print(f"Error: {file_path} not found.")
        return None
    except Exception as e:
        print(f"Error loading CSV: {e}")
        return None

def read_fasta(file_path):
    """Read protein sequence from FASTA file."""
    try:
        for record in SeqIO.parse(file_path, "fasta"):
            return str(record.seq)
    except FileNotFoundError:
        print(f"Error: {file_path} not found.")
        return None
    except Exception as e:
        print(f"Error reading FASTA: {e}")
        return None

def parse_mutation(attribute):
    """Parse mutation details from attributes (e.g., A123T, A123del)."""
    match = re.search(r'([A-Za-z]+)(\d+)([A-Za-z]+|del|ins.*)', attribute)
    if match:
        ref, pos, alt = match.groups()
        return f"{ref}{pos}{alt}"
    return None

def score_peptide(peptide):
    """Score peptide for physicochemical properties."""
    try:
        analysis = ProteinAnalysis(peptide)
        return {
            "hydrophobicity": analysis.gravy(),
            "instability": analysis.instability_index(),
            "molecular_weight": analysis.molecular_weight()
        }
    except Exception as e:
        print(f"Error scoring peptide {peptide}: {e}")
        return None

# --- Main Processing Functions ---

def filter_ecd_mutations(mutation_data):
    """Filter mutations in HER2 extracellular domain (ECD, 1-645)."""
    ecd_mutations = mutation_data[(mutation_data['start'] >= 1) & (mutation_data['end'] <= 645)].copy()
    ecd_mutations['mutation_detail'] = ecd_mutations['attributes'].apply(parse_mutation)
    ecd_mutations = ecd_mutations.dropna(subset=['mutation_detail'])
    print(f"Filtered {len(ecd_mutations)} ECD mutations.")
    return ecd_mutations

def map_mutations_to_sequence(mutations, sequence, verbose=True):
    """Map mutations to the protein sequence and validate."""
    results = []
    for _, row in mutations.iterrows():
        pos = row['start'] - 1  # 0-based indexing
        if not (0 <= pos < len(sequence)):
            print(f"Error: Position {pos+1} out of bounds (sequence length: {len(sequence)}).")
            continue
        ref, alt = row['mutation_detail'][0], row['mutation_detail'][-1]
        seq_aa = sequence[pos]
        if seq_aa != ref and verbose:
            print(f"Mismatch at {pos+1}: Expected {ref}, Found {seq_aa}")
        results.append({'pos': pos+1, 'ref': ref, 'alt': alt, 'seq_aa': seq_aa})
    return pd.DataFrame(results)

def assign_functional_region(pos, regions):
    """Assign mutation to a functional region."""
    region_names = {
        (1, 195): "Subdomain I",
        (196, 320): "Subdomain II",
        (321, 488): "Subdomain III",
        (489, 645): "Subdomain IV",
        (712, 987): "Kinase Domain"
    }
    for start, end in regions:
        if start <= pos <= end:
            return region_names.get((start, end), "Unknown")
    return "Outside Defined Regions"

def is_in_binding_site(pos, binding_range=(529, 627)):
    """Check if mutation is in trastuzumab binding site."""
    return binding_range[0] <= pos <= binding_range[1]

def extract_peptides_with_mutation(sequence, mutations, window_size=10):
    """Extract wild-type and mutant peptides around mutation sites."""
    peptides = []
    for _, row in mutations.iterrows():
        pos = row['start'] - 1
        ref, alt = row['mutation_detail'][0], row['mutation_detail'][-1]
        start = max(0, pos - window_size)
        end = min(len(sequence), pos + window_size + 1)
        wt_peptide = sequence[start:end]
        mut_peptide = wt_peptide[:pos-start] + alt + wt_peptide[pos-start+1:]
        scores = score_peptide(wt_peptide)
        peptides.append((row['mutation_detail'], wt_peptide, mut_peptide, scores))
    return peptides

def save_peptides_to_fasta(peptides, filename="peptides.fasta"):
    """Save peptides to FASTA format for downstream analysis."""
    with open(filename, "w") as f:
        for i, (mutation, wt_peptide, mut_peptide, _) in enumerate(peptides):
            f.write(f">{mutation}_wildtype_{i}\n{wt_peptide}\n>{mutation}_mutant_{i}\n{mut_peptide}\n")

def plot_mutation_distribution(mutations, output_file="mutation_distribution.png"):
    """Visualize mutation distribution."""
    plt.hist(mutations['start'], bins=20, edgecolor='black')
    plt.xlabel("Position in HER2")
    plt.ylabel("Number of Mutations")
    plt.title("Mutation Distribution in HER2")
    plt.savefig(output_file)
    plt.close()

# --- Main Workflow ---

  # Define the helper function
def process_row(row, sequence, window_size):
    return extract_peptides_with_mutation(sequence, pd.DataFrame([row]), window_size)[0]


def main(csv_file, fasta_file, window_size=10):
    """Main function to process HER2 mutation data."""
    # Load data
    mutation_data = load_mutation_data(csv_file)
    if mutation_data is None:
        raise SystemExit("Stopping due to CSV loading error.")
    
    her2_sequence = read_fasta(fasta_file)
    if her2_sequence is None:
        raise SystemExit("Stopping due to FASTA loading error.")
    print(f"Loaded HER2 sequence, length: {len(her2_sequence)}")

    # Filter ECD mutations
    ecd_mutations = filter_ecd_mutations(mutation_data)
    ecd_mutations.to_csv("filtered_ecd_mutations.csv", index=False)

    # Map mutations to sequence
    mapped_results = map_mutations_to_sequence(ecd_mutations, her2_sequence)
    print("Mapped Mutations:\n", mapped_results.head())

    # Define functional regions
    functional_regions = [(1, 195), (196, 320), (321, 488), (489, 645), (712, 987)]
    ecd_mutations['region'] = ecd_mutations['start'].apply(lambda x: assign_functional_region(x, functional_regions))
    ecd_mutations['in_trastuzumab_site'] = ecd_mutations['start'].apply(is_in_binding_site)

    # Filter key mutations in functional regions
    key_mutations = ecd_mutations[
        ecd_mutations.apply(lambda row: any(start <= row['start'] <= end for start, end in functional_regions), axis=1)
    ]
    key_mutations.to_csv("final_filtered_ecd_mutations.csv", index=False)
    print("Key Mutations in Functional Regions:\n", key_mutations[['start', 'mutation_detail', 'region', 'in_trastuzumab_site']])

    window_size = 10
    # Use multiprocessing with the helper function
    with Pool() as pool:
        peptide_candidates = pool.starmap(
            process_row,
            [(row, her2_sequence, window_size) for _, row in key_mutations.iterrows()]
        )

    # Save peptide data
    peptide_df = pd.DataFrame(peptide_candidates, columns=["Mutation", "Wild-Type Peptide", "Mutant Peptide", "Scores"])
    peptide_df.to_csv("peptide_candidates.csv", index=False)
    save_peptides_to_fasta(peptide_candidates)
    print("Peptide Candidates (Top 5):")
    for mutation, wt, mut, scores in peptide_candidates[:5]:
        print(f"Mutation: {mutation}, WT: {wt}, Mutant: {mut}, Scores: {scores}")

    # Visualize mutation distribution
    plot_mutation_distribution(key_mutations)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze HER2 mutations for peptide design.")
    parser.add_argument("--csv", default="her2_mutations.csv", help="Mutation CSV file")
    parser.add_argument("--fasta", default="P04626.fasta", help="HER2 FASTA file")
    parser.add_argument("--window", type=int, default=10, help="Peptide window size")
    args = parser.parse_args()
    main(args.csv, args.fasta, args.window)