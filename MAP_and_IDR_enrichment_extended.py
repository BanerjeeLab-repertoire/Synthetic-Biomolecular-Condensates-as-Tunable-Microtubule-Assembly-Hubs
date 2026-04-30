# This script processes an XLSX file containing protein sequences and their disordered regions, calculates the relative frequencies of amino acids
    # in both the full protein sequences (MAPs) and the disordered regions (IDRs), and generates histograms and enrichment plots.
# It also saves the results in CSV files and plots the data for further analysis.

import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import os
import numpy as np

def process_sequences(file_path):
    # Read the data from the XLSX file
    df = pd.read_excel(file_path)

    # Extract the relevant columns for MAPs and IDRs
    map_sequences = df.iloc[:, 2]  # Full protein sequence in column 3 (index 2)
    idr_sequences = df.iloc[:, 3:32]  # Disordered regions in columns 4 to 32 (index 3 to 31)
    
    # Initialize lists to store amino acids for both MAPs and IDRs
    amino_acids_all_proteins = []
    amino_acids_all_idrs = []
    
    # Counters for occurrences of "X"
    count_x_in_maps = 0
    count_x_in_idrs = 0
    
    proteins_with_disordered_regions = 0  # For tracking the number of proteins with IDRs

    # Process the full protein sequence (MAPs)
    for seq in map_sequences:
        if isinstance(seq, str):
            count_x_in_maps += seq.count("X")  # Count occurrences of "X"
            amino_acids_all_proteins.extend([aa for aa in seq if aa != "X"])  # Exclude "X"

    # Process the disordered regions (IDRs)
    for _, row in idr_sequences.iterrows():
        has_idr = False
        for disordered_region in row:
            if pd.notna(disordered_region) and disordered_region != '':
                count_x_in_idrs += disordered_region.count("X")  # Count occurrences of "X"
                has_idr = True
                amino_acids_all_idrs.extend([aa for aa in disordered_region if aa != "X"])  # Exclude "X"
        
        if has_idr:
            proteins_with_disordered_regions += 1

    # Count occurrences of each amino acid
    count_all_proteins = Counter(amino_acids_all_proteins)
    count_all_idrs = Counter(amino_acids_all_idrs)

    # Calculate total amino acids in MAPs and IDRs
    total_aa_proteins = sum(count_all_proteins.values())
    total_aa_idrs = sum(count_all_idrs.values())

    # Print the total number of proteins and number of proteins with IDRs
    print(f"Proteins with disordered regions: {proteins_with_disordered_regions}")
    print(f"Total amino acids in MAPs: {total_aa_proteins}")
    print(f"Total amino acids in IDRs: {total_aa_idrs}")
    print(f"Number of 'X' in MAPs: {count_x_in_maps}")
    print(f"Number of 'X' in IDRs: {count_x_in_idrs}")

    # Normalize the counts to get relative frequencies for both MAPs and IDRs
    normalized_proteins = {aa: count / total_aa_proteins for aa, count in count_all_proteins.items()}
    normalized_idrs = {aa: count / total_aa_idrs for aa, count in count_all_idrs.items()}

    # Create tables for relative frequencies (MAPs and IDRs)
    df_proteins = pd.DataFrame(list(normalized_proteins.items()), columns=["Amino Acid", "Relative Frequency (MAPs)"])
    df_idrs = pd.DataFrame(list(normalized_idrs.items()), columns=["Amino Acid", "Relative Frequency (IDRs)"])

    # Save these tables to CSV files
    df_proteins.to_csv("Amino_Acid_Frequencies_in_MAPs.csv", index=False)
    df_idrs.to_csv("Amino_Acid_Frequencies_in_IDRs.csv", index=False)
    print("CSV files for relative frequencies saved.")

    # Plotting histograms for MAPs
    sorted_amino_acids_proteins = sorted(normalized_proteins.items())
    amino_acids = [aa for aa, _ in sorted_amino_acids_proteins]
    relative_frequencies_proteins = [freq for _, freq in sorted_amino_acids_proteins]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(amino_acids, relative_frequencies_proteins, color='blue', edgecolor='black')
    ax.set_title('Relative Frequency of Amino Acid Content in MAPs')
    ax.set_xlabel('Amino Acid')
    ax.set_ylabel('Relative Frequency')
    plt.tight_layout()
    output_file_name = "Amino_Acid_Content_MAPs_Relative_Frequency.png"
    plt.savefig(output_file_name)
    plt.close()
    print(f"Histogram for MAPs saved as: {output_file_name}")

    # Plotting histograms for IDRs
    sorted_amino_acids_idrs = sorted(normalized_idrs.items())
    amino_acids = [aa for aa, _ in sorted_amino_acids_idrs]
    relative_frequencies_idrs = [freq for _, freq in sorted_amino_acids_idrs]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(amino_acids, relative_frequencies_idrs, color='green', edgecolor='black')
    ax.set_title('Relative Frequency of Amino Acid Content in IDRs')
    ax.set_xlabel('Amino Acid')
    ax.set_ylabel('Relative Frequency')
    plt.tight_layout()
    output_file_name = "Amino_Acid_Content_IDRs_Relative_Frequency.png"
    plt.savefig(output_file_name)
    plt.close()
    print(f"Histogram for IDRs saved as: {output_file_name}")

    # Amino acids present in both MAPs and IDRs
    amino_acids = sorted(set(count_all_proteins.keys()).intersection(set(count_all_idrs.keys())))

    # Compute log10 enrichment
    enrichment = {}
    for aa in amino_acids:
        f_map = count_all_proteins[aa] / total_aa_proteins
        f_idr = count_all_idrs[aa] / total_aa_idrs
        enrichment[aa] = np.log10(f_idr / f_map)  # Log10 enrichment ratio

    # Sorting enrichment values for plotting
    sorted_enrichment = sorted(enrichment.items())
    aa_labels = [aa for aa, _ in sorted_enrichment]
    enrichment_values = [enrich for _, enrich in sorted_enrichment]

    # Create table for enrichment values
    df_enrichment = pd.DataFrame(list(enrichment.items()), columns=["Amino Acid", "Log10 Enrichment"])
    df_enrichment.to_csv("Amino_Acid_Enrichment_in_IDRs.csv", index=False)
    print("CSV file for enrichment values saved.")

    # Plotting enrichment plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(aa_labels, enrichment_values, color='purple', edgecolor='black')
    ax.axhline(0, color='black', linewidth=1)
    ax.set_title('Log10 Enrichment of Amino Acids in IDRs Relative to MAPs')
    ax.set_xlabel('Amino Acid')
    ax.set_ylabel('log10( F_IDR / F_MAP )')
    plt.tight_layout()
    output_file_name = "Amino_Acid_Enrichment_in_IDRs.png"
    plt.savefig(output_file_name)
    plt.close()
    print(f"Enrichment plot saved as: {output_file_name}")


    # Define entire proteome frequencies (convert percentage to fraction)
    proteome_freqs = {
        'A': 7.07, 'B': 0.00, 'C': 1.99, 'D': 4.75, 'E': 6.75, 'F': 3.68, 'G': 6.65,
        'H': 2.53, 'I': 4.38, 'K': 5.03, 'L': 9.90, 'M': 2.28, 'N': 3.50, 'P': 5.84,
        'Q': 4.82, 'R': 5.90, 'S': 8.15, 'T': 6.00, 'U': 0.00, 'V': 5.92, 'W': 1.58,
        'X': 0.12, 'Y': 3.17, 'Z': 0.00
    }
    proteome_freqs = {aa: freq / 100 for aa, freq in proteome_freqs.items()}  # Normalize

    # MAPs vs Proteome
    enrichment_map_vs_proteome = {}
    for aa in normalized_proteins:
        if aa in proteome_freqs and proteome_freqs[aa] > 0:
            enrichment_map_vs_proteome[aa] = np.log10(normalized_proteins[aa] / proteome_freqs[aa])

    df_enrichment_map = pd.DataFrame(list(enrichment_map_vs_proteome.items()),
                                     columns=["Amino Acid", "Log10 Enrichment (MAPs vs Proteome)"])
    df_enrichment_map.to_csv("Enrichment_MAPs_vs_Proteome.csv", index=False)

    fig, ax = plt.subplots(figsize=(8, 6))
    aa_labels = sorted(enrichment_map_vs_proteome.keys())
    enrichment_vals = [enrichment_map_vs_proteome[aa] for aa in aa_labels]
    ax.bar(aa_labels, enrichment_vals, color='blue', edgecolor='black')
    ax.axhline(0, color='black', linewidth=1)
    ax.set_title("Log10 Enrichment of MAPs Relative to Proteome")
    ax.set_xlabel("Amino Acid")
    ax.set_ylabel("log10( F_MAP / F_PROTEOME )")
    plt.tight_layout()
    plt.savefig("MAPs_vs_Proteome_Enrichment.png")
    plt.close()
    print("Enrichment plot for MAPs vs Proteome saved.")

    # IDRs vs Proteome
    enrichment_idr_vs_proteome = {}
    for aa in normalized_idrs:
        if aa in proteome_freqs and proteome_freqs[aa] > 0:
            enrichment_idr_vs_proteome[aa] = np.log10(normalized_idrs[aa] / proteome_freqs[aa])

    df_enrichment_idr = pd.DataFrame(list(enrichment_idr_vs_proteome.items()),
                                     columns=["Amino Acid", "Log10 Enrichment (IDRs vs Proteome)"])
    df_enrichment_idr.to_csv("Enrichment_IDRs_vs_Proteome.csv", index=False)

    fig, ax = plt.subplots(figsize=(8, 6))
    aa_labels = sorted(enrichment_idr_vs_proteome.keys())
    enrichment_vals = [enrichment_idr_vs_proteome[aa] for aa in aa_labels]
    ax.bar(aa_labels, enrichment_vals, color='purple', edgecolor='black')
    ax.axhline(0, color='black', linewidth=1)
    ax.set_title("Log10 Enrichment of IDRs Relative to Proteome")
    ax.set_xlabel("Amino Acid")
    ax.set_ylabel("log10( F_IDR / F_PROTEOME )")
    plt.tight_layout()
    plt.savefig("IDRs_vs_Proteome_Enrichment.png")
    plt.close()
    print("Enrichment plot for IDRs vs Proteome saved.")


    # Combined Enrichment Plot: MAPs vs Proteome and IDRs vs Proteome

    all_aa = sorted(set(enrichment_map_vs_proteome.keys()) | set(enrichment_idr_vs_proteome.keys()))

    map_enrichment_vals = [enrichment_map_vs_proteome.get(aa, 0) for aa in all_aa]
    idr_enrichment_vals = [enrichment_idr_vs_proteome.get(aa, 0) for aa in all_aa]

    x = np.arange(len(all_aa))
    width = 0.35

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(x - width / 2, map_enrichment_vals, width, label='MAPs vs Proteome', color='blue', edgecolor='black')
    ax.bar(x + width / 2, idr_enrichment_vals, width, label='IDRs vs Proteome', color='purple', edgecolor='black')

    # Horizontal baseline
    ax.axhline(0, color='black', linewidth=1)

    # Tick labels and spacing
    ax.set_xticks(x)
    ax.set_xticklabels(all_aa, rotation=0)
    ax.set_xlabel("Amino Acid")
    ax.set_ylabel("log10( Frequency / Proteome )")
    ax.set_title("Amino Acid Enrichment: MAPs & IDRs vs Proteome")
    ax.legend()

    # Add dashed vertical lines between amino acids
    for i in range(len(all_aa)):
        ax.axvline(x=i, color='gray', linestyle='--', linewidth=0.5, zorder=0)

    plt.tight_layout()
    plt.savefig("MAPs_IDRs_vs_Proteome_Combined_Enrichment.png")
    plt.close()

    print("Combined enrichment plot (MAPs & IDRs vs Proteome) saved.")

# file path
file_path = r"QuickGO-annotations-1742935775768-03 25 25 - CLEANED.xlsx"
process_sequences(file_path)
