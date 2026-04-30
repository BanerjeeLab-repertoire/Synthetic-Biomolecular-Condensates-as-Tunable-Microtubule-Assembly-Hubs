# This script fetches UniProt data for proteins listed in an Excel file, extracts disordered regions, and saves the results back to a new Excel file.

import xml.etree.ElementTree as ET
import requests
import pandas as pd
import time

# Record start time
start_time = time.time()

# Load data from the Excel file
# file path
file_path = r"QuickGO-annotations-1742935775768-03 25 25.xlsx"
df = pd.read_excel(file_path)

# Define the namespace (from the XML file)
namespaces = {'uniprot': 'http://uniprot.org/uniprot'}

# Counters for statistics
total_proteins = len(df)
found_in_uniprot = 0
proteins_with_disordered_regions = 0

# Function to fetch UniProt data with retries
def fetch_uniprot_data(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
    retries = 2  # Max number of retries
    for i in range(retries):
        try:
            response = requests.get(url)
            if response.status_code == 200:
                return response.text
            else:
                print(f"Error fetching data for UniProt ID {uniprot_id}. HTTP Status code: {response.status_code}")
                return None
        except requests.exceptions.RequestException as e:
            print(f"Request error for UniProt ID {uniprot_id}: {e}")
            if i < retries - 1:
                print(f"Retrying... {i + 1}/{retries}")
                time.sleep(5)  # Wait before retrying
            else:
                print(f"Max retries reached for {uniprot_id}. Skipping.")
                return None

# Iterate over each row in the DataFrame
for index, row in df.iterrows():
    # Only proceed if the 'GENE PRODUCT DB' column value is 'UniProtKB'
    if row['GENE PRODUCT DB'] == 'UniProtKB':
        # Get the UniProt ID from the 'GENE PRODUCT ID' column (assuming it contains the UniProt ID directly)
        uniprot_id = row['GENE PRODUCT ID']
        
        # Fetch the XML data from UniProt with retries
        xml_data = fetch_uniprot_data(uniprot_id)
        if xml_data is None:
            continue  # Skip to next row if fetching data failed

        found_in_uniprot += 1  # Increment for proteins found in UniProt

        # Parse the XML data
        try:
            root = ET.fromstring(xml_data)
        except ET.ParseError as e:
            print(f"Error parsing XML for UniProt ID {uniprot_id}: {e}")
            continue

        # Initialize variables
        full_protein_sequence = ""
        disordered_sequences = []

        # Iterate through the entries in the XML
        for entry in root.findall('.//uniprot:entry', namespaces=namespaces):
            # Extract the full protein sequence
            sequence_tag = entry.find('uniprot:sequence', namespaces=namespaces)
            if sequence_tag is not None:
                full_protein_sequence = sequence_tag.text

            # Look for <feature> tags with type="region of interest" and description="Disordered"
            for feature in entry.findall('.//uniprot:feature', namespaces=namespaces):
                if feature.get('type') == 'region of interest' and feature.get('description') == 'Disordered':
                    # Extract the location (start and end positions)
                    location = feature.find('uniprot:location', namespaces=namespaces)
                    begin = location.find('uniprot:begin', namespaces=namespaces)
                    end = location.find('uniprot:end', namespaces=namespaces)
                    start = int(begin.get('position'))  # Extract the start position
                    end_position = int(end.get('position'))  # Extract the end position

                    # Ensure the end position does not exceed sequence length
                    if end_position > len(full_protein_sequence):
                        end_position = len(full_protein_sequence)

                    # Extract the disordered region from the sequence
                    disordered_region = full_protein_sequence[start-1:end_position]  # Adjusting for 0-based indexing in Python
                    disordered_sequences.append(disordered_region)

            # If disordered regions are found, increment the count
            if disordered_sequences:
                proteins_with_disordered_regions += 1  # Increment for proteins with disordered regions

        # Add the full protein sequence to the DataFrame
        df.at[index, 'Full Protein Sequence'] = full_protein_sequence

        # Add the disordered regions to the DataFrame
        if disordered_sequences:
            for i, disordered_sequence in enumerate(disordered_sequences):
                df.at[index, f'Disordered Region {i + 1}'] = disordered_sequence

# Print statistics
print(f"Total proteins in the xlsx file: {total_proteins}")
print(f"Proteins found in UniProt: {found_in_uniprot}")
print(f"Proteins with disordered regions: {proteins_with_disordered_regions}")

# Save the updated DataFrame to a new Excel file
output_file_path = r"C:\Users\shaml\Box\Shamli M Banerjee Lab Research Folder\IDR sequence analysis\Data files\Sukanya\QuickGO-annotations-1742935775768-03 25 25 - NEW.xlsx"

df.to_excel(output_file_path, index=False)

print(f"Updated Excel file saved at: {output_file_path}")

# Record end time
end_time = time.time()

# Calculate and print elapsed time
elapsed_time = end_time - start_time
print(f"Time taken to complete the request: {elapsed_time:.2f} seconds")
