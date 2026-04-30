import pandas as pd
import numpy as np
import requests
import json
import time
import os

# file paths
INPUT_FILE = r"QuickGO-annotations-1742935775768-03 25 25 - NEW.xlsx"
OUTPUT_FILE = r"IUPRED_VALIDATED.xlsx"
SUMMARY_FILE = r"IUPRED_SUMMARY.csv"

IUPRED_MODE = "long"   # "long" or "short"
THRESHOLD = 0.5

CACHE_DIR = "iupred_cache"
os.makedirs(CACHE_DIR, exist_ok=True)


def get_iupred_scores(uniprot_id):
    cache_file = os.path.join(CACHE_DIR, f"{uniprot_id}_{IUPRED_MODE}.json")

    if os.path.exists(cache_file):
        with open(cache_file, "r") as f:
            data = json.load(f)
    else:
        url = f"http://iupred2a.elte.hu/iupred2a/{IUPRED_MODE}/{uniprot_id}.json"
        try:
            response = requests.get(url)
            if response.status_code != 200:
                print(f"[FAIL] {uniprot_id}")
                return None
            data = response.json()

            with open(cache_file, "w") as f:
                json.dump(data, f)

            time.sleep(0.5)

        except Exception as e:
            print(f"[ERROR] {uniprot_id}: {e}")
            return None

    scores = data["iupred2"]
    return scores


df = pd.read_excel(INPUT_FILE)

# Identify IDR columns dynamically
idr_columns = [col for col in df.columns if "Disordered Region" in col]

results = []


for idx, row in df.iterrows():

    uniprot_id = row.get("Protein accession")
    full_seq = row.get("Full Protein Sequence")

    if not isinstance(full_seq, str) or not isinstance(uniprot_id, str):
        continue

    print(f"Processing {idx+1}/{len(df)}: {uniprot_id}")

    scores = get_iupred_scores(uniprot_id)
    if scores is None:
        continue

    # sanity check
    if len(scores) != len(full_seq):
        print(f"[WARNING] Length mismatch for {uniprot_id}")
        continue

    for col in idr_columns:
        idr_seq = row.get(col)

        if not isinstance(idr_seq, str) or len(idr_seq) == 0:
            continue

        # Find IDR in full sequence
        start = full_seq.find(idr_seq)

        if start == -1:
            print(f"[MISS] Could not map IDR in {uniprot_id}")
            continue

        end = start + len(idr_seq)

        segment_scores = scores[start:end]

        if len(segment_scores) == 0:
            continue

        mean_score = np.mean(segment_scores)
        frac_disordered = np.sum(np.array(segment_scores) > THRESHOLD) / len(segment_scores)

        results.append({
            "UniProt_ID": uniprot_id,
            "IDR_column": col,
            "IDR_length": len(idr_seq),
            "Mean_IUPred": mean_score,
            "Fraction_Disordered": frac_disordered
        })


df_results = pd.DataFrame(results)
df_results.to_csv(SUMMARY_FILE, index=False)

print(f"\nSaved summary: {SUMMARY_FILE}")


if len(df_results) > 0:
    print("\nSUMMARY")
    print(f"Total IDRs analyzed: {len(df_results)}")
    print(f"Mean disorder score: {df_results['Mean_IUPred'].mean():.3f}")
    print(f"% IDRs > 0.5 fraction disordered: {(df_results['Fraction_Disordered'] > 0.5).mean()*100:.1f}%")
