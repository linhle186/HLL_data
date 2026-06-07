# %%
import argparse
import os
import sys
import glob
import scvelo as scv
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import anndata as ad

scv.logging.print_version()
scv.settings.verbosity = 3  # show errors and warnings
scv.settings.presenter_view = True  # sets matplotlib parameters for nice plots

# Add the script's folder directory to the Python path to import local modules
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_dir)

# Import local function
from get_qc_metrics import get_qc_metrics


# Helper function modified to load and parse explicit loom files provided by Snakemake
def load_and_label_loom(loom_path):
    # 1. Extract the sample/GSM ID from the filename
    filename = os.path.basename(loom_path)
    sample_id = filename.split('_Aligned')[0]

    print(f"Loading {sample_id} from: {loom_path}")

    # 2. Use Scanpy's reader for optimized matrix construction
    adata_subset = sc.read_loom(loom_path)

    # COMPATIBILITY CHECK: Decode byte strings to standard strings if loaded as bytes
    if adata_subset.obs_names.dtype == object and len(adata_subset.obs_names) > 0:
        if isinstance(adata_subset.obs_names[0], bytes):
            adata_subset.obs_names = adata_subset.obs_names.str.decode('utf-8')
            
    if adata_subset.var_names.dtype == object and len(adata_subset.var_names) > 0:
        if isinstance(adata_subset.var_names[0], bytes):
            adata_subset.var_names = adata_subset.var_names.str.decode('utf-8')

    # 3. Add the sample ID metadata
    adata_subset.obs['sample'] = sample_id

    # 4. Swap to readable Gene Names (Check your adata.var to see if it's 'Gene' or 'Accession')
    if 'Gene' in adata_subset.var.columns:
        adata_subset.var_names = adata_subset.var['Gene']

    # 5. Make GENE names unique
    adata_subset.var_names_make_unique()

    # 6. Make CELL barcodes unique across batches
    adata_subset.obs_names = sample_id + "_" + adata_subset.obs_names
    adata_subset.obs_names_make_unique()

    return adata_subset


def main():
    # 1. Parse Snakemake inputs, outputs, and dataset IDs
    parser = argparse.ArgumentParser(description="Downstream preprocessing and QC of loom files.")
    parser.add_argument(
        "--inputs", nargs="+", required=True, help="List of resolved loom files from Snakemake."
    )
    parser.add_argument(
        "--output", required=True, help="Destination path for final processed .h5ad output."
    )
    parser.add_argument(
        "--srp", required=True, help="SRP dataset ID (e.g., SRP250304)."
    )
    # FIX BUG-1: Added explicit --qc_dir parameter
    parser.add_argument(
        "--qc_dir", default=None, help="Directory for raw/light QC backups. Defaults to output directory's folder."
    )
    args = parser.parse_args()

    # 2. Create dictionary and load loom data dynamically
    loom_data = {}
    for loom_file in args.inputs:
        filename = os.path.basename(loom_file)
        sample_id = filename.split('_Aligned')[0]
        loom_data[sample_id] = load_and_label_loom(loom_file)

    # 3. Concatenate dataset
    adata = ad.concat(loom_data, join='outer')
    print("Merge complete! Total cells and genes:", adata.shape)

    # Preprocessing for RNA velocity
    print("Raw merged data shape:", adata.shape)
    adata_raw = adata.copy()

    # Print original proportions
    scv.utils.show_proportions(adata_raw)

    # --- MINIMAL QC FILTERING ---
    sc.pp.filter_cells(adata, min_genes=100)

    print("Minimal QC complete! Remaining cells and genes:", adata.shape)
    scv.utils.show_proportions(adata)

    # --- SAVE METRICS AND BACKUPS ---
    # FIX BUG-1: Dynamically assign backup folder, defaulting to parent folder of output
    qc_dir = args.qc_dir or os.path.dirname(os.path.abspath(args.output)) or "."
    os.makedirs(qc_dir, exist_ok=True)

    raw_h5ad_path = os.path.join(qc_dir, f"raw_{args.srp}.h5ad")
    light_h5ad_path = os.path.join(qc_dir, f"light_QC_{args.srp}.h5ad")

    adata_raw.write(raw_h5ad_path)
    adata.write(light_h5ad_path)

    # Save to the target output expected by Snakemake's dependency tracker
    adata.write(args.output)
    print(f"Processed file successfully saved to target output: {args.output}")


if __name__ == "__main__":
    main()