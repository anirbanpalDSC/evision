import pandas as pd
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
experiment_file = os.path.join(BASE_DIR, "..", "data", "VESICLEPEDIA_EXPERIMENT_DETAILS_5.1.txt")
data_file = os.path.join(BASE_DIR, "..", "data", "VESICLEPEDIA_PROTEIN_MRNA_DETAILS_5.1.txt")
output_gmt = os.path.join(BASE_DIR, "..", "data", "vesiclepedia_db.gmt")

def convert_vesiclepedia_to_gmt(
    experiment_file=experiment_file,
    data_file=data_file,
    output_gmt=output_gmt
):
    """
    Converts raw Vesiclepedia text files into a GMT file for FunEnrich.
    """
    print(f"Reading {experiment_file}...")
    try:
        # Attempt to read Experiment Details
        # Note: Vesiclepedia files are often tab-delimited
        exp_df = pd.read_csv(experiment_file, sep='\t', encoding='latin-1', on_bad_lines='skip')
        exp_df.columns = exp_df.columns.str.upper()
    except Exception as e:
        print(f"Error reading experiment file: {e}")
        return

    print(f"Reading {data_file}...")
    try:
        # Attempt to read Protein/mRNA Details
        data_df = pd.read_csv(data_file, sep='\t', encoding='latin-1', on_bad_lines='skip')
        data_df.columns = data_df.columns.str.upper()
    except Exception as e:
        print(f"Error reading data file: {e}")
        return

    # --- Standardize Column Names ---
    # Strip whitespace from column headers
    exp_df.columns = exp_df.columns.str.strip()
    data_df.columns = data_df.columns.str.strip()

    print("Merging datasets...")
    # The common link is typically 'EXPERIMENT ID'
    if 'EXPERIMENT ID' not in exp_df.columns or 'EXPERIMENT ID' not in data_df.columns:
        print("Error: Could not find 'EXPERIMENT ID' column in one of the files.")
        print(f"Exp Cols: {exp_df.columns}")
        print(f"Data Cols: {data_df.columns}")
        return

    # Merge metadata with gene data
    merged = pd.merge(data_df, exp_df, on='EXPERIMENT ID', how='inner')

    print(f"Merged dataset has {len(merged)} rows.")

    # --- Define Categories to Create Sets For ---
    # We want to create gene sets for different metadata categories
    categories_to_group_by = [
        'VESICLE TYPE', 
        'SAMPLE NAME', 
        'ISOLATION METHOD', 
        'SPECIES'
    ]

    gene_col = 'GENE SYMBOL' # Adjust if your file uses 'Gene Name' or 'Entrez ID'
    
    if gene_col not in merged.columns:
        print(f"Warning: '{gene_col}' not found. Searching for alternatives...")
        possible_cols = [c for c in merged.columns if 'GENE' in c or 'SYMBOL' in c]
        if possible_cols:
            gene_col = possible_cols[0]
            print(f"Using '{gene_col}' as gene identifier.")
        else:
            print("Error: Could not identify a Gene column.")
            return

    # --- Build GMT Data ---
    gmt_lines = []
    
    print("Building Gene Sets...")
    for category in categories_to_group_by:
        if category not in merged.columns:
            continue
            
        # Group by the category (e.g., group by "Vesicle Type")
        grouped = merged.groupby(category)[gene_col].apply(lambda x: set(x.dropna()))
        
        for group_name, genes in grouped.items():
            # Clean up the group name
            clean_name = str(group_name).strip().replace(" ", "_").upper()
            clean_category = category.upper().replace(" ", "_")
            
            # Create a unique ID: e.g., VESICLE_TYPE:EXOSOMES
            term_id = f"{clean_category}:{clean_name}"
            
            # Description
            description = f"Genes found in {category}: {group_name}"
            
            # Remove empty strings or invalid genes
            valid_genes = [g for g in genes if isinstance(g, str) and len(g) > 1]
            
            if len(valid_genes) > 5: # Only keep sets with at least 5 genes
                # GMT Format: ID \t Description \t Gene1 \t Gene2 ...
                line = f"{term_id}\t{description}\t" + "\t".join(valid_genes)
                gmt_lines.append(line)

    # --- Write to File ---
    with open(output_gmt, 'w', encoding='utf-8') as f:
        f.write("\n".join(gmt_lines))
    
    print(f"Success! Created {output_gmt} with {len(gmt_lines)} gene sets.")

if __name__ == "__main__":
    convert_vesiclepedia_to_gmt()