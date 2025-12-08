import pandas as pd
import numpy as np
import plotly.express as px
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import warnings
import os
import networkx as nx
from itertools import combinations

class VesiclePediaPipeline:
    """
    A Python tool for functional enrichment analysis.
    
    Features:   
    - Fisher's Exact Test for enrichment (Hypergeometric).
    - Benjamini-Hochberg FDR correction.
    - Support for custom databases (GMT format).
    - Interactive visualization using Plotly.
    """

    def __init__(self):
        """Initialize the VesiclePediaPipeline."""
        self.database = {}  # Stores term -> set of genes
        self.term_descriptions = {} # Stores term -> description (optional)
        self.universe = set() # All genes in the loaded database
        self.base_dir = os.path.dirname(os.path.abspath(__file__))
        self.experiment_file = os.path.join(self.base_dir, "..", "data", "VESICLEPEDIA_EXPERIMENT_DETAILS_5.1.txt")
        self.data_file = os.path.join(self.base_dir, "..", "data", "VESICLEPEDIA_PROTEIN_MRNA_DETAILS_5.1.txt")
        self.output_gmt = os.path.join(self.base_dir, "..", "data", "vesiclepedia_db.gmt")

    def load_gmt(self, filepath):
        """
        Loads a GMT (Gene Matrix Transposed) file.
        Format: TERM_ID \t DESCRIPTION \t GENE_1 \t GENE_2 ...
        """
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < 3:
                        continue
                    term_id = parts[0]
                    description = parts[1]
                    genes = set(parts[2:])
                    
                    self.database[term_id] = genes
                    self.term_descriptions[term_id] = description
                    self.universe.update(genes)
            print(f"Successfully loaded {len(self.database)} terms from {filepath}.")
        except Exception as e:
            print(f"Error loading GMT file: {e}")

    def load_custom_dict(self, data_dict):
        """
        Load a dictionary directly: {'Term_Name': [gene1, gene2, ...]}
        """
        for term, genes in data_dict.items():
            gene_set = set(genes)
            self.database[term] = gene_set
            self.term_descriptions[term] = term
            self.universe.update(gene_set)

    def run_enrichment(self, gene_list, background_list=None, correction_method='fdr_bh'):
        """
        Performs the enrichment analysis (Fisher's Exact Test).
        Returns a DataFrame sorted by p-value.
        """
        # 1. Standardize inputs
        genes_of_interest = set(gene_list)
        
        if background_list is None:
            background = self.universe
            warnings.warn("No background list provided. Using all genes in the database as background. This may inflate significance.")
        else:
            background = set(background_list)

        # Filter genes of interest to ensure they are in the background
        genes_of_interest = genes_of_interest.intersection(background)
        
        # 2. Prepare stats
        results = []
        total_background_count = len(background)
        total_interest_count = len(genes_of_interest)

        # 3. Iterate through each term in the database
        for term, term_genes in self.database.items():
            term_genes_in_bg = term_genes.intersection(background)
            
            if not term_genes_in_bg:
                continue

            # Calculate Contingency Table for Fisher's Exact Test
            hits = genes_of_interest.intersection(term_genes_in_bg)
            a = len(hits)               # In Interest AND In Term
            b = total_interest_count - a # In Interest AND NOT In Term
            c = len(term_genes_in_bg) - a # In Background AND In Term
            d = (total_background_count - total_interest_count) - c # In Background AND NOT In Term

            if a > 0:
                odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
                
                results.append({
                    'Term_ID': term,
                    'Description': self.term_descriptions.get(term, term),
                    'p_value': p_value,
                    'Odds_Ratio': odds_ratio,
                    'Count': a, 
                    'Term_Size': len(term_genes_in_bg),
                    'Genes': ", ".join(list(hits)) # Added space for better tooltip reading
                })

        # 4. Create DataFrame and Correct for Multiple Testing
        if not results:
            return pd.DataFrame()

        df = pd.DataFrame(results)
        
        reject, pvals_corrected, _, _ = multipletests(df['p_value'], method=correction_method)
        df['adj_p_value'] = pvals_corrected
        
        df = df.sort_values(by='p_value')
        
        return df

    def plot_results(self, df, top_n=10, title="Enrichment Analysis", save_path=None):
        """
        Visualizes the top enriched terms using an interactive Plotly bar chart.
        """
        if df.empty:
            print("No results to plot.")
            return

        # Select top N terms
        plot_data = df.head(top_n).copy()
        
        # Calculate -log10(p-value) for visualization
        plot_data['-log10(p)'] = -np.log10(plot_data['p_value'])
        
        # Create Interactive Bar Chart
        fig = px.bar(plot_data, 
                     x='-log10(p)', 
                     y='Term_ID',
                     orientation='h',
                     title=title,
                     color='-log10(p)', # Color bars by significance strength
                     color_continuous_scale='Viridis',
                     hover_data={
                         'Term_ID': False, # Already on axis
                         '-log10(p)': ':.2f',
                         'p_value': ':.2e',
                         'adj_p_value': ':.2e',
                         'Count': True,
                         'Genes': True # Shows gene names on hover
                     },
                     labels={'Term_ID': 'Pathway / Term', '-log10(p)': '-log10(p-value)'}
                    )

        # Ensure the most significant (top) item is at the top of the chart
        fig.update_layout(yaxis={'categoryorder':'total ascending'})
        
        if save_path:
            # If saving to HTML, the interactivity is preserved
            fig.write_html(save_path)
            print(f"Interactive plot saved to {save_path}")
        
        return fig

    def convert_vesiclepedia_to_gmt(
        self,
        experiment_file,
        data_file,
        output_gmt
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
        self.merged_dataset = merged

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

    def convert_to_gmt(
        self, 
        experiment_file, 
        data_file, 
        mapping, 
        output_gmt,
        sep='\t',
        allowed_species=None
    ):
        """
        Converts linked files into GMT, optionally filtering by species.
        """
        print(f"Reading files with separator: {repr(sep)}")
        try:
            exp_df = pd.read_csv(experiment_file, sep=sep, encoding='latin-1', on_bad_lines='skip')
            data_df = pd.read_csv(data_file, sep=sep, encoding='latin-1', on_bad_lines='skip')
        except Exception as e:
            print(f"Error reading files: {e}")
            return "Error reading files"

        # 1. Standardize Merge Keys
        exp_id = mapping['id_col_exp']
        data_id = mapping['id_col_data']
        
        # Standardize column names for easier matching
        exp_df.columns = exp_df.columns.str.strip()
        data_df.columns = data_df.columns.str.strip()
        
        if exp_id not in exp_df.columns or data_id not in data_df.columns:
            return f"Error: ID columns not found. Looked for '{exp_id}' and '{data_id}'"

        exp_df = exp_df.rename(columns={exp_id: 'JOIN_ID'})
        data_df = data_df.rename(columns={data_id: 'JOIN_ID'})

        # 2. Merge
        print("Merging datasets...")
        merged = pd.merge(data_df, exp_df, on='JOIN_ID', how='inner')
        
        if allowed_species:
            # We try to find the species column automatically to filter
            species_col = next((c for c in merged.columns if 'SPECIES' in c.upper()), None)
            if species_col:
                print(f"Filtering for species: {allowed_species}")
                # Filter rows where the species column matches our allowed list
                # If allowed_species is a single string, wrap it in list
                if isinstance(allowed_species, str):
                    allowed_species = [allowed_species]
                    
                merged = merged[merged[species_col].isin(allowed_species)]
                
                if merged.empty:
                    return f"Error: Filtering for '{allowed_species}' resulted in empty data."
            else:
                print("Warning: allowed_species provided but no 'SPECIES' column found in data.")
        # ---------------------------------

        self.merged_dataset = merged 
        
        # 3. Standardize Gene Column
        gene_col_user = mapping['gene_col']
        if gene_col_user not in merged.columns:
             return f"Error: Gene column '{gene_col_user}' not found in data."
             
        # 4. Build GMT
        gmt_lines = []
        categories = mapping['categories']
        
        print(f"Building sets for categories: {categories}")
        
        for category in categories:
            if category not in merged.columns:
                continue
                
            grouped = merged.groupby(category)[gene_col_user].apply(lambda x: set(x.dropna()))
            
            for group_name, genes in grouped.items():
                clean_name = str(group_name).strip().replace(" ", "_").upper()
                clean_category = category.upper().replace(" ", "_")
                
                term_id = f"{clean_category}:{clean_name}"
                description = f"Genes found in {category}: {group_name}"
                
                # Filter invalid genes
                valid_genes = [str(g).upper() for g in genes if len(str(g)) > 1]
                
                if len(valid_genes) >= 5: 
                    line = f"{term_id}\t{description}\t" + "\t".join(valid_genes)
                    gmt_lines.append(line)

        # 5. Write to File
        with open(output_gmt, 'w', encoding='utf-8') as f:
            f.write("\n".join(gmt_lines))
        
        return f"Success! Created {output_gmt} with {len(gmt_lines)} gene sets."

    # Utility to create a mock Vesiclepedia database for testing
    def create_mock_vesiclepedia_gmt(self, filename="mock_vesiclepedia.gmt"):
        with open(filename, 'w') as f:
            f.write("EV_SURFACE\tExtracellular Vesicle Surface\tCD9\tCD63\tCD81\tTSG101\tALIX\n")
            f.write("EXOSOME_CORE\tExosome Core Proteins\tHSP70\tHSP90\tGAPDH\tACTB\n")
            f.write("MITOCHONDRIA\tMitochondrial Components\tMT-CO1\tMT-ND1\tCYCS\n")
            f.write("NUCLEUS\tNuclear Components\tTP53\tBRCA1\tHIST1H1\n")

    def build_enrichment_network(self, df, similarity_threshold=0.2):
        """
        Creates a NetworkX graph where:
        - Nodes = Enriched Terms
        - Edges = Significant gene overlap (Jaccard Index > threshold)
        """
        G = nx.Graph()
        
        # 1. Create Nodes
        # We use the DataFrame index or Term_ID as the node identifier
        for _, row in df.iterrows():
            # Scale node size by significance (-log10 p-value)
            # We assume -log10(p) is already calculated or we calc it here
            nlogp = -np.log10(row['p_value']) if row['p_value'] > 0 else 50
            
            G.add_node(
                row['Term_ID'], 
                title=f"{row['Term_ID']}\nGenes: {row['Count']}\nP: {row['p_value']:.2e}", # Tooltip
                label=row['Term_ID'],
                value=nlogp, # Size of node
                group='Enriched Term'
            )

        # 2. Create Edges (based on gene overlap)
        # We need the actual gene sets for these terms. 
        # The 'Genes' column in your df is a string "GeneA, GeneB". Let's split it back.
        term_genes = {}
        for term in df['Term_ID']:
            # Retrieve from the main database to ensure we have the full set
            # OR use the 'Genes' col from the df (which is just the intersecting genes)
            # Using intersecting genes (Genes col) is usually better for 'Contextual Similarity'
            genes_str = df.loc[df['Term_ID'] == term, 'Genes'].values[0]
            term_genes[term] = set([g.strip() for g in genes_str.split(',') if g.strip()])

        # Calculate Jaccard Index for every pair
        terms = list(term_genes.keys())
        for term_a, term_b in combinations(terms, 2):
            set_a = term_genes[term_a]
            set_b = term_genes[term_b]
            
            intersection = len(set_a.intersection(set_b))
            union = len(set_a.union(set_b))
            
            if union > 0:
                jaccard = intersection / union
                
                # If overlap is significant, draw an edge
                if jaccard > similarity_threshold:
                    # Width of edge depends on how similar they are
                    G.add_edge(term_a, term_b, weight=jaccard, title=f"Overlap: {jaccard:.2f}")

        return G