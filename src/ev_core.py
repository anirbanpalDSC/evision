import pandas as pd
import numpy as np
import plotly.express as px
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import warnings

class EVision:
    """
    A Python tool for functional enrichment analysis.
    
    Features:
    - Fisher's Exact Test for enrichment (Hypergeometric).
    - Benjamini-Hochberg FDR correction.
    - Support for custom databases (GMT format).
    - Interactive visualization using Plotly.
    """

    def __init__(self):
        """Initialize the EVision tool."""
        self.database = {}  # Stores term -> set of genes
        self.term_descriptions = {} # Stores term -> description (optional)
        self.universe = set() # All genes in the loaded database

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
        else:
            fig.show()

# Utility to create a mock Vesiclepedia database for testing
def create_mock_vesiclepedia_gmt(filename="mock_vesiclepedia.gmt"):
    with open(filename, 'w') as f:
        f.write("EV_SURFACE\tExtracellular Vesicle Surface\tCD9\tCD63\tCD81\tTSG101\tALIX\n")
        f.write("EXOSOME_CORE\tExosome Core Proteins\tHSP70\tHSP90\tGAPDH\tACTB\n")
        f.write("MITOCHONDRIA\tMitochondrial Components\tMT-CO1\tMT-ND1\tCYCS\n")
        f.write("NUCLEUS\tNuclear Components\tTP53\tBRCA1\tHIST1H1\n")