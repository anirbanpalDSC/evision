import pandas as pd
from src.vp_pipeline import EVision, create_mock_vesiclepedia_gmt

def run_demo():
    print("--- Starting EVision Demo ---")
    
    # 1. Setup a mock database (simulating Vesiclepedia or GO)
    # In a real scenario, you would download a .gmt file from GSEA msigdb or Vesiclepedia
    gmt_file = "mock_vesiclepedia.gmt"
    create_mock_vesiclepedia_gmt(gmt_file)
    
    # 2. Initialize the tool
    enricher = EVision()
    enricher.load_gmt(gmt_file)
    
    # 3. Create synthetic data
    # Let's say we did Proteomics on an EV sample and found these proteins:
    # (Intentionally selecting EV markers to show enrichment)
    my_gene_list = ["CD9", "CD63", "CD81", "HSP70", "GAPDH", "TP53"]
    
    # Define a background (Simulating the total proteome we could have detected)
    # It must contain the gene list + others
    my_background = ["CD9", "CD63", "CD81", "TSG101", "ALIX", 
                     "HSP70", "HSP90", "GAPDH", "ACTB", 
                     "MT-CO1", "MT-ND1", "CYCS", 
                     "TP53", "BRCA1", "HIST1H1", 
                     "RANDOM1", "RANDOM2", "RANDOM3"]
    
    print(f"\nInput List: {my_gene_list}")
    print(f"Background Size: {len(my_background)}")
    
    # 4. Run Analysis
    print("\nRunning Fisher's Exact Test...")
    results = enricher.run_enrichment(gene_list=my_gene_list, 
                                      background_list=my_background)
    
    # 5. Display Results
    if not results.empty:
        print("\n--- Enrichment Results (Top 5) ---")
        # Display specific columns for readability
        cols = ['Term_ID', 'p_value', 'adj_p_value', 'Count', 'Term_Size', 'Genes']
        print(results[cols].head(5).to_string(index=False))
        
        # 6. Plot
        print("\nGenerating plot...")
        enricher.plot_results(results, title="Enrichment of Mock EV Data")
    else:
        print("No significant enrichment found.")

if __name__ == "__main__":
    run_demo()