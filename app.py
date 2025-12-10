import streamlit as st
import os
import pandas as pd
import numpy as np
from PIL import Image
import tempfile
from io import StringIO
from pyvis.network import Network
import streamlit.components.v1 as components
from src.vp_pipeline import VesiclePediaPipeline

# ---------------------------------------------------------------------
# App Config
# ---------------------------------------------------------------------
st.set_page_config(
    page_title="EVision: Extracellular Vesicles Data Analysis Suite",
    layout="wide"
)

logo_path = "images/evision.png"
# st.logo(
#     image=logo_path,
#     link="https://github.com/anirbanpalDSC/evision",
#     # Adjust the size (small, medium, large)
#     size="large"
# )

if os.path.exists(logo_path):
    st.image(logo_path, width=100)
else:
    st.warning("Logo not found. Please add 'images/evision.png' to your project.")

st.markdown("[GitHub Repository](https://github.com/anirbanpalDSC/evision)")

# ---------------------------------------------------------------------
# Initialize Pipeline in Session
# ---------------------------------------------------------------------
if "pipeline" not in st.session_state:
    st.session_state.pipeline = VesiclePediaPipeline()

pipeline = st.session_state.pipeline

# ---------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------
def parse_gene_list(text: str):
    """Parse a comma- or newline-separated gene list into a cleaned list."""
    if not text:
        return []
    genes = [
        g.strip().upper()
        for g in text.replace(",", "\n").split("\n")
        if g.strip()
    ]
    return genes

# ---------------------------------------------------------------------
# Sidebar Navigation
# ---------------------------------------------------------------------
st.sidebar.title("🧬 Analysis Suite")

with st.sidebar.expander("🧪 Vesiclepedia Tools", expanded=True):
    vesicle_option = st.radio(
        "Vesiclepedia Modules",
        [
            "🏠 Home",
            "📥 Upload & Convert Vesiclepedia Files",
            "📚 Load GMT Database",
            "🔍 Enrichment Analysis",
            "⚖️ Compare Two Gene Lists",
            "🧪 EV QC (MISEV Marker Check)",
            "🔧 Build Custom Gene Sets",
            "🗂 Explore Metadata",
            "📜 About",
        ],
        key="vesiclepedia_option"
    )

# PRIDE section is informational only for now
with st.sidebar.expander("🧬 PRIDE Tools (Coming Soon)", expanded=False):
    st.markdown(
        """
        PRIDE integration is planned for a future release.

        Planned capabilities:
        - Search PRIDE datasets related to EVs
        - Pull processed protein lists from PRIDE
        - Run enrichment directly on PRIDE-derived gene lists
        - Explore project metadata and QC
        """
    )

# =========================================================
# SEARCHABLE DROPDOWN SPECIES FILTER
# =========================================================

st.sidebar.divider()
st.sidebar.header("🌍 Global Species Filter")

# Helper to load data if uploaded directly in sidebar
def load_sidebar_data(uploaded_file):
    try:
        df = pd.read_csv(uploaded_file, sep=None, engine='python', encoding='latin-1', on_bad_lines='skip')
        df.columns = df.columns.str.upper().str.strip()
        return df
    except Exception as e:
        st.sidebar.error(f"Error reading file: {e}")
        return None

# Check if metadata is loaded
if hasattr(pipeline, "merged_dataset") and pipeline.merged_dataset is not None:
    df = pipeline.merged_dataset
    
    # 1. Identify Species Column
    species_col = next((c for c in df.columns if 'SPECIES' in c), None)
    
    if species_col:
        # Get unique species list
        all_species = sorted(df[species_col].dropna().unique().tolist())
        
        # Add an "All Species" option at the very top
        search_options = ["All Species"] + all_species
        
        # 2. Searchable Dropdown Widget
        with st.sidebar.expander("⚙️ Filter Settings", expanded=True):
            # st.selectbox allows typing to filter the list automatically
            selected_species = st.selectbox(
                "Select Species (Type to Search):",
                options=search_options,
                index=0, # Default to "All Species"
                key="global_species_selection",
                help="Select the specific species of your experiment. This refines the background for accurate statistics."
            )
         
            if st.button("Apply Filter", type="primary"):
                if selected_species == "All Species":
                    # Use everything
                    subset = df
                    st.toast("Using full database (All Species).")
                else:
                    # Filter for the ONE specific species
                    subset = df[df[species_col] == selected_species]
                
                # Find the gene column to rebuild the universe
                gene_col = next((c for c in subset.columns if 'GENE' in c or 'SYMBOL' in c), None)
                
                if gene_col:
                    valid_genes = set(subset[gene_col].dropna().astype(str).str.upper())
                    
                    # UPDATE THE PIPELINE UNIVERSE
                    pipeline.universe = valid_genes
                    st.session_state['custom_universe'] = valid_genes
                    
                    st.success(f"✅ Filter Applied! Analysis restricted to: {selected_species}")
                    st.caption(f"Universe Size: {len(valid_genes)} genes")
                else:
                    st.error("Could not auto-detect gene column.")
    else:
        st.sidebar.warning("Loaded data has no 'SPECIES' column.")

else:
    # Fallback Uploader
    st.sidebar.info("Upload **Protein/mRNA Details** to enable filtering.")
    sidebar_file = st.sidebar.file_uploader("📂 Load Metadata", type=["txt", "csv", "tsv"], key="sidebar_uploader")
    
    if sidebar_file:
        with st.spinner("Loading metadata..."):
            loaded_df = load_sidebar_data(sidebar_file)
            if loaded_df is not None:
                pipeline.merged_dataset = loaded_df
                st.rerun()

# Main mode is driven entirely by Vesiclepedia options for now
app_mode = vesicle_option

# ---------------------------------------------------------------------
# Home
# ---------------------------------------------------------------------
if app_mode == "🏠 Home":
    st.title("EVision: Extracellular Vesicles Data Analysis Suite")
    st.markdown(
        """
        This Streamlit app provides a full-scale analysis environment for exploring
        extracellular vesicle (EV) datasets from Vesiclepedia, generating gene sets,
        performing functional enrichment, building custom pathways, and checking EV sample quality.

        **Use Cases Enabled:**
        - Convert raw Vesiclepedia text files → GMT
        - Load Vesiclepedia or custom gene sets
        - Run enrichment analysis
        - Compare healthy vs diseased lists
        - Check EV QC markers (CD9, CD63, CD81, TSG101)
        - Browse metadata
        - Build custom pathways/dictionaries
        """
    )

# ---------------------------------------------------------------------
# Upload & Convert Vesiclepedia Raw Files
# ---------------------------------------------------------------------

elif app_mode == "📥 Upload & Convert Vesiclepedia Files":
    st.title("📥 Convert Raw Files → GMT")
    st.info("Upload your Experiment Metadata and Protein/mRNA data files.")

    col1, col2 = st.columns(2)
    with col1:
        exp_file = st.file_uploader("1. Experiment/Metadata File", type=["txt", "tsv", "csv"])
    with col2:
        data_file = st.file_uploader("2. Protein/Gene Data File", type=["txt", "tsv", "csv"])

    separator = st.selectbox("File Separator", ["Tab (\\t)", "Comma (,)"], index=0)
    sep_char = '\t' if "Tab" in separator else ','

    # --- THE COLUMN MAPPER UI ---
    if exp_file and data_file:
        st.divider()
        st.subheader("🔧 Map Your Columns")
        
        # 1. Read headers only (efficient)
        exp_file.seek(0)
        data_file.seek(0)
        exp_headers = pd.read_csv(exp_file, sep=sep_char, nrows=0, encoding='latin-1').columns.tolist()
        data_headers = pd.read_csv(data_file, sep=sep_char, nrows=0, encoding='latin-1').columns.tolist()
        
        # Reset pointers for actual processing later
        exp_file.seek(0)
        data_file.seek(0)

        c1, c2, c3 = st.columns(3)
        
        # User selects the linking ID in Experiment File
        id_col_exp = c1.selectbox("Experiment ID Column (Metadata File)", exp_headers, 
                                index=exp_headers.index('EXPERIMENT ID') if 'EXPERIMENT ID' in exp_headers else 0)
        
        # User selects the linking ID in Data File
        id_col_data = c2.selectbox("Experiment ID Column (Data File)", data_headers,
                                index=data_headers.index('EXPERIMENT ID') if 'EXPERIMENT ID' in data_headers else 0)
        
        # User selects the Gene Column
        gene_col = c3.selectbox("Gene Symbol Column", data_headers,
                                index=data_headers.index('GENE SYMBOL') if 'GENE SYMBOL' in data_headers else 0)

        # User selects which metadata columns to create sets from
        st.write("Select categories to create gene sets from:")
        categories = st.multiselect("Metadata Categories", exp_headers, default=[h for h in exp_headers if 'TYPE' in h or 'SPECIES' in h])

        st.divider()
        st.subheader("🌍 Species Filter (Optional)")
        st.info("Filter the GMT to a specific species before creating gene sets. Leave as 'All Species' to include everything.")
        
        # Quick preview to get species list
        exp_file.seek(0)
        exp_preview = pd.read_csv(exp_file, sep=sep_char, encoding='latin-1', on_bad_lines='skip')
        exp_file.seek(0)
        
        species_col = next((c for c in exp_preview.columns if 'SPECIES' in c.upper()), None)
        
        if species_col:
            available_species = sorted(exp_preview[species_col].dropna().unique().tolist())
            species_options = ["All Species"] + available_species
            
            selected_species = st.selectbox(
                "Select Species:",
                options=species_options,
                index=0,
                help="This will filter the entire dataset before creating gene sets. Recommended for species-specific analysis."
            )
        else:
            st.warning("No 'SPECIES' column found. Filter not available.")
            selected_species = "All Species"

        # --- CONVERSION BUTTON ---
        if st.button("Convert to GMT"):
            if not categories:
                st.error("Please select at least one metadata category.")
            else:
                # Save to temp
                with tempfile.NamedTemporaryFile(delete=False, suffix=".txt", mode="wb") as f_exp:
                    f_exp.write(exp_file.getvalue())
                    exp_path = f_exp.name
                
                with tempfile.NamedTemporaryFile(delete=False, suffix=".txt", mode="wb") as f_data:
                    f_data.write(data_file.getvalue())
                    data_path = f_data.name

                # Determine output filename based on species
                if selected_species == "All Species":
                    output_gmt = "vesiclepedia_all_species.gmt"
                else:
                    safe_name = selected_species.replace(" ", "_").lower()
                    output_gmt = f"vesiclepedia_{safe_name}.gmt"
                
                # Create the mapping dict
                mapping_config = {
                    'id_col_exp': id_col_exp,
                    'id_col_data': id_col_data,
                    'gene_col': gene_col,
                    'categories': categories,
                    'species_filter': None if selected_species == "All Species" else selected_species  # ⭐ NEW
                }

                # Run the pipeline function
                with st.spinner("Converting... This may take a moment."):
                    status = pipeline.convert_to_gmt(
                        experiment_file=exp_path,
                        data_file=data_path,
                        mapping=mapping_config,
                        output_gmt=output_gmt,
                        sep=sep_char
                    )

                if "Success" in status:
                    st.success(status)
                    
                    # Show what was filtered
                    if selected_species != "All Species":
                        st.info(f"✅ GMT filtered to: **{selected_species}**")
                    
                    with open(output_gmt, "rb") as f_out:
                        st.download_button(
                            "⬇ Download GMT", f_out, output_gmt
                        )
                else:
                    st.error(status)

# ---------------------------------------------------------------------
# Load GMT Database (UPDATED WITH TERM INSPECTOR)
# ---------------------------------------------------------------------
elif app_mode == "📚 Load GMT Database":
    st.title("📚 Load Gene Sets (GMT)")

    # --- PART A: PERSISTENT STATUS DASHBOARD ---
    num_sets = len(pipeline.database)
    
    st.subheader("Current Database Status")
    m1, m2 = st.columns(2)
    m1.metric("Loaded Gene Sets", num_sets)
    m2.metric("Total Unique Genes", len(pipeline.universe))

    # --- NEW FEATURE: TERM INSPECTOR ---
    if num_sets > 0:
        st.divider()
        st.subheader("🔍 Inspect Specific Term & Copy Genes")
        st.info("Select a term below to see the genes inside it. Useful for verifying data or creating test lists.")

        # 1. Select a term
        all_terms = sorted(list(pipeline.database.keys()))
        selected_term = st.selectbox("Select a Pathway / Term:", all_terms)

        if selected_term:
            genes_in_term = list(pipeline.database[selected_term])
            
            # 2. Show details
            c1, c2 = st.columns([1, 3])
            with c1:
                st.write(f"**ID:** `{selected_term}`")
                st.write(f"**Count:** {len(genes_in_term)} genes")
            with c2:
                desc = pipeline.term_descriptions.get(selected_term, "No description")
                st.write(f"**Description:** {desc}")

            # 3. Copyable Text Area
            st.caption("Genes in this set (Copy this list to test the Enrichment Analysis):")
            st.text_area(
                "Gene List", 
                value=", ".join(genes_in_term), 
                height=150,
                help="Click inside, press Ctrl+A (Select All) then Ctrl+C (Copy)."
            )

    st.divider()

    # --- PART B: THE UPLOADER (Keep existing logic) ---
    st.subheader("Add More Gene Sets")
    gmt_file = st.file_uploader("Upload GMT File", type=["gmt"], key="gmt_uploader")

    if gmt_file:
        file_hash = hash(gmt_file.getvalue())
        if "last_loaded_file" not in st.session_state or st.session_state.last_loaded_file != file_hash:
            path = tempfile.mktemp(suffix=".gmt")
            with open(path, "w", encoding="utf-8") as f:
                f.write(gmt_file.getvalue().decode())
            pipeline.load_gmt(path)
            st.session_state.last_loaded_file = file_hash
            st.rerun()
            
    # Clear Database Button
    if num_sets > 0:
        if st.button("🗑️ Clear Database", type="secondary"):
            pipeline.database = {}
            pipeline.term_descriptions = {}
            pipeline.universe = set()
            st.rerun()

# ---------------------------------------------------------------------
# Enrichment Analysis
# ---------------------------------------------------------------------
elif app_mode == "🔍 Enrichment Analysis":
    st.title("🔍 Enrichment Analysis")
    st.markdown("Run a Fisher's Exact Test to find over-represented terms in your gene list.")

    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("1. Your Gene List")
        gene_text = st.text_area(
            "Paste genes (Interests)", 
            height=200, 
            help="These are the 'significant' genes from your experiment."
        )

    with col2:
        st.subheader("2. Background Universe")
        st.markdown(
            """
            <small><b>Why upload a background?</b> To calculate accurate p-values, 
            the statistical model needs to know which genes were <i>screened</i> 
            but <i>not selected</i>.</small>
            """, 
            unsafe_allow_html=True
        )
        
        bg_paste = st.text_area("Option A: Paste Background", height=100)
        bg_file = st.file_uploader("Option B: Upload Background File (.txt/.csv)", type=["txt", "csv"])

    # Settings
    with st.expander("⚙️ Advanced Statistics Settings"):
        correction_method = st.selectbox(
            "Multiple Testing Correction", 
            ['fdr_bh', 'bonferroni', 'holm', 'none'],
            index=0
        )
        top_n = st.number_input("Top N terms to plot", 5, 50, 15)

    # --- ACTION: RUN ANALYSIS ---
    if st.button("🚀 Run Analysis", type="primary"):
        # 1. Parse Inputs
        genes = parse_gene_list(gene_text)
        background = []
        
        if bg_file is not None:
            stringio = StringIO(bg_file.getvalue().decode("utf-8"))
            background = parse_gene_list(stringio.read())
            st.success(f"Loaded Custom Universe: {len(background)} genes")
        elif bg_paste:
            background = parse_gene_list(bg_paste)
            st.info(f"Loaded Custom Universe (from text): {len(background)} genes")
        else:
            background = None 

        if not genes:
            st.error("Please provide at least one gene in the gene list.")
        else:
            # 2. Run Pipeline
            df_results = pipeline.run_enrichment(
                genes, 
                background_list=background, 
                correction_method=correction_method
            )
            
            # 3. SAVE TO SESSION STATE (The Critical Fix)
            st.session_state['enrichment_results'] = df_results
            # We also save the 'top_n' setting used at the time of running, 
            # though the user can adjust the slider later if we wanted.
            st.session_state['last_top_n'] = top_n

    # --- VIEW: DISPLAY RESULTS ---
    # This block runs if results exist in memory, regardless of button press
    if 'enrichment_results' in st.session_state:
        df = st.session_state['enrichment_results']
        
        # Ensure top_n is consistent
        current_top_n = st.session_state.get('last_top_n', 15)

        if df.empty:
            st.warning("No enriched terms found.")
        else:
            st.divider()
            st.subheader("Results")
            
            # Formatting for display
            display_df = df.copy()
            if 'p_value' in display_df.columns:
                display_df['p_value'] = display_df['p_value'].map('{:.2e}'.format)
            if 'adj_p_value' in display_df.columns:
                display_df['adj_p_value'] = display_df['adj_p_value'].map('{:.2e}'.format)
            
            st.dataframe(display_df, use_container_width=True)

            # --- VISUALIZATION TAB ---
            st.subheader("Visualization")
            tab1, tab2 = st.tabs(["📊 Bar Chart", "🕸️ Enrichment Map"])
            
            with tab1:
                try:
                    # Pass top_n dynamically
                    fig = pipeline.plot_results(df, top_n=current_top_n)
                    if fig is not None:
                        st.plotly_chart(fig, use_container_width=True)
                except TypeError:
                     st.info("Visuals handled internally.")

            with tab2:
                col_net1, col_net2 = st.columns(2)
                with col_net1:
                    # This slider now works because it is OUTSIDE the button block
                    sim_threshold = st.slider(
                        "Edge Similarity Threshold", 
                        0.0, 1.0, 0.2, 0.05,
                        key="slider_threshold" # Unique key is good practice
                    )
                with col_net2:
                    physics = st.checkbox("Keep Physics Active", value=False)

                top_df = df.head(current_top_n) 
                nx_graph = pipeline.build_enrichment_network(top_df, similarity_threshold=sim_threshold)
                
                if nx_graph.number_of_nodes() > 0:
                    net = Network(height="750px", width="100%", bgcolor="#ffffff", font_color="black")
                    net.from_nx(nx_graph)
                    
                    physics_bool = "true" if physics else "false"
                    
                    net.set_options(f"""
                    var options = {{
                      "physics": {{
                        "enabled": {physics_bool},
                        "forceAtlas2Based": {{
                          "gravitationalConstant": -50,
                          "centralGravity": 0.01,
                          "springLength": 100,
                          "springConstant": 0.08,
                          "damping": 0.4
                        }},
                        "minVelocity": 0.75,
                        "solver": "forceAtlas2Based"
                      }}
                    }}
                    """)

                    try:
                        path = tempfile.mktemp(suffix=".html")
                        net.save_graph(path)
                        with open(path, 'r', encoding='utf-8') as f:
                            html_string = f.read()
                        
                        components.html(html_string, height=800, scrolling=True)
                        
                    except Exception as e:
                        st.error(f"Error generating network: {e}")
                else:
                    st.warning("No nodes to display at this threshold.")
            
            # Download
            csv = df.to_csv(index=False)
            st.download_button("⬇ Download Results CSV", csv, "enrichment_results.csv", "text/csv")

# ---------------------------------------------------------------------
# Compare Two Gene Lists
# ---------------------------------------------------------------------
elif app_mode == "⚖️ Compare Two Gene Lists":
    st.title("⚖️ Compare Two Gene Sets")

    col1, col2 = st.columns(2)
    with col1:
        genes_A = st.text_area("Gene List A")
    with col2:
        genes_B = st.text_area("Gene List B")

    if st.button("Compare Enrichment"):
        list_A = parse_gene_list(genes_A)
        list_B = parse_gene_list(genes_B)

        if not list_A or not list_B:
            st.error("Please provide both Gene List A and Gene List B.")
        else:
            dfA = pipeline.run_enrichment(list_A)
            dfB = pipeline.run_enrichment(list_B)

            st.subheader("Enrichment — List A")
            st.dataframe(dfA)
            try:
                figA = pipeline.plot_results(dfA, top_n=10, title="Enrichment — List A")
                if figA is not None:
                    st.plotly_chart(figA, use_container_width=True)
            except TypeError:
                pass

            st.subheader("Enrichment — List B")
            st.dataframe(dfB)
            try:
                figB = pipeline.plot_results(dfB, top_n=10, title="Enrichment — List B")
                if figB is not None:
                    st.plotly_chart(figB, use_container_width=True)
            except TypeError:
                pass

            # Overlap of enriched terms
            if not dfA.empty and not dfB.empty:
                overlap = set(dfA["Term_ID"]).intersection(dfB["Term_ID"])
                st.subheader("Shared Enriched Terms")
                st.write(list(overlap))
            else:
                st.info("One or both lists had no enriched terms.")

# ---------------------------------------------------------------------
# EV QC (MISEV-inspired)
# ---------------------------------------------------------------------
elif app_mode == "🧪 EV QC (MISEV Marker Check)":
    st.title("🧪 EV QC — MISEV Marker Screening")

    qc_text = st.text_area("Paste your EV protein/mRNA list")

    core_markers = {"CD9", "CD63", "CD81", "TSG101", "ALIX"}
    contaminants = {
        "Mitochondrial": {"MT-CO1", "MT-ND1", "CYCS"},
        "Nuclear": {"TP53", "BRCA1", "HIST1H1"},
        "Cytoskeletal": {"ACTB", "ACTG1", "TUBB"},
    }

    if st.button("Run QC Check"):
        genes = set(parse_gene_list(qc_text))

        if not genes:
            st.error("Please paste at least one gene.")
        else:
            st.subheader("Core EV Markers")
            found = genes.intersection(core_markers)
            if found:
                st.success(f"Detected markers: {', '.join(sorted(found))}")
            else:
                st.warning("No core EV markers detected.")

            st.subheader("Contamination Markers")
            for label, markers in contaminants.items():
                hits = genes.intersection(markers)
                if hits:
                    st.warning(f"{label}: {', '.join(sorted(hits))}")
                else:
                    st.info(f"{label}: None detected")

# ---------------------------------------------------------------------
# Build Custom Gene Sets
# ---------------------------------------------------------------------
elif app_mode == "🔧 Build Custom Gene Sets":
    st.title("🔧 Add Custom Gene Sets")

    term_name = st.text_input("Term name")
    term_genes = st.text_area("Genes (newline or comma separated)")

    if st.button("Add Term"):
        if not term_name.strip():
            st.error("Please provide a term name.")
        else:
            genes = parse_gene_list(term_genes)
            if not genes:
                st.error("Please provide at least one gene for the term.")
            else:
                pipeline.load_custom_dict({term_name: genes})
                st.success(f"Added custom gene set '{term_name}' (n={len(genes)})")
                st.json({term_name: genes})

# ---------------------------------------------------------------------
# Explore Metadata
# ---------------------------------------------------------------------
elif app_mode == "🗂 Explore Metadata":
    st.title("🗂 Vesiclepedia Metadata Explorer")

    if hasattr(pipeline, "merged_dataset") and pipeline.merged_dataset is not None:
        df = pipeline.merged_dataset
        if df.empty:
            st.warning("Merged metadata is empty. Try re-running the conversion step.")
        else:
            st.subheader("Raw Metadata Table")
            st.dataframe(df)

            st.subheader("Summary Statistics")
            st.write(df.describe(include="all"))

            # Use only object-type (categorical/string) columns for group by
            cat_cols = [c for c in df.columns if df[c].dtype == "object"]
            if cat_cols:
                category = st.selectbox("Group By", cat_cols)
                grouped = df.groupby(category).size()
                st.subheader(f"Counts by {category}")
                st.bar_chart(grouped)
            else:
                st.info("No categorical columns available for grouping.")
    else:
        st.info(
            "Metadata becomes available after converting Vesiclepedia raw files "
            "using the 'Upload & Convert Vesiclepedia Files' module."
        )

# ---------------------------------------------------------------------
# About
# ---------------------------------------------------------------------
elif app_mode == "📜 About":
    st.title("📜 About This App")
    st.markdown(
        """
        **EVision** is a demonstration and research tool built on top of the
        `VesiclePediaPipeline` for:

        - Parsing Vesiclepedia experimental data
        - Constructing EV-specific gene sets (GMT)
        - Running enrichment analysis for EV-related gene lists
        - Performing basic QC checks inspired by MISEV guidelines
        - Exploring metadata derived from Vesiclepedia experiments

        The PRIDE section in the sidebar is reserved for future integration
        with PRIDE proteomics datasets to extend EV analysis using independent
        experimental sources.
        """
    )