import streamlit as st
import pandas as pd
import numpy as np
import json
import tempfile
from io import StringIO
from src.vp_pipeline import VesiclePediaPipeline

# ---------------------------------------------------------------------
# App Config
# ---------------------------------------------------------------------
st.set_page_config(
    page_title="EVision: Extracellular Vesicles Data Analysis Suite",
    layout="wide"
)

logo_path = "images/evision.png"
st.logo(
    image=logo_path,
    link="https://github.com/anirbanpalDSC/evision",
    # Adjust the size (small, medium, large)
    size="large"
)

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

# Main mode is driven entirely by Vesiclepedia options for now
app_mode = vesicle_option

# ---------------------------------------------------------------------
# 1 — Home
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
# 2 — Upload & Convert Vesiclepedia Raw Files
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

                output_gmt = "custom_db.gmt"
                
                # Create the mapping dict
                mapping_config = {
                    'id_col_exp': id_col_exp,
                    'id_col_data': id_col_data,
                    'gene_col': gene_col,
                    'categories': categories
                }

                # Run the new generic pipeline function
                status = pipeline.convert_to_gmt(
                    experiment_file=exp_path,
                    data_file=data_path,
                    mapping=mapping_config,
                    output_gmt=output_gmt,
                    sep=sep_char
                )

                if "Success" in status:
                    st.success(status)
                    with open(output_gmt, "rb") as f_out:
                        st.download_button(
                            "⬇ Download GMT", f_out, "custom_db.gmt"
                        )
                else:
                    st.error(status)

# ---------------------------------------------------------------------
# 3 — Load GMT Database
# ---------------------------------------------------------------------
elif app_mode == "📚 Load GMT Database":
    st.title("📚 Load Vesiclepedia or Custom GMT File")

    gmt_file = st.file_uploader("Upload GMT File", type=["gmt"])

    if gmt_file:
        path = "uploaded.gmt"
        with open(path, "w", encoding="utf-8") as f:
            f.write(gmt_file.getvalue().decode())

        pipeline.load_gmt(path)

        st.success(f"Loaded {len(pipeline.database)} gene sets.")
        st.markdown("**Example of first few sets (up to 5 sets, 5 genes each):**")
        preview = {
            k: list(v)[:5] for k, v in list(pipeline.database.items())[:5]
        }
        st.json(preview)

# ---------------------------------------------------------------------
# 4 — Enrichment Analysis
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# 4 — Enrichment Analysis (UPDATED WITH RIGOROUS BACKGROUND)
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
            but <i>not selected</i>. Usually, this is the full list of all proteins/genes 
            detected by your instrument.</small>
            """, 
            unsafe_allow_html=True
        )
        
        # Option A: Paste (Good for small lists)
        bg_paste = st.text_area("Option A: Paste Background", height=100)
        
        # Option B: Upload (Essential for full genomes/proteomes)
        bg_file = st.file_uploader("Option B: Upload Background File (.txt/.csv)", type=["txt", "csv"])

    # Settings
    with st.expander("⚙️ Advanced Statistics Settings"):
        correction_method = st.selectbox(
            "Multiple Testing Correction", 
            ['fdr_bh', 'bonferroni', 'holm', 'none'],
            index=0,
            help="FDR (Benjamini-Hochberg) is standard. Bonferroni is very strict."
        )
        top_n = st.number_input("Top N terms to plot", 5, 50, 15)

    run_btn = st.button("🚀 Run Analysis", type="primary")

    if run_btn:
        # 1. Parse Interest List
        genes = parse_gene_list(gene_text)
        
        # 2. Parse Background
        background = []
        
        # Priority to File Upload (handles 20k+ lines better)
        if bg_file is not None:
            # Read file as string
            stringio = StringIO(bg_file.getvalue().decode("utf-8"))
            background = parse_gene_list(stringio.read())
            st.success(f"Loaded Custom Universe: {len(background)} genes")
        elif bg_paste:
            background = parse_gene_list(bg_paste)
            st.info(f"Loaded Custom Universe (from text): {len(background)} genes")
        else:
            background = None # Pipeline will default to GMT universe
            st.warning("⚠️ No background provided. Using default GMT universe (Statistical rigor reduced).")

        if not genes:
            st.error("Please provide at least one gene in the gene list.")
        else:
            # Pass the correction method and background to the pipeline
            df = pipeline.run_enrichment(
                genes, 
                background_list=background, 
                correction_method=correction_method
            )

            if df.empty:
                st.warning("No enriched terms found (or p-values were not significant).")
            else:
                st.subheader("Results Table")
                # Format p-values for readability
                display_df = df.copy()
                display_df['p_value'] = display_df['p_value'].map('{:.2e}'.format)
                display_df['adj_p_value'] = display_df['adj_p_value'].map('{:.2e}'.format)
                st.dataframe(display_df, use_container_width=True)

                st.subheader("Visualization")
                try:
                    fig = pipeline.plot_results(df, top_n=top_n)
                    if fig is not None:
                        st.plotly_chart(fig, use_container_width=True)
                except TypeError:
                     st.info("Visuals handled internally.")

                # Download Logic
                csv = df.to_csv(index=False)
                st.download_button(
                    "⬇ Download Full Results",
                    csv,
                    "enrichment_results.csv",
                    "text/csv",
                )

# ---------------------------------------------------------------------
# 5 — Compare Two Gene Lists
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
# 6 — EV QC (MISEV-inspired)
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
# 7 — Build Custom Gene Sets
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
# 8 — Explore Metadata
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
# 9 — About
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
