<p align="left">
  <img src="images/evision.png" width="200" alt="EVision Logo">
</p>

# 📦 EVision: An Automated Quality Control and Enrichment Engine for EV Omics

*A lightweight, modular, open-source pipeline for processing, analyzing, and visualizing extracellular vesicle (EV) cargo datasets.*

## 🌟 Overview

EVision is a Python and Streamlit based pipeline designed to bring reproducibility, rigorous quality control, dynamic statistical analysis to extracellular vesicle (EV) research.

## 🔬 Current Features

### 🧬 Dynamic ETL Pipeline

* **Schema-Agnostic Ingestion:** Instantly map heterogeneous CSV/Excel/Text files into standardized Gene Matrix Transposed (GMT) formats, or load predefined GMT files for analysis and comparison.

* **Species-Specific Filtering:** Optionally, segregate data by organism (e.g., $Homo sapiens$ vs. $Mus musculus$) to eliminate cross-species noise.

* **Live Updates:** Regenerate your background "gene universe" instantly as new public data becomes available; no more waiting for software patches.

### 🛡️ Automated MISEV Quality Control

* **The "Gatekeeper" Module:** Validates purity of the sample before analysis. Uses a MISEV compliant curated dictionary of positive (e.g., CD63, TSG101) and negative (e.g., Albumin, Histones) markers.
  
* **Purity Scoring:** Calculates a quantitative "Purity Score" ($0.0 - 1.0$) for every dataset, flagging contaminated samples that could skew enrichment results.

### 📊 Robust Statistical Inference

* **Engine:** Performs Over Representation Analysis (ORA) using Fisher’s Exact Test.

* **Correction:** Applies Benjamini–Hochberg False Discovery Rate (FDR) correction by default to control false positives in high-dimensional omics data.

* **Context-Aware:** Calculates enrichment against a rigorous, domain-specific and customizable background (the EV proteome) rather than the entire genome.

### 🕸️ Network Visualization

* **Enrichment Maps:** Moves beyond static bar charts by projecting results as force-directed network graphs.

* **Redundancy Reduction:** Uses Jaccard Similarity to cluster related terms (e.g., "Angiogenesis" and "Blood Vessel Development"), revealing high-level biological themes and pathway crosstalk.

## ⚙️ How It Works

* **Extract:** Upload your raw protein/gene lists or public repository exports.

* **Transform:** EVision maps metadata, filters by species, and applies QC metrics.

* **Load:** The cleaned data defines a custom "Gene Universe" ($N$) for the session.

* **Analyze:** Input your query list to identify statistically over-represented biological functions.
  
* **Visualize:** Explore relationships between pathways in an interactive network graph.

## 🧬 Future Roadmap



## 🧠 Why This Pipeline Matters

Extracellular vesicles carry molecular signatures that reflect disease states, cell-type identity, and intercellular communication pathways.
However:

* Open databases, e.g. **Vesiclepedia** provides curated cargo, not raw or standardized datasets.

~* PRIDE provides raw datasets, but they are large and hard to parse.

EVision bridges this gap by providing a fast, reproducible, scriptable pipeline for building EV gene sets and enabling downstream enrichment, comparison, and publication-ready analysis — without requiring heavy proteomics tools.
~
## 🤝 Contributions

Pull requests, feature requests, and issue reports are welcome!
This project is evolving fast, and community feedback is appreciated.

## 📄 License

This project is licensed under the Creative Commons Attribution 4.0
International (CC BY 4.0) License.

© 2025 Anirban Pal

~PRIDE Notes:
https://www.ebi.ac.uk/pride/markdownpage/prideapi~

## 📜 How to Cite

If you use EVision in academic research, publications, or derivative tools,
please cite:

> Pal, A. (2025). *EVision: An Automated Quality Control and Enrichment Engine for EV Omics*.
> GitHub repository. https://github.com/anirbanpalDSC/evision/
