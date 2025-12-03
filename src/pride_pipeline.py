import os
import requests
import subprocess
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pyteomics import mztab


class PrideEVPipeline:
    def __init__(self, output_dir="./pride_data"):
        self.output_dir = output_dir
        self.api_url = "https://www.ebi.ac.uk/pride/ws/archive/v2"
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        self.MISEV_POSITIVE_MARKERS = [
            "CD9",
            "CD63",
            "CD81",
            "TSG101",
            "ALIX",
            "PDCD6IP",
        ]  # Tetraspanins + Cytosolic EV markers

        self.MISEV_NEGATIVE_MARKERS = [
            "CANX",
            "CALR",
            "LAMP1",
            "GM130",
            "HSPA5",
            "ALB",
            "ACTB",
        ]

    def search_datasets(self, disease_terms, ev_terms):
        """
        Phase 1: Discovery
        Queries PRIDE API for datasets matching Disease AND EV terms.
        """
        print("=== Searching PRIDE with multi-keyword intersection (Disease AND EV) ===")

        # --------------------------------------------------------------
        # Helper function (inline) to handle a *single keyword* search
        # --------------------------------------------------------------
        def single_keyword_search(keyword):
            params = {
                "keyword": keyword,
                "pageSize": 100,
                "sortDirection": "DESC",
                "sortFields": "submissionDate"
            }

            try:
                response = requests.get(f"{self.api_url}/search/projects", params=params)
                response.raise_for_status()
                results = response.json()

                # PRIDE returns inconsistent formats:
                # HAL JSON → dict
                # Simple fallback → list
                if isinstance(results, dict):
                    projects = results.get("_embedded", {}).get("compactprojects", [])
                    # Some endpoints also return {"list": [...]}
                    if not projects and "list" in results:
                        projects = results["list"]
                elif isinstance(results, list):
                    projects = results
                else:
                    print(f"Unexpected response type for keyword '{keyword}':", type(results))
                    return pd.DataFrame()

                df = pd.DataFrame(projects)
                return df

            except Exception as e:
                print(f"Search failed for keyword '{keyword}': {e}")
                return pd.DataFrame()

        # --------------------------------------------------------------
        # 1. Search disease keywords (OR within disease group)
        # --------------------------------------------------------------
        disease_results = []
        for term in disease_terms:
            df = single_keyword_search(term)
            if not df.empty:
                disease_results.append(df)

        if not disease_results:
            print("No hits found for disease-related terms.")
            return pd.DataFrame()

        disease_pxds = set().union(*[set(df["accession"]) for df in disease_results])

        # --------------------------------------------------------------
        # 2. Search EV keywords (OR within EV group)
        # --------------------------------------------------------------
        ev_results = []
        for term in ev_terms:
            df = single_keyword_search(term)
            if not df.empty:
                ev_results.append(df)

        if not ev_results:
            print("No hits found for EV-related terms.")
            return pd.DataFrame()

        ev_pxds = set().union(*[set(df["accession"]) for df in ev_results])

        # --------------------------------------------------------------
        # 3. Boolean AND across groups: (Disease) AND (EV)
        # --------------------------------------------------------------
        final_pxds = disease_pxds.intersection(ev_pxds)

        if not final_pxds:
            print("No datasets match BOTH disease terms AND EV terms.")
            return pd.DataFrame()

        print(f"Found {len(final_pxds)} datasets matching BOTH conditions.")

        # Merge all results into one dataframe
        combined_df = pd.concat(disease_results + ev_results, ignore_index=True)
        clean_cols = ["accession", "title", "submissionDate", "projectDescription"]
        final_df = (
            combined_df[combined_df["accession"].isin(final_pxds)][clean_cols]
            .drop_duplicates(subset=["accession"])
        )

        # Sort by submission date when available
        if "submissionDate" in final_df.columns:
            final_df = final_df.sort_values(by="submissionDate", ascending=False)

        return final_df

    def download_mztab(self, pxd_accession):
        """
        Phase 2 (fixed): Retrieve RESULT files using the correct pridepy CLI.
        1. Query file list
        2. Filter for category == RESULT
        3. Download each RESULT file individually
        """

        target_dir = os.path.join(self.output_dir, pxd_accession)
        os.makedirs(target_dir, exist_ok=True)

        print(f"--- Listing files for {pxd_accession} using pridepy ---")

        # STEP 1 — Get file list
        list_cmd = [
            "pridepy", "file", "list",
            "--project", pxd_accession
        ]

        try:
            list_result = subprocess.run(
                list_cmd,
                capture_output=True,
                text=True,
                check=True
            )
        except Exception as e:
            print(f"Error running pridepy file list: {e}")
            return []

        # pridepy returns JSON on stdout
        try:
            file_list = json.loads(list_result.stdout)
        except Exception as e:
            print("Could not decode pridepy file list output:", e)
            print("Raw output:", list_result.stdout)
            return []

        # STEP 2 — Filter for RESULT category
        result_files = [
            f for f in file_list
            if f.get("category", "").upper() == "RESULT"
        ]

        if not result_files:
            print("No RESULT files found for this project.")
            return []

        print(f"Found {len(result_files)} RESULT files. Downloading...")

        downloaded_files = []

        # STEP 3 — Download each RESULT file
        for f in result_files:
            file_id = f["fileId"]
            filename = f["fileName"]

            print(f"  Downloading {filename} ...")

            dl_cmd = [
                "pridepy", "file", "download",
                "--project", pxd_accession,
                "--fileId", str(file_id),
                "--output", target_dir
            ]

            try:
                subprocess.run(dl_cmd, check=True)
                downloaded_files.append(os.path.join(target_dir, filename))
            except Exception as e:
                print(f"Failed to download {filename}: {e}")

        if not downloaded_files:
            print("Finished, but no files were downloaded.")
        else:
            print(f"Download complete. {len(downloaded_files)} files saved.")

        return downloaded_files


    def load_and_clean_data(self, mztab_file):
        """
        Phase 3: Cleaning & Standardization
        Parses mzTab, extracts abundance data, and handles missing values.
        """
        print(f"--- Parsing {os.path.basename(mztab_file)} ---")

        try:
            # Pyteomics mzTab parser
            tables = mztab.MzTab(mztab_file)

            # Extract Protein Quantification Table
            prot_df = tables.protein_table

            if prot_df is None or prot_df.empty:
                print("No protein table found in mzTab.")
                return None

            # 1. Clean Column Names: Identify abundance columns
            # mzTab standard uses 'protein_abundance_assay[1]', etc.
            abundance_cols = [c for c in prot_df.columns if "abundance" in c]
            meta_cols = ["accession", "description", "gene_name"]

            # Keep only relevant columns if they exist
            cols_to_keep = [c for c in meta_cols if c in prot_df.columns] + abundance_cols
            clean_df = prot_df[cols_to_keep].copy()

            # 2. Filter Contaminants (Keratin, Trypsin)
            # Standard PRIDE contaminant flag is often "CONT_" or specific accessions
            if "accession" in clean_df.columns:
                clean_df = clean_df[~clean_df["accession"].str.contains("CONT_", na=False)]

            # 3. Handle Missing Values
            if abundance_cols:
                clean_df.dropna(subset=abundance_cols, how="all", inplace=True)
            else:
                print("Warning: No abundance columns found. Returning identification data only.")

            return clean_df

        except Exception as e:
            print(f"Parsing failed: {e}")
            return None

    def perform_ev_eda(self, df, pxd_id):
        """
        Phase 4: Exploratory Data Analysis (EDA) & MISEV Quality Check
        """
        if df is None or df.empty:
            print("No data available for EDA.")
            return

        print(f"--- Performing EDA on {pxd_id} ---")

        # 4.1 Define MISEV (Minimal Information for Studies of Extracellular Vesicles) Markers
        # Based on MISEV2018 guidelines (example subset)
        positive_markers = [
            "CD9",
            "CD63",
            "CD81",
            "TSG101",
            "ALIX",
            "PDCD6IP",
        ]  # Tetraspanins + Cytosolic EV markers

        negative_markers = [
            "CANX",
            "CALR",
            "LAMP1",
            "GM130",
            "HSPA5",
            "ALB",
            "ACTB",
        ]  # ER, Golgi, Nucleus/other contaminants

        # 4.2 Score the Dataset
        # Check if markers are present in 'gene_name'; if not, fall back to 'description'
        if "gene_name" in df.columns:
            df["is_EV_marker"] = df["gene_name"].apply(
                lambda x: x in positive_markers if pd.notnull(x) else False
            )
            df["is_contaminant"] = df["gene_name"].apply(
                lambda x: x in negative_markers if pd.notnull(x) else False
            )
        elif "description" in df.columns:
            df["is_EV_marker"] = df["description"].apply(
                lambda x: any(m in str(x) for m in positive_markers) if pd.notnull(x) else False
            )
            df["is_contaminant"] = df["description"].apply(
                lambda x: any(m in str(x) for m in negative_markers) if pd.notnull(x) else False
            )
        else:
            # If neither gene_name nor description exist, create flags but they will be all False
            df["is_EV_marker"] = False
            df["is_contaminant"] = False
            print("Warning: no 'gene_name' or 'description' columns to check markers.")

        ev_hits = df[df["is_EV_marker"]]
        contam_hits = df[df["is_contaminant"]]

        print("MISEV Quality Check:")
        if "gene_name" in df.columns:
            print(f"  - Positive Markers Found: {ev_hits['gene_name'].dropna().unique()}")
            print(f"  - Negative Markers Found: {contam_hits['gene_name'].dropna().unique()}")
        else:
            print(f"  - Positive Marker Rows: {len(ev_hits)}")
            print(f"  - Negative Marker Rows: {len(contam_hits)}")

        # 4.3 Visualization: Abundance Distribution
        abundance_cols = [c for c in df.columns if "abundance" in c]
        if abundance_cols:
            melted = df.melt(value_vars=abundance_cols, var_name="Sample", value_name="Abundance")

            plt.figure(figsize=(10, 6))
            sns.boxplot(x="Sample", y="Abundance", data=melted)
            plt.yscale("log")
            plt.title(f"Protein Abundance Distribution - {pxd_id}")
            plt.xticks(rotation=45, ha="right")
            plt.tight_layout()
            plt.show()
        else:
            print("No quantitative columns found for plotting.")

    def assess_misev_markers(self, df):
        """
        Input: Clean protein table from mzTab
        Output: Dict containing marker presence and purity score
        """

        if df is None or df.empty:
            return {
                "positive_found": [],
                "negative_found": [],
                "ev_purity_score": None
            }

        # identify gene name column
        name_col = "gene_name" if "gene_name" in df.columns else "description"

        df[name_col] = df[name_col].astype(str)

        # positive markers
        positive_hits = [
            m for m in self.MISEV_POSITIVE_MARKERS
            if df[name_col].str.contains(m, case=False, na=False).any()
        ]

        # negative markers
        negative_hits = [
            m for m in self.MISEV_NEGATIVE_MARKERS
            if df[name_col].str.contains(m, case=False, na=False).any()
        ]

        purity_score = len(positive_hits) - len(negative_hits)

        return {
            "positive_found": positive_hits,
            "negative_found": negative_hits,
            "ev_purity_score": purity_score
        }
    

    def build_misev_scorecard(self, pipeline, pxds):
        score_rows = []

        for pxd in pxds:
            files = pipeline.download_mztab(pxd)
            if not files:
                continue

            df = pipeline.load_and_clean_data(files[0])
            q = self.assess_misev_markers(df)

            score_rows.append({
                "pxd": pxd,
                "pos_markers": len(q["positive_found"]),
                "neg_markers": len(q["negative_found"]),
                "ev_purity_score": q["ev_purity_score"],
                "positive_found_list": q["positive_found"],
                "negative_found_list": q["negative_found"],
            })

        return pd.DataFrame(score_rows)
    
    def plot_marker_heatmap(self, score_df):
        marker_set = self.MISEV_POSITIVE_MARKERS + self.MISEV_NEGATIVE_MARKERS
        mat = []

        for _, row in score_df.iterrows():
            presence = [
                1 if m in row["positive_found_list"] or m in row["negative_found_list"] else 0
                for m in marker_set
            ]
            mat.append(presence)

        marker_df = pd.DataFrame(mat, columns=marker_set, index=score_df["pxd"])

        plt.figure(figsize=(10,6))
        sns.heatmap(marker_df, cmap="YlGnBu", annot=True)
        plt.title("MISEV Marker Presence Across Datasets")
        plt.show()

    def plot_purity_ranking(self, score_df):
        sorted_df = score_df.sort_values("ev_purity_score", ascending=False)

        plt.figure(figsize=(10,4))
        sns.barplot(
            x="pxd",
            y="ev_purity_score",
            data=sorted_df,
            palette="viridis"
        )
        plt.title("EV Purity Score Ranking")
        plt.xticks(rotation=45)
        plt.show()


# --- Execution ---

if __name__ == "__main__":
    pipeline = PrideEVPipeline()

    # 1. Discovery
    # Example disease terms for AD – adjust to your needs
    datasets = pipeline.search_datasets(
        disease_terms=["Alzheimer’s", "PXD015578"],
        ev_terms = [
            "extracellular vesicle",
            "extracellular vesicles",
            "exosome",
            "exosomes",
            "microvesicle",
            "microvesicles",
            "exosomal"
        ]
    )

    if not datasets.empty:
        print("\nTop 3 Discovered Datasets:")
        print(datasets[["accession", "title", "submissionDate"]].head(3))

        # 2. Select a dataset for the pipeline
        # PXD015578 is a known AD EV dataset (Guo et al.), if present
        target_pxd = "PXD015578"

        if target_pxd not in datasets["accession"].values:
            first_accession = datasets.iloc[0]["accession"]
            print(
                f"\nTarget {target_pxd} not in search results, "
                f"defaulting to first result: {first_accession}"
            )
            target_pxd = first_accession

        # 3. Download
        downloaded_files = pipeline.download_mztab(target_pxd)

        # 4. Process & EDA
        if downloaded_files:
            # Process the first mzTab found
            first_mztab = downloaded_files[0]
            df = pipeline.load_and_clean_data(first_mztab)
            pipeline.perform_ev_eda(df, target_pxd)
        else:
            print("No mzTab files found. The dataset might be RAW-only or use older formats.")
    else:
        print("No datasets discovered with the given search terms.")

 

