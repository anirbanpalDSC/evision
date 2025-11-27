import requests
import pandas as pd
import re
import gzip
import io
import os

# --- Configuration ---
PRIDE_API_PROJECTS = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects"
PRIDE_API_FILES = "https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject"

# Regex for UniProt Accession Numbers (Standard Pattern)
# Matches: P12345, Q9XYZ1, A0A024RBG1, etc.
UNIPROT_PATTERN = re.compile(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')

def search_pride_projects(keyword, page_size=5):
    """
    Searches PRIDE for projects matching a keyword.
    Returns a list of PXD Accessions.
    """
    print(f"Searching PRIDE for: '{keyword}'...")
    params = {
        "query": keyword,
        "pageSize": page_size,
        "sortDirection": "DESC",
        "sortConditions": "submission_date"
    }
    
    try:
        response = requests.get(PRIDE_API_PROJECTS, params=params)
        response.raise_for_status()
        data = response.json()
        
        projects = []
        if "_embedded" in data and "projects" in data["_embedded"]:
            for proj in data["_embedded"]["projects"]:
                accession = proj["accession"]
                title = proj["title"]
                print(f"Found: {accession} - {title[:60]}...")
                projects.append((accession, title))
        
        return projects
    except Exception as e:
        print(f"API Error: {e}")
        return []

def get_result_files(pxd_accession):
    """
    Fetches file list for a project and filters for RESULT/SEARCH files.
    We want to avoid RAW files (huge).
    """
    url = f"{PRIDE_API_FILES}?accession={pxd_accession}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        files = response.json()
        
        # Filter for interesting small files
        # 'RESULT' often contains mzIdentML (.mzid)
        # 'OTHER' often contains summary .txt/.xlsx tables provided by authors
        target_files = []
        for f in files:
            fname = f["fileName"].lower()
            ftype = f["fileCategory"]["value"]
            fsize_mb = f["fileSizeBytes"] / (1024 * 1024)
            
            # Logic: look for result files, mzid, or small text tables
            # Skip RAW files immediately
            if ftype == 'RAW':
                continue
                
            if (ftype == 'RESULT' or ftype == 'SEARCH' or 'peptide' in fname or 'protein' in fname) and fsize_mb < 500:
                target_files.append({
                    'name': f["fileName"],
                    'url': f["publicFileLocations"][0]["value"], # usually FTP/HTTP
                    'size_mb': fsize_mb
                })
                
        return target_files
    except Exception as e:
        print(f"Error fetching files for {pxd_accession}: {e}")
        return []

def extract_uniprot_ids_from_url(file_url):
    """
    Streams a file from a URL and uses Regex to find UniProt IDs.
    This works on .xml, .mzid, .txt, .csv without needing complex parsers.
    """
    print(f"  Streaming & Scanning: {file_url} ...")
    found_proteins = set()
    
    try:
        # Stream the download
        with requests.get(file_url, stream=True) as r:
            r.raise_for_status()
            
            # Handle GZIP if necessary
            if file_url.endswith('.gz'):
                f_stream = gzip.open(r.raw, mode='rt', encoding='latin-1')
            else:
                r.raw.decode_content = True
                f_stream = io.TextIOWrapper(r.raw, encoding='latin-1')
            
            # Read line by line to keep memory low
            for line in f_stream:
                matches = UNIPROT_PATTERN.findall(line)
                found_proteins.update(matches)
                
    except Exception as e:
        print(f"  Error processing file: {e}")
        
    return found_proteins

def uniprot_to_gene_symbol(uniprot_ids):
    """
    Simple batch converter using MyGene.info API.
    Needed because your Vesiclepedia GMT uses Gene Symbols.
    """
    print(f"  Converting {len(uniprot_ids)} UniProt IDs to Gene Symbols...")
    
    # Split into chunks of 1000 for API
    uniprot_list = list(uniprot_ids)
    chunk_size = 1000
    mapping = {}
    
    for i in range(0, len(uniprot_list), chunk_size):
        chunk = uniprot_list[i:i+chunk_size]
        try:
            res = requests.post(
                'https://mygene.info/v3/query',
                json={
                    "q": chunk,
                    "scopes": "uniprot",
                    "fields": "symbol",
                    "species": "human" # Defaulting to human for this example
                }
            )
            for item in res.json():
                if 'symbol' in item and 'query' in item:
                    mapping[item['query']] = item['symbol']
        except Exception as e:
            print(f"  Mapping error: {e}")
            
    return set(mapping.values())

def run_pride_pipeline(keyword="exosomes", output_gmt="pride_db.gmt"):
    """
    Main Orchestrator
    """
    # 1. Find Projects
    projects = search_pride_projects(keyword)
    
    gmt_lines = []
    
    # 2. Process each project (Limiting to first 2 for demo)
    for pxd, title in projects[:2]:
        print(f"\nProcessing Project: {pxd}")
        
        # 3. Find suitable result files
        files = get_result_files(pxd)
        
        if not files:
            print("  No suitable result files found.")
            continue
            
        # Heuristic: Pick the file with 'protein' in name, or just the largest result file
        # that isn't huge. For now, we take the first one that looks promising.
        target_file = files[0] 
        for f in files:
            if 'protein' in f['name'].lower() and f['name'].endswith('.txt'):
                target_file = f
                break
        
        print(f"  Selected file: {target_file['name']} ({target_file['size_mb']:.2f} MB)")
        
        # 4. Extract UniProt IDs
        uniprot_ids = extract_uniprot_ids_from_url(target_file['url'])
        print(f"  Found {len(uniprot_ids)} unique UniProt IDs.")
        
        if len(uniprot_ids) < 5:
            print("  Too few proteins found, skipping.")
            continue
            
        # 5. Convert to Gene Symbols
        gene_symbols = uniprot_to_gene_symbol(uniprot_ids)
        print(f"  Mapped to {len(gene_symbols)} Gene Symbols.")
        
        if not gene_symbols:
            continue
            
        # 6. Format for GMT
        # ID: PRIDE:PXDxxxxxx
        term_id = f"PRIDE:{pxd}"
        description = f"{title} ({target_file['name']})"
        
        # Remove empty strings
        valid_genes = [g for g in gene_symbols if g]
        
        line = f"{term_id}\t{description}\t" + "\t".join(valid_genes)
        gmt_lines.append(line)
        print(f"  Added {pxd} to database.")

    # 7. Save
    if gmt_lines:
        with open(output_gmt, 'w', encoding='utf-8') as f:
            f.write("\n".join(gmt_lines))
        print(f"\nSuccess! Created {output_gmt} with {len(gmt_lines)} datasets.")
    else:
        print("\nNo data extracted.")

if __name__ == "__main__":
    run_pride_pipeline(keyword="exosomes")