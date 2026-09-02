#!/usr/bin/env python3
import os
import sys
import json
import urllib.request
import xml.etree.ElementTree as ET
import time
from collections import Counter

ABUNDANCE_FILE = "output_AAy/abundance_table_species.tsv"
ASSIGNMENTS_FILE = "output_AAy/reads_assignments/AmbarAyunda_minimap2_16S_lineages.minimap2.assignments.tsv"
CACHE_FILE = "output_AAy/taxonomy_cache.json"

def main():
    print("[taxonomy] Starting NCBI taxonomy resolution...")
    
    # 1. Read abundance table to get all unique paths
    if not os.path.exists(ABUNDANCE_FILE):
        print(f"[taxonomy] Error: Abundance file not found: {ABUNDANCE_FILE}")
        sys.exit(1)
        
    ab_paths = []
    with open(ABUNDANCE_FILE, "r") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            if parts:
                ab_paths.append(parts[0])

    print(f"[taxonomy] Found {len(ab_paths)} lineages in abundance table.")
    
    # 2. Read assignments to map species names to TaxIDs
    if not os.path.exists(ASSIGNMENTS_FILE):
        print(f"[taxonomy] Error: Assignments file not found: {ASSIGNMENTS_FILE}")
        sys.exit(1)
        
    print("[taxonomy] Parsing assignments file for species-to-taxid mapping...")
    species_to_taxids = {}
    with open(ASSIGNMENTS_FILE, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                taxid = parts[2].strip()
                lineage = parts[4].strip()
                if taxid and taxid != "0" and lineage:
                    # Leaf species name is the last part
                    leaf_name = lineage.split("|")[-1].strip()
                    if leaf_name:
                        if leaf_name not in species_to_taxids:
                            species_to_taxids[leaf_name] = []
                        species_to_taxids[leaf_name].append(taxid)
                        
    # For species with multiple taxids, resolve to the most common one
    resolved_species_taxid = {}
    for species, tids in species_to_taxids.items():
        counter = Counter(tids)
        resolved_species_taxid[species] = counter.most_common(1)[0][0]
        
    print(f"[taxonomy] Mapped {len(resolved_species_taxid)} unique species to NCBI TaxIDs from assignments.")

    # 3. Identify paths that need taxids
    needed_taxids = set()
    path_to_species = {}
    
    for path in ab_paths:
        if path == "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown":
            continue
            
        parts = path.split(";")
        if len(parts) >= 8:
            species_name = parts[7].strip()
            path_to_species[path] = species_name
            taxid = resolved_species_taxid.get(species_name)
            if taxid:
                needed_taxids.add(taxid)
            else:
                # Fallback: check clean name
                clean_name = species_name.replace("[", "").replace("]", "")
                taxid = resolved_species_taxid.get(clean_name)
                if taxid:
                    needed_taxids.add(taxid)
                    path_to_species[path] = clean_name

    print(f"[taxonomy] Need to fetch NCBI taxonomy for {len(needed_taxids)} unique species taxids.")

    # 4. Fetch taxonomy from NCBI in batches of 300
    taxid_list = list(needed_taxids)
    batch_size = 300
    new_mappings = {}
    
    for i in range(0, len(taxid_list), batch_size):
        batch = taxid_list[i:i+batch_size]
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={','.join(batch)}&retmode=xml"
        print(f"[taxonomy] Querying NCBI batch {i//batch_size + 1}/{-(-len(taxid_list)//batch_size)} ({len(batch)} IDs)...")
        
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req) as response:
                xml_data = response.read()
                
            root = ET.fromstring(xml_data)
            for taxon in root.findall('Taxon'):
                taxid = int(taxon.find('TaxId').text)
                name = taxon.find('ScientificName').text
                
                # Retrieve structured lineage
                lineage_ex = taxon.find('LineageEx')
                full_lineage = []
                if lineage_ex is not None:
                    for parent in lineage_ex.findall('Taxon'):
                        p_taxid = int(parent.find('TaxId').text)
                        p_name = parent.find('ScientificName').text
                        # Skip root/cellular organisms
                        if p_name in ["cellular organisms", "root"]:
                            continue
                        full_lineage.append((p_name, p_taxid))
                            
                # Add the taxon itself to the lineage
                full_lineage.append((name, taxid))
                
                # Build path prefixes and map them to taxids
                path_prefix = []
                for idx, (p_name, p_taxid) in enumerate(full_lineage):
                    path_prefix.append(p_name)
                    constructed_path = ";".join(path_prefix)
                    new_mappings[constructed_path] = p_taxid
                    new_mappings[f"name:{p_name}"] = p_taxid
                    
            # Be polite to NCBI API
            time.sleep(1.0)
            
        except Exception as e:
            print(f"[taxonomy] Error fetching batch: {e}")
            time.sleep(2.0)

    # 5. Resolve abundance table paths to taxids
    cache = {}
    for path in ab_paths:
        if path == "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown":
            cache[path] = 0
            continue
            
        parts = path.split(";")
        current_prefix = []
        for j in range(len(parts)):
            node_name = parts[j]
            current_prefix.append(node_name)
            prefix_path = ";".join(current_prefix)
            
            if prefix_path in cache:
                continue
                
            taxid = None
            
            # Check direct prefix path match
            if prefix_path in new_mappings:
                taxid = new_mappings[prefix_path]
            # Try 7-rank path match (by removing the second rank "Bacillati")
            elif j >= 1:
                prefix_7 = ";".join([parts[0]] + parts[2:j+1])
                if prefix_7 in new_mappings:
                    taxid = new_mappings[prefix_7]
                    
            # Fallback to simple name lookup
            if not taxid:
                taxid = new_mappings.get(f"name:{node_name}")
                
            # Final fallback: default to 0 if not found
            cache[prefix_path] = taxid if taxid is not None else 0

    # 6. Write cache file
    try:
        with open(CACHE_FILE, "w") as f:
            json.dump(cache, f, indent=2)
        print(f"[taxonomy] Saved {len(cache)} path-to-taxid mappings to cache.")
    except Exception as e:
        print(f"[taxonomy] Error saving cache: {e}")

if __name__ == "__main__":
    main()
