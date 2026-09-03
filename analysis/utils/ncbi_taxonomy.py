#!/usr/bin/env python3
# =============================================================================
# NCBI Taxonomy Resolver & Cache Manager
# =============================================================================

import os
import sys
import json
import argparse
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
import tempfile
import hashlib
from collections import Counter

TOOL_NAME = "ont_wf16s_postprocess"

def normalize_abundance_path_to_7(path):
    """Normalize 8-rank abundance path (omits kingdom) to match 7-rank minimap2 lineage."""
    parts = path.split(";")
    if len(parts) >= 8:
        # omit rank index 1 (kingdom)
        return "|".join([parts[0]] + parts[2:8])
    return "|".join(parts)

def compute_sha256(filepath):
    if not filepath or not os.path.exists(filepath):
        return None
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description="Resolve NCBI TaxIDs from wf-16s data")
    parser.add_argument("--abundance", required=True, help="Path to abundance_table_species.tsv")
    parser.add_argument("--assignments", default=None, help="Path to read assignments TSV")
    parser.add_argument("--cache", required=True, help="Path to taxonomy_cache.json")
    parser.add_argument("--mode", choices=["cache_only", "refresh"], default="cache_only")
    parser.add_argument("--email", default=None, help="NCBI Contact Email")
    parser.add_argument("--api-key", default=None, help="NCBI API Key")
    parser.add_argument("--unresolved-policy", choices=["warn", "error"], default="warn")
    parser.add_argument("--unresolved-tsv", default=None, help="Path to write unresolved TaxIDs")
    parser.add_argument("--provenance", default=None, help="Path to write resolver provenance JSON")
    args = parser.parse_args()

    print(f"[taxonomy] Running in '{args.mode}' mode...")

    # 1. Load abundance table lineages
    if not os.path.exists(args.abundance):
        print(f"[taxonomy] ERROR: Abundance file not found: {args.abundance}", file=sys.stderr)
        sys.exit(1)

    ab_paths = []
    with open(args.abundance, "r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            if parts and parts[0]:
                ab_paths.append(parts[0])

    print(f"[taxonomy] Loaded {len(ab_paths)} lineages from abundance table.")

    # 2. Parse assignments for normalized lineage -> TaxID mapping
    lineage_to_taxids = {}
    if args.assignments and os.path.exists(args.assignments):
        print(f"[taxonomy] Parsing assignments file: {args.assignments}")
        with open(args.assignments, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 5:
                    taxid = parts[2].strip()
                    lineage = parts[4].strip()
                    if taxid and taxid != "0" and lineage:
                        if lineage not in lineage_to_taxids:
                            lineage_to_taxids[lineage] = []
                        lineage_to_taxids[lineage].append(taxid)

    # Resolve conflicts via majority vote
    resolved_lineage_taxid = {}
    conflict_records = []
    for lin, tids in lineage_to_taxids.items():
        counts = Counter(tids)
        most_common = counts.most_common()
        winner = most_common[0][0]
        resolved_lineage_taxid[lin] = int(winner)
        if len(most_common) > 1:
            conflict_records.append({
                "lineage": lin,
                "winner_taxid": int(winner),
                "counts": dict(counts)
            })

    if conflict_records:
        print(f"[taxonomy] Resolved {len(conflict_records)} multi-TaxID conflicts via majority voting.")

    # 3. Load existing cache
    cache = {}
    if os.path.exists(args.cache):
        try:
            with open(args.cache, "r", encoding="utf-8") as f:
                cache = json.load(f)
            print(f"[taxonomy] Loaded existing cache with {len(cache)} entries.")
        except Exception as e:
            print(f"[taxonomy] WARNING: Could not read cache: {e}", file=sys.stderr)

    # 4. In offline mode, identify unresolved nodes
    unresolved_nodes = []
    for path in ab_paths:
        if path.startswith("Unclassified"):
            continue
        parts = path.split(";")
        for j in range(1, len(parts) + 1):
            subpath = ";".join(parts[:j])
            current_taxid = cache.get(subpath, 0)
            if current_taxid == 0:
                # Try matching from assignments if leaf
                if j == 8:
                    norm_7 = normalize_abundance_path_to_7(subpath)
                    asgn_taxid = resolved_lineage_taxid.get(norm_7)
                    if asgn_taxid:
                        cache[subpath] = asgn_taxid
                        current_taxid = asgn_taxid
                
            if current_taxid == 0:
                unresolved_nodes.append({
                    "path": subpath,
                    "depth": j,
                    "name": parts[j - 1]
                })

    # De-duplicate unresolved nodes
    seen_unresolved = set()
    dedup_unresolved = []
    for u in unresolved_nodes:
        if u["path"] not in seen_unresolved:
            seen_unresolved.add(u["path"])
            dedup_unresolved.append(u)

    print(f"[taxonomy] Unresolved nodes in cache: {len(dedup_unresolved)}")

    # 5. Refresh mode (if requested and authorized)
    if args.mode == "refresh" and dedup_unresolved:
        email = args.email or os.environ.get("NCBI_EMAIL")
        if not email:
            print("[taxonomy] ERROR: NCBI_EMAIL environment variable or --email is required for refresh mode.", file=sys.stderr)
            sys.exit(1)
        api_key = args.api_key or os.environ.get("NCBI_API_KEY")

        # Collect unique species names to query
        needed_names = list({u["name"] for u in dedup_unresolved if u["depth"] == 8})
        print(f"[taxonomy] Querying NCBI for {len(needed_names)} species...")

        # Batch querying NCBI esearch / efetch with bounded retries
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        delay = 0.35 if api_key else 1.0

        for sp_name in needed_names:
            params = {
                "db": "taxonomy",
                "term": sp_name,
                "retmode": "json",
                "tool": TOOL_NAME,
                "email": email
            }
            if api_key:
                params["api_key"] = api_key

            url = f"{base_url}?{urllib.parse.urlencode(params)}"
            retries = 3
            found_taxid = 0
            while retries > 0:
                try:
                    req = urllib.request.Request(url, headers={"User-Agent": f"{TOOL_NAME}/1.0"})
                    with urllib.request.urlopen(req, timeout=10) as resp:
                        res_json = json.loads(resp.read().decode("utf-8"))
                        id_list = res_json.get("esearchresult", {}).get("idlist", [])
                        if id_list:
                            found_taxid = int(id_list[0])
                    break
                except Exception as ex:
                    retries -= 1
                    time.sleep(2.0)

            if found_taxid > 0:
                for path in ab_paths:
                    if path.endswith(f";{sp_name}"):
                        cache[path] = found_taxid

            time.sleep(delay)

        # Atomic write back to cache
        cache_dir = os.path.dirname(os.path.abspath(args.cache))
        with tempfile.NamedTemporaryFile("w", dir=cache_dir, delete=False, encoding="utf-8") as tmp_f:
            json.dump(cache, tmp_f, indent=2)
            tmp_path = tmp_f.name
        os.replace(tmp_path, args.cache)
        print(f"[taxonomy] Successfully updated cache atomically: {args.cache}")

    # 6. Write unresolved TSV
    if args.unresolved_tsv:
        os.makedirs(os.path.dirname(os.path.abspath(args.unresolved_tsv)), exist_ok=True)
        with open(args.unresolved_tsv, "w", encoding="utf-8") as f:
            f.write("Depth\tNodeName\tTaxonPath\n")
            for u in dedup_unresolved:
                f.write(f"{u['depth']}\t{u['name']}\t{u['path']}\n")
        print(f"[taxonomy] Wrote unresolved nodes to: {args.unresolved_tsv}")

    # 7. Write provenance sidecar
    if args.provenance:
        os.makedirs(os.path.dirname(os.path.abspath(args.provenance)), exist_ok=True)
        prov = {
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "mode": args.mode,
            "tool": TOOL_NAME,
            "abundance_sha256": compute_sha256(args.abundance),
            "assignments_sha256": compute_sha256(args.assignments) if args.assignments else None,
            "cache_sha256": compute_sha256(args.cache),
            "total_lineages": len(ab_paths),
            "unresolved_count": len(dedup_unresolved),
            "conflicts_count": len(conflict_records),
            "conflicts": conflict_records[:50]
        }
        with open(args.provenance, "w", encoding="utf-8") as f:
            json.dump(prov, f, indent=2)

    # 8. Check policy
    if args.unresolved_policy == "error" and len(dedup_unresolved) > 0:
        print(f"[taxonomy] ERROR: unresolved_policy is 'error' and {len(dedup_unresolved)} nodes unresolved.", file=sys.stderr)
        sys.exit(1)

    print("[taxonomy] Resolution complete.")
    sys.exit(0)

if __name__ == "__main__":
    main()
