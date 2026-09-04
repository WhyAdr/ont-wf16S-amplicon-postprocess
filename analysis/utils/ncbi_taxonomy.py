#!/usr/bin/env python3
"""Resolve NCBI TaxIDs without mutating the source cache in offline mode."""

import argparse
import csv
import hashlib
import json
import os
import sys
import tempfile
import time
import urllib.parse
import urllib.request
from collections import Counter

TOOL_NAME = "ont_wf16s_postprocess"


def normalize_abundance_path_to_7(path):
    """Omit the eight-rank abundance schema's kingdom field."""
    parts = path.split(";")
    return "|".join([parts[0]] + parts[2:8]) if len(parts) == 8 else "|".join(parts)


def compute_sha256(filepath):
    if not filepath or not os.path.exists(filepath):
        return None
    digest = hashlib.sha256()
    with open(filepath, "rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_json(path, payload):
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    temp_path = None
    try:
        with tempfile.NamedTemporaryFile("w", dir=directory, delete=False, encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            temp_path = handle.name
        os.replace(temp_path, path)
    finally:
        if temp_path and os.path.exists(temp_path):
            os.unlink(temp_path)


def load_cache(path):
    if not os.path.exists(path):
        return {}
    try:
        with open(path, "r", encoding="utf-8") as handle:
            cache = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not read taxonomy cache safely: {exc}") from exc
    if not isinstance(cache, dict):
        raise ValueError("Taxonomy cache root must be a JSON object.")
    for taxon_path, taxid in cache.items():
        if not isinstance(taxon_path, str) or not isinstance(taxid, int) or isinstance(taxid, bool) or taxid < 0:
            raise ValueError(f"Invalid cache entry for {taxon_path!r}: expected a non-negative integer TaxID.")
    return cache


def read_abundance_paths(path, tax_column):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or tax_column not in reader.fieldnames:
            raise ValueError(f"Abundance table does not contain tax column {tax_column!r}.")
        paths = [row[tax_column].strip() for row in reader if row.get(tax_column, "").strip()]
    return list(dict.fromkeys(paths))


def read_assignment_taxids(paths):
    lineage_to_taxids = {}
    for path in paths:
        with open(path, "r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.rstrip("\r\n"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 5:
                    raise ValueError(f"{path}:{line_number}: expected exactly 5 assignment fields.")
                taxid_text = fields[2].strip()
                lineage = fields[4].strip()
                if not taxid_text.isdigit():
                    raise ValueError(f"{path}:{line_number}: invalid TaxID {taxid_text!r}.")
                taxid = int(taxid_text)
                if taxid > 0 and lineage:
                    lineage_to_taxids.setdefault(lineage, []).append(taxid)

    resolved = {}
    conflicts = []
    for lineage in sorted(lineage_to_taxids):
        counts = Counter(lineage_to_taxids[lineage])
        maximum = max(counts.values())
        winner = min(taxid for taxid, count in counts.items() if count == maximum)
        resolved[lineage] = winner
        if len(counts) > 1:
            conflicts.append({
                "lineage": lineage,
                "winner_taxid": winner,
                "counts": {str(key): counts[key] for key in sorted(counts)},
            })
    return resolved, conflicts


def find_unresolved(abundance_paths, cache):
    unresolved = {}
    for path in abundance_paths:
        if path.startswith("Unclassified;"):
            continue
        parts = path.split(";")
        for depth in range(1, len(parts) + 1):
            subpath = ";".join(parts[:depth])
            if cache.get(subpath, 0) <= 0:
                unresolved[subpath] = {"path": subpath, "depth": depth, "name": parts[depth - 1]}
    return [unresolved[path] for path in sorted(unresolved)]


def query_exact_scientific_name(name, email, api_key, attempts=3):
    params = {
        "db": "taxonomy",
        "term": f'"{name}"[Scientific Name]',
        "retmode": "json",
        "tool": TOOL_NAME,
        "email": email,
    }
    if api_key:
        params["api_key"] = api_key
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(params)
    last_error = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": f"{TOOL_NAME}/1.0"})
            with urllib.request.urlopen(request, timeout=15) as response:
                result = json.loads(response.read().decode("utf-8"))
            ids = result.get("esearchresult", {}).get("idlist", [])
            return (min(map(int, ids)) if ids else 0), None
        except Exception as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            if attempt + 1 < attempts:
                time.sleep(2 ** attempt)
    return 0, last_error


def write_unresolved_tsv(path, unresolved):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Depth", "NodeName", "TaxonPath"])
        for item in unresolved:
            writer.writerow([item["depth"], item["name"], item["path"]])


def write_conflicts_tsv(path, conflicts):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Lineage", "WinnerTaxID", "TaxIDCountsJSON"])
        for item in conflicts:
            writer.writerow([item["lineage"], item["winner_taxid"], json.dumps(item["counts"], sort_keys=True)])


def main():
    parser = argparse.ArgumentParser(description="Resolve NCBI TaxIDs from wf-16s data")
    parser.add_argument("--abundance", required=True)
    parser.add_argument("--tax-column", default="tax")
    parser.add_argument("--assignments", action="append", default=[])
    parser.add_argument("--cache", required=True, help="Read-only source cache in cache_only mode")
    parser.add_argument("--resolved-cache", required=True, help="Run-local resolved cache output")
    parser.add_argument("--mode", choices=["cache_only", "refresh"], default="cache_only")
    parser.add_argument("--email-env", default="NCBI_EMAIL")
    parser.add_argument("--api-key-env", default="NCBI_API_KEY")
    parser.add_argument("--unresolved-policy", choices=["warn", "error"], default="warn")
    parser.add_argument("--unresolved-tsv", required=True)
    parser.add_argument("--conflicts-tsv", required=True)
    parser.add_argument("--provenance", required=True)
    args = parser.parse_args()

    try:
        if not os.path.exists(args.abundance):
            raise ValueError(f"Abundance file not found: {args.abundance}")
        abundance_paths = read_abundance_paths(args.abundance, args.tax_column)
        cache_sha_before = compute_sha256(args.cache)
        cache = load_cache(args.cache)
        assignment_map, conflicts = read_assignment_taxids(args.assignments)

        for path in abundance_paths:
            parts = path.split(";")
            if len(parts) == 8 and cache.get(path, 0) <= 0:
                assignment_taxid = assignment_map.get(normalize_abundance_path_to_7(path), 0)
                if assignment_taxid > 0:
                    cache[path] = assignment_taxid

        unresolved = find_unresolved(abundance_paths, cache)
        query_failures = []
        cache_updated = False

        if args.mode == "refresh" and unresolved:
            email = os.environ.get(args.email_env, "").strip()
            if not email:
                raise ValueError(f"Environment variable {args.email_env!r} is required for refresh mode.")
            api_key = os.environ.get(args.api_key_env, "").strip() or None
            delay = 0.12 if api_key else 0.35
            name_results = {}
            for item in unresolved:
                name = item["name"]
                if name not in name_results:
                    name_results[name] = query_exact_scientific_name(name, email, api_key)
                    time.sleep(delay)
                taxid, error = name_results[name]
                if error:
                    query_failures.append({"path": item["path"], "name": name, "error": error})
                elif taxid > 0:
                    cache[item["path"]] = taxid

            unresolved = find_unresolved(abundance_paths, cache)
            if not query_failures:
                atomic_write_json(args.cache, cache)
                cache_updated = True

        atomic_write_json(args.resolved_cache, cache)
        write_unresolved_tsv(args.unresolved_tsv, unresolved)
        write_conflicts_tsv(args.conflicts_tsv, conflicts)

        provenance = {
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "mode": args.mode,
            "tool": TOOL_NAME,
            "abundance_sha256": compute_sha256(args.abundance),
            "assignments_sha256": {path: compute_sha256(path) for path in args.assignments},
            "source_cache_sha256_before": cache_sha_before,
            "source_cache_sha256_after": compute_sha256(args.cache),
            "resolved_cache_sha256": compute_sha256(args.resolved_cache),
            "source_cache_updated": cache_updated,
            "total_lineages": len(abundance_paths),
            "unresolved_count": len(unresolved),
            "conflicts_count": len(conflicts),
            "conflicts": conflicts,
            "query_failures": query_failures,
        }
        atomic_write_json(args.provenance, provenance)

        if query_failures:
            print(f"[taxonomy] ERROR: {len(query_failures)} NCBI query failure(s); source cache preserved.", file=sys.stderr)
            return 1
        if args.unresolved_policy == "error" and unresolved:
            print(f"[taxonomy] ERROR: {len(unresolved)} taxonomy nodes remain unresolved.", file=sys.stderr)
            return 1
        print(f"[taxonomy] Resolution complete: {len(unresolved)} unresolved, {len(conflicts)} conflicts.")
        return 0
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"[taxonomy] ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
