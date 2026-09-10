#!/usr/bin/env python3
"""Resolve NCBI TaxIDs without mutating the source cache in offline mode."""

import argparse
import contextlib
import csv
import gzip
import hashlib
import json
import os
import re
import socket
import sys
import tempfile
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter

TOOL_NAME = "ont_wf16s_postprocess"
MAX_SAFE_INTEGER = 9007199254740991
MAX_READ_LENGTH = 2147483647
PLACEHOLDER_NAMES = {"unknown", "unclassified", "uncultured", "unidentified"}
READ_LENGTH_RE = re.compile(r"^[0-9]+$|^[0-9]+\|[1-9][0-9]*$")
EXPECTED_RANKS = ("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")


@contextlib.contextmanager
def acquire_cache_lock(cache_path, timeout=10.0, poll_interval=0.05):
    lock_path = cache_path + ".lock"
    os.makedirs(os.path.dirname(os.path.abspath(lock_path)), exist_ok=True)
    deadline = time.time() + timeout
    fd = None
    acquired = False

    while time.time() < deadline:
        try:
            flags = os.O_RDWR | os.O_CREAT
            if hasattr(os, "O_BINARY"):
                flags |= os.O_BINARY
            candidate_fd = os.open(lock_path, flags, 0o666)
            if sys.platform == "win32":
                import msvcrt
                try:
                    msvcrt.locking(candidate_fd, msvcrt.LK_NBLCK, 1)
                    fd = candidate_fd
                    acquired = True
                    break
                except OSError:
                    os.close(candidate_fd)
            else:
                import fcntl
                try:
                    # R's filelock package uses POSIX record locks on Unix.
                    # lockf() uses that same fcntl lock family; flock() does
                    # not contend with it on Linux and would allow concurrent
                    # R/Python cache writers.
                    fcntl.lockf(candidate_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    fd = candidate_fd
                    acquired = True
                    break
                except (OSError, IOError):
                    os.close(candidate_fd)
        except OSError:
            pass
        time.sleep(poll_interval)

    if not acquired:
        owner_info = "unknown"
        if os.path.exists(lock_path):
            try:
                with open(lock_path, "r", encoding="utf-8") as h:
                    owner_info = h.read().strip()
            except OSError:
                pass
        raise SystemExit(
            f"[taxonomy] ERROR: E_TAXONOMY_CACHE_BUSY: cache lock '{lock_path}' is held by another process: {owner_info}"
        )

    try:
        os.lseek(fd, 0, os.SEEK_SET)
        os.ftruncate(fd, 0)
        owner_payload = json.dumps({
            "pid": os.getpid(),
            "hostname": socket.gethostname(),
            "start_time": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        })
        os.write(fd, owner_payload.encode("utf-8"))
        yield lock_path
    finally:
        if fd is not None:
            try:
                os.lseek(fd, 0, os.SEEK_SET)
                os.ftruncate(fd, 0)
                if sys.platform == "win32":
                    import msvcrt
                    msvcrt.locking(fd, msvcrt.LK_UNLCK, 1)
                else:
                    import fcntl
                    fcntl.lockf(fd, fcntl.LOCK_UN)
            except OSError:
                pass
            try:
                os.close(fd)
            except OSError:
                pass
            # Keep a stable lock inode. Unlinking after unlock is unsafe on POSIX:
            # a waiter can acquire the old inode while a third process creates and
            # locks a new file at the same pathname.


def parse_taxid(value, context="TaxID"):
    """Parse the shared canonical decimal-string TaxID contract."""
    text = str(value)
    if not text.isascii() or not text.isdigit() or (len(text) > 1 and text.startswith("0")):
        raise ValueError(f"{context}: expected canonical unsigned decimal TaxID, found {text!r}.")
    parsed = int(text)
    if parsed > MAX_SAFE_INTEGER:
        raise ValueError(f"{context}: TaxID exceeds {MAX_SAFE_INTEGER}.")
    return parsed


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


def compute_json_sha256(payload):
    encoded = json.dumps(payload, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    return hashlib.sha256(encoded).hexdigest()


def verify_expected_file(path, expected_sha256, context):
    if expected_sha256 is None:
        return
    actual = compute_sha256(path)
    if actual != expected_sha256:
        raise ValueError(f"{context}: input changed before or after read for {path!r}.")


def load_cache(path, expected_sha256=None):
    verify_expected_file(path, expected_sha256, "taxonomy cache")
    if not os.path.exists(path):
        return {}
    try:
        with open(path, "r", encoding="utf-8") as handle:
            cache = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Could not read taxonomy cache safely: {exc}") from exc
    if not isinstance(cache, dict):
        raise ValueError("Taxonomy cache root must be a JSON object.")
    normalized = {}
    for taxon_path, taxid in cache.items():
        if not isinstance(taxon_path, str) or isinstance(taxid, bool) or not isinstance(taxid, (int, str)):
            raise ValueError(f"Invalid cache entry for {taxon_path!r}: expected a canonical TaxID.")
        normalized[taxon_path] = parse_taxid(taxid, f"cache entry {taxon_path!r}")
    verify_expected_file(path, expected_sha256, "taxonomy cache")
    return normalized


def read_abundance_paths(path, tax_column, expected_sha256=None):
    verify_expected_file(path, expected_sha256, "abundance")
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or tax_column not in reader.fieldnames:
            raise ValueError(f"Abundance table does not contain tax column {tax_column!r}.")
        paths = [row[tax_column].strip() for row in reader if row.get(tax_column, "").strip()]
    verify_expected_file(path, expected_sha256, "abundance")
    return list(dict.fromkeys(paths))


def read_assignment_taxids(paths, expected_hashes=None):
    lineage_to_taxids = {}
    expected_hashes = expected_hashes or {}
    for path in paths:
        verify_expected_file(path, expected_hashes.get(path), "assignment")
        seen_read_ids = set()
        opener = gzip.open if str(path).lower().endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.rstrip("\r\n"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 5:
                    raise ValueError(f"{path}:{line_number}: expected exactly 5 assignment fields.")
                status, read_id, length_field, lineage = fields[0], fields[1], fields[3], fields[4]
                if status not in {"C", "U"}:
                    raise ValueError(f"{path}:{line_number}: invalid status {status!r}.")
                if not read_id or read_id != read_id.strip() or any(ord(ch) < 32 for ch in read_id):
                    raise ValueError(f"{path}:{line_number}: unsafe or empty read ID.")
                if read_id in seen_read_ids:
                    raise ValueError(f"{path}:{line_number}: duplicate read ID {read_id!r}.")
                seen_read_ids.add(read_id)
                if not READ_LENGTH_RE.fullmatch(length_field):
                    raise ValueError(f"{path}:{line_number}: malformed read length field.")
                read_length = int(length_field.rsplit("|", 1)[-1])
                if read_length <= 0 or read_length > MAX_READ_LENGTH:
                    raise ValueError(f"{path}:{line_number}: read length must satisfy 0 < length <= {MAX_READ_LENGTH}.")
                taxid_text = fields[2]
                if any(field != field.strip() or any(ord(ch) < 32 for ch in field)
                       for field in (fields[1], taxid_text, lineage)):
                    raise ValueError(f"{path}:{line_number}: unsafe boundary whitespace or control character.")
                taxid = parse_taxid(taxid_text, f"{path}:{line_number}")
                if status == "U" and taxid != 0:
                    raise ValueError(f"{path}:{line_number}: status U requires TaxID 0.")
                if taxid > 0 and (not lineage or any(not part for part in lineage.split("|"))):
                    raise ValueError(f"{path}:{line_number}: positive TaxID requires a complete lineage.")
                if taxid > 0:
                    lineage_to_taxids.setdefault(lineage, []).append(taxid)
        verify_expected_file(path, expected_hashes.get(path), "assignment")

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
            if parts[depth - 1].strip().lower() in PLACEHOLDER_NAMES:
                continue
            if cache.get(subpath, 0) <= 0:
                unresolved[subpath] = {"path": subpath, "depth": depth, "name": parts[depth - 1]}
    return [unresolved[path] for path in sorted(unresolved)]


def iter_taxon_nodes(abundance_paths):
    """Return every classified lineage node in deterministic path order."""
    nodes = {}
    for path in abundance_paths:
        if path.startswith("Unclassified;"):
            continue
        parts = path.split(";")
        for depth in range(1, len(parts) + 1):
            subpath = ";".join(parts[:depth])
            nodes[subpath] = {"path": subpath, "depth": depth, "name": parts[depth - 1]}
    return [nodes[path] for path in sorted(nodes)]


def fetch_taxonomy_context(taxid, email, api_key):
    params = {
        "db": "taxonomy",
        "id": str(taxid),
        "retmode": "xml",
        "tool": TOOL_NAME,
        "email": email,
    }
    if api_key:
        params["api_key"] = api_key
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(url, headers={"User-Agent": f"{TOOL_NAME}/1.0"})
    with urllib.request.urlopen(request, timeout=15) as response:
        root = ET.fromstring(response.read())
    taxon = root.find(".//Taxon")
    if taxon is None:
        raise ValueError(f"NCBI returned no taxonomy record for TaxID {taxid}.")
    returned_taxid = taxon.findtext("TaxId")
    if returned_taxid is None or parse_taxid(returned_taxid, "NCBI returned TaxID") != int(taxid):
        raise ValueError(f"NCBI taxonomy response did not identify requested TaxID {taxid}.")
    scientific_name = taxon.findtext("ScientificName")
    rank = taxon.findtext("Rank")
    lineage = [node.findtext("ScientificName") for node in taxon.findall("./LineageEx/Taxon")]
    lineage = [value for value in lineage if value]
    if not scientific_name or not rank:
        raise ValueError(f"NCBI taxonomy record for TaxID {taxid} lacks name or rank.")
    return {"scientific_name": scientific_name, "rank": rank, "lineage": lineage}


def _normalized_taxon_name(value):
    return " ".join(str(value).split()).casefold()


def validate_taxonomy_context(record, name, expected_rank=None, ancestor_names=None):
    if _normalized_taxon_name(record["scientific_name"]) != _normalized_taxon_name(name):
        return "NCBI TaxID scientific name does not exactly match the requested name."
    if expected_rank and record["rank"].casefold() != expected_rank.casefold():
        return (f"NCBI TaxID rank mismatch for {name!r}: expected {expected_rank!r}, "
                f"found {record['rank']!r}.")
    if not ancestor_names:
        return None
    normalized_lineage = [_normalized_taxon_name(val) for val in record.get("lineage", [])]
    expected_ancestors = [_normalized_taxon_name(val) for val in ancestor_names]

    missing = [anc for anc in expected_ancestors if anc not in normalized_lineage]
    if missing:
        return (f"NCBI TaxID ancestry mismatch for {name!r}; missing ancestor context: "
                + ", ".join(missing))

    last_idx = -1
    for anc in expected_ancestors:
        indices = [i for i, val in enumerate(normalized_lineage) if val == anc]
        if len(indices) > 1:
            return (f"NCBI TaxID ancestry ambiguity for {name!r}; ancestor {anc!r} "
                    f"appears {len(indices)} times in NCBI lineage.")
        idx = indices[0]
        if idx <= last_idx:
            return (f"NCBI TaxID ancestry order mismatch for {name!r}; ancestor {anc!r} "
                    f"appears out of order or reversed in NCBI lineage.")
        last_idx = idx
    return None


def query_exact_scientific_name(name, email, api_key, attempts=3,
                                expected_rank=None, ancestor_names=None):
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
            ids = sorted({int(taxid) for taxid in result.get("esearchresult", {}).get("idlist", [])})
            if len(ids) > 1:
                return 0, None, ids
            if not ids:
                return 0, None, []
            record = fetch_taxonomy_context(ids[0], email, api_key)
            context_error = validate_taxonomy_context(
                record, name, expected_rank=expected_rank, ancestor_names=ancestor_names
            )
            if context_error:
                return 0, context_error, []
            return ids[0], None, []
        except Exception as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            if attempt + 1 < attempts:
                time.sleep(2 ** attempt)
    return 0, last_error, []


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


def write_resolution_sources_tsv(path, nodes, cache, resolution_sources):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["TaxonPath", "TaxID", "ResolutionSource"])
        for item in nodes:
            taxon_path = item["path"]
            taxid = cache.get(taxon_path, 0)
            source = resolution_sources.get(taxon_path, "unresolved") if taxid > 0 else "unresolved"
            writer.writerow([taxon_path, taxid, source])


def main():
    parser = argparse.ArgumentParser(description="Resolve NCBI TaxIDs from wf-16s data")
    parser.add_argument("--abundance", required=True)
    parser.add_argument("--tax-column", default="tax")
    parser.add_argument("--assignments", action="append", default=[])
    parser.add_argument("--expected-input", action="append", default=[],
                        help="Expected input fingerprint as PATH<TAB>SHA256; may be repeated")
    parser.add_argument("--cache", required=True, help="Read-only source cache in cache_only mode")
    parser.add_argument("--resolved-cache", help="Run-local resolved cache output")
    parser.add_argument("--mode", choices=["cache_only", "refresh"], default="cache_only")
    parser.add_argument("--email-env", default="NCBI_EMAIL")
    parser.add_argument("--api-key-env", default="NCBI_API_KEY")
    parser.add_argument("--unresolved-policy", choices=["warn", "error"], default="warn")
    parser.add_argument("--unresolved-tsv")
    parser.add_argument("--conflicts-tsv")
    parser.add_argument("--resolution-sources-tsv")
    parser.add_argument("--provenance")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--online-preflight", action="store_true")
    parser.add_argument("--defer-cache-commit", action="store_true")
    parser.add_argument("--cache-lock-held", action="store_true")
    parser.add_argument("--cache-lock-owner-pid", type=int)
    parser.add_argument("--transaction-id")
    args = parser.parse_args()

    try:
        if args.cache_lock_held:
            if args.mode != "refresh" or not args.defer_cache_commit:
                raise ValueError("--cache-lock-held requires deferred refresh mode.")
            if args.cache_lock_owner_pid not in {os.getpid(), os.getppid()}:
                raise ValueError(
                    "--cache-lock-held requires --cache-lock-owner-pid matching this process or its parent."
                )
        elif args.cache_lock_owner_pid is not None:
            raise ValueError("--cache-lock-owner-pid requires --cache-lock-held.")
        if args.transaction_id is not None and not re.fullmatch(
                r"tx-[0-9a-f]{64}", args.transaction_id):
            raise ValueError(
                "--transaction-id must be tx- followed by 64 lowercase hex characters."
            )
        expected_inputs = {}
        for specification in args.expected_input:
            if "\t" in specification:
                path, _, expected_sha256 = specification.partition("\t")
            elif "=" in specification:
                path, _, expected_sha256 = specification.partition("=")
            else:
                raise ValueError("--expected-input must be PATH<TAB>64-hex-SHA256 or PATH=64-hex-SHA256.")
            if not path or not re.fullmatch(r"[0-9a-fA-F]{64}", expected_sha256):
                raise ValueError("--expected-input must specify a non-empty path and 64-hex-SHA256.")
            normalized_sha = expected_sha256.lower()
            if path in expected_inputs and expected_inputs[path] != normalized_sha:
                raise ValueError(f"Conflicting --expected-input for {path!r}: {expected_inputs[path]} vs {normalized_sha}.")
            expected_inputs[path] = normalized_sha

        if not os.path.exists(args.abundance):
            raise ValueError(f"Abundance file not found: {args.abundance}")
        abundance_paths = read_abundance_paths(
            args.abundance, args.tax_column, expected_inputs.get(args.abundance)
        )

        lock_ctx = (
            acquire_cache_lock(args.cache)
            if (args.mode == "refresh" and not args.validate_only and not args.cache_lock_held)
            else contextlib.nullcontext()
        )
        with lock_ctx:
            cache_sha_before = compute_sha256(args.cache)
            cache = load_cache(args.cache, expected_inputs.get(args.cache))
            resolution_sources = {
                taxon_path: "source_cache" for taxon_path, taxid in cache.items() if taxid > 0
            }
            assignment_hashes = {path: expected_inputs[path] for path in args.assignments if path in expected_inputs}
            assignment_map, conflicts = read_assignment_taxids(args.assignments, assignment_hashes)
            conflict_lineages = {item["lineage"] for item in conflicts}

            for path in abundance_paths:
                parts = path.split(";")
                if len(parts) == 8 and cache.get(path, 0) <= 0:
                    assignment_taxid = assignment_map.get(normalize_abundance_path_to_7(path), 0)
                    if assignment_taxid > 0:
                        cache[path] = assignment_taxid
                        resolution_sources[path] = (
                            "assignment_conflict"
                            if normalize_abundance_path_to_7(path) in conflict_lineages
                            else "assignment"
                        )

            unresolved = find_unresolved(abundance_paths, cache)
            query_failures = []
            ambiguous_queries = []
            cache_updated = False

            if args.mode == "refresh" and unresolved:
                email = os.environ.get(args.email_env, "").strip()
                if not email:
                    raise ValueError(f"Environment variable {args.email_env!r} is required for refresh mode.")
                if args.validate_only and not args.online_preflight:
                    print(
                        f"[taxonomy] ERROR: E_ONLINE_PREFLIGHT_REQUIRED: {len(unresolved)} node(s) require online refresh, but --online-preflight was not specified.",
                        file=sys.stderr,
                    )
                    return 1
                api_key = os.environ.get(args.api_key_env, "").strip() or None
                delay = 0.12 if api_key else 0.35
                name_results = {}
                for item in unresolved:
                    name = item["name"]
                    expected_rank = EXPECTED_RANKS[item["depth"] - 1] if item["depth"] <= len(EXPECTED_RANKS) else None
                    ancestor_names = [part for part in item["path"].split(";")[:-1]
                                      if part.strip().lower() not in PLACEHOLDER_NAMES]
                    query_key = (name, expected_rank, tuple(ancestor_names))
                    if query_key not in name_results:
                        name_results[query_key] = query_exact_scientific_name(
                            name, email, api_key, expected_rank=expected_rank,
                            ancestor_names=ancestor_names
                        )
                        time.sleep(delay)
                    taxid, error, ambiguous_taxids = name_results[query_key]
                    if error:
                        query_failures.append({"path": item["path"], "name": name, "error": error})
                    elif ambiguous_taxids:
                        ambiguous_queries.append({
                            "path": item["path"], "name": name, "taxids": ambiguous_taxids
                        })
                    elif taxid > 0:
                        cache[item["path"]] = taxid
                        resolution_sources[item["path"]] = "ncbi_refresh"

                unresolved = find_unresolved(abundance_paths, cache)
            candidate_cache = {key: str(value) for key, value in cache.items()}
            cache_sha_candidate = compute_json_sha256(candidate_cache)
            if args.validate_only:
                if args.unresolved_policy == "error" and unresolved:
                    raise ValueError(f"{len(unresolved)} taxonomy nodes remain unresolved.")
                print(f"[taxonomy] Preflight complete: {len(unresolved)} unresolved, {len(conflicts)} conflicts.")
                return 0

            required_outputs = (args.resolved_cache, args.unresolved_tsv, args.conflicts_tsv,
                                args.resolution_sources_tsv, args.provenance)
            if any(not value for value in required_outputs):
                raise ValueError("Resolver output paths are required unless --validate-only is used.")

            atomic_write_json(args.resolved_cache, {key: str(value) for key, value in cache.items()})
            write_unresolved_tsv(args.unresolved_tsv, unresolved)
            write_conflicts_tsv(args.conflicts_tsv, conflicts)
            taxon_nodes = iter_taxon_nodes(abundance_paths)
            write_resolution_sources_tsv(
                args.resolution_sources_tsv, taxon_nodes, cache, resolution_sources
            )
            source_counts = Counter(
                resolution_sources.get(item["path"], "unresolved")
                if cache.get(item["path"], 0) > 0 else "unresolved"
                for item in taxon_nodes
            )
            source_counts = {
                label: source_counts.get(label, 0)
                for label in (
                    "source_cache", "assignment", "assignment_conflict",
                    "ncbi_refresh", "unresolved"
                )
            }

            provenance = {
                "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "transaction_id": args.transaction_id,
                "mode": args.mode,
                "tool": TOOL_NAME,
                "abundance_sha256": compute_sha256(args.abundance),
                "assignments_sha256": {path: compute_sha256(path) for path in args.assignments},
                "source_cache_sha256_before": cache_sha_before,
                "source_cache_sha256_after": compute_sha256(args.cache),
                "source_cache_sha256_candidate": cache_sha_candidate,
                "source_cache_sha256_committed": cache_sha_before
                if (args.mode == "cache_only" or args.defer_cache_commit) else None,
                "resolved_cache_sha256": compute_sha256(args.resolved_cache),
                "source_cache_updated": cache_updated,
                "source_cache_commit_deferred": bool(
                    args.mode == "refresh" and args.defer_cache_commit
                ),
                "total_lineages": len(abundance_paths),
                "unresolved_count": len(unresolved),
                "conflicts_count": len(conflicts),
                "conflicts": conflicts,
                "query_failures": query_failures,
                "ambiguous_queries": ambiguous_queries,
                "resolution_source_counts": source_counts,
            }
            atomic_write_json(args.provenance, provenance)

            if query_failures:
                print(f"[taxonomy] ERROR: {len(query_failures)} NCBI query failure(s); source cache preserved.", file=sys.stderr)
                return 1
            if args.unresolved_policy == "error" and unresolved:
                print(f"[taxonomy] ERROR: {len(unresolved)} taxonomy nodes remain unresolved.", file=sys.stderr)
                return 1
            if args.mode == "refresh":
                current_cache_sha = compute_sha256(args.cache)
                if current_cache_sha != cache_sha_before:
                    print(
                        f"[taxonomy] ERROR: E_TAXONOMY_CACHE_CHANGED: source cache '{args.cache}' "
                        f"was mutated during query execution (expected {cache_sha_before}, found {current_cache_sha}).",
                        file=sys.stderr,
                    )
                    return 1
                if not args.defer_cache_commit:
                    atomic_write_json(args.cache, candidate_cache)
                    cache_updated = True
                    provenance["source_cache_updated"] = True
                    provenance["source_cache_sha256_after"] = compute_sha256(args.cache)
                    provenance["source_cache_sha256_committed"] = provenance["source_cache_sha256_after"]
                    atomic_write_json(args.provenance, provenance)
            print(f"[taxonomy] Resolution complete: {len(unresolved)} unresolved, {len(conflicts)} conflicts.")
            return 0
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"[taxonomy] ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
