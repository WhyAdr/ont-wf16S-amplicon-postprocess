#!/usr/bin/env python3
"""Resolve NCBI TaxIDs without mutating the source cache in offline mode."""

import argparse
import contextlib
import csv
import dataclasses
import datetime
import email.utils
import gzip
import hashlib
import http.client
import json
import os
import re
import socket
import sys
import tempfile
import threading
import time
import urllib.error
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
CELLULAR_DOMAINS = {"bacteria", "archaea", "eukaryota"}
FATAL_OUTCOMES = {"request_failed", "response_invalid"}
SEMANTIC_OUTCOMES = {"not_found", "ambiguous", "context_mismatch"}
MAX_REQUEST_ATTEMPTS = 3
MAX_RETRY_AFTER_SECONDS = 30.0
_ACTIVE_CACHE_LOCKS = set()
_ACTIVE_CACHE_LOCKS_GUARD = threading.Lock()


@dataclasses.dataclass(frozen=True)
class LookupOutcome:
    status: str
    taxid: int = 0
    candidates: tuple = ()
    code: str = None
    message: str = None
    expected_rank: str = None
    returned_rank: str = None
    rank_rule: str = None


class SafeRequestFailure(Exception):
    """A credential-safe, classified request failure."""

    def __init__(self, code, message):
        super().__init__(message)
        self.code = code
        self.safe_message = message


class InvalidResponse(Exception):
    """A credential-safe response validation failure."""

    def __init__(self, message):
        super().__init__(message)
        self.safe_message = message


class RequestLimiter:
    def __init__(self, interval, monotonic=time.monotonic, sleep=time.sleep):
        self.interval = float(interval)
        self.monotonic = monotonic
        self.sleep = sleep
        self._last_request = None

    def wait_for_request_slot(self):
        now = self.monotonic()
        if self._last_request is not None:
            remaining = self.interval - (now - self._last_request)
            if remaining > 0:
                self.sleep(remaining)
                now = self.monotonic()
        self._last_request = now


class EventWriter:
    """Flush bounded, credential-safe JSONL diagnostics for execution failures."""

    MAX_MESSAGE_LENGTH = 500

    def __init__(self, diagnostics_dir=None, transaction_id=None):
        self.path = None
        self.summary_path = None
        self._handle = None
        self.started = time.monotonic()
        if diagnostics_dir:
            if not transaction_id:
                raise ValueError("--diagnostics-dir requires --transaction-id.")
            directory = os.path.abspath(diagnostics_dir)
            if not os.path.isdir(directory) or os.path.islink(directory):
                raise ValueError("Diagnostics directory must be a pre-created, non-symlink directory.")
            self.path = os.path.join(directory, "taxonomy_events.jsonl")
            self.summary_path = os.path.join(directory, "taxonomy_failure.json")
            flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
            fd = os.open(self.path, flags, 0o600)
            self._handle = os.fdopen(fd, "w", encoding="utf-8", newline="\n")

    @staticmethod
    def _safe_text(value):
        text = " ".join(str(value).split())
        text = re.sub(r"(?i)(api[_-]?key|email)=([^&\s]+)", r"\1=[REDACTED]", text)
        text = re.sub(r"https?://\S+", "[REDACTED_URL]", text)
        return text[:EventWriter.MAX_MESSAGE_LENGTH]

    def emit(self, event, **fields):
        record = {
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "elapsed_seconds": round(time.monotonic() - self.started, 3),
            "event": self._safe_text(event),
        }
        for key, value in fields.items():
            if value is not None:
                record[key] = self._safe_text(value) if isinstance(value, str) else value
        line = json.dumps(record, sort_keys=True, ensure_ascii=True)
        print(line, file=sys.stderr, flush=True)
        if self._handle is not None:
            self._handle.write(line + "\n")
            self._handle.flush()

    def write_failure_summary(self, outcome, completed, total):
        if not self.summary_path:
            return
        payload = {
            "status": "failed",
            "code": outcome.code,
            "outcome": outcome.status,
            "message": self._safe_text(outcome.message or "taxonomy resolution failed"),
            "completed": completed,
            "total": total,
            "events_file": os.path.basename(self.path),
        }
        atomic_write_json(self.summary_path, payload)

    def close(self):
        if self._handle is not None:
            self._handle.close()
            self._handle = None


def _try_lock_fd(fd):
    """Acquire one non-blocking record lock using the cross-runtime contract."""
    if sys.platform == "win32":
        import msvcrt
        msvcrt.locking(fd, msvcrt.LK_NBLCK, 1)
    else:
        import fcntl
        # R's filelock package uses POSIX record locks on Unix.  lockf() uses
        # that same fcntl lock family; flock() does not contend with it.
        fcntl.lockf(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)


def _unlock_fd(fd):
    if sys.platform == "win32":
        import msvcrt
        msvcrt.locking(fd, msvcrt.LK_UNLCK, 1)
    else:
        import fcntl
        fcntl.lockf(fd, fcntl.LOCK_UN)


def _write_all(fd, payload):
    """Write complete owner metadata even when os.write() is short."""
    view = memoryview(payload)
    while view:
        written = os.write(fd, view)
        if written <= 0:
            raise OSError("lock-owner metadata write made no progress")
        view = view[written:]


def _read_lock_owner(lock_path):
    try:
        with open(lock_path, "r", encoding="utf-8") as handle:
            return handle.read().strip() or "unknown"
    except OSError:
        return "unknown"


@contextlib.contextmanager
def acquire_cache_lock(cache_path, timeout=10.0, poll_interval=0.05):
    cache_identity = os.path.realpath(os.path.abspath(os.fspath(cache_path)))
    lock_path = cache_identity + ".lock"
    os.makedirs(os.path.dirname(os.path.abspath(lock_path)), exist_ok=True)
    lock_identity = os.path.normcase(os.path.realpath(lock_path))
    registered = False
    fd = None
    try:
        with _ACTIVE_CACHE_LOCKS_GUARD:
            if lock_identity in _ACTIVE_CACHE_LOCKS:
                raise SystemExit(
                    f"[taxonomy] ERROR: E_TAXONOMY_CACHE_BUSY: cache lock "
                    f"'{lock_path}' is already held by this process"
                )
            _ACTIVE_CACHE_LOCKS.add(lock_identity)
            registered = True

        deadline = time.monotonic() + timeout
        while fd is None and time.monotonic() < deadline:
            candidate_fd = None
            try:
                flags = os.O_RDWR | os.O_CREAT
                if hasattr(os, "O_BINARY"):
                    flags |= os.O_BINARY
                candidate_fd = os.open(lock_path, flags, 0o666)
                try:
                    _try_lock_fd(candidate_fd)
                except (OSError, IOError):
                    pass
                else:
                    fd = candidate_fd
                    candidate_fd = None
            except OSError:
                pass
            finally:
                if candidate_fd is not None:
                    try:
                        os.close(candidate_fd)
                    except OSError:
                        pass

            if fd is None:
                remaining = deadline - time.monotonic()
                if remaining > 0:
                    time.sleep(min(poll_interval, remaining))

        if fd is None:
            owner_info = _read_lock_owner(lock_path)
            raise SystemExit(
                f"[taxonomy] ERROR: E_TAXONOMY_CACHE_BUSY: cache lock '{lock_path}' "
                f"is held by another process: {owner_info}"
            )

        os.lseek(fd, 0, os.SEEK_SET)
        os.ftruncate(fd, 0)
        owner_payload = json.dumps({
            "pid": os.getpid(),
            "hostname": socket.gethostname(),
            "start_time": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        })
        _write_all(fd, owner_payload.encode("utf-8"))
        yield lock_path
    finally:
        if fd is not None:
            try:
                os.lseek(fd, 0, os.SEEK_SET)
                os.ftruncate(fd, 0)
                _unlock_fd(fd)
            except OSError:
                pass
            try:
                os.close(fd)
            except OSError:
                pass
        if registered:
            with _ACTIVE_CACHE_LOCKS_GUARD:
                _ACTIVE_CACHE_LOCKS.discard(lock_identity)
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


def _retry_after_seconds(headers):
    value = headers.get("Retry-After") if headers is not None else None
    try:
        return min(MAX_RETRY_AFTER_SECONDS, max(0.0, float(value)))
    except (TypeError, ValueError):
        try:
            parsed = email.utils.parsedate_to_datetime(value)
            if parsed.tzinfo is None:
                parsed = parsed.replace(tzinfo=datetime.timezone.utc)
            remaining = parsed.timestamp() - time.time()
            return min(MAX_RETRY_AFTER_SECONDS, max(0.0, remaining))
        except (TypeError, ValueError, OverflowError):
            return None


def _safe_http_message(endpoint, status, retryable):
    disposition = "transient" if retryable else "permanent"
    return f"{endpoint} {disposition} HTTP failure (status {status})"


def request_payload(endpoint, params, limiter, events, attempts=MAX_REQUEST_ATTEMPTS,
                    opener=None, sleep=None):
    """Fetch one endpoint with request-level limiting and bounded classified retries."""
    opener = opener or urllib.request.urlopen
    sleep = sleep or time.sleep
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}.fcgi?" +
           urllib.parse.urlencode(params))
    for attempt in range(1, attempts + 1):
        limiter.wait_for_request_slot()
        events.emit("request_start", phase=endpoint, attempt=attempt)
        try:
            request = urllib.request.Request(url, headers={"User-Agent": f"{TOOL_NAME}/1.0"})
            with opener(request, timeout=15) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            retryable = exc.code == 429 or 500 <= exc.code <= 599
            message = _safe_http_message(endpoint, exc.code, retryable)
            if not retryable or attempt == attempts:
                raise SafeRequestFailure("E_NCBI_REQUEST", message) from None
            retry_after = _retry_after_seconds(exc.headers)
            delay = retry_after if retry_after is not None else min(2 ** (attempt - 1), MAX_RETRY_AFTER_SECONDS)
            events.emit("request_retry", phase=endpoint, attempt=attempt,
                        code="E_NCBI_REQUEST", message=message, retry_after_seconds=delay)
            sleep(delay)
        except http.client.IncompleteRead:
            message = f"{endpoint} transient HTTP body failure (IncompleteRead)"
            if attempt == attempts:
                raise SafeRequestFailure("E_NCBI_REQUEST", message) from None
            delay = min(2 ** (attempt - 1), MAX_RETRY_AFTER_SECONDS)
            events.emit("request_retry", phase=endpoint, attempt=attempt,
                        code="E_NCBI_REQUEST", message=message, retry_after_seconds=delay)
            sleep(delay)
        except http.client.HTTPException as exc:
            message = f"{endpoint} transient HTTP protocol failure ({type(exc).__name__})"
            if attempt == attempts:
                raise SafeRequestFailure("E_NCBI_REQUEST", message) from None
            delay = min(2 ** (attempt - 1), MAX_RETRY_AFTER_SECONDS)
            events.emit("request_retry", phase=endpoint, attempt=attempt,
                        code="E_NCBI_REQUEST", message=message, retry_after_seconds=delay)
            sleep(delay)
        except (urllib.error.URLError, TimeoutError, socket.timeout, ConnectionError) as exc:
            message = f"{endpoint} transient request failure ({type(exc).__name__})"
            if attempt == attempts:
                raise SafeRequestFailure("E_NCBI_REQUEST", message) from None
            delay = min(2 ** (attempt - 1), MAX_RETRY_AFTER_SECONDS)
            events.emit("request_retry", phase=endpoint, attempt=attempt,
                        code="E_NCBI_REQUEST", message=message, retry_after_seconds=delay)
            sleep(delay)
        except (KeyboardInterrupt, SystemExit):
            raise
        except OSError as exc:
            message = f"{endpoint} request failure ({type(exc).__name__})"
            raise SafeRequestFailure("E_NCBI_REQUEST", message) from None
    raise AssertionError("unreachable")


def _parse_search_payload(payload):
    try:
        result = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        raise InvalidResponse("esearch returned malformed JSON") from None
    if not isinstance(result, dict) or any(str(key).casefold() == "error" for key in result):
        raise InvalidResponse("esearch returned an explicit API error or non-object payload")
    search = result.get("esearchresult")
    if not isinstance(search, dict) or "idlist" not in search or "count" not in search:
        raise InvalidResponse("esearch response is missing esearchresult.idlist or count")
    if any(str(key).casefold() == "error" for key in search):
        raise InvalidResponse("esearch returned an explicit API error")
    idlist = search["idlist"]
    if not isinstance(idlist, list):
        raise InvalidResponse("esearchresult.idlist is not an array")
    try:
        count = parse_taxid(search["count"], "esearch result count")
        parsed_ids = [parse_taxid(value, "esearch TaxID") for value in idlist]
    except (TypeError, ValueError) as exc:
        raise InvalidResponse(f"esearch response contains an invalid count or TaxID: {exc}") from None
    unique_ids = tuple(sorted(set(parsed_ids)))
    if any(taxid <= 0 for taxid in unique_ids):
        raise InvalidResponse("esearch returned a non-positive candidate TaxID")
    if count != len(idlist):
        raise InvalidResponse(
            f"esearch result count {count} does not match returned ID count {len(idlist)}"
        )
    if len(unique_ids) != len(idlist):
        raise InvalidResponse("esearch returned duplicate candidate TaxIDs")
    return unique_ids


def _parse_taxonomy_context(payload, taxid):
    try:
        root = ET.fromstring(payload)
    except (ET.ParseError, LookupError):
        raise InvalidResponse("efetch returned malformed XML") from None
    error_node = root.find(".//ERROR")
    if error_node is not None:
        raise InvalidResponse("efetch returned an explicit API error")
    taxa = root.findall("./Taxon") if root.tag == "TaxaSet" else root.findall(".//Taxon")
    if len(taxa) != 1:
        raise InvalidResponse(f"efetch returned {len(taxa)} top-level taxonomy records")
    taxon = taxa[0]
    returned_taxid = taxon.findtext("TaxId")
    try:
        returned = parse_taxid(returned_taxid, "NCBI returned TaxID")
    except (TypeError, ValueError):
        raise InvalidResponse("efetch record has a missing or invalid TaxID") from None
    if returned != int(taxid):
        raise InvalidResponse(f"efetch TaxID {returned} does not match requested TaxID {taxid}")
    scientific_name = taxon.findtext("ScientificName")
    rank = taxon.findtext("Rank")
    lineage = [node.findtext("ScientificName") for node in taxon.findall("./LineageEx/Taxon")]
    if not scientific_name or not rank or any(not value for value in lineage):
        raise InvalidResponse(f"efetch record for TaxID {taxid} lacks required taxonomy fields")
    return {"scientific_name": scientific_name, "rank": rank, "lineage": lineage}


def fetch_taxonomy_context(taxid, email, api_key, limiter=None, events=None,
                           attempts=MAX_REQUEST_ATTEMPTS, opener=None, sleep=None):
    params = {"db": "taxonomy", "id": str(taxid), "retmode": "xml",
              "tool": TOOL_NAME, "email": email}
    if api_key:
        params["api_key"] = api_key
    limiter = limiter or RequestLimiter(0.12 if api_key else 0.35)
    events = events or EventWriter()
    payload = request_payload("efetch", params, limiter, events, attempts, opener, sleep)
    return _parse_taxonomy_context(payload, taxid)


def _normalized_taxon_name(value):
    return " ".join(str(value).split()).casefold()


def compatible_rank(expected, returned, name):
    expected_normalized = str(expected).strip().casefold()
    returned_normalized = str(returned).strip().casefold()
    if expected_normalized == returned_normalized:
        return True, "exact"
    if (_normalized_taxon_name(name) in CELLULAR_DOMAINS and
            {expected_normalized, returned_normalized} == {"superkingdom", "domain"}):
        return True, "cellular_domain_legacy_alias"
    return False, None


def validate_taxonomy_context(record, name, expected_rank=None, ancestor_names=None):
    if _normalized_taxon_name(record["scientific_name"]) != _normalized_taxon_name(name):
        return ("E_TAXONOMY_NAME", "NCBI TaxID scientific name does not exactly match the requested name.", None)
    rank_rule = "exact"
    if expected_rank:
        accepted, rank_rule = compatible_rank(expected_rank, record["rank"], name)
        if not accepted:
            return ("E_TAXONOMY_CONTEXT",
                    f"NCBI TaxID rank mismatch for {name!r}: expected {expected_rank!r}, found {record['rank']!r}.",
                    None)
    if not ancestor_names:
        return None, None, rank_rule
    normalized_lineage = [_normalized_taxon_name(val) for val in record.get("lineage", [])]
    expected_ancestors = [_normalized_taxon_name(val) for val in ancestor_names]

    missing = [anc for anc in expected_ancestors if anc not in normalized_lineage]
    if missing:
        return ("E_TAXONOMY_CONTEXT",
                f"NCBI TaxID ancestry mismatch for {name!r}; missing ancestor context: " + ", ".join(missing),
                None)

    last_idx = -1
    for anc in expected_ancestors:
        indices = [i for i, val in enumerate(normalized_lineage) if val == anc]
        if len(indices) > 1:
            return ("E_TAXONOMY_CONTEXT",
                    f"NCBI TaxID ancestry ambiguity for {name!r}; ancestor {anc!r} appears {len(indices)} times in NCBI lineage.",
                    None)
        idx = indices[0]
        if idx <= last_idx:
            return ("E_TAXONOMY_CONTEXT",
                    f"NCBI TaxID ancestry order mismatch for {name!r}; ancestor {anc!r} appears out of order or reversed in NCBI lineage.",
                    None)
        last_idx = idx
    return None, None, rank_rule


def query_exact_scientific_name(name, email, api_key, attempts=MAX_REQUEST_ATTEMPTS,
                                expected_rank=None, ancestor_names=None, limiter=None,
                                events=None, opener=None, sleep=None):
    params = {
        "db": "taxonomy",
        "term": f'"{name}"[Scientific Name]',
        "retmode": "json",
        "tool": TOOL_NAME,
        "email": email,
    }
    if api_key:
        params["api_key"] = api_key
    limiter = limiter or RequestLimiter(0.12 if api_key else 0.35)
    events = events or EventWriter()
    try:
        search_payload = request_payload("esearch", params, limiter, events, attempts, opener, sleep)
        ids = _parse_search_payload(search_payload)
        if not ids:
            return LookupOutcome("not_found", code="E_TAXONOMY_NOT_FOUND",
                                 message=f"No exact NCBI scientific-name match for {name!r}")
        if len(ids) > 1:
            return LookupOutcome("ambiguous", candidates=ids, code="E_TAXONOMY_AMBIGUOUS",
                                 message=f"Multiple exact NCBI scientific-name matches for {name!r}")
        record = fetch_taxonomy_context(ids[0], email, api_key, limiter=limiter,
                                        events=events, attempts=attempts,
                                        opener=opener, sleep=sleep)
        code, context_error, rank_rule = validate_taxonomy_context(
            record, name, expected_rank=expected_rank, ancestor_names=ancestor_names
        )
        if context_error:
            return LookupOutcome("context_mismatch", code=code, message=context_error,
                                 expected_rank=expected_rank, returned_rank=record["rank"])
        return LookupOutcome("resolved", taxid=ids[0], expected_rank=expected_rank,
                             returned_rank=record["rank"], rank_rule=rank_rule)
    except SafeRequestFailure as exc:
        return LookupOutcome("request_failed", code=exc.code, message=exc.safe_message)
    except InvalidResponse as exc:
        return LookupOutcome("response_invalid", code="E_NCBI_RESPONSE", message=exc.safe_message)


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
    parser.add_argument("--diagnostics-dir",
                        help="Pre-created execution-only directory for durable operational diagnostics")
    args = parser.parse_args()

    events = None
    try:
        if args.online_preflight and not args.validate_only:
            raise ValueError(
                "E_ONLINE_PREFLIGHT_MODE: --online-preflight is only valid with --validate-only; "
                "remove it for normal one-pass refresh execution."
            )
        if args.validate_only and args.diagnostics_dir:
            raise ValueError("--diagnostics-dir is not permitted with --validate-only.")
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
        events = EventWriter(args.diagnostics_dir, args.transaction_id)
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
        if not os.path.isfile(args.cache):
            raise ValueError(
                f"Taxonomy cache must already exist as a regular file: {args.cache}"
            )
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
            taxonomy_mismatches = []
            lookup_outcomes = []
            rank_compatibility = []
            cache_updated = False

            if args.mode == "refresh" and unresolved:
                email = os.environ.get(args.email_env, "").strip()
                if not email:
                    raise ValueError(f"Environment variable {args.email_env!r} is required for refresh mode.")
                if args.validate_only and not args.online_preflight:
                    print(
                        json.dumps({"taxonomy_resolution": "pending_online",
                                    "pending_count": len(unresolved)}, sort_keys=True),
                        flush=True,
                    )
                    return 0
                api_key = os.environ.get(args.api_key_env, "").strip() or None
                limiter = RequestLimiter(0.12 if api_key else 0.35)
                name_results = {}
                fatal_outcome = None
                total_unresolved = len(unresolved)
                for completed, item in enumerate(unresolved, start=1):
                    name = item["name"]
                    expected_rank = EXPECTED_RANKS[item["depth"] - 1] if item["depth"] <= len(EXPECTED_RANKS) else None
                    ancestor_names = [part for part in item["path"].split(";")[:-1]
                                      if part.strip().lower() not in PLACEHOLDER_NAMES]
                    query_key = (name, expected_rank, tuple(ancestor_names))
                    if query_key not in name_results:
                        name_results[query_key] = query_exact_scientific_name(
                            name, email, api_key, expected_rank=expected_rank,
                            ancestor_names=ancestor_names, limiter=limiter, events=events
                        )
                    outcome = name_results[query_key]
                    outcome_record = {
                        "path": item["path"], "name": name, "status": outcome.status,
                        "code": outcome.code, "message": outcome.message,
                    }
                    if outcome.expected_rank is not None:
                        outcome_record["expected_rank"] = outcome.expected_rank
                    if outcome.returned_rank is not None:
                        outcome_record["returned_rank"] = outcome.returned_rank
                    if outcome.rank_rule is not None:
                        outcome_record["rank_rule"] = outcome.rank_rule
                    lookup_outcomes.append(outcome_record)
                    events.emit("lookup_progress", phase="taxonomy_resolution", name=name,
                                code=outcome.code, outcome=outcome.status,
                                message=outcome.message,
                                expected_rank=outcome.expected_rank,
                                returned_rank=outcome.returned_rank,
                                rank_rule=outcome.rank_rule,
                                completed=completed, total=total_unresolved)
                    if outcome.status in FATAL_OUTCOMES:
                        query_failures.append(outcome_record)
                        fatal_outcome = outcome
                        events.emit("lookup_failure", phase="taxonomy_resolution", name=name,
                                    code=outcome.code, message=outcome.message,
                                    completed=completed, total=total_unresolved)
                        events.write_failure_summary(outcome, completed, total_unresolved)
                        break
                    if outcome.status == "context_mismatch":
                        taxonomy_mismatches.append(outcome_record)
                    elif outcome.status == "ambiguous":
                        ambiguous_queries.append({
                            "path": item["path"], "name": name, "taxids": list(outcome.candidates)
                        })
                    elif outcome.status == "resolved":
                        cache[item["path"]] = outcome.taxid
                        resolution_sources[item["path"]] = "ncbi_refresh"
                        if outcome.rank_rule != "exact":
                            rank_compatibility.append({
                                "path": item["path"], "name": name,
                                "expected_rank": outcome.expected_rank,
                                "returned_rank": outcome.returned_rank,
                                "rule": outcome.rank_rule,
                            })

                unresolved = find_unresolved(abundance_paths, cache)
                if fatal_outcome is not None:
                    print(
                        f"[taxonomy] ERROR: {fatal_outcome.code}: "
                        f"{EventWriter._safe_text(fatal_outcome.message)}; "
                        "source cache preserved.", file=sys.stderr, flush=True
                    )
                    return 1
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
                "taxonomy_mismatches": taxonomy_mismatches,
                "lookup_outcomes": lookup_outcomes,
                "rank_compatibility": rank_compatibility,
                "resolution_source_counts": source_counts,
            }
            atomic_write_json(args.provenance, provenance)

            if query_failures:
                print(f"[taxonomy] ERROR: {len(query_failures)} NCBI query failure(s); source cache preserved.", file=sys.stderr)
                return 1
            if args.unresolved_policy == "error" and unresolved:
                events.write_failure_summary(
                    LookupOutcome("context_mismatch", code="E_TAXONOMY_UNRESOLVED",
                                  message=f"{len(unresolved)} taxonomy nodes remain unresolved"),
                    len(lookup_outcomes), len(lookup_outcomes)
                )
                print(f"[taxonomy] ERROR: {len(unresolved)} taxonomy nodes remain unresolved.", file=sys.stderr)
                return 1
            if args.mode == "refresh":
                current_cache_sha = compute_sha256(args.cache)
                if current_cache_sha != cache_sha_before:
                    events.write_failure_summary(
                        LookupOutcome("request_failed", code="E_TAXONOMY_CACHE_CHANGED",
                                      message="source taxonomy cache changed during query execution"),
                        len(lookup_outcomes), len(lookup_outcomes)
                    )
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
        if events is not None:
            outcome = LookupOutcome("response_invalid", code="E_TAXONOMY_EXECUTION",
                                    message=EventWriter._safe_text(exc))
            events.emit("execution_failure", code=outcome.code, message=outcome.message)
            events.write_failure_summary(outcome, 0, 0)
        print(f"[taxonomy] ERROR: {EventWriter._safe_text(exc)}", file=sys.stderr)
        return 1
    finally:
        if events is not None:
            events.close()


if __name__ == "__main__":
    sys.exit(main())
