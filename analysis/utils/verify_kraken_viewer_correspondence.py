#!/usr/bin/env python3
"""Independently verify kreport -> sidecar -> JSON -> HTML correspondence.

This verifier deliberately has its own report and TSV parsers. It does not
import the producer module, so a producer/parser bug cannot certify itself.
"""

from __future__ import annotations

import argparse
import base64
import json
import math
import re
import sys
import unicodedata
from decimal import Decimal, InvalidOperation
from pathlib import Path


MAX_SAFE_INTEGER = 9007199254740991
RANK_BY_DEPTH = {1: "D", 2: "K", 3: "P", 4: "C", 5: "O", 6: "F", 7: "G", 8: "S"}
INTEGER_RE = re.compile(r"(?:0|[1-9][0-9]*)\Z")
DECIMAL_RE = re.compile(r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?\Z")
CONTROL_RE = re.compile(r"[\x00-\x1f\x7f]")
HEADER = [
    "SampleID", "Depth", "RankCode", "NodeName", "TaxonPath", "TaxID",
    "Status", "ResolutionSource",
]


class CorrespondenceError(ValueError):
    """A correspondence or security-contract failure."""


def fail(message: str) -> None:
    raise CorrespondenceError(message)


def read_file(path: Path, label: str) -> str:
    if path.is_symlink() or not path.is_file():
        fail(f"{label} is not an existing regular file: {path}")
    try:
        return path.read_bytes().decode("utf-8", errors="strict")
    except (OSError, UnicodeError) as exc:
        fail(f"could not read {label} '{path}': {exc}")
    raise AssertionError("unreachable")


def integer(value: str, label: str) -> int:
    if not isinstance(value, str) or not INTEGER_RE.fullmatch(value):
        fail(f"{label} is not a canonical non-negative integer")
    result = int(value)
    if result > MAX_SAFE_INTEGER:
        fail(f"{label} exceeds the JavaScript-safe integer range")
    return result


def percentage(value: str, label: str) -> float:
    if not isinstance(value, str) or not DECIMAL_RE.fullmatch(value):
        fail(f"{label} is not a finite decimal percentage")
    try:
        result = float(Decimal(value))
    except (InvalidOperation, ValueError) as exc:
        fail(f"{label} is not a finite decimal percentage")
        raise AssertionError("unreachable") from exc
    if not math.isfinite(result) or result < 0 or result > 100:
        fail(f"{label} is outside 0..100")
    return result


def label(value: str, label_name: str) -> str:
    if not value or CONTROL_RE.search(value) or any(unicodedata.category(c) == "Cc" for c in value):
        fail(f"{label_name} is empty or contains a control character")
    if value.startswith(" "):
        fail(f"{label_name} begins with ambiguous indentation")
    return value


def check_percent(value: float, clade: int, total: int, label_name: str) -> None:
    if abs(value - (100.0 * clade / total)) > 0.005000001:
        fail(f"{label_name} percentage does not agree with its count")


def parse_report_independently(path: Path, expected_total: int) -> dict[str, object]:
    if expected_total <= 0:
        fail("expected total must be positive")
    text = read_file(path, "kreport")
    raw_lines = text.split("\n")
    if raw_lines and raw_lines[-1] == "":
        raw_lines.pop()
    if not raw_lines:
        fail("kreport is empty")

    rows: list[dict[str, object]] = []
    stack: list[str] = []
    for line_no, line in enumerate(raw_lines, start=1):
        if line.endswith("\r"):
            line = line[:-1]
        if not line:
            fail(f"kreport line {line_no} is blank")
        columns = line.split("\t")
        if len(columns) != 6:
            fail(f"kreport line {line_no} must have six columns")
        pct = percentage(columns[0], f"kreport line {line_no} percentage")
        clade = integer(columns[1], f"kreport line {line_no} clade")
        direct = integer(columns[2], f"kreport line {line_no} direct")
        rank, taxid, raw_name = columns[3], columns[4], columns[5]
        if line_no == 1:
            if (rank, taxid, raw_name) != ("U", "0", "unclassified") or clade != direct:
                fail("invalid unclassified record")
            rows.append({"path": "", "parent": "", "name": raw_name, "depth": 0,
                         "rank": rank, "taxid": taxid, "direct": direct, "clade": clade})
            continue
        if line_no == 2:
            if (rank, taxid, raw_name) != ("R", "1", "root") or direct != 0:
                fail("invalid classified-root record")
            rows.append({"path": "", "parent": "", "name": raw_name, "depth": 0,
                         "rank": rank, "taxid": taxid, "direct": direct, "clade": clade})
            continue

        if rank not in RANK_BY_DEPTH.values() or not INTEGER_RE.fullmatch(taxid):
            fail(f"invalid rank or TaxID at kreport line {line_no}")
        if int(taxid) > MAX_SAFE_INTEGER:
            fail(f"TaxID exceeds safe range at kreport line {line_no}")
        spaces = len(raw_name) - len(raw_name.lstrip(" "))
        if raw_name.startswith("\t") or spaces % 2:
            fail(f"invalid indentation at kreport line {line_no}")
        depth = spaces // 2
        name = label(raw_name[spaces:], f"kreport line {line_no} name")
        if depth not in RANK_BY_DEPTH or RANK_BY_DEPTH[depth] != rank:
            fail(f"rank/depth mismatch at kreport line {line_no}")
        if depth > len(stack) + 1:
            fail(f"depth jump at kreport line {line_no}")
        check_percent(pct, clade, expected_total, f"kreport line {line_no}")
        stack = stack[: depth - 1]
        parent = stack[-1] if stack else ""
        path_key = f"{parent};{name}" if parent else name
        if any(str(row["path"]) == path_key for row in rows[2:]):
            fail(f"duplicate path at kreport line {line_no}")
        rows.append({"path": path_key, "parent": parent, "name": name, "depth": depth,
                     "rank": rank, "taxid": str(int(taxid)), "direct": direct, "clade": clade})
        stack.append(path_key)

    if len(rows) < 2:
        fail("kreport must include U and R records")
    unclassified = int(rows[0]["clade"])
    classified = int(rows[1]["clade"])
    total = unclassified + classified
    if total <= 0 or total != expected_total:
        fail("kreport total does not match expected total")
    check_percent(percentage(raw_lines[0].rstrip("\r").split("\t")[0], "U percentage"),
                  unclassified, total, "U")
    check_percent(percentage(raw_lines[1].rstrip("\r").split("\t")[0], "R percentage"),
                  classified, total, "R")
    nodes = rows[2:]
    if sum(int(row["direct"]) for row in nodes) != classified:
        fail("classified direct total mismatch")
    children: dict[str, list[dict[str, object]]] = {}
    for row in nodes:
        children.setdefault(str(row["parent"]), []).append(row)
    if sum(int(row["clade"]) for row in children.get("", [])) != classified:
        fail("top-level clade total mismatch")
    for row in nodes:
        if int(row["clade"]) != int(row["direct"]) + sum(
                int(child["clade"]) for child in children.get(str(row["path"]), [])):
            fail(f"clade arithmetic mismatch at {row['path']}")
    order: list[str] = []

    def walk(parent: str) -> None:
        for child in sorted(children.get(parent, []),
                            key=lambda item: (-int(item["clade"]), str(item["name"]))):
            order.append(str(child["path"]))
            walk(str(child["path"]))

    walk("")
    if order != [str(row["path"]) for row in nodes]:
        fail("kreport order is not deterministic DFS order")
    return {
        "total": total, "classified": classified, "unclassified": unclassified,
        "nodes": [
            {"path": str(row["path"]), "parent_path": str(row["parent"]),
             "name": str(row["name"]), "depth": int(row["depth"]),
             "rank_code": str(row["rank"]), "taxid": str(row["taxid"]),
             "direct": int(row["direct"]), "clade": int(row["clade"])}
            for row in nodes
        ],
    }


def parse_resolution_independently(path: Path, sample_id: str) -> dict[str, dict[str, object]]:
    text = read_file(path, "taxonomy resolution TSV")
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines or lines[0].rstrip("\r").split("\t") != HEADER:
        fail("taxonomy resolution TSV header mismatch")
    result: dict[str, dict[str, object]] = {}
    seen: set[tuple[str, str]] = set()
    for line_no, line in enumerate(lines[1:], start=2):
        if line.endswith("\r"):
            line = line[:-1]
        if not line:
            fail(f"taxonomy resolution TSV line {line_no} is blank")
        fields = line.split("\t")
        if len(fields) != 8 or any(CONTROL_RE.search(field) for field in fields):
            fail(f"taxonomy resolution TSV line {line_no} is malformed")
        row_sample = label(fields[0], f"resolution line {line_no} SampleID")
        depth = integer(fields[1], f"resolution line {line_no} Depth")
        if depth not in RANK_BY_DEPTH or fields[2] != RANK_BY_DEPTH[depth]:
            fail(f"resolution line {line_no} rank/depth mismatch")
        name = label(fields[3], f"resolution line {line_no} NodeName")
        taxon_path = label(fields[4], f"resolution line {line_no} TaxonPath")
        taxid = str(integer(fields[5], f"resolution line {line_no} TaxID"))
        if fields[6] not in {"Resolved", "Unresolved", "Conflicted"}:
            fail(f"resolution line {line_no} has an invalid status")
        source = label(fields[7], f"resolution line {line_no} ResolutionSource")
        key = (row_sample, taxon_path)
        if key in seen:
            fail(f"resolution line {line_no} duplicates SampleID/TaxonPath")
        seen.add(key)
        if row_sample == sample_id:
            result[taxon_path] = {"depth": depth, "rank_code": fields[2], "name": name,
                                  "path": taxon_path, "taxid": taxid, "status": fields[6],
                                  "resolution_source": source}
    return result


def expected_payload(report: dict[str, object], resolution: dict[str, dict[str, object]],
                     sample_id: str, renderer_version: str) -> dict[str, object]:
    nodes = report["nodes"]
    assert isinstance(nodes, list)
    if {str(node["path"]) for node in nodes} != set(resolution):
        fail("resolution sidecar/report path sets differ")
    status_names = {"Resolved": "resolved", "Unresolved": "unresolved", "Conflicted": "conflicted"}
    merged = []
    for node in nodes:
        row = resolution[str(node["path"])]
        for field in ("depth", "rank_code", "name", "path", "taxid"):
            if row[field] != node[field]:
                fail(f"resolution sidecar disagrees at {node['path']} field {field}")
        merged.append({**node, "status": status_names[str(row["status"])],
                       "resolution_source": str(row["resolution_source"])})
    return {
        "schema_version": 1,
        "sample_id": sample_id,
        "renderer": "builtin_kraken_report_explorer",
        "renderer_version": renderer_version,
        "official_pavian_compatibility": "kraken_report_input_contract_only",
        "count_model": "direct abundance-table taxon counts plus canonical unclassified count",
        "denominator": "TotalReads",
        "classified_definition": "sum of direct positive-count classified taxonomy rows",
        "totals": {"total": int(report["total"]), "classified": int(report["classified"]),
                   "unclassified": int(report["unclassified"])},
        "nodes": merged,
    }


def canonical_bytes(value: dict[str, object]) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True,
                       separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")


def embedded_payload_bytes(html_text: str) -> bytes:
    matches = re.findall(
        r'<script\s+id="wf16s-payload"\s+type="application/octet-stream">([A-Za-z0-9+/=]*)</script>',
        html_text,
    )
    if len(matches) != 1:
        fail("HTML must contain exactly one base64 payload element")
    try:
        return base64.b64decode(matches[0], validate=True)
    except (ValueError, TypeError) as exc:
        fail(f"HTML payload is not strict base64: {exc}")
        raise AssertionError("unreachable") from exc


def reject_network_capabilities(html_text: str) -> None:
    lowered = html_text.lower()
    required = ["default-src 'none'", "connect-src 'none'", "script-src 'unsafe-inline'"]
    if any(item not in lowered for item in required):
        fail("HTML CSP does not enforce the offline viewer contract")
    forbidden = ["fetch(", "xmlhttprequest", "websocket", "eventsource", "sendbeacon",
                 "innerhtml", "<script src=", "<link href=", "<img src=", "http://", "https://"]
    for token in forbidden:
        if token in lowered:
            fail(f"HTML contains a forbidden network or unsafe DOM capability: {token}")


def read_renderer_version() -> str:
    version_file = Path(__file__).resolve().parents[2] / "VERSION"
    value = read_file(version_file, "pipeline VERSION").rstrip("\n")
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", value):
        fail("pipeline VERSION is not a SemVer")
    return value


def verify(kreport: Path, resolution_tsv: Path, sample_id: str, json_path: Path,
           html_path: Path | None, expected_total: int) -> None:
    report = parse_report_independently(kreport, expected_total)
    resolution = parse_resolution_independently(resolution_tsv, sample_id)
    expected = expected_payload(report, resolution, sample_id, read_renderer_version())
    json_bytes = json_path.read_bytes()
    try:
        actual = json.loads(json_bytes.decode("utf-8"), parse_constant=lambda value: fail(
            f"JSON contains non-finite value {value}"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"could not parse JSON output: {exc}")
    if actual != expected or json_bytes != canonical_bytes(expected):
        fail("normalized JSON is not the canonical payload for the supplied inputs")
    if html_path is not None:
        html_text = read_file(html_path, "HTML output")
        if embedded_payload_bytes(html_text) != json_bytes:
            fail("embedded HTML payload differs byte-for-byte from JSON output")
        reject_network_capabilities(html_text)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Verify Kraken viewer correspondence")
    parser.add_argument("--kreport", required=True, type=Path)
    parser.add_argument("--resolution-tsv", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--json", required=True, type=Path, dest="json_path")
    parser.add_argument("--html", type=Path)
    parser.add_argument("--expected-total", required=True)
    args = parser.parse_args(argv)
    try:
        expected_total = integer(args.expected_total, "expected-total")
        verify(args.kreport, args.resolution_tsv, args.sample_id, args.json_path,
               args.html, expected_total)
        return 0
    except Exception as exc:
        print(f"verify_kraken_viewer_correspondence: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
