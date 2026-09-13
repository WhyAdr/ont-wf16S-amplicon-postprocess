#!/usr/bin/env python3
"""Build deterministic offline Kraken-report explorer artifacts.

The ``.kreport`` remains the official upstream Pavian interoperability
artifact. This file implements an original static viewer and does not vendor
or relabel the upstream Pavian Shiny application.
"""

from __future__ import annotations

import argparse
import base64
import html
import json
import math
import os
import re
import sys
import tempfile
import unicodedata
from decimal import Decimal, InvalidOperation
from pathlib import Path


MAX_SAFE_INTEGER = 9007199254740991
RANK_BY_DEPTH = {1: "D", 2: "K", 3: "P", 4: "C", 5: "O", 6: "F", 7: "G", 8: "S"}
INTEGER_RE = re.compile(r"(?:0|[1-9][0-9]*)\Z")
DECIMAL_RE = re.compile(r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?\Z")
CONTROL_RE = re.compile(r"[\x00-\x1f\x7f]")
RESOLUTION_HEADER = [
    "SampleID", "Depth", "RankCode", "NodeName", "TaxonPath", "TaxID",
    "Status", "ResolutionSource",
]


class ViewerError(ValueError):
    """A concise, user-actionable viewer input or rendering error."""


def _fail(message: str) -> None:
    raise ViewerError(message)


def _read_utf8(path: Path, label: str) -> str:
    if path.is_symlink() or not path.is_file():
        _fail(f"{label} is not an existing regular file: {path}")
    try:
        return path.read_bytes().decode("utf-8", errors="strict")
    except (OSError, UnicodeError) as exc:
        _fail(f"could not read {label} '{path}': {exc}")
    raise AssertionError("unreachable")


def parse_nonnegative_integer(value: str, label: str) -> int:
    if not isinstance(value, str) or not INTEGER_RE.fullmatch(value):
        _fail(f"{label} must be a canonical non-negative integer")
    parsed = int(value)
    if parsed > MAX_SAFE_INTEGER:
        _fail(f"{label} exceeds the JavaScript-safe integer range")
    return parsed


def parse_percentage(value: str, label: str) -> float:
    if not isinstance(value, str) or not DECIMAL_RE.fullmatch(value):
        _fail(f"{label} must be a finite decimal percentage")
    try:
        decimal = Decimal(value)
    except InvalidOperation as exc:
        _fail(f"{label} must be a finite decimal percentage")
        raise AssertionError("unreachable") from exc
    parsed = float(decimal)
    if not math.isfinite(parsed) or parsed < 0.0 or parsed > 100.0:
        _fail(f"{label} must be between 0 and 100")
    return parsed


def validate_label(value: str, label: str) -> str:
    if not isinstance(value, str) or not value:
        _fail(f"{label} must be non-empty")
    if CONTROL_RE.search(value) or any(unicodedata.category(char) == "Cc" for char in value):
        _fail(f"{label} may not contain control characters")
    if value.startswith(" "):
        _fail(f"{label} may not begin with ambiguous indentation")
    return value


def validate_percentage(value: float, clade: int, total: int, label: str) -> None:
    expected = 100.0 * clade / total
    if abs(value - expected) > 0.005000001:
        _fail(f"{label} percentage is inconsistent with its clade count")


def parse_kreport(path: Path, expected_total: int) -> dict[str, object]:
    """Parse and independently validate the fixed six-column report contract."""

    if expected_total <= 0:
        _fail("expected total reads must be greater than zero")
    text = _read_utf8(path, "kreport")
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines:
        _fail("kreport contains no records")

    parsed: list[dict[str, object]] = []
    active_stack: list[str] = []
    for row_number, raw_line in enumerate(lines, start=1):
        if raw_line.endswith("\r"):
            raw_line = raw_line[:-1]
        if not raw_line:
            _fail(f"kreport line {row_number}: blank lines are not allowed")
        fields = raw_line.split("\t")
        if len(fields) != 6:
            _fail(f"kreport line {row_number}: expected exactly six tab-separated fields")
        percentage = parse_percentage(fields[0], f"kreport line {row_number} percentage")
        clade = parse_nonnegative_integer(fields[1], f"kreport line {row_number} clade")
        direct = parse_nonnegative_integer(fields[2], f"kreport line {row_number} direct")
        rank = fields[3]
        taxid = fields[4]
        raw_name = fields[5]

        if row_number == 1:
            if (rank, taxid, raw_name) != ("U", "0", "unclassified"):
                _fail("kreport first record must be U, TaxID 0, named unclassified")
            if clade != direct:
                _fail("kreport unclassified clade and direct counts must match")
            validate_percentage(percentage, clade, expected_total, "kreport unclassified")
            parsed.append({
                "percentage": percentage, "clade": clade, "direct": direct,
                "rank": rank, "taxid": taxid, "name": raw_name, "depth": 0,
                "path": "", "parent_path": "",
            })
            continue

        if row_number == 2:
            if (rank, taxid, raw_name) != ("R", "1", "root"):
                _fail("kreport second record must be R, TaxID 1, named root")
            if direct != 0:
                _fail("kreport classified root direct count must be zero")
            validate_percentage(percentage, clade, expected_total, "kreport root")
            parsed.append({
                "percentage": percentage, "clade": clade, "direct": direct,
                "rank": rank, "taxid": taxid, "name": raw_name, "depth": 0,
                "path": "", "parent_path": "",
            })
            continue

        if rank not in set(RANK_BY_DEPTH.values()):
            _fail(f"kreport line {row_number}: invalid rank code '{rank}'")
        if raw_name.startswith("\t"):
            _fail(f"kreport line {row_number}: taxonomy name begins with a tab")
        leading_spaces = len(raw_name) - len(raw_name.lstrip(" "))
        if leading_spaces % 2:
            _fail(f"kreport line {row_number}: indentation must use two ASCII spaces per depth")
        depth = leading_spaces // 2
        name = raw_name[leading_spaces:]
        validate_label(name, f"kreport line {row_number} taxonomy name")
        if depth < 1 or depth > 8:
            _fail(f"kreport line {row_number}: taxonomy depth is outside 1..8")
        if RANK_BY_DEPTH[depth] != rank:
            _fail(f"kreport line {row_number}: rank code does not match depth")
        if active_stack and depth > len(active_stack) + 1:
            _fail(f"kreport line {row_number}: taxonomy depth jumps by more than one level")
        if not INTEGER_RE.fullmatch(taxid) or int(taxid) > MAX_SAFE_INTEGER:
            _fail(f"kreport line {row_number}: TaxID must be a canonical safe integer")

        if depth > len(active_stack) + 1:
            _fail(f"kreport line {row_number}: classified node has no active parent")
        active_stack = active_stack[: depth - 1]
        parent_path = active_stack[-1] if active_stack else ""
        path_key = f"{parent_path};{name}" if parent_path else name
        if any(str(item["path"]) == path_key for item in parsed[2:]):
            _fail(f"kreport line {row_number}: duplicate taxonomy path '{path_key}'")
        validate_percentage(percentage, clade, expected_total, f"kreport line {row_number}")
        parsed.append({
            "percentage": percentage, "clade": clade, "direct": direct,
            "rank": rank, "taxid": taxid, "name": name, "depth": depth,
            "path": path_key, "parent_path": parent_path,
        })
        active_stack.append(path_key)

    if len(parsed) < 2:
        _fail("kreport must contain both U and R records")
    unclassified = int(parsed[0]["clade"])
    classified = int(parsed[1]["clade"])
    total = unclassified + classified
    if total <= 0:
        _fail("kreport total reads must be greater than zero")
    if total != expected_total:
        _fail(f"kreport total {total} does not match expected total {expected_total}")

    nodes = parsed[2:]
    if sum(int(node["direct"]) for node in nodes) != classified:
        _fail("kreport classified direct counts do not equal root clade count")
    if sum(int(node["clade"]) for node in nodes if int(node["depth"]) == 1) != classified:
        _fail("kreport top-level clade counts do not equal root clade count")
    by_parent: dict[str, list[dict[str, object]]] = {}
    for node in nodes:
        by_parent.setdefault(str(node["parent_path"]), []).append(node)
    for node in nodes:
        child_sum = sum(int(child["clade"]) for child in by_parent.get(str(node["path"]), []))
        if int(node["clade"]) != int(node["direct"]) + child_sum:
            _fail(f"kreport clade arithmetic failed at '{node['path']}'")

    expected_order: list[str] = []

    def visit(parent_path: str) -> None:
        children = sorted(by_parent.get(parent_path, []),
                          key=lambda node: (-int(node["clade"]), str(node["name"])))
        for child in children:
            expected_order.append(str(child["path"]))
            visit(str(child["path"]))

    visit("")
    actual_order = [str(node["path"]) for node in nodes]
    if actual_order != expected_order:
        _fail("kreport nodes are not in the required deterministic DFS order")

    return {
        "total": total,
        "classified": classified,
        "unclassified": unclassified,
        "nodes": [
            {
                "path": str(node["path"]),
                "parent_path": str(node["parent_path"]),
                "name": str(node["name"]),
                "depth": int(node["depth"]),
                "rank_code": str(node["rank"]),
                "taxid": str(node["taxid"]),
                "direct": int(node["direct"]),
                "clade": int(node["clade"]),
            }
            for node in nodes
        ],
    }


def load_resolution_rows(path: Path, sample_id: str) -> dict[str, dict[str, object]]:
    text = _read_utf8(path, "taxonomy resolution TSV")
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines or lines[0].rstrip("\r").split("\t") != RESOLUTION_HEADER:
        _fail("taxonomy resolution TSV has an unexpected header")
    selected: dict[str, dict[str, object]] = {}
    seen: set[tuple[str, str]] = set()
    for row_number, raw_line in enumerate(lines[1:], start=2):
        if raw_line.endswith("\r"):
            raw_line = raw_line[:-1]
        if not raw_line:
            _fail(f"taxonomy resolution TSV line {row_number}: blank lines are not allowed")
        fields = raw_line.split("\t")
        if len(fields) != len(RESOLUTION_HEADER):
            _fail(f"taxonomy resolution TSV line {row_number}: expected eight fields")
        if any(CONTROL_RE.search(field) for field in fields):
            _fail(f"taxonomy resolution TSV line {row_number}: control characters are not allowed")
        row_sample = validate_label(fields[0], f"taxonomy resolution TSV line {row_number} SampleID")
        depth = parse_nonnegative_integer(fields[1], f"taxonomy resolution TSV line {row_number} Depth")
        if depth < 1 or depth > 8:
            _fail(f"taxonomy resolution TSV line {row_number}: Depth is outside 1..8")
        rank_code = fields[2]
        if rank_code != RANK_BY_DEPTH[depth]:
            _fail(f"taxonomy resolution TSV line {row_number}: RankCode does not match Depth")
        node_name = validate_label(fields[3], f"taxonomy resolution TSV line {row_number} NodeName")
        taxon_path = validate_label(fields[4], f"taxonomy resolution TSV line {row_number} TaxonPath")
        taxid_int = parse_nonnegative_integer(fields[5], f"taxonomy resolution TSV line {row_number} TaxID")
        status = fields[6]
        if status not in {"Resolved", "Unresolved", "Conflicted"}:
            _fail(f"taxonomy resolution TSV line {row_number}: invalid Status")
        resolution_source = validate_label(
            fields[7], f"taxonomy resolution TSV line {row_number} ResolutionSource"
        )
        key = (row_sample, taxon_path)
        if key in seen:
            _fail(f"taxonomy resolution TSV line {row_number}: duplicate SampleID/TaxonPath")
        seen.add(key)
        if row_sample == sample_id:
            selected[taxon_path] = {
                "sample_id": row_sample,
                "depth": depth,
                "rank_code": rank_code,
                "name": node_name,
                "path": taxon_path,
                "taxid": str(taxid_int),
                "status": status,
                "resolution_source": resolution_source,
            }
    return selected


def merge_resolution(report: dict[str, object], resolution: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    nodes = report["nodes"]
    assert isinstance(nodes, list)
    report_paths = {str(node["path"]) for node in nodes}
    resolution_paths = set(resolution)
    if report_paths != resolution_paths:
        missing = sorted(report_paths - resolution_paths)
        extra = sorted(resolution_paths - report_paths)
        _fail(f"resolution sidecar/report path mismatch; missing={missing}, extra={extra}")
    status_map = {"Resolved": "resolved", "Unresolved": "unresolved", "Conflicted": "conflicted"}
    merged: list[dict[str, object]] = []
    for node in nodes:
        assert isinstance(node, dict)
        record = resolution[str(node["path"])]
        for field in ("depth", "rank_code", "name", "path", "taxid"):
            if record[field] != node[field]:
                _fail(f"resolution sidecar disagrees with kreport at '{node['path']}' field '{field}'")
        merged.append({
            "path": str(node["path"]),
            "parent_path": str(node["parent_path"]),
            "name": str(node["name"]),
            "depth": int(node["depth"]),
            "rank_code": str(node["rank_code"]),
            "taxid": str(node["taxid"]),
            "direct": int(node["direct"]),
            "clade": int(node["clade"]),
            "status": status_map[str(record["status"])],
            "resolution_source": str(record["resolution_source"]),
        })
    return merged


def canonical_json_bytes(payload: dict[str, object]) -> bytes:
    text = json.dumps(payload, ensure_ascii=False, sort_keys=True,
                      separators=(",", ":"), allow_nan=False)
    return text.encode("utf-8") + b"\n"


def build_html(payload_bytes: bytes, sample_id: str) -> bytes:
    payload_b64 = base64.b64encode(payload_bytes).decode("ascii")
    safe_title = html.escape(sample_id, quote=True)
    csp = ("default-src 'none'; script-src 'unsafe-inline'; style-src 'unsafe-inline'; "
           "img-src data:; connect-src 'none'; font-src 'none'; object-src 'none'; "
           "base-uri 'none'; form-action 'none'")
    script = r"""
(function () {
  const encoded = document.getElementById("wf16s-payload").textContent.trim();
  const binary = atob(encoded);
  const bytes = Uint8Array.from(binary, function (character) { return character.charCodeAt(0); });
  const payload = JSON.parse(new TextDecoder("utf-8").decode(bytes));
  const nodes = payload.nodes;
  const statusLabel = {resolved: "Resolved", unresolved: "Unresolved", conflicted: "Conflicted"};
  document.getElementById("total").textContent = String(payload.totals.total);
  document.getElementById("classified").textContent = String(payload.totals.classified);
  document.getElementById("unclassified").textContent = String(payload.totals.unclassified);
  const rank = document.getElementById("rank");
  ["all", "D", "K", "P", "C", "O", "F", "G", "S"].forEach(function (value) {
    const option = document.createElement("option");
    option.value = value;
    option.textContent = value === "all" ? "All ranks" : value;
    rank.appendChild(option);
  });
  const search = document.getElementById("search");
  const useClade = document.getElementById("use-clade");
  const tree = document.getElementById("tree");
  const body = document.getElementById("rows");
  function visible(node) {
    const needle = search.value.toLocaleLowerCase();
    return (!needle || node.name.toLocaleLowerCase().indexOf(needle) >= 0 || node.path.toLocaleLowerCase().indexOf(needle) >= 0)
      && (rank.value === "all" || node.rank_code === rank.value);
  }
  function render() {
    while (tree.firstChild) tree.removeChild(tree.firstChild);
    while (body.firstChild) body.removeChild(body.firstChild);
    nodes.filter(visible).forEach(function (node) {
      const amount = useClade.checked ? node.clade : node.direct;
      const li = document.createElement("li");
      li.style.paddingLeft = String(node.depth * 1.25) + "em";
      li.textContent = node.name + " [" + statusLabel[node.status] + "] — " + String(amount) + " (" + String(node.taxid) + ")";
      tree.appendChild(li);
      const row = document.createElement("tr");
      [node.path, node.rank_code, node.taxid, String(node.direct), String(node.clade), statusLabel[node.status]].forEach(function (value) {
        const cell = document.createElement("td");
        cell.textContent = value;
        row.appendChild(cell);
      });
      body.appendChild(row);
    });
  }
  [search, rank, useClade].forEach(function (control) { control.addEventListener("input", render); });
  const download = document.getElementById("download");
  download.href = URL.createObjectURL(new Blob([bytes], {type: "application/json"}));
  render();
}());
"""
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta http-equiv="Content-Security-Policy" content="{csp}">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Kraken report explorer — {safe_title}</title>
<style>
body {{ font-family: system-ui, sans-serif; margin: 2rem; color: #202124; }}
header {{ border-bottom: 1px solid #d8dbe0; margin-bottom: 1rem; }}
.summary {{ display: flex; gap: 1.5rem; flex-wrap: wrap; }}
label {{ margin-right: 1rem; }}
table {{ border-collapse: collapse; width: 100%; margin-top: 1rem; }}
th, td {{ border: 1px solid #d8dbe0; padding: .35rem .5rem; text-align: left; }}
th {{ background: #f4f6f8; }}
ul {{ list-style: none; padding-left: 0; }}
</style>
</head>
<body>
<header><h1>Kraken report explorer — {safe_title}</h1>
<p>Original offline viewer; the sibling <code>.kreport</code> remains the official Pavian interoperability artifact.</p></header>
<section class="summary" aria-label="Read totals">
<span>Total: <strong id="total"></strong></span>
<span>Classified: <strong id="classified"></strong></span>
<span>Unclassified: <strong id="unclassified"></strong></span>
</section>
<p><label>Search <input id="search" type="search"></label>
<label>Rank <select id="rank"></select></label>
<label><input id="use-clade" type="checkbox" checked> use clade counts</label>
<a id="download" download="{safe_title}.pavian.json">Download normalized JSON</a></p>
<details open><summary>Taxonomy tree</summary><ul id="tree"></ul></details>
<details open><summary>Node table</summary><table><thead><tr><th>Taxon path</th><th>Rank</th><th>TaxID</th><th>Direct</th><th>Clade</th><th>Status</th></tr></thead><tbody id="rows"></tbody></table></details>
<script id="wf16s-payload" type="application/octet-stream">{payload_b64}</script>
<script>{script}</script>
</body>
</html>
"""
    return document.encode("utf-8")


def atomic_write_bytes(path: Path, data: bytes) -> None:
    path = Path(path)
    if path.exists() or path.is_symlink():
        _fail(f"refusing to overwrite existing builder output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.tmp-", dir=str(path.parent))
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def read_renderer_version() -> str:
    version_path = Path(__file__).resolve().parents[2] / "VERSION"
    version = _read_utf8(version_path, "pipeline VERSION").rstrip("\n")
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", version):
        _fail("pipeline VERSION must contain one SemVer value")
    return version


def build_payload(kreport: Path, resolution_tsv: Path, sample_id: str,
                  expected_total: int) -> dict[str, object]:
    if not sample_id or CONTROL_RE.search(sample_id):
        _fail("sample_id must be non-empty and free of control characters")
    report = parse_kreport(kreport, expected_total)
    resolution = load_resolution_rows(resolution_tsv, sample_id)
    nodes = merge_resolution(report, resolution)
    return {
        "schema_version": 1,
        "sample_id": sample_id,
        "renderer": "builtin_kraken_report_explorer",
        "renderer_version": read_renderer_version(),
        "official_pavian_compatibility": "kraken_report_input_contract_only",
        "count_model": "direct abundance-table taxon counts plus canonical unclassified count",
        "denominator": "TotalReads",
        "classified_definition": "sum of direct positive-count classified taxonomy rows",
        "totals": {
            "total": int(report["total"]),
            "classified": int(report["classified"]),
            "unclassified": int(report["unclassified"]),
        },
        "nodes": nodes,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build the offline Kraken-report explorer")
    parser.add_argument("--kreport", required=True, type=Path)
    parser.add_argument("--resolution-tsv", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--json-out", required=True, type=Path)
    parser.add_argument("--html-out", type=Path)
    parser.add_argument("--expected-total", required=True)
    args = parser.parse_args(argv)
    try:
        expected_total = parse_nonnegative_integer(args.expected_total, "expected-total")
        payload = build_payload(args.kreport, args.resolution_tsv, args.sample_id, expected_total)
        payload_bytes = canonical_json_bytes(payload)
        atomic_write_bytes(args.json_out, payload_bytes)
        if args.html_out is not None:
            atomic_write_bytes(args.html_out, build_html(payload_bytes, args.sample_id))
        return 0
    except Exception as exc:
        print(f"kraken_report_viewer: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
