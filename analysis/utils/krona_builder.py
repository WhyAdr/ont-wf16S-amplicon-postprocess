#!/usr/bin/env python3
"""Build deterministic, self-contained Krona-compatible HTML.

The pipeline writes direct-count lineage TSV files.  This module turns those
files into the small XML data model consumed by the vendored Krona 2.0
renderer and embeds the renderer and image assets in one offline HTML file.
Only the Python standard library is used so the builtin renderer is available
in both the locked R environment and a plain Python installation.
"""

from __future__ import annotations

import argparse
import base64
import binascii
import hashlib
import html
import json
import mimetypes
import os
from pathlib import Path, PurePosixPath
import re
import sys
import tempfile
import unicodedata
import xml.etree.ElementTree as ET


REPO_ROOT = Path(__file__).resolve().parents[2]
SEMVER_RE = re.compile(r"[0-9]+\.[0-9]+\.[0-9]+\n\Z")


def read_strict_semver(path: os.PathLike[str] | str) -> str:
    """Read the byte-strict pipeline version from the repository VERSION file."""

    try:
        text = Path(path).read_bytes().decode("utf-8")
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"could not read pipeline VERSION '{path}': {exc}") from exc
    if not SEMVER_RE.fullmatch(text):
        raise ValueError("VERSION must contain exactly one newline-terminated SemVer value")
    return text[:-1]


BUILDER_VERSION = read_strict_semver(REPO_ROOT / "VERSION")
VENDOR_TAG = "v2.8.1"
EXPECTED_VENDOR_FILES = (
    "LICENSE.txt",
    "src/krona-2.0.js",
    "img/favicon.ico",
    "img/hidden.png",
    "img/loading.gif",
    "img/logo-med.png",
)
INTEGER_RE = re.compile(r"(?:0|[1-9][0-9]*)\Z")
CONTROL_RE = re.compile(r"[\x00-\x1f\x7f]")
ATTRIBUTION_URL = "https://github.com/marbl/Krona/wiki"


class BuilderError(ValueError):
    """A concise, user-actionable input or renderer error."""


def _fail(message: str) -> None:
    raise BuilderError(message)


def _require_text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or CONTROL_RE.search(value):
        _fail(f"{label} must be a non-empty string without control characters")
    return value


def parse_nonnegative_integer(value: object, label: str) -> int:
    """Parse the canonical non-negative integer spelling used by Krona TSV."""

    if isinstance(value, bool):
        _fail(f"{label} must be a canonical non-negative integer")
    text = str(value)
    if not INTEGER_RE.fullmatch(text):
        _fail(f"{label} must be a canonical non-negative integer")
    return int(text)


def _validate_label(label: str, row_number: int) -> str:
    if not label or not label.strip():
        _fail(f"line {row_number}: taxonomy labels must be non-empty")
    if CONTROL_RE.search(label):
        _fail(f"line {row_number}: taxonomy labels may not contain control characters")
    if any(unicodedata.category(char) == "Cc" for char in label):
        _fail(f"line {row_number}: taxonomy labels may not contain control characters")
    return label


def parse_krona_input(input_path: os.PathLike[str] | str) -> list[tuple[int, tuple[str, ...]]]:
    """Read and validate a direct-count Krona TSV.

    Each non-empty record is ``magnitude<TAB>label<TAB>...``.  Empty records,
    fractional/negative magnitudes, malformed UTF-8, and control characters
    are rejected rather than silently normalized.
    """

    path = Path(input_path)
    if not path.is_file():
        _fail(f"input is not an existing regular file: {path}")
    try:
        text = path.read_text(encoding="utf-8", errors="strict")
    except (OSError, UnicodeError) as exc:
        _fail(f"could not read input '{path}': {exc}")

    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines:
        _fail("input contains no data lines")

    records: list[tuple[int, tuple[str, ...]]] = []
    for row_number, raw_line in enumerate(lines, start=1):
        if raw_line.endswith("\r"):
            raw_line = raw_line[:-1]
        if not raw_line:
            _fail(f"line {row_number}: blank lines are not allowed")
        fields = raw_line.split("\t")
        if len(fields) < 2:
            _fail(f"line {row_number}: expected magnitude and at least one taxonomy label")
        magnitude = parse_nonnegative_integer(fields[0], f"line {row_number} magnitude")
        labels = tuple(_validate_label(label, row_number) for label in fields[1:])
        records.append((magnitude, labels))
    return records


def _new_node() -> dict[str, object]:
    return {"direct": 0, "clade": 0, "children": {}}


def build_tree(records: list[tuple[int, tuple[str, ...]]]) -> dict[str, object]:
    """Aggregate duplicate paths while retaining direct and clade totals."""

    if not records:
        _fail("input contains no data lines")
    root = _new_node()
    for magnitude, labels in records:
        node = root
        node["clade"] = int(node["clade"]) + magnitude
        for label in labels:
            children = node["children"]
            assert isinstance(children, dict)
            child = children.setdefault(label, _new_node())
            assert isinstance(child, dict)
            child["clade"] = int(child["clade"]) + magnitude
            node = child
        node["direct"] = int(node["direct"]) + magnitude
    return root


def validate_tree(node: dict[str, object]) -> int:
    """Validate and return a node's clade total from direct and child totals."""

    direct = int(node["direct"])
    clade = int(node["clade"])
    if direct < 0:
        _fail("Krona tree direct magnitude cannot be negative")
    children = node["children"]
    assert isinstance(children, dict)
    child_total = 0
    for child in children.values():
        if not isinstance(child, dict):
            _fail("Krona tree contains a malformed child node")
        child_total += validate_tree(child)
    if clade != direct + child_total or clade < direct:
        _fail("Krona tree direct/clade invariant failed")
    return clade


def _element_text(parent: ET.Element, tag: str, text: str) -> ET.Element:
    child = ET.SubElement(parent, tag)
    child.text = text
    return child


def _append_node(parent: ET.Element, name: str, node: dict[str, object]) -> None:
    xml_node = ET.SubElement(parent, "node", {"name": name})
    magnitude = ET.SubElement(xml_node, "magnitude")
    _element_text(magnitude, "val", str(int(node["clade"])))
    unassigned = ET.SubElement(xml_node, "magnitudeUnassigned")
    _element_text(unassigned, "val", str(int(node["direct"])))
    children = node["children"]
    assert isinstance(children, dict)
    for child_name in sorted(children):
        child = children[child_name]
        assert isinstance(child, dict)
        _append_node(xml_node, child_name, child)


def build_krona_xml(
    records: list[tuple[int, tuple[str, ...]]], dataset_name: str
) -> bytes:
    """Return deterministic Krona 2.0-compatible XML bytes."""

    dataset = _require_text(dataset_name, "dataset name")
    tree = build_tree(records)
    validate_tree(tree)
    root = ET.Element("krona", {"collapse": "false", "key": "false"})
    attributes = ET.SubElement(root, "attributes", {"magnitude": "magnitude"})
    attribute = ET.SubElement(attributes, "attribute", {"display": "Total"})
    attribute.text = "magnitude"
    unassigned_attribute = ET.SubElement(
        attributes, "attribute", {"display": "Unassigned"}
    )
    unassigned_attribute.text = "magnitudeUnassigned"
    datasets = ET.SubElement(root, "datasets")
    _element_text(datasets, "dataset", dataset)
    _append_node(root, dataset, tree)
    return ET.tostring(root, encoding="utf-8", xml_declaration=True) + b"\n"


def _load_manifest(vendor_dir: os.PathLike[str] | str) -> tuple[dict[str, object], Path]:
    root = Path(vendor_dir)
    if root.is_symlink() or not root.is_dir():
        _fail(f"Krona vendor directory is missing: {root}")
    manifest_path = root / "SOURCE.json"
    if not manifest_path.is_file():
        _fail(f"Krona vendor manifest is missing: {manifest_path}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        _fail(f"could not read Krona vendor manifest '{manifest_path}': {exc}")
    if not isinstance(manifest, dict):
        _fail("Krona vendor manifest must be a JSON object")
    if manifest.get("upstream") != "marbl/Krona" or manifest.get("tag") != VENDOR_TAG:
        _fail("Krona vendor manifest does not identify marbl/Krona v2.8.1")
    if not isinstance(manifest.get("commit"), str) or not re.fullmatch(r"[0-9a-f]{40}", manifest["commit"]):
        _fail("Krona vendor manifest has an invalid upstream commit")
    entries = manifest.get("files")
    if not isinstance(entries, list):
        _fail("Krona vendor manifest files must be a list")
    declared: set[str] = set()
    declared_casefolded: set[str] = set()
    for entry in entries:
        if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
            _fail("Krona vendor manifest contains an invalid file entry")
        relative = entry["path"]
        if (
            not relative
            or "\\" in relative
            or relative.startswith("/")
            or re.match(r"^[A-Za-z]:", relative)
        ):
            _fail(f"Krona vendor manifest path is unsafe: {relative}")
        parts = relative.split("/")
        pure = PurePosixPath(relative)
        if (
            pure.is_absolute()
            or any(part in {"", ".", ".."} for part in parts)
            or pure.as_posix() != relative
        ):
            _fail(f"Krona vendor manifest path is unsafe: {relative}")
        folded = relative.casefold()
        if relative in declared or folded in declared_casefolded:
            _fail(f"Krona vendor manifest contains duplicate or case-colliding path: {relative}")
        declared.add(relative)
        declared_casefolded.add(folded)
        expected_hash = entry.get("sha256")
        if not isinstance(expected_hash, str) or not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
            _fail(f"Krona vendor manifest has an invalid SHA-256 for {relative}")
        file_path = root.joinpath(*pure.parts)
        if file_path.is_symlink() or not file_path.is_file():
            _fail(f"Krona vendor file is missing: {file_path}")
        actual_hash = hashlib.sha256(file_path.read_bytes()).hexdigest()
        if actual_hash != expected_hash:
            _fail(f"Krona vendor SHA-256 mismatch for {relative}")
    for candidate in root.rglob("*"):
        if candidate.is_symlink():
            _fail(f"Krona vendor tree contains a symlink: {candidate}")
    actual = {
        candidate.relative_to(root).as_posix()
        for candidate in root.rglob("*")
        if candidate.is_file() and candidate != manifest_path
    }
    actual_casefolded = {path.casefold() for path in actual}
    if len(actual_casefolded) != len(actual):
        _fail("Krona vendor tree contains case-colliding files")
    if actual != declared:
        _fail("declared and actual Krona vendor inventories differ")
    if declared != set(EXPECTED_VENDOR_FILES):
        _fail("Krona vendor inventory differs from the renderer contract")
    return manifest, root


def validate_vendor(vendor_dir: os.PathLike[str] | str) -> dict[str, object]:
    """Validate pinned renderer assets and return their manifest."""

    manifest, _ = _load_manifest(vendor_dir)
    return manifest


def _read_vendor_asset(root: Path, relative: str) -> bytes:
    path = root.joinpath(*PurePosixPath(relative).parts)
    try:
        return path.read_bytes()
    except OSError as exc:
        _fail(f"could not read Krona vendor asset '{relative}': {exc}")


def _data_uri(relative: str, content: bytes) -> str:
    mime = {
        "img/favicon.ico": "image/x-icon",
        "img/hidden.png": "image/png",
        "img/loading.gif": "image/gif",
        "img/logo-med.png": "image/png",
    }.get(relative) or mimetypes.guess_type(relative)[0]
    if not mime:
        _fail(f"could not determine MIME type for Krona asset '{relative}'")
    return f"data:{mime};base64,{base64.b64encode(content).decode('ascii')}"


def build_html(xml_bytes: bytes, dataset_name: str, vendor_dir: os.PathLike[str] | str) -> bytes:
    """Embed Krona JavaScript, images, and XML into a standalone HTML document."""

    dataset = _require_text(dataset_name, "dataset name")
    manifest, root = _load_manifest(vendor_dir)
    vendor_version = str(manifest["tag"]).removeprefix("v")
    javascript = _read_vendor_asset(root, "src/krona-2.0.js").decode("utf-8")
    hidden_uri = _data_uri("img/hidden.png", _read_vendor_asset(root, "img/hidden.png"))
    loading_uri = _data_uri("img/loading.gif", _read_vendor_asset(root, "img/loading.gif"))
    favicon_uri = _data_uri("img/favicon.ico", _read_vendor_asset(root, "img/favicon.ico"))
    logo_uri = _data_uri("img/logo-med.png", _read_vendor_asset(root, "img/logo-med.png"))
    try:
        xml_text = xml_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        _fail(f"Krona XML is not valid UTF-8: {exc}")
    if "</script" in javascript.lower():
        _fail("vendored Krona JavaScript contains a closing script tag")

    title = html.escape(f"Krona - {dataset}", quote=True)
    document = (
        "<!doctype html>\n"
        '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">\n'
        " <head>\n"
        '  <meta charset="utf-8"/>\n'
        f'  <meta name="generator" content="AAy Amplicon builtin {BUILDER_VERSION} + Krona {vendor_version}"/>\n'
        f"  <title>{title}</title>\n"
        f'  <link rel="shortcut icon" href="{favicon_uri}"/>\n'
        '  <script type="text/javascript">\n'
        f"{javascript}\n"
        "  </script>\n"
        " </head>\n"
        " <body>\n"
        f'  <img id="hiddenImage" src="{hidden_uri}" style="display:none" alt="Hidden Image"/>\n'
        f'  <img id="loadingImage" src="{loading_uri}" style="display:none" alt="Loading Indicator"/>\n'
        f'  <img id="logo" src="{logo_uri}" style="display:none" alt="Logo of Krona"/>\n'
        "  <noscript>Javascript must be enabled to view this page.</noscript>\n"
        '  <div style="display:none">\n'
        f"{xml_text}"
        "  </div>\n"
        " </body>\n"
        "</html>\n"
    )
    return document.encode("utf-8")


def atomic_write(path: os.PathLike[str] | str, content: bytes) -> None:
    """Atomically replace *path* using a sibling temporary file."""

    output = Path(path)
    parent = output.parent
    temporary: Path | None = None
    try:
        parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="wb", prefix=output.name + ".tmp-", dir=parent, delete=False
        ) as handle:
            temporary = Path(handle.name)
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, output)
        temporary = None
    except Exception as exc:
        if temporary is not None:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass
        if isinstance(exc, BuilderError):
            raise
        _fail(f"could not atomically write '{output}': {exc}")


def render(
    input_path: os.PathLike[str] | str,
    output_path: os.PathLike[str] | str,
    dataset_name: str,
    expected_total: object,
    vendor_dir: os.PathLike[str] | str,
) -> None:
    records = parse_krona_input(input_path)
    total = parse_nonnegative_integer(expected_total, "expected total")
    observed = sum(magnitude for magnitude, _labels in records)
    if observed != total:
        _fail(f"input magnitude sum ({observed}) does not equal expected total ({total})")
    xml_bytes = build_krona_xml(records, dataset_name)
    atomic_write(output_path, build_html(xml_bytes, dataset_name, vendor_dir))


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_pos", nargs="?", help="Krona TSV input")
    parser.add_argument("output_pos", nargs="?", help="standalone HTML output")
    parser.add_argument("--input", dest="input_path", help="Krona TSV input")
    parser.add_argument("--output", dest="output_path", help="standalone HTML output")
    parser.add_argument("--dataset-name", "--dataset", dest="dataset_name")
    parser.add_argument("--expected-total", dest="expected_total")
    parser.add_argument("--vendor-dir", required=True, help="pinned Krona vendor directory")
    parser.add_argument("--validate-only", action="store_true", help="validate assets without writing output")
    parser.add_argument("--version", action="version", version=f"{Path(__file__).name} {BUILDER_VERSION}")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        validate_vendor(args.vendor_dir)
        if args.validate_only:
            if args.input_path or args.input_pos:
                input_path = args.input_path or args.input_pos
                parse_krona_input(input_path)
            return 0
        input_path = args.input_path or args.input_pos
        output_path = args.output_path or args.output_pos
        if not input_path or not output_path:
            _fail("--input and --output are required unless --validate-only is used")
        if args.dataset_name is None:
            _fail("--dataset-name is required")
        if args.expected_total is None:
            _fail("--expected-total is required")
        render(input_path, output_path, args.dataset_name, args.expected_total, args.vendor_dir)
        return 0
    except (BuilderError, OSError, UnicodeError, binascii.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
