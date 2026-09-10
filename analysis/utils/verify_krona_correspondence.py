#!/usr/bin/env python3
"""Verify exact taxonomy/count correspondence between Krona TSV and HTML."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from xml.etree import ElementTree as ET


def _canonical_integer(value: str | None, location: str) -> int:
    if value is None or not value.isdigit() or (len(value) > 1 and value[0] == "0"):
        raise ValueError(f"non-canonical integer at {location}")
    return int(value)


def verify_krona_correspondence(
    html_path: str | Path,
    tsv_path: str | Path,
    sample_id: str,
    expected_total: int,
) -> int:
    html_path = Path(html_path)
    tsv_path = Path(tsv_path)
    direct_expected: dict[tuple[str, ...], int] = {}
    for line_number, line in enumerate(
        tsv_path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        fields = line.split("\t")
        if len(fields) < 2 or any(not label for label in fields[1:]):
            raise ValueError(f"invalid Krona TSV row {line_number}")
        magnitude = _canonical_integer(fields[0], f"TSV row {line_number}")
        key = tuple(fields[1:])
        direct_expected[key] = direct_expected.get(key, 0) + magnitude

    clade_expected = {(): sum(direct_expected.values())}
    for key, value in direct_expected.items():
        for depth in range(1, len(key) + 1):
            prefix = key[:depth]
            clade_expected[prefix] = clade_expected.get(prefix, 0) + value

    text = html_path.read_text(encoding="utf-8")
    starts = list(re.finditer(r"<krona(?:\s[^>]*)?>", text))
    if len(starts) != 1 or text.count("</krona>") != 1:
        raise ValueError("expected exactly one Krona XML fragment")
    start = starts[0].start()
    end = text.index("</krona>", start) + len("</krona>")
    document = ET.fromstring(text[start:end])

    datasets = [node.text for node in document.findall("./datasets/dataset")]
    if datasets != [tsv_path.stem]:
        raise ValueError("Krona dataset identity does not match TSV basename")
    dataset_nodes = document.findall("./node")
    if len(dataset_nodes) != 1:
        raise ValueError("expected exactly one dataset node")
    if dataset_nodes[0].get("name") != sample_id:
        raise ValueError("Krona root dataset name does not match sample identity")

    observed: dict[tuple[str, ...], tuple[int, int]] = {}

    def check(node: ET.Element, path: tuple[str, ...] = ()) -> int:
        clade = _canonical_integer(node.findtext("./magnitude/val"), "magnitude")
        direct = _canonical_integer(
            node.findtext("./magnitudeUnassigned/val"), "magnitudeUnassigned"
        )
        children = node.findall("./node")
        names = [child.get("name") for child in children]
        if any(name is None or name == "" for name in names) or len(names) != len(
            set(names)
        ):
            raise ValueError("missing or duplicate Krona sibling label")
        child_clades = sum(
            check(child, path + (child.get("name"),)) for child in children
        )
        if clade < direct or clade != direct + child_clades:
            raise ValueError("Krona clade/direct arithmetic mismatch")
        observed[path] = (direct, clade)
        return clade

    root_total = check(dataset_nodes[0])
    expected_keys = set(clade_expected)
    if set(observed) != expected_keys:
        raise ValueError("Krona HTML and TSV taxonomy paths differ")
    for key in expected_keys:
        expected_direct = direct_expected.get(key, 0) if key else 0
        if observed[key] != (expected_direct, clade_expected[key]):
            raise ValueError("Krona HTML and TSV direct/clade counts differ")
    if root_total != expected_total:
        raise ValueError("Krona root total differs from expected total")
    return root_total


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", required=True)
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--expected-total", required=True, type=int)
    args = parser.parse_args()
    print(
        verify_krona_correspondence(
            args.html, args.tsv, args.sample_id, args.expected_total
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
