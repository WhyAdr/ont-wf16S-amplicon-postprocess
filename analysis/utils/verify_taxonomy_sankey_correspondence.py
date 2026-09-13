"""Independently verify a taxonomy Sankey artifact against its source payload.

This module deliberately does not import the producer.  It contains a small
reference model with separate data flow and validation code so release checks
can detect a producer and verifier making the same mistake.
"""

from __future__ import annotations

import argparse
import base64
import bisect
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any


RANKS = ("D", "K", "P", "C", "O", "F", "G", "S")
RANK_INDEX = {rank: index for index, rank in enumerate(RANKS)}
SAFE_MAX = 9007199254740991
ENTRY_ID = "synthetic:entry:classified"
COUNT_MODEL = "classified clade-read flow with persistent explicit residual lanes"
NODE_ORDER = {"entry": 0, "taxon": 1, "residual": 2, "residual_carry": 3}
LINK_ORDER = {"biological": 0, "other_hidden": 1, "assigned_above": 2, "carry": 3}
SUBTYPE_ORDER = {"other_hidden": 0, "assigned_above": 1}
CONTROL = re.compile(r"[\x00-\x1f\x7f]")
INTEGER_TEXT = re.compile(r"(?:0|[1-9][0-9]*)\Z")


class CorrespondenceError(ValueError):
    """A fail-closed correspondence or artifact error."""


def fail(message: str) -> None:
    raise CorrespondenceError(message)


def read_json(path: Path, label: str) -> tuple[dict[str, Any], bytes]:
    if path.is_symlink() or not path.is_file():
        fail(f"{label} is not a regular file: {path}")
    raw = path.read_bytes()
    if raw.startswith(b"\xef\xbb\xbf"):
        fail(f"{label} contains a UTF-8 BOM")
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        fail(f"{label} is not valid UTF-8 JSON: {exc}")
    if not isinstance(value, dict):
        fail(f"{label} root must be an object")
    return value, raw


def text(value: Any, label: str, *, empty: bool = False) -> str:
    if not isinstance(value, str) or (not empty and not value):
        fail(f"{label} must be a non-empty string")
    if CONTROL.search(value):
        fail(f"{label} contains a control character")
    return value


def integer(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        fail(f"{label} must be a non-negative integer")
    if value > SAFE_MAX:
        fail(f"{label} exceeds the JavaScript-safe integer range")
    return value


def sha(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def utf8(value: str) -> bytes:
    return value.encode("utf-8")


def canonical(value: Any) -> bytes:
    """Apply the contract's byte rules independently of the producer."""

    def check(item: Any, path: str) -> None:
        if item is None or isinstance(item, (bool, str)):
            return
        if isinstance(item, int):
            integer(item, path)
            return
        if isinstance(item, float):
            fail(f"{path} contains a floating-point number")
        if isinstance(item, list):
            for index, child in enumerate(item):
                check(child, f"{path}[{index}]")
            return
        if isinstance(item, dict):
            for key, child in item.items():
                text(key, f"{path} key")
                check(child, f"{path}.{key}")
            return
        fail(f"{path} contains an unsupported JSON value")

    check(value, "$")
    return (json.dumps(value, ensure_ascii=False, sort_keys=False,
                       separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")


def validate_source(payload: dict[str, Any]) -> dict[str, Any]:
    needed = {"schema_version", "sample_id", "renderer", "renderer_version", "totals", "nodes"}
    absent = sorted(needed - set(payload))
    if absent:
        fail(f"source payload missing: {', '.join(absent)}")
    if payload["schema_version"] != 1 or payload["renderer"] != "builtin_kraken_report_explorer":
        fail("source payload identity is not the v1 builtin explorer contract")
    sample = text(payload["sample_id"], "source sample_id")
    version = text(payload["renderer_version"], "source renderer_version")
    totals = payload["totals"]
    if not isinstance(totals, dict):
        fail("source totals must be an object")
    total = integer(totals.get("total"), "source total")
    classified = integer(totals.get("classified"), "source classified")
    unclassified = integer(totals.get("unclassified"), "source unclassified")
    if total != classified + unclassified:
        fail("source totals do not conserve reads")
    raw_nodes = payload["nodes"]
    if not isinstance(raw_nodes, list):
        fail("source nodes must be an array")

    nodes: list[dict[str, Any]] = []
    by_path: dict[str, dict[str, Any]] = {}
    by_rank: dict[str, list[dict[str, Any]]] = {rank: [] for rank in RANKS}
    for index, raw in enumerate(raw_nodes):
        if not isinstance(raw, dict):
            fail(f"source nodes[{index}] must be an object")
        path = text(raw.get("path"), f"source nodes[{index}].path")
        parent_path = text(raw.get("parent_path"), f"source nodes[{index}].parent_path", empty=True)
        name = text(raw.get("name"), f"source nodes[{index}].name")
        rank = text(raw.get("rank_code"), f"source nodes[{index}].rank_code")
        if rank not in RANK_INDEX:
            fail(f"source nodes[{index}] has an invalid rank")
        depth = integer(raw.get("depth"), f"source nodes[{index}].depth")
        if depth != RANK_INDEX[rank] + 1:
            fail(f"source nodes[{index}] has an invalid depth")
        taxid = raw.get("taxid")
        if not isinstance(taxid, str) or not INTEGER_TEXT.fullmatch(taxid):
            fail(f"source nodes[{index}].taxid is not canonical")
        integer(int(taxid), f"source nodes[{index}].taxid")
        direct = integer(raw.get("direct"), f"source nodes[{index}].direct")
        clade = integer(raw.get("clade"), f"source nodes[{index}].clade")
        status = text(raw.get("status"), f"source nodes[{index}].status")
        if status not in {"resolved", "unresolved", "conflicted"}:
            fail(f"source nodes[{index}] has an invalid status")
        if direct > clade:
            fail(f"source nodes[{index}] direct exceeds clade")
        expected = f"{parent_path};{name}" if parent_path else name
        if path != expected:
            fail(f"source path/name mismatch at '{path}'")
        if path in by_path:
            fail(f"duplicate source path '{path}'")
        record = {
            "id": f"taxon:{rank}:{sha(path)}", "taxid": taxid, "rank": rank,
            "name": name, "path": path, "parent_path": parent_path,
            "parent_id": None, "clade": clade, "direct": direct, "status": status,
            "source_order_index": index,
        }
        nodes.append(record)
        by_path[path] = record
        by_rank[rank].append(record)

    for record in nodes:
        parent_path = record["parent_path"]
        if parent_path:
            parent = by_path.get(parent_path)
            if parent is None:
                fail(f"missing source ancestor for '{record['path']}'")
            if RANK_INDEX[parent["rank"]] + 1 != RANK_INDEX[record["rank"]]:
                fail(f"non-adjacent source ancestor for '{record['path']}'")
            record["parent_id"] = parent["id"]
        elif RANK_INDEX[record["rank"]] != 0:
            fail(f"non-Domain source node '{record['path']}' has no parent")

    children: dict[str, list[dict[str, Any]]] = {}
    for record in nodes:
        children.setdefault(record["parent_path"], []).append(record)
    for record in nodes:
        child_sum = sum(child["clade"] for child in children.get(record["path"], []))
        if record["clade"] != record["direct"] + child_sum:
            fail(f"source clade arithmetic failed at '{record['path']}'")
    if nodes and sum(node["clade"] for node in by_rank["D"]) != classified:
        fail("source Domain clades do not equal classified reads")
    if sum(node["direct"] for node in nodes) != classified:
        fail("source direct counts do not equal classified reads")
    return {"sample_id": sample, "renderer_version": version, "total": total,
            "classified": classified, "unclassified": unclassified, "nodes": nodes,
            "by_path": by_path, "by_rank": by_rank}


def ranks_arg(value: str) -> list[str]:
    ranks = value.split(",")
    if len(ranks) < 2 or len(ranks) > 8 or len(set(ranks)) != len(ranks):
        fail("invalid selected ranks")
    if any(rank not in RANK_INDEX for rank in ranks):
        fail("selected ranks contain an invalid code")
    if ranks != sorted(ranks, key=RANK_INDEX.get):
        fail("selected ranks are not a canonical subsequence")
    return ranks


def select(source: dict[str, Any], ranks: list[str], limit: int) -> tuple[dict[str, set[str]], dict[str, str]]:
    if not isinstance(limit, int) or isinstance(limit, bool) or not 1 <= limit <= 100:
        fail("invalid top-N limit")
    retained = {rank: set() for rank in ranks}
    reason: dict[str, str] = {}
    for rank in ranks:
        candidates = [node for node in source["by_rank"][rank] if node["clade"] > 0]
        candidates.sort(key=lambda n: (-n["clade"], utf8(n["path"]), n["source_order_index"], n["id"]))
        for node in candidates[:limit]:
            retained[rank].add(node["path"])
            reason[node["path"]] = "top_n"
    for later_index, later_rank in enumerate(ranks[1:], 1):
        for path in list(retained[later_rank]):
            selected = source["by_path"][path]
            for earlier_rank in ranks[:later_index]:
                current = selected
                while current["rank"] != earlier_rank:
                    if not current["parent_path"]:
                        fail(f"cannot close ancestor for '{path}'")
                    current = source["by_path"][current["parent_path"]]
                retained[earlier_rank].add(current["path"])
                reason.setdefault(current["path"], "ancestor_closure")
    return retained, reason


def digest(paths: list[str]) -> str:
    return sha("\n".join(sorted(paths, key=utf8)))


def node(**fields: Any) -> dict[str, Any]:
    order = ("id", "kind", "subtype", "biological", "carried", "lane_id", "taxid",
             "rank", "name", "path", "parent_id", "clade", "direct", "status",
             "selection_reason", "source_order_index", "visual_order", "origin_source_id",
             "origin_target_rank", "column_rank", "value", "member_count", "member_paths_sha256")
    return {key: fields.get(key) for key in order}


def lane(subtype: str, origin: str, origin_rank: str) -> str:
    return sha("\0".join((subtype, origin, origin_rank)))


def residual(subtype: str, origin: str, origin_index: int | None, origin_rank: str,
             column: str, amount: int, members: list[str] | None, carried: bool) -> dict[str, Any]:
    lane_id = lane(subtype, origin, origin_rank)
    return node(
        id=f"synthetic:residual:{lane_id}:{column}",
        kind="residual_carry" if carried else "residual", subtype=subtype,
        biological=False, carried=carried, lane_id=lane_id, taxid=None, rank=None,
        name=f"Other {column}" if subtype == "other_hidden" else f"Assigned above {column}",
        path=None, parent_id=None, clade=None, direct=None, status=None,
        selection_reason="carry" if carried else "residual", source_order_index=origin_index,
        visual_order=None, origin_source_id=origin, origin_target_rank=origin_rank,
        column_rank=column, value=amount,
        member_count=len(members) if members is not None else None,
        member_paths_sha256=digest(members) if members is not None else None,
    )


def biological(source_node: dict[str, Any], column: str, reason: str) -> dict[str, Any]:
    return node(
        id=source_node["id"], kind="taxon", subtype=None, biological=True, carried=False,
        lane_id=None, taxid=source_node["taxid"], rank=source_node["rank"],
        name=source_node["name"], path=source_node["path"], parent_id=source_node["parent_id"],
        clade=source_node["clade"], direct=source_node["direct"], status=source_node["status"],
        selection_reason=reason, source_order_index=source_node["source_order_index"],
        visual_order=None, origin_source_id=None, origin_target_rank=None, column_rank=column,
        value=source_node["clade"], member_count=None, member_paths_sha256=None,
    )


def edge(source: str, target: str, kind: str, amount: int, transition: int) -> dict[str, Any]:
    return {"id": "link:" + sha("\0".join((source, target, kind))), "source": source,
            "target": target, "kind": kind, "value": amount, "transition_index": transition}


def descendant(path: str, ancestor: str) -> bool:
    return path == ancestor or path.startswith(ancestor + ";")


def visual_order(items: list[dict[str, Any]], source: dict[str, Any], ranks: list[str]) -> None:
    for rank in ranks:
        bio = sorted((item for item in items if item["kind"] == "taxon" and item["column_rank"] == rank),
                     key=lambda item: item["source_order_index"])
        position = {item["id"]: index for index, item in enumerate(bio)}
        paths = sorted((utf8(item["path"]) for item in bio))
        ranked: list[tuple[tuple[Any, ...], dict[str, Any]]] = [
            ((position[item["id"]], 1, 0, "", item["id"]), item) for item in bio
        ]
        for item in items:
            if item["column_rank"] != rank or item["kind"] not in {"residual", "residual_carry"}:
                continue
            origin_id = item["origin_source_id"]
            if origin_id == ENTRY_ID:
                slot = len(bio)
            else:
                origin = next((candidate for candidate in source["nodes"] if candidate["id"] == origin_id), None)
                if origin is None:
                    fail(f"unknown residual origin '{origin_id}'")
                descendants = [position[b["id"]] for b in bio if descendant(b["path"], origin["path"])]
                slot = max(descendants) + 1 if descendants else bisect.bisect_right(paths, utf8(origin["path"]))
            ranked.append(((slot, 0, SUBTYPE_ORDER[item["subtype"]], item["lane_id"], item["id"]), item))
        ranked.sort(key=lambda pair: pair[0])
        for index, (_, item) in enumerate(ranked):
            item["visual_order"] = index


def build_view(source: dict[str, Any], ranks: list[str], limit: int) -> dict[str, Any]:
    classified = source["classified"]
    if classified == 0:
        return {"nodes": [], "links": [], "conservation": {
            "classified": 0,
            "column_totals": [{"rank": rank, "value": 0} for rank in ranks],
            "transition_totals": [{"from": "ENTRY" if i == 0 else ranks[i - 1],
                                   "to": rank, "value": 0} for i, rank in enumerate(ranks)],
            "rightmost_flow": 0,
        }}
    retained, reasons = select(source, ranks, limit)
    nodes: list[dict[str, Any]] = [node(
        id=ENTRY_ID, kind="entry", subtype=None, biological=False, carried=False, lane_id=None,
        taxid=None, rank=None, name="Classified reads", path=None, parent_id=None, clade=None,
        direct=None, status=None, selection_reason="entry", source_order_index=None,
        visual_order=0, origin_source_id=None, origin_target_rank=None, column_rank="ENTRY",
        value=classified, member_count=None, member_paths_sha256=None,
    )]
    links: list[dict[str, Any]] = []
    columns: dict[str, list[dict[str, Any]]] = {}
    residuals: list[dict[str, Any]] = []
    first = ranks[0]
    first_bio = [biological(source["by_path"][path], first, reasons[path]) for path in retained[first]]
    first_bio.sort(key=lambda item: item["source_order_index"])
    columns[first] = first_bio
    nodes.extend(first_bio)
    all_first = source["by_rank"][first]
    hidden = [item for item in all_first if item["path"] not in retained[first]]
    hidden_amount = sum(item["clade"] for item in hidden)
    above_amount = classified - sum(item["clade"] for item in all_first)
    if above_amount < 0:
        fail("first selected rank exceeds classified flow")
    for subtype, amount, members in (("other_hidden", hidden_amount, [item["path"] for item in hidden]),
                                     ("assigned_above", above_amount, None)):
        if amount:
            item = residual(subtype, ENTRY_ID, None, first, first, amount, members, False)
            residuals.append(item)
            nodes.append(item)
            links.append(edge(ENTRY_ID, item["id"], subtype, amount, 0))
    for item in first_bio:
        links.append(edge(ENTRY_ID, item["id"], "biological", item["value"], 0))

    active_residuals = residuals
    for transition, (source_rank, target_rank) in enumerate(zip(ranks, ranks[1:]), 1):
        target_bio = [biological(source["by_path"][path], target_rank, reasons[path])
                      for path in retained[target_rank]]
        target_bio.sort(key=lambda item: item["source_order_index"])
        columns[target_rank] = target_bio
        nodes.extend(target_bio)
        target_by_path = {item["path"]: item for item in target_bio}
        next_residuals: list[dict[str, Any]] = []
        for source_item in columns[source_rank]:
            source_record = source["by_path"][source_item["path"]]
            descendants = [item for item in source["by_rank"][target_rank]
                           if descendant(item["path"], source_item["path"])]
            visible = [item for item in descendants if item["path"] in target_by_path]
            hidden_targets = [item for item in descendants if item["path"] not in target_by_path]
            assigned = source_record["clade"] - sum(item["clade"] for item in descendants)
            if assigned < 0:
                fail(f"negative assigned-above flow at '{source_item['path']}'")
            for target in visible:
                links.append(edge(source_item["id"], target_by_path[target["path"]]["id"],
                                  "biological", target["clade"], transition))
            for subtype, amount, members in (("other_hidden", sum(item["clade"] for item in hidden_targets),
                                              [item["path"] for item in hidden_targets]),
                                             ("assigned_above", assigned, None)):
                if amount:
                    item = residual(subtype, source_item["id"], source_item["source_order_index"],
                                    target_rank, target_rank, amount, members, False)
                    next_residuals.append(item)
                    nodes.append(item)
                    links.append(edge(source_item["id"], item["id"], subtype, amount, transition))
        for prior in active_residuals:
            carry = residual(prior["subtype"], prior["origin_source_id"], prior["source_order_index"],
                             prior["origin_target_rank"], target_rank, prior["value"],
                             None if prior["subtype"] == "assigned_above" else [], True)
            if prior["subtype"] == "other_hidden":
                carry["member_count"] = prior["member_count"]
                carry["member_paths_sha256"] = prior["member_paths_sha256"]
            next_residuals.append(carry)
            nodes.append(carry)
            links.append(edge(prior["id"], carry["id"], "carry", prior["value"], transition))
        active_residuals = next_residuals

    visual_order(nodes, source, ranks)
    nodes.sort(key=lambda item: (-1 if item["kind"] == "entry" else RANK_INDEX[item["column_rank"]],
                                 NODE_ORDER[item["kind"]],
                                 -1 if item["source_order_index"] is None else item["source_order_index"],
                                 item["id"]))
    node_order = {item["id"]: index for index, item in enumerate(nodes)}
    links.sort(key=lambda item: (item["transition_index"], node_order[item["source"]],
                                 node_order[item["target"]], LINK_ORDER[item["kind"]], item["id"]))
    conservation = {"classified": classified,
                    "column_totals": [{"rank": rank, "value": classified} for rank in ranks],
                    "transition_totals": [{"from": "ENTRY" if i == 0 else ranks[i - 1],
                                           "to": rank, "value": classified}
                                          for i, rank in enumerate(ranks)],
                    "rightmost_flow": classified}
    # Independently check the compact reference result before comparing it.
    if any(item["value"] <= 0 for item in links):
        fail("reference model emitted a zero-valued link")
    for rank in ranks:
        if sum(item["value"] for item in nodes if item["column_rank"] == rank) != classified:
            fail(f"reference column {rank} does not conserve")
    if any(sum(item["value"] for item in links if item["transition_index"] == i) != classified
           for i in range(len(ranks))):
        fail("reference transition does not conserve")
    return {"nodes": nodes, "links": links, "conservation": conservation}


def expected_document(payload: dict[str, Any], source_bytes: bytes, ranks: list[str],
                      limit: int, renderer_version: str) -> dict[str, Any]:
    source = validate_source(payload)
    source_nodes = [{"id": item["id"], "taxid": item["taxid"], "rank": item["rank"],
                     "name": item["name"], "path": item["path"], "parent_id": item["parent_id"],
                     "clade": item["clade"], "direct": item["direct"], "status": item["status"],
                     "source_order_index": item["source_order_index"]} for item in source["nodes"]]
    return {
        "schema_version": 1,
        "sample_id": source["sample_id"],
        "renderer": "builtin_taxonomy_sankey",
        "renderer_version": renderer_version,
        "source_renderer": "builtin_kraken_report_explorer",
        "source_schema_version": 1,
        "source_payload_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "count_model": COUNT_MODEL,
        "denominator": "TotalReads",
        "totals": {"total": source["total"], "classified": source["classified"],
                   "unclassified": source["unclassified"]},
        "defaults": {"ranks": ranks, "max_taxa_per_rank": limit},
        "rank_order": ranks,
        "source_nodes": source_nodes,
        "entry": {"id": ENTRY_ID, "value": source["classified"]},
        "default_view": build_view(source, ranks, limit),
    }


def verify_html(path: Path, json_bytes: bytes) -> None:
    if path.is_symlink() or not path.is_file():
        fail(f"Sankey HTML is not a regular file: {path}")
    html = path.read_bytes().decode("utf-8")
    if "\ufeff" in html or re.search(r"(?:src|href)=[\"']https?://", html, re.I):
        fail("Sankey HTML contains a BOM or external resource")
    if "connect-src 'none'" not in html or "default-src 'none'" not in html:
        fail("Sankey HTML does not have the required restrictive CSP")
    if re.search(r"<script\b[^>]*\bsrc=", html, re.I):
        fail("Sankey HTML has an external script")
    forbidden = ("fetch(", "XMLHttpRequest", "WebSocket", "navigator.sendBeacon",
                 "EventSource", "importScripts(", "navigator.serviceWorker",
                 "document.cookie", "localStorage", "sessionStorage")
    if any(token in html for token in forbidden):
        fail("Sankey HTML contains a forbidden capability")
    match = re.search(r'<script id="wf16s-payload"[^>]*>([^<]*)</script>', html)
    if match is None:
        fail("Sankey HTML does not embed the canonical payload")
    try:
        embedded = base64.b64decode(match.group(1), validate=True)
    except Exception as exc:
        fail(f"embedded Sankey payload is not base64: {exc}")
    if embedded != json_bytes:
        fail("embedded Sankey payload differs from Sankey JSON bytes")


def verify(source_path: Path, sankey_path: Path, html_path: Path | None,
           ranks: list[str], limit: int, renderer_version: str) -> None:
    payload, source_bytes = read_json(source_path, "source Pavian JSON")
    observed, observed_bytes = read_json(sankey_path, "Sankey JSON")
    expected = expected_document(payload, source_bytes, ranks, limit, renderer_version)
    if observed != expected:
        fail("Sankey JSON does not match the independent reference model")
    if observed_bytes != canonical(observed):
        fail("Sankey JSON is not canonical UTF-8 JSON")
    if html_path is not None:
        verify_html(html_path, observed_bytes)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", required=True, type=Path)
    parser.add_argument("--sankey-json", required=True, type=Path)
    parser.add_argument("--sankey-html", type=Path)
    parser.add_argument("--ranks", required=True)
    parser.add_argument("--max-taxa-per-rank", required=True, type=int)
    parser.add_argument("--renderer-version", required=True)
    args = parser.parse_args(argv)
    try:
        verify(args.payload, args.sankey_json, args.sankey_html, ranks_arg(args.ranks),
               args.max_taxa_per_rank, text(args.renderer_version, "renderer version"))
        return 0
    except Exception as exc:
        print(f"verify_taxonomy_sankey_correspondence: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
