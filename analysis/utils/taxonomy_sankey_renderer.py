"""Build a deterministic, offline taxonomy-flow Sankey artifact.

The input is the already verified builtin Kraken-report explorer payload.  The
renderer is intentionally independent of the upstream Pavian application and
does not perform taxonomy lookups or network access.
"""

from __future__ import annotations

import argparse
import base64
import bisect
import hashlib
import html
import json
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Any


RANKS = ("D", "K", "P", "C", "O", "F", "G", "S")
RANK_INDEX = {rank: index for index, rank in enumerate(RANKS)}
MAX_SAFE_INTEGER = 9007199254740991
ENTRY_ID = "synthetic:entry:classified"
COUNT_MODEL = "classified clade-read flow with persistent explicit residual lanes"
NODE_KINDS = {"entry": 0, "taxon": 1, "residual": 2, "residual_carry": 3}
LINK_KINDS = {"biological": 0, "other_hidden": 1, "assigned_above": 2, "carry": 3}
SUBTYPE_ORDER = {"other_hidden": 0, "assigned_above": 1}
CONTROL_RE = re.compile(r"[\x00-\x1f\x7f]")
SAFE_INTEGER_TEXT = re.compile(r"(?:0|[1-9][0-9]*)\Z")
RANK_TEXT = re.compile(r"^[DKPCOFGS]\Z")


class SankeyError(ValueError):
    """A concise, fail-closed renderer or verifier error."""


def fail(message: str) -> None:
    raise SankeyError(message)


def _read_bytes(path: Path, label: str) -> bytes:
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        fail(f"{label} is not an existing regular file: {path}")
    try:
        return path.read_bytes()
    except OSError as exc:
        fail(f"could not read {label} '{path}': {exc}")


def _read_json(path: Path, label: str) -> tuple[dict[str, Any], bytes]:
    raw = _read_bytes(path, label)
    if raw.startswith(b"\xef\xbb\xbf"):
        fail(f"{label} must not contain a UTF-8 BOM")
    try:
        text = raw.decode("utf-8", errors="strict")
        value = json.loads(text)
    except (UnicodeError, json.JSONDecodeError) as exc:
        fail(f"could not parse {label} '{path}': {exc}")
    if not isinstance(value, dict):
        fail(f"{label} root must be a JSON object")
    return value, raw


def _safe_int(value: Any, label: str, *, allow_negative: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        fail(f"{label} must be an integer")
    if not allow_negative and value < 0:
        fail(f"{label} must be non-negative")
    if value < -MAX_SAFE_INTEGER or value > MAX_SAFE_INTEGER:
        fail(f"{label} exceeds the JavaScript-safe integer range")
    return value


def _safe_integer_text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not SAFE_INTEGER_TEXT.fullmatch(value):
        fail(f"{label} must be a canonical non-negative integer string")
    parsed = int(value)
    _safe_int(parsed, label)
    return value


def _text(value: Any, label: str, *, allow_empty: bool = False) -> str:
    if not isinstance(value, str) or (not allow_empty and not value):
        fail(f"{label} must be a non-empty string")
    if CONTROL_RE.search(value):
        fail(f"{label} contains a control character")
    return value


def _utf8_key(value: str) -> bytes:
    return value.encode("utf-8", errors="strict")


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def canonical_json_bytes(value: Any) -> bytes:
    """Serialize the v0.4.8 contract without platform-dependent formatting."""

    def validate(item: Any, path: str) -> None:
        if isinstance(item, bool) or item is None or isinstance(item, str):
            return
        if isinstance(item, int):
            _safe_int(item, path, allow_negative=False)
            return
        if isinstance(item, float):
            fail(f"{path} must not contain floating-point JSON numbers")
        if isinstance(item, list):
            for index, child in enumerate(item):
                validate(child, f"{path}[{index}]")
            return
        if isinstance(item, dict):
            for key, child in item.items():
                _text(key, f"{path} key")
                validate(child, f"{path}.{key}")
            return
        fail(f"{path} has an unsupported JSON value")

    validate(value, "$")
    try:
        text = json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=False,
            separators=(",", ":"),
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        fail(f"could not serialize canonical JSON: {exc}")
    return text.encode("utf-8") + b"\n"


def _atomic_write_new(path: Path, data: bytes) -> None:
    path = Path(path)
    if path.exists() or path.is_symlink():
        fail(f"refusing to overwrite existing renderer output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.is_symlink():
        fail(f"renderer output appeared unexpectedly: {path}")
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.tmp-", dir=str(path.parent))
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        if path.exists() or path.is_symlink():
            fail(f"renderer output appeared during publication: {path}")
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def _version() -> str:
    path = Path(__file__).resolve().parents[2] / "VERSION"
    value = _read_bytes(path, "pipeline VERSION").decode("utf-8")
    if not value.endswith("\n"):
        fail("pipeline VERSION must be LF newline-terminated")
    value = value[:-1]
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", value):
        fail("pipeline VERSION must contain one SemVer value")
    return value


def _node_id(rank: str, path: str) -> str:
    return f"taxon:{rank}:{_sha256_text(path)}"


def _lane_id(subtype: str, origin_source_id: str, origin_target_rank: str) -> str:
    return _sha256_text("\0".join((subtype, origin_source_id, origin_target_rank)))


def _residual_id(lane_id: str, column_rank: str) -> str:
    return f"synthetic:residual:{lane_id}:{column_rank}"


def _link_id(source_id: str, target_id: str, kind: str) -> str:
    return "link:" + _sha256_text("\0".join((source_id, target_id, kind)))


def _validate_payload(payload: dict[str, Any]) -> dict[str, Any]:
    required = {
        "schema_version", "sample_id", "renderer", "renderer_version",
        "totals", "nodes",
    }
    missing = sorted(required - set(payload))
    if missing:
        fail(f"source payload missing key(s): {', '.join(missing)}")
    if payload["schema_version"] != 1:
        fail("source payload schema_version must be 1")
    _text(payload["sample_id"], "source sample_id")
    if payload["renderer"] != "builtin_kraken_report_explorer":
        fail("source payload has an unexpected renderer identity")
    _text(payload["renderer_version"], "source renderer_version")
    totals = payload["totals"]
    if not isinstance(totals, dict):
        fail("source totals must be an object")
    total = _safe_int(totals.get("total"), "source totals.total")
    classified = _safe_int(totals.get("classified"), "source totals.classified")
    unclassified = _safe_int(totals.get("unclassified"), "source totals.unclassified")
    if total != classified + unclassified:
        fail("source totals do not conserve total reads")
    nodes = payload["nodes"]
    if not isinstance(nodes, list):
        fail("source nodes must be an array")

    normalized: list[dict[str, Any]] = []
    by_path: dict[str, dict[str, Any]] = {}
    by_rank: dict[str, list[dict[str, Any]]] = {rank: [] for rank in RANKS}
    for index, raw in enumerate(nodes):
        if not isinstance(raw, dict):
            fail(f"source nodes[{index}] must be an object")
        for key in ("path", "parent_path", "name", "rank_code", "taxid", "status"):
            if key not in raw:
                fail(f"source nodes[{index}] is missing '{key}'")
        path = _text(raw["path"], f"source nodes[{index}].path")
        parent_path = _text(raw["parent_path"], f"source nodes[{index}].parent_path", allow_empty=True)
        name = _text(raw["name"], f"source nodes[{index}].name")
        rank = _text(raw["rank_code"], f"source nodes[{index}].rank_code")
        if not RANK_TEXT.fullmatch(rank):
            fail(f"source nodes[{index}] has an invalid rank code")
        depth = _safe_int(raw.get("depth"), f"source nodes[{index}].depth")
        if depth != RANK_INDEX[rank] + 1:
            fail(f"source nodes[{index}] depth does not match rank")
        taxid = _safe_integer_text(raw["taxid"], f"source nodes[{index}].taxid")
        direct = _safe_int(raw.get("direct"), f"source nodes[{index}].direct")
        clade = _safe_int(raw.get("clade"), f"source nodes[{index}].clade")
        if direct > clade:
            fail(f"source nodes[{index}] direct count exceeds clade count")
        status = _text(raw["status"], f"source nodes[{index}].status")
        if status not in {"resolved", "unresolved", "conflicted"}:
            fail(f"source nodes[{index}] has an invalid status")
        expected_path = f"{parent_path};{name}" if parent_path else name
        if path != expected_path:
            fail(f"source nodes[{index}] path does not match parent_path/name")
        if path in by_path:
            fail(f"source payload contains duplicate path '{path}'")
        record: dict[str, Any] = {
            "id": _node_id(rank, path),
            "taxid": taxid,
            "rank": rank,
            "name": name,
            "path": path,
            "parent_path": parent_path,
            "parent_id": None,
            "clade": clade,
            "direct": direct,
            "status": status,
            "source_order_index": index,
        }
        by_path[path] = record
        by_rank[rank].append(record)
        normalized.append(record)

    for record in normalized:
        parent_path = record["parent_path"]
        if parent_path:
            parent = by_path.get(parent_path)
            if parent is None:
                fail(f"source node '{record['path']}' has a missing ancestor")
            if RANK_INDEX[parent["rank"]] + 1 != RANK_INDEX[record["rank"]]:
                fail(f"source node '{record['path']}' has a non-adjacent ancestor")
            record["parent_id"] = parent["id"]
        elif RANK_INDEX[record["rank"]] != 0:
            fail(f"non-Domain source node '{record['path']}' has no parent")

    children: dict[str, list[dict[str, Any]]] = {}
    for record in normalized:
        children.setdefault(record["parent_path"], []).append(record)
    for record in normalized:
        child_sum = sum(child["clade"] for child in children.get(record["path"], []))
        if record["clade"] != record["direct"] + child_sum:
            fail(f"source clade arithmetic failed at '{record['path']}'")
    top_sum = sum(record["clade"] for record in by_rank["D"])
    if top_sum != classified and normalized:
        fail("source Domain clades do not equal classified reads")
    if sum(record["direct"] for record in normalized) != classified:
        fail("source direct counts do not equal classified reads")
    return {
        "sample_id": payload["sample_id"],
        "renderer": payload["renderer"],
        "schema_version": payload["schema_version"],
        "renderer_version": payload["renderer_version"],
        "total": total,
        "classified": classified,
        "unclassified": unclassified,
        "nodes": normalized,
        "by_path": by_path,
        "by_rank": by_rank,
    }


def parse_ranks(value: str | list[str]) -> list[str]:
    ranks = value.split(",") if isinstance(value, str) else list(value)
    if len(ranks) < 2 or len(ranks) > len(RANKS):
        fail("ranks must contain at least two and at most eight rank codes")
    if any(not isinstance(rank, str) or not RANK_TEXT.fullmatch(rank) for rank in ranks):
        fail("ranks must use only D,K,P,C,O,F,G,S")
    if len(set(ranks)) != len(ranks):
        fail("ranks must be unique")
    if ranks != sorted(ranks, key=RANK_INDEX.get):
        fail("ranks must be a canonical subsequence of D,K,P,C,O,F,G,S")
    return ranks


def _selection(source: dict[str, Any], ranks: list[str], max_n: int) -> tuple[dict[str, set[str]], dict[str, str]]:
    if isinstance(max_n, bool) or not isinstance(max_n, int) or not 1 <= max_n <= 100:
        fail("max taxa per rank must be an integer from 1 through 100")
    retained: dict[str, set[str]] = {rank: set() for rank in ranks}
    reasons: dict[str, str] = {}
    for rank in ranks:
        candidates = [node for node in source["by_rank"][rank] if node["clade"] > 0]
        candidates.sort(key=lambda node: (-node["clade"], _utf8_key(node["path"]), node["source_order_index"], node["id"]))
        for node in candidates[:max_n]:
            retained[rank].add(node["path"])
            reasons[node["path"]] = "top_n"

    for later_index in range(1, len(ranks)):
        later_rank = ranks[later_index]
        for path in list(retained[later_rank]):
            current = source["by_path"][path]
            for earlier_rank in ranks[:later_index]:
                while current["rank"] != earlier_rank:
                    parent_path = current["parent_path"]
                    if not parent_path:
                        fail(f"retained node '{path}' has no ancestor at rank {earlier_rank}")
                    current = source["by_path"][parent_path]
                retained[earlier_rank].add(current["path"])
                reasons.setdefault(current["path"], "ancestor_closure")
                current = source["by_path"][path]
    return retained, reasons


def _node_record(*, node_id: str, kind: str, subtype: str | None, biological: bool,
                 carried: bool, lane_id: str | None, taxid: str | None, rank: str | None,
                 name: str | None, path: str | None, parent_id: str | None,
                 clade: int | None, direct: int | None, status: str | None,
                 selection_reason: str | None, source_order_index: int | None,
                 visual_order: int | None, origin_source_id: str | None,
                 origin_target_rank: str | None, column_rank: str,
                 value: int, member_count: int | None,
                 member_paths_sha256: str | None) -> dict[str, Any]:
    return {
        "id": node_id,
        "kind": kind,
        "subtype": subtype,
        "biological": biological,
        "carried": carried,
        "lane_id": lane_id,
        "taxid": taxid,
        "rank": rank,
        "name": name,
        "path": path,
        "parent_id": parent_id,
        "clade": clade,
        "direct": direct,
        "status": status,
        "selection_reason": selection_reason,
        "source_order_index": source_order_index,
        "visual_order": visual_order,
        "origin_source_id": origin_source_id,
        "origin_target_rank": origin_target_rank,
        "column_rank": column_rank,
        "value": value,
        "member_count": member_count,
        "member_paths_sha256": member_paths_sha256,
    }


def _entry_node(value: int) -> dict[str, Any]:
    return _node_record(
        node_id=ENTRY_ID, kind="entry", subtype=None, biological=False, carried=False,
        lane_id=None, taxid=None, rank=None, name="Classified reads", path=None,
        parent_id=None, clade=None, direct=None, status=None, selection_reason="entry",
        source_order_index=None, visual_order=0, origin_source_id=None,
        origin_target_rank=None, column_rank="ENTRY", value=value,
        member_count=None, member_paths_sha256=None,
    )


def _biological_node(source_node: dict[str, Any], column_rank: str,
                     reason: str) -> dict[str, Any]:
    return _node_record(
        node_id=source_node["id"], kind="taxon", subtype=None, biological=True,
        carried=False, lane_id=None, taxid=source_node["taxid"],
        rank=source_node["rank"], name=source_node["name"], path=source_node["path"],
        parent_id=source_node["parent_id"], clade=source_node["clade"],
        direct=source_node["direct"], status=source_node["status"],
        selection_reason=reason, source_order_index=source_node["source_order_index"],
        visual_order=None, origin_source_id=None, origin_target_rank=None,
        column_rank=column_rank, value=source_node["clade"], member_count=None,
        member_paths_sha256=None,
    )


def _member_digest(paths: list[str]) -> str:
    blob = "\n".join(sorted(paths, key=_utf8_key))
    return _sha256_text(blob)


def _residual_node(subtype: str, origin_source_id: str, origin_source_index: int | None,
                   origin_target_rank: str, column_rank: str, value: int,
                   member_paths: list[str] | None, carried: bool) -> dict[str, Any]:
    lane = _lane_id(subtype, origin_source_id, origin_target_rank)
    return _node_record(
        node_id=_residual_id(lane, column_rank),
        kind="residual_carry" if carried else "residual",
        subtype=subtype, biological=False, carried=carried, lane_id=lane,
        taxid=None, rank=None, name=(f"Other {column_rank}" if subtype == "other_hidden"
                                     else f"Assigned above {column_rank}"),
        path=None, parent_id=None, clade=None, direct=None, status=None,
        selection_reason="carry" if carried else "residual",
        source_order_index=origin_source_index, visual_order=None,
        origin_source_id=origin_source_id, origin_target_rank=origin_target_rank,
        column_rank=column_rank, value=value,
        member_count=(len(member_paths) if member_paths is not None else None),
        member_paths_sha256=(_member_digest(member_paths) if member_paths is not None else None),
    )


def _link(source_id: str, target_id: str, kind: str, value: int,
          transition_index: int) -> dict[str, Any]:
    return {
        "id": _link_id(source_id, target_id, kind),
        "source": source_id,
        "target": target_id,
        "kind": kind,
        "value": value,
        "transition_index": transition_index,
    }


def _is_descendant(node: dict[str, Any], ancestor_path: str) -> bool:
    return node["path"] == ancestor_path or node["path"].startswith(ancestor_path + ";")


def _assign_visual_order(nodes: list[dict[str, Any]], source: dict[str, Any], ranks: list[str]) -> None:
    by_rank: dict[str, list[dict[str, Any]]] = {rank: [] for rank in ranks}
    for node in nodes:
        if node["kind"] == "taxon":
            by_rank[node["column_rank"]].append(node)
    for rank in ranks:
        biological = sorted(by_rank[rank], key=lambda node: node["source_order_index"])
        bio_index = {node["id"]: index for index, node in enumerate(biological)}
        path_order = sorted((_utf8_key(node["path"]) for node in biological))
        items: list[tuple[tuple[Any, ...], dict[str, Any]]] = []
        for node in biological:
            items.append(((bio_index[node["id"]], 1, 0, "", node["id"]), node))
        for node in nodes:
            if node["column_rank"] != rank or node["kind"] not in {"residual", "residual_carry"}:
                continue
            origin_id = node["origin_source_id"]
            if origin_id == ENTRY_ID:
                slot = len(biological)
            else:
                origin = next((item for item in source["nodes"] if item["id"] == origin_id), None)
                if origin is None:
                    fail(f"residual '{node['id']}' has an unknown origin source")
                descendants = [
                    bio_index[item["id"]]
                    for item in biological
                    if _is_descendant(item, origin["path"])
                ]
                if descendants:
                    slot = max(descendants) + 1
                else:
                    slot = bisect.bisect_right(path_order, _utf8_key(origin["path"]))
            subtype_order = SUBTYPE_ORDER[node["subtype"]]
            items.append(((slot, 0, subtype_order, node["lane_id"], node["id"]), node))
        items.sort(key=lambda item: item[0])
        for visual_order, (_, node) in enumerate(items):
            node["visual_order"] = visual_order


def _canonical_nodes(nodes: list[dict[str, Any]], ranks: list[str]) -> list[dict[str, Any]]:
    def key(node: dict[str, Any]) -> tuple[Any, ...]:
        column = -1 if node["kind"] == "entry" else RANK_INDEX[node["column_rank"]]
        kind = NODE_KINDS[node["kind"]]
        anchor = node["source_order_index"]
        if anchor is None:
            anchor = -1
        return column, kind, anchor, node["id"]
    return sorted(nodes, key=key)


def _canonical_links(links: list[dict[str, Any]], nodes: list[dict[str, Any]],
                    ranks: list[str]) -> list[dict[str, Any]]:
    node_order = {node["id"]: index for index, node in enumerate(nodes)}
    seen: set[tuple[str, str, str]] = set()
    for link in links:
        triple = (link["source"], link["target"], link["kind"])
        if triple in seen:
            fail(f"duplicate link source/target/kind: {triple}")
        seen.add(triple)
        if link["value"] <= 0:
            fail(f"zero-valued link is not allowed: {link['id']}")
        if link["source"] not in node_order or link["target"] not in node_order:
            fail(f"link has an unknown endpoint: {link['id']}")
    return sorted(links, key=lambda link: (
        link["transition_index"], node_order[link["source"]],
        node_order[link["target"]], LINK_KINDS[link["kind"]], link["id"],
    ))


def _assert_conservation(source: dict[str, Any], ranks: list[str], nodes: list[dict[str, Any]],
                         links: list[dict[str, Any]], entry: dict[str, Any]) -> dict[str, Any]:
    classified = source["classified"]
    by_id = {node["id"]: node for node in nodes}
    if classified == 0:
        if nodes or links:
            fail("zero-classified views must not emit graph nodes or links")
        return {
            "classified": 0,
            "column_totals": [{"rank": rank, "value": 0} for rank in ranks],
            "transition_totals": [{
                "from": "ENTRY" if i == 0 else ranks[i - 1],
                "to": rank,
                "value": 0,
            } for i, rank in enumerate(ranks)],
            "rightmost_flow": 0,
        }
    if entry["value"] != classified:
        fail("entry value does not equal classified reads")
    if sum(link["value"] for link in links if link["transition_index"] == 0) != classified:
        fail("entry transition does not conserve classified reads")
    column_totals: list[dict[str, Any]] = []
    for rank in ranks:
        value = sum(node["value"] for node in nodes if node["column_rank"] == rank)
        if value != classified:
            fail(f"column {rank} total {value} does not equal classified reads")
        column_totals.append({"rank": rank, "value": value})

    transition_totals: list[dict[str, Any]] = []
    for index, rank in enumerate(ranks):
        value = sum(link["value"] for link in links if link["transition_index"] == index)
        if value != classified:
            fail(f"transition {index} total {value} does not equal classified reads")
        transition_totals.append({
            "from": "ENTRY" if index == 0 else ranks[index - 1],
            "to": rank,
            "value": value,
        })

    outgoing: dict[str, list[dict[str, Any]]] = {}
    incoming: dict[str, list[dict[str, Any]]] = {}
    for link in links:
        outgoing.setdefault(link["source"], []).append(link)
        incoming.setdefault(link["target"], []).append(link)
    for node in nodes:
        if node["kind"] == "taxon" and node["column_rank"] != ranks[-1]:
            actual = sum(link["value"] for link in outgoing.get(node["id"], []))
            if actual != node["clade"]:
                fail(f"biological source '{node['id']}' does not conserve its clade")
        elif node["kind"] in {"residual", "residual_carry"} and node["column_rank"] != ranks[-1]:
            outgoing_links = outgoing.get(node["id"], [])
            if len(outgoing_links) != 1 or outgoing_links[0]["value"] != node["value"]:
                fail(f"residual lane '{node['id']}' does not have one equal carry")
    for node in nodes:
        if node["kind"] == "taxon" and node["column_rank"] != ranks[0]:
            biological_in = sum(
                link["value"] for link in incoming.get(node["id"], [])
                if link["kind"] == "biological"
            )
            if biological_in != node["clade"]:
                fail(f"biological target '{node['id']}' has incorrect incoming flow")
    rightmost = sum(node["value"] for node in nodes if node["column_rank"] == ranks[-1])
    if rightmost != classified:
        fail("rightmost flow does not equal classified reads")
    return {
        "classified": classified,
        "column_totals": column_totals,
        "transition_totals": transition_totals,
        "rightmost_flow": rightmost,
    }


def build_view(source: dict[str, Any], ranks: list[str], max_n: int) -> dict[str, Any]:
    classified = source["classified"]
    if classified == 0:
        entry = {"id": ENTRY_ID, "value": 0}
        conservation = _assert_conservation(source, ranks, [], [], entry)
        return {"nodes": [], "links": [], "conservation": conservation}

    retained, reasons = _selection(source, ranks, max_n)
    nodes: list[dict[str, Any]] = []
    links: list[dict[str, Any]] = []
    entry_node = _entry_node(classified)
    nodes.append(entry_node)
    bio_by_rank: dict[str, list[dict[str, Any]]] = {}
    residuals: list[dict[str, Any]] = []

    first_rank = ranks[0]
    first_bio = [
        _biological_node(source["by_path"][path], first_rank, reasons[path])
        for path in retained[first_rank]
    ]
    first_bio.sort(key=lambda node: node["source_order_index"])
    bio_by_rank[first_rank] = first_bio
    nodes.extend(first_bio)
    all_first = source["by_rank"][first_rank]
    visible_first = set(retained[first_rank])
    hidden_first = [node for node in all_first if node["path"] not in visible_first]
    hidden_mass = sum(node["clade"] for node in hidden_first)
    entry_above = classified - sum(node["clade"] for node in all_first)
    if entry_above < 0:
        fail("first selected rank accounts for more than classified reads")
    for subtype, value, member_paths in (
        ("other_hidden", hidden_mass, [node["path"] for node in hidden_first]),
        ("assigned_above", entry_above, None),
    ):
        if value <= 0:
            continue
        residual = _residual_node(
            subtype, ENTRY_ID, None, first_rank, first_rank, value,
            member_paths, False,
        )
        residuals.append(residual)
        nodes.append(residual)
        links.append(_link(ENTRY_ID, residual["id"], subtype, value, 0))
    for target in first_bio:
        links.append(_link(ENTRY_ID, target["id"], "biological", target["value"], 0))

    current_residuals = residuals
    for transition_index, (source_rank, target_rank) in enumerate(zip(ranks, ranks[1:]), start=1):
        target_bio = [
            _biological_node(source["by_path"][path], target_rank, reasons[path])
            for path in retained[target_rank]
        ]
        target_bio.sort(key=lambda node: node["source_order_index"])
        bio_by_rank[target_rank] = target_bio
        nodes.extend(target_bio)
        target_by_path = {node["path"]: node for node in target_bio}
        new_residuals: list[dict[str, Any]] = []
        for source_node in bio_by_rank[source_rank]:
            source_path = source_node["path"]
            source_record = source["by_path"][source_path]
            all_targets = [
                node for node in source["by_rank"][target_rank]
                if _is_descendant(node, source_path)
            ]
            visible_targets = [node for node in all_targets if node["path"] in target_by_path]
            hidden_targets = [node for node in all_targets if node["path"] not in target_by_path]
            visible_mass = sum(node["clade"] for node in visible_targets)
            all_target_mass = sum(node["clade"] for node in all_targets)
            assigned_above = source_record["clade"] - all_target_mass
            if assigned_above < 0:
                fail(f"source '{source_path}' has negative assigned-above mass")
            for target in visible_targets:
                links.append(_link(source_node["id"], target_by_path[target["path"]]["id"],
                                   "biological", target["clade"], transition_index))
            for subtype, value, member_paths in (
                ("other_hidden", sum(node["clade"] for node in hidden_targets),
                 [node["path"] for node in hidden_targets]),
                ("assigned_above", assigned_above, None),
            ):
                if value <= 0:
                    continue
                residual = _residual_node(
                    subtype, source_node["id"], source_node["source_order_index"],
                    target_rank, target_rank, value, member_paths, False,
                )
                new_residuals.append(residual)
                nodes.append(residual)
                links.append(_link(source_node["id"], residual["id"], subtype,
                                   value, transition_index))
        for prior in current_residuals:
            carry = _residual_node(
                prior["subtype"], prior["origin_source_id"], prior["source_order_index"],
                prior["origin_target_rank"], target_rank, prior["value"],
                None if prior["subtype"] == "assigned_above" else [], True,
            )
            if prior["subtype"] == "other_hidden":
                carry["member_count"] = prior["member_count"]
                carry["member_paths_sha256"] = prior["member_paths_sha256"]
            new_residuals.append(carry)
            nodes.append(carry)
            links.append(_link(prior["id"], carry["id"], "carry", prior["value"], transition_index))
        current_residuals = new_residuals

    _assign_visual_order(nodes, source, ranks)
    canonical_nodes = _canonical_nodes(nodes, ranks)
    canonical_links = _canonical_links(links, canonical_nodes, ranks)
    entry = {"id": ENTRY_ID, "value": classified}
    conservation = _assert_conservation(source, ranks, canonical_nodes, canonical_links, entry)
    return {"nodes": canonical_nodes, "links": canonical_links, "conservation": conservation}


def build_document(payload: dict[str, Any], source_bytes: bytes, ranks: list[str],
                   max_n: int, renderer_version: str) -> dict[str, Any]:
    source = _validate_payload(payload)
    view = build_view(source, ranks, max_n)
    source_nodes: list[dict[str, Any]] = []
    for node in source["nodes"]:
        source_nodes.append({
            "id": node["id"],
            "taxid": node["taxid"],
            "rank": node["rank"],
            "name": node["name"],
            "path": node["path"],
            "parent_id": node["parent_id"],
            "clade": node["clade"],
            "direct": node["direct"],
            "status": node["status"],
            "source_order_index": node["source_order_index"],
        })
    return {
        "schema_version": 1,
        "sample_id": source["sample_id"],
        "renderer": "builtin_taxonomy_sankey",
        "renderer_version": renderer_version,
        "source_renderer": source["renderer"],
        "source_schema_version": source["schema_version"],
        "source_payload_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "count_model": COUNT_MODEL,
        "denominator": "TotalReads",
        "totals": {
            "total": source["total"],
            "classified": source["classified"],
            "unclassified": source["unclassified"],
        },
        "defaults": {"ranks": ranks, "max_taxa_per_rank": max_n},
        "rank_order": ranks,
        "source_nodes": source_nodes,
        "entry": {"id": ENTRY_ID, "value": source["classified"]},
        "default_view": view,
    }


def _client_source(path: Path) -> str:
    return _read_bytes(path, "Sankey client").decode("utf-8")


def build_html(payload_bytes: bytes, sample_id: str, client_text: str) -> bytes:
    encoded = base64.b64encode(payload_bytes).decode("ascii")
    title = html.escape(sample_id, quote=True)
    csp = ("default-src 'none'; script-src 'unsafe-inline'; style-src 'unsafe-inline'; "
           "img-src data:; connect-src 'none'; font-src 'none'; object-src 'none'; "
           "base-uri 'none'; form-action 'none'")
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta http-equiv="Content-Security-Policy" content="{csp}">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Taxonomy Sankey — {title}</title>
<style>
body {{ font-family: system-ui, sans-serif; margin: 1rem; color: #202124; }}
header {{ border-bottom: 1px solid #d8dbe0; margin-bottom: .75rem; }}
.controls {{ display:flex; flex-wrap:wrap; gap:.75rem; align-items:center; }}
.summary {{ display:flex; flex-wrap:wrap; gap:1rem; margin:.75rem 0; }}
#status {{ font-weight:600; }}
svg {{ width:100%; min-height:30rem; border:1px solid #d8dbe0; background:#fff; }}
.node rect {{ stroke:#263238; stroke-width:1; }}
.node text {{ font-size:12px; dominant-baseline:middle; }}
.link {{ fill:none; stroke-opacity:.35; }}
.residual rect {{ fill:#eceff1; }}
.carry rect {{ fill:#f5f5f5; stroke-dasharray:3 2; }}
.unresolved rect {{ fill:#fff3cd; }}
.conflicted rect {{ fill:#f8d7da; }}
button, select, input {{ font:inherit; }}
</style>
</head>
<body>
<header><h1>Taxonomy Sankey — {title}</h1>
<p>Original offline taxonomy-flow viewer. The sibling <code>.kreport</code> remains the official Pavian interoperability artifact.</p></header>
<section class="summary" aria-label="Read totals">
<span>Total: <strong id="total"></strong></span>
<span>Classified: <strong id="classified"></strong></span>
<span>Unclassified: <strong id="unclassified"></strong></span>
</section>
<section class="controls" aria-label="Sankey controls">
<label>Ranks <select id="ranks" multiple size="1"></select></label>
<label>Top N <input id="max-n" type="number" min="1" max="100" step="1"></label>
<label><input id="percent" type="checkbox"> percentages</label>
<input id="search" type="search" placeholder="Highlight taxa">
<button id="reset" type="button">Reset verified default</button>
<a id="download-json" download="{title}.sankey.json">Download canonical JSON</a>
<button id="download-svg" type="button">Download current SVG</button>
</section>
<p id="status" role="status"></p>
<svg id="sankey" role="img" aria-label="Taxonomy flow diagram" viewBox="0 0 1200 640"></svg>
<script id="wf16s-payload" type="application/octet-stream">{encoded}</script>
<script>{client_text}</script>
</body>
</html>
"""
    return document.encode("utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload", required=True, type=Path)
    parser.add_argument("--json-out", required=True, type=Path)
    parser.add_argument("--html-out", type=Path)
    parser.add_argument("--ranks", required=True)
    parser.add_argument("--max-taxa-per-rank", required=True, type=int)
    parser.add_argument("--client-js", type=Path,
                        default=Path(__file__).with_name("taxonomy_sankey_client.js"))
    parser.add_argument("--renderer-version", default=None)
    args = parser.parse_args(argv)
    try:
        payload, source_bytes = _read_json(args.payload, "source Pavian JSON")
        ranks = parse_ranks(args.ranks)
        max_n = args.max_taxa_per_rank
        renderer_version = args.renderer_version or _version()
        _text(renderer_version, "renderer version")
        document = build_document(payload, source_bytes, ranks, max_n, renderer_version)
        json_bytes = canonical_json_bytes(document)
        client_text = _client_source(args.client_js)
        forbidden = (
            "fetch(", "XMLHttpRequest", "WebSocket", "navigator.sendBeacon",
            "EventSource", "importScripts(", "navigator.serviceWorker",
            "document.cookie", "localStorage", "sessionStorage",
        )
        if any(token in client_text for token in forbidden):
            fail("Sankey client contains a forbidden network capability")
        _atomic_write_new(args.json_out, json_bytes)
        if args.html_out is not None:
            _atomic_write_new(args.html_out, build_html(json_bytes, str(payload["sample_id"]), client_text))
        return 0
    except Exception as exc:
        print(f"taxonomy_sankey_renderer: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
