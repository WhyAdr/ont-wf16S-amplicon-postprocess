from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from analysis.utils.taxonomy_sankey_renderer import (
    ENTRY_ID,
    SankeyError,
    build_document,
    canonical_json_bytes,
)


RANK_DEPTH = {"D": 1, "K": 2, "P": 3, "C": 4, "O": 5, "F": 6, "G": 7, "S": 8}


def make_payload(*, zero: bool = False, unicode_name: bool = False) -> dict:
    if zero:
        return {
            "schema_version": 1,
            "sample_id": "S1",
            "renderer": "builtin_kraken_report_explorer",
            "renderer_version": "0.4.8",
            "totals": {"total": 100, "classified": 0, "unclassified": 100},
            "nodes": [],
        }

    names = {"A": "Álpha" if unicode_name else "Alpha", "B": "Beta", "C": "Gamma"}
    specs = [
        ("D", names["A"], "", 60, 0),
        ("K", "A1", names["A"], 40, 0),
        ("P", "A1a", f"{names['A']};A1", 20, 20),
        ("P", "A1b", f"{names['A']};A1", 20, 20),
        ("K", "A2", names["A"], 20, 0),
        ("P", "A2a", f"{names['A']};A2", 20, 20),
        ("D", names["B"], "", 25, 0),
        ("K", "B1", names["B"], 25, 0),
        ("P", "B1a", f"{names['B']};B1", 25, 25),
        ("D", names["C"], "", 15, 15),
    ]
    # The path is supplied separately so names can contain Unicode without
    # changing the deterministic source order fixture.
    nodes = []
    for rank, name, parent_path, clade, direct in specs:
        path = f"{parent_path};{name}" if parent_path else name
        nodes.append({
            "path": path,
            "parent_path": parent_path,
            "name": name,
            "depth": RANK_DEPTH[rank],
            "rank_code": rank,
            "taxid": str(100 + len(nodes)),
            "direct": direct,
            "clade": clade,
            "status": "resolved",
        })
    return {
        "schema_version": 1,
        "sample_id": "S1",
        "renderer": "builtin_kraken_report_explorer",
        "renderer_version": "0.4.8",
        "official_pavian_compatibility": "kraken_report_input_contract_only",
        "totals": {"total": 100, "classified": 100, "unclassified": 0},
        "nodes": nodes,
    }


def build_fixture(**kwargs):
    payload = make_payload(**kwargs)
    source_bytes = canonical_json_bytes(payload)
    return payload, source_bytes


def build_stress_fixture():
    ranks = ["D", "K", "P", "C", "O", "F", "G", "S"]
    nodes = []
    parents = [""]
    for rank_index, rank in enumerate(ranks):
        next_parents = []
        clade = 3 ** (len(ranks) - rank_index - 1)
        for parent_path in parents:
            for child_index in range(3):
                name = f"{rank}{child_index}"
                path = f"{parent_path};{name}" if parent_path else name
                nodes.append({
                    "path": path,
                    "parent_path": parent_path,
                    "name": name,
                    "depth": rank_index + 1,
                    "rank_code": rank,
                    "taxid": str(100000 + len(nodes)),
                    "direct": 1 if rank == "S" else 0,
                    "clade": 1 if rank == "S" else clade,
                    "status": "resolved",
                })
                next_parents.append(path)
        parents = next_parents
    assert len(nodes) == 9840
    payload = {
        "schema_version": 1,
        "sample_id": "stress",
        "renderer": "builtin_kraken_report_explorer",
        "renderer_version": "0.4.8",
        "totals": {"total": len(parents), "classified": len(parents), "unclassified": 0},
        "nodes": nodes,
    }
    return payload, canonical_json_bytes(payload)


class TaxonomySankeyRendererTests(unittest.TestCase):
    def test_first_rank_top_n_and_persistent_residual_carries(self):
        payload, source_bytes = build_fixture()
        document = build_document(payload, source_bytes, ["D", "K", "P"], 2, "0.4.8")
        view = document["default_view"]
        self.assertEqual(document["entry"], {"id": ENTRY_ID, "value": 100})
        self.assertEqual(view["conservation"]["column_totals"], [
            {"rank": "D", "value": 100},
            {"rank": "K", "value": 100},
            {"rank": "P", "value": 100},
        ])
        residuals = [node for node in view["nodes"] if node["kind"] == "residual"]
        self.assertEqual({node["subtype"] for node in residuals}, {"other_hidden"})
        carry_nodes = [node for node in view["nodes"] if node["kind"] == "residual_carry"]
        self.assertGreaterEqual(len(carry_nodes), 3)
        self.assertTrue(all(node["carried"] for node in carry_nodes))
        self.assertEqual(view["conservation"]["rightmost_flow"], 100)
        self.assertTrue(all(link["value"] > 0 for link in view["links"]))

    def test_direct_counts_and_ancestral_closure(self):
        payload, source_bytes = build_fixture()
        document = build_document(payload, source_bytes, ["D", "P"], 1, "0.4.8")
        taxa = [node for node in document["default_view"]["nodes"] if node["kind"] == "taxon"]
        self.assertEqual({node["rank"] for node in taxa}, {"D", "P"})
        self.assertTrue(any(node["selection_reason"] == "ancestor_closure" for node in taxa))
        self.assertEqual(document["default_view"]["conservation"]["rightmost_flow"], 100)

    def test_assigned_above_is_a_persistent_lane(self):
        payload, source_bytes = build_fixture()
        document = build_document(payload, source_bytes, ["D", "K", "P"], 3, "0.4.8")
        view = document["default_view"]
        assigned = [node for node in view["nodes"]
                    if node["subtype"] == "assigned_above" and node["kind"] == "residual"]
        self.assertEqual(len(assigned), 1)
        self.assertEqual(assigned[0]["value"], 15)
        self.assertTrue(any(link["kind"] == "carry" for link in view["links"]))
        self.assertEqual(view["conservation"]["rightmost_flow"], 100)

    def test_unicode_and_tie_serialization_are_repeatable(self):
        payload, source_bytes = build_fixture(unicode_name=True)
        first = canonical_json_bytes(build_document(payload, source_bytes, ["D", "K", "P"], 2, "0.4.8"))
        second = canonical_json_bytes(build_document(payload, source_bytes, ["D", "K", "P"], 2, "0.4.8"))
        self.assertEqual(first, second)
        self.assertIn("Álpha".encode("utf-8"), first)
        self.assertNotIn(b"\\u00", first)

    def test_zero_classified_is_empty_but_keeps_accounting(self):
        payload, source_bytes = build_fixture(zero=True)
        document = build_document(payload, source_bytes, ["D", "K"], 10, "0.4.8")
        self.assertEqual(document["entry"]["value"], 0)
        self.assertEqual(document["default_view"]["nodes"], [])
        self.assertEqual(document["default_view"]["links"], [])
        self.assertEqual(document["default_view"]["conservation"]["rightmost_flow"], 0)

    def test_invalid_missing_ancestor_and_unsafe_count_fail_closed(self):
        payload, source_bytes = build_fixture()
        payload["nodes"][1]["parent_path"] = "missing"
        with self.assertRaisesRegex(SankeyError, "path does not match|missing ancestor"):
            build_document(payload, canonical_json_bytes(payload), ["D", "K"], 2, "0.4.8")

        payload, source_bytes = build_fixture()
        payload["nodes"][0]["clade"] = 9007199254740992
        with self.assertRaisesRegex(SankeyError, "classified reads|direct counts|clade"):
            build_document(payload, canonical_json_bytes(payload), ["D", "K"], 2, "0.4.8")

    def test_destination_is_not_overwritten(self):
        payload, source_bytes = build_fixture()
        document = build_document(payload, source_bytes, ["D", "K"], 2, "0.4.8")
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "payload.json"
            path.write_bytes(b"existing")
            with self.assertRaisesRegex(SankeyError, "overwrite"):
                # Exercise the public writer through the CLI-like output path
                # without depending on a browser runtime.
                from analysis.utils.taxonomy_sankey_renderer import _atomic_write_new
                _atomic_write_new(path, canonical_json_bytes(document))

    def test_generated_9840_node_taxonomy_conserves_across_rank_and_top_n_views(self):
        payload, source_bytes = build_stress_fixture()
        self.assertEqual(len(payload["nodes"]), 9840)
        rank_sets = [
            ["D", "K"], ["D", "P"], ["D", "C"], ["D", "G"],
            ["K", "P"], ["K", "F"], ["P", "G"], ["C", "S"],
            ["D", "K", "P"], ["D", "P", "G"],
        ]
        for ranks in rank_sets:
            for max_n in (1, 2):
                document = build_document(payload, source_bytes, ranks, max_n, "0.4.8")
                self.assertEqual(document["default_view"]["conservation"]["rightmost_flow"], len(payload["nodes"][-6561:]))
                self.assertTrue(all(link["value"] > 0 for link in document["default_view"]["links"]))


if __name__ == "__main__":
    unittest.main()
