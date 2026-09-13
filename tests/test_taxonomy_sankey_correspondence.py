from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from analysis.utils.taxonomy_sankey_renderer import (
    build_document,
    build_html,
    canonical_json_bytes,
)
from analysis.utils.verify_taxonomy_sankey_correspondence import (
    CorrespondenceError,
    verify,
)
from tests.test_taxonomy_sankey_renderer import build_fixture


class TaxonomySankeyCorrespondenceTests(unittest.TestCase):
    def test_independent_reference_model_accepts_json_and_html(self):
        payload, source_bytes = build_fixture(unicode_name=True)
        document = build_document(payload, source_bytes, ["D", "K", "P"], 2, "0.4.8")
        client_path = Path("analysis/utils/taxonomy_sankey_client.js")
        client_text = client_path.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            payload_path = root / "sample.pavian.json"
            json_path = root / "sample.sankey.json"
            html_path = root / "sample.sankey.html"
            payload_path.write_bytes(canonical_json_bytes(payload))
            json_path.write_bytes(canonical_json_bytes(document))
            html_path.write_bytes(build_html(json_path.read_bytes(), "S1", client_text))
            verify(payload_path, json_path, html_path, ["D", "K", "P"], 2, "0.4.8")

    def test_reference_rejects_changed_value_and_missing_ancestor(self):
        payload, source_bytes = build_fixture()
        document = build_document(payload, source_bytes, ["D", "K", "P"], 2, "0.4.8")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            payload_path = root / "sample.pavian.json"
            json_path = root / "sample.sankey.json"
            payload_path.write_bytes(canonical_json_bytes(payload))
            changed = dict(document)
            changed["default_view"] = dict(document["default_view"])
            changed["default_view"]["nodes"] = list(document["default_view"]["nodes"])
            changed["default_view"]["nodes"][0] = dict(changed["default_view"]["nodes"][0])
            changed["default_view"]["nodes"][0]["value"] = 99
            json_path.write_bytes(canonical_json_bytes(changed))
            with self.assertRaisesRegex(CorrespondenceError, "does not match"):
                verify(payload_path, json_path, None, ["D", "K", "P"], 2, "0.4.8")

            payload["nodes"][1]["parent_path"] = "missing"
            payload_path.write_bytes(canonical_json_bytes(payload))
            json_path.write_bytes(canonical_json_bytes(document))
            with self.assertRaisesRegex(CorrespondenceError, "path/name mismatch|missing source ancestor"):
                verify(payload_path, json_path, None, ["D", "K", "P"], 2, "0.4.8")


if __name__ == "__main__":
    unittest.main()
