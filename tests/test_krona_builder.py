"""Unit tests for the stdlib-only offline Krona renderer."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET


REPO_ROOT = Path(__file__).resolve().parents[1]
BUILDER_PATH = REPO_ROOT / "analysis" / "utils" / "krona_builder.py"
VENDOR_DIR = REPO_ROOT / "analysis" / "vendor" / "krona-2.8.1"

spec = importlib.util.spec_from_file_location("krona_builder", BUILDER_PATH)
assert spec is not None and spec.loader is not None
krona_builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(krona_builder)


class KronaBuilderTests(unittest.TestCase):
    def test_pinned_vendor_manifest_and_assets_validate(self) -> None:
        manifest = krona_builder.validate_vendor(VENDOR_DIR)
        self.assertEqual(manifest["upstream"], "marbl/Krona")
        self.assertEqual(manifest["tag"], "v2.8.1")

    def test_duplicate_paths_are_aggregated_and_xml_is_escaped(self) -> None:
        records = [
            (2, ("Bacteria & friends", "Firmicutes")),
            (3, ("Bacteria & friends", "Firmicutes")),
            (1, ("Bacteria & friends", "Proteobacteria")),
        ]
        xml = krona_builder.build_krona_xml(records, "sample <one>")
        root = ET.fromstring(xml)
        self.assertEqual(root.tag, "krona")
        dataset_node = root.find("node")
        self.assertIsNotNone(dataset_node)
        self.assertEqual(dataset_node.attrib["name"], "sample <one>")
        firmicutes = root.find("node/node/node")
        self.assertIsNotNone(firmicutes)
        self.assertEqual(firmicutes.findtext("magnitude/val"), "5")
        self.assertIn(b"Bacteria &amp; friends", xml)

    def test_parser_rejects_malformed_counts_and_controls(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "bad.tsv"
            path.write_text("01\tBacteria\n", encoding="utf-8")
            with self.assertRaises(krona_builder.BuilderError):
                krona_builder.parse_krona_input(path)
            path.write_text("1\tBacteria\tbad\x00label\n", encoding="utf-8")
            with self.assertRaises(krona_builder.BuilderError):
                krona_builder.parse_krona_input(path)

    def test_render_is_deterministic_self_contained_and_atomic(self) -> None:
        with tempfile.TemporaryDirectory(prefix="krona path ") as temporary:
            root = Path(temporary)
            input_path = root / "sample input.tsv"
            output_path = root / "nested output" / "sample chart.html"
            input_path.write_text(
                "2\tBacteria\tFirmicutes\n"
                "3\tBacteria\tFirmicutes\n"
                "1\tBacteria\tProteobacteria\n",
                encoding="utf-8",
            )
            krona_builder.render(input_path, output_path, "sample <one>", 6, VENDOR_DIR)
            first = output_path.read_bytes()
            krona_builder.render(input_path, output_path, "sample <one>", 6, VENDOR_DIR)
            self.assertEqual(first, output_path.read_bytes())
            text = first.decode("utf-8")
            self.assertIn('id="hiddenImage" src="data:image/png;base64,', text)
            self.assertIn('id="loadingImage" src="data:image/gif;base64,', text)
            self.assertIn('id="logo" src="data:image/png;base64,', text)
            self.assertIn("<krona", text)
            self.assertNotIn('<script src="http', text)
            self.assertNotIn('<link rel="shortcut icon" href="http', text)
            self.assertFalse(list(output_path.parent.glob(output_path.name + ".tmp-*")))

    def test_cli_validate_only_does_not_create_output(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "never-created.html"
            result = subprocess.run(
                [
                    sys.executable,
                    str(BUILDER_PATH),
                    "--validate-only",
                    "--vendor-dir",
                    str(VENDOR_DIR),
                    "--output",
                    str(output),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
