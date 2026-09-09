"""Unit tests for the stdlib-only offline Krona renderer."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
BUILDER_PATH = REPO_ROOT / "analysis" / "utils" / "krona_builder.py"
VENDOR_DIR = REPO_ROOT / "analysis" / "vendor" / "krona-2.8.1"
KRONA_FIXTURE = REPO_ROOT / "tests" / "fixtures" / "krona" / "direct_clade.tsv"
KRONA_EXPECTED = REPO_ROOT / "tests" / "fixtures" / "krona" / "direct_clade.expected.json"

spec = importlib.util.spec_from_file_location("krona_builder", BUILDER_PATH)
assert spec is not None and spec.loader is not None
krona_builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(krona_builder)


class KronaBuilderTests(unittest.TestCase):
    def test_pinned_vendor_manifest_and_assets_validate(self) -> None:
        manifest = krona_builder.validate_vendor(VENDOR_DIR)
        self.assertEqual(manifest["upstream"], "marbl/Krona")
        self.assertEqual(manifest["tag"], "v2.8.1")

    def test_builder_version_comes_from_repository_version(self) -> None:
        expected = (REPO_ROOT / "VERSION").read_bytes().decode("utf-8").removesuffix("\n")
        self.assertEqual(krona_builder.BUILDER_VERSION, expected)

    def test_xml_preserves_direct_and_clade_magnitudes(self) -> None:
        records = krona_builder.parse_krona_input(KRONA_FIXTURE)
        expected = json.loads(KRONA_EXPECTED.read_text(encoding="utf-8"))
        root = ET.fromstring(
            krona_builder.build_krona_xml(records, "dataset", "direct_clade")
        )

        self.assertEqual(root.attrib, {
            "collapse": expected["document"]["collapse"],
            "key": expected["document"]["key"],
        })
        self.assertEqual(
            root.find("attributes").attrib["magnitude"],
            expected["document"]["magnitudeAttribute"],
        )

        self.assertEqual(
            [[attribute.attrib["display"], attribute.text]
             for attribute in root.findall("attributes/attribute")],
            expected["document"]["attributes"],
        )
        self.assertEqual(root.findtext("datasets/dataset"), expected["document"]["datasetLabel"])
        dataset = root.find("node")
        self.assertIsNotNone(dataset)
        self.assertEqual(dataset.attrib["name"], expected["document"]["rootName"])
        self.assertEqual(dataset.findtext("magnitude/val"), str(expected["dataset"]["magnitude"]))
        self.assertEqual(
            dataset.findtext("magnitudeUnassigned/val"),
            str(expected["dataset"]["magnitudeUnassigned"]),
        )
        node_a = dataset.find("node")
        self.assertIsNotNone(node_a)
        self.assertEqual(node_a.attrib["name"], "A")
        self.assertEqual(node_a.findtext("magnitude/val"), str(expected["A"]["magnitude"]))
        self.assertEqual(
            node_a.findtext("magnitudeUnassigned/val"),
            str(expected["A"]["magnitudeUnassigned"]),
        )
        self.assertEqual(
            node_a.find("node").findtext("magnitudeUnassigned/val"),
            str(expected["A/B"]["magnitudeUnassigned"]),
        )
        self.assertEqual(
            node_a.findall("node")[1].findtext("magnitudeUnassigned/val"),
            str(expected["A/C"]["magnitudeUnassigned"]),
        )

    def test_xml_is_byte_deterministic_under_input_permutation(self) -> None:
        records = [(3, ("A",)), (2, ("A", "B")), (4, ("A", "C")), (2, ("A", "B"))]
        self.assertEqual(
            krona_builder.build_krona_xml(records, "dataset"),
            krona_builder.build_krona_xml(list(reversed(records)), "dataset"),
        )

    def test_vendor_inventory_rejects_extra_and_case_colliding_files(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            vendor = Path(temporary) / "vendor"
            shutil.copytree(VENDOR_DIR, vendor)
            (vendor / "extra.txt").write_text("unexpected\n", encoding="utf-8")
            with self.assertRaises(krona_builder.BuilderError):
                krona_builder.validate_vendor(vendor)

            (vendor / "extra.txt").unlink()
            manifest_path = vendor / "SOURCE.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["files"].append({
                "path": "license.txt",
                "sha256": manifest["files"][0]["sha256"],
            })
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaises(krona_builder.BuilderError):
                krona_builder.validate_vendor(vendor)

    def test_vendor_inventory_rejects_duplicate_manifest_paths(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            vendor = Path(temporary) / "vendor"
            shutil.copytree(VENDOR_DIR, vendor)
            manifest_path = vendor / "SOURCE.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["files"].append(dict(manifest["files"][0]))
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaises(krona_builder.BuilderError):
                krona_builder.validate_vendor(vendor)

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
            self.assertIn('<krona collapse="true" key="true">', text)
            self.assertIn('<dataset>sample input</dataset>', text)
            self.assertNotIn('<script src="http', text)
            self.assertNotIn('<link rel="shortcut icon" href="http', text)
            self.assertFalse(list(output_path.parent.glob(output_path.name + ".tmp-*")))

    def test_atomic_write_failure_preserves_prior_and_leaves_no_temp(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "result.html"
            output.write_bytes(b"old")
            with mock.patch.object(krona_builder.os, "replace", side_effect=OSError("boom")):
                with self.assertRaises(krona_builder.BuilderError):
                    krona_builder.atomic_write(output, b"new")
            self.assertEqual(output.read_bytes(), b"old")
            self.assertEqual(list(root.glob(output.name + ".tmp-*")), [])

            missing = root / "missing.html"
            with mock.patch.object(krona_builder.os, "replace", side_effect=OSError("boom")):
                with self.assertRaises(krona_builder.BuilderError):
                    krona_builder.atomic_write(missing, b"new")
            self.assertFalse(missing.exists())
            self.assertEqual(list(root.glob(missing.name + ".tmp-*")), [])

    def test_atomic_write_flush_failure_leaves_no_temp(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "result.html"
            with mock.patch.object(krona_builder.os, "fsync", side_effect=OSError("flush boom")):
                with self.assertRaises(krona_builder.BuilderError):
                    krona_builder.atomic_write(output, b"new")
            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob(output.name + ".tmp-*")), [])

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
