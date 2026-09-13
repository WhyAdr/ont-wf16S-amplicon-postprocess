import base64
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


VIEWER = load_module("viewer_for_correspondence_fixture", "analysis/utils/kraken_report_viewer.py")
VERIFY = load_module(
    "independent_correspondence_verifier",
    "analysis/utils/verify_kraken_viewer_correspondence.py",
)

HEADER = "SampleID\tDepth\tRankCode\tNodeName\tTaxonPath\tTaxID\tStatus\tResolutionSource"


class CorrespondenceMutationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.report = self.root / "S1.kreport"
        self.sidecar = self.root / "taxonomy_resolution.tsv"
        self.json = self.root / "S1.pavian.json"
        self.html = self.root / "S1.pavian.html"
        self.report.write_text(
            "50.00\t50\t50\tU\t0\tunclassified\n"
            "50.00\t50\t0\tR\t1\troot\n"
            "50.00\t50\t50\tD\t2\t  Bacteria\n",
            encoding="utf-8",
        )
        self.sidecar.write_text(
            HEADER + "\nS1\t1\tD\tBacteria\tBacteria\t2\tResolved\tcache\n",
            encoding="utf-8",
        )
        self.assertEqual(
            VIEWER.main([
                "--kreport", str(self.report), "--resolution-tsv", str(self.sidecar),
                "--sample-id", "S1", "--json-out", str(self.json),
                "--html-out", str(self.html), "--expected-total", "100",
            ]),
            0,
        )

    def tearDown(self):
        self.temp.cleanup()

    def verify(self):
        VERIFY.verify(self.report, self.sidecar, "S1", self.json, self.html, 100)

    def test_valid_correspondence(self):
        self.verify()

    def test_json_mutation_is_rejected(self):
        payload = json.loads(self.json.read_text(encoding="utf-8"))
        payload["nodes"][0]["taxid"] = "9"
        self.json.write_bytes(VIEWER.canonical_json_bytes(payload))
        with self.assertRaises(VERIFY.CorrespondenceError):
            self.verify()

    def test_sidecar_mutations_are_rejected(self):
        self.sidecar.write_text(
            HEADER + "\nS1\t1\tD\tBacteria\tBacteria\t2\tConflicted\tcache\n",
            encoding="utf-8",
        )
        with self.assertRaises(VERIFY.CorrespondenceError):
            self.verify()

    def test_embedded_payload_mutation_is_rejected(self):
        html_text = self.html.read_text(encoding="utf-8")
        marker = '<script id="wf16s-payload" type="application/octet-stream">'
        start = html_text.index(marker) + len(marker)
        end = html_text.index("</script>", start)
        payload = bytearray(base64.b64decode(html_text[start:end]))
        payload[-2] = ord("0") if payload[-2] != ord("0") else ord("1")
        mutated = base64.b64encode(bytes(payload)).decode("ascii")
        self.html.write_text(html_text[:start] + mutated + html_text[end:], encoding="utf-8")
        with self.assertRaises(VERIFY.CorrespondenceError):
            self.verify()

    def test_network_capability_mutation_is_rejected(self):
        self.html.write_text(
            self.html.read_text(encoding="utf-8").replace("connect-src 'none'", "connect-src *"),
            encoding="utf-8",
        )
        with self.assertRaises(VERIFY.CorrespondenceError):
            self.verify()

    def test_verifier_has_no_producer_import(self):
        source = (ROOT / "analysis/utils/verify_kraken_viewer_correspondence.py").read_text(encoding="utf-8")
        self.assertNotRegex(source, r"(?:from|import)\s+.*kraken_report_viewer")


if __name__ == "__main__":
    unittest.main()
