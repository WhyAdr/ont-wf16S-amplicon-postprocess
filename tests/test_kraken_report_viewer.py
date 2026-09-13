import importlib.util
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


VIEWER = load_module("wf16s_kraken_report_viewer", "analysis/utils/kraken_report_viewer.py")
VERIFY = load_module(
    "wf16s_kraken_viewer_correspondence",
    "analysis/utils/verify_kraken_viewer_correspondence.py",
)


HEADER = "SampleID\tDepth\tRankCode\tNodeName\tTaxonPath\tTaxID\tStatus\tResolutionSource"


class KrakenViewerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def write_inputs(self, report, resolution):
        report_path = self.root / "sample.kreport"
        resolution_path = self.root / "taxonomy_resolution.tsv"
        report_path.write_text(report, encoding="utf-8", newline="")
        resolution_path.write_text(resolution, encoding="utf-8", newline="")
        return report_path, resolution_path

    def test_build_and_independently_verify_json_and_html(self):
        report = (
            "10.00\t10\t10\tU\t0\tunclassified\n"
            "90.00\t90\t0\tR\t1\troot\n"
            "90.00\t90\t0\tD\t2\t  Bacteria\n"
            "90.00\t90\t90\tK\t3\t    Bacilli\n"
        )
        resolution = (
            HEADER + "\n"
            "S1\t1\tD\tBacteria\tBacteria\t2\tConflicted\tassignment_conflict\n"
            "S1\t2\tK\tBacilli\tBacteria;Bacilli\t3\tResolved\tcache\n"
        )
        report_path, resolution_path = self.write_inputs(report, resolution)
        json_path = self.root / "S1.pavian.json"
        html_path = self.root / "S1.pavian.html"
        self.assertEqual(
            VIEWER.main([
                "--kreport", str(report_path),
                "--resolution-tsv", str(resolution_path),
                "--sample-id", "S1",
                "--json-out", str(json_path),
                "--html-out", str(html_path),
                "--expected-total", "100",
            ]),
            0,
        )
        VERIFY.verify(report_path, resolution_path, "S1", json_path, html_path, 100)
        payload = json_path.read_text(encoding="utf-8")
        self.assertIn('"status":"conflicted"', payload)
        html_text = html_path.read_text(encoding="utf-8")
        self.assertIn("connect-src 'none'", html_text)
        self.assertNotIn("innerHTML", html_text)
        self.assertNotIn("<img", html_text.lower())

    def test_all_unclassified_accepts_header_only_resolution_sidecar(self):
        report, resolution = self.write_inputs(
            "100.00\t100\t100\tU\t0\tunclassified\n"
            "0.00\t0\t0\tR\t1\troot\n",
            HEADER + "\n",
        )
        json_path = self.root / "S1.pavian.json"
        self.assertEqual(
            VIEWER.main([
                "--kreport", str(report), "--resolution-tsv", str(resolution),
                "--sample-id", "S1", "--json-out", str(json_path),
                "--expected-total", "100",
            ]),
            0,
        )
        VERIFY.verify(report, resolution, "S1", json_path, None, 100)
        self.assertIn('"nodes":[]', json_path.read_text(encoding="utf-8"))

    def test_rejects_malformed_rank_and_sidecar_mismatch(self):
        report, resolution = self.write_inputs(
            "100.00\t100\t100\tU\t0\tunclassified\n"
            "0.00\t0\t0\tR\t1\troot\n"
            "0.00\t0\t0\tK\t2\t  BadRank\n",
            HEADER + "\nS1\t2\tK\tBadRank\tBadRank\t2\tUnresolved\tcache\n",
        )
        result = VIEWER.main([
            "--kreport", str(report), "--resolution-tsv", str(resolution),
            "--sample-id", "S1", "--json-out", str(self.root / "bad.json"),
            "--expected-total", "100",
        ])
        self.assertEqual(result, 2)

    def test_sample_and_taxon_labels_are_encoded_not_executed(self):
        report, resolution = self.write_inputs(
            "100.00\t100\t100\tU\t0\tunclassified\n"
            "0.00\t0\t0\tR\t1\troot\n"
            "0.00\t0\t0\tD\t2\t  <img src=x>\n",
            HEADER + "\nS<&\t1\tD\t<img src=x>\t<img src=x>\t2\tUnresolved\tcache\n",
        )
        json_path = self.root / "x.json"
        html_path = self.root / "x.html"
        self.assertEqual(
            VIEWER.main([
                "--kreport", str(report), "--resolution-tsv", str(resolution),
                "--sample-id", "S<&", "--json-out", str(json_path),
                "--html-out", str(html_path), "--expected-total", "100",
            ]),
            0,
        )
        html_text = html_path.read_text(encoding="utf-8")
        self.assertIn("S&lt;&amp;", html_text)
        self.assertNotIn("<img src=x>", html_text)
        VERIFY.verify(report, resolution, "S<&", json_path, html_path, 100)


if __name__ == "__main__":
    unittest.main()
