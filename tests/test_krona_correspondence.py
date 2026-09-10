import importlib.util
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "analysis"
    / "utils"
    / "verify_krona_correspondence.py"
)
SPEC = importlib.util.spec_from_file_location("verify_krona_correspondence", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


BASE_HTML = """<!doctype html><html><body>
<krona><datasets><dataset>sample.krona</dataset></datasets>
<node name="S1"><magnitude><val>6</val></magnitude><magnitudeUnassigned><val>0</val></magnitudeUnassigned>
<node name="Unclassified"><magnitude><val>3</val></magnitude><magnitudeUnassigned><val>3</val></magnitudeUnassigned></node>
<node name="Bacteria"><magnitude><val>3</val></magnitude><magnitudeUnassigned><val>0</val></magnitudeUnassigned>
<node name="Firmicutes"><magnitude><val>2</val></magnitude><magnitudeUnassigned><val>2</val></magnitudeUnassigned></node>
<node name="Proteobacteria"><magnitude><val>1</val></magnitude><magnitudeUnassigned><val>1</val></magnitudeUnassigned></node>
</node></node></krona></body></html>
"""


class KronaCorrespondenceTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root = Path(self.temp_dir.name)
        self.tsv = self.root / "sample.krona.tsv"
        self.html = self.root / "sample.krona.html"
        self.tsv.write_text(
            "3\tUnclassified\n"
            "2\tBacteria\tFirmicutes\n"
            "1\tBacteria\tProteobacteria\n",
            encoding="utf-8",
        )

    def tearDown(self):
        self.temp_dir.cleanup()

    def verify(self, html=BASE_HTML):
        self.html.write_text(html, encoding="utf-8")
        return MODULE.verify_krona_correspondence(self.html, self.tsv, "S1", 6)

    def test_exact_correspondence_passes(self):
        self.assertEqual(self.verify(), 6)

    def test_semantic_mutations_fail_closed(self):
        mutations = {
            "changed label": BASE_HTML.replace("Firmicutes", "Actinomycetota"),
            "changed direct count": BASE_HTML.replace(
                "<magnitudeUnassigned><val>2</val></magnitudeUnassigned>",
                "<magnitudeUnassigned><val>1</val></magnitudeUnassigned>",
                1,
            ),
            "missing path": BASE_HTML.replace(
                '<node name="Proteobacteria"><magnitude><val>1</val></magnitude><magnitudeUnassigned><val>1</val></magnitudeUnassigned></node>',
                "",
            ),
            "duplicate sibling": BASE_HTML.replace(
                "</node></node></krona>",
                '<node name="Firmicutes"><magnitude><val>0</val></magnitude><magnitudeUnassigned><val>0</val></magnitudeUnassigned></node></node></node></krona>',
            ),
            "dataset identity": BASE_HTML.replace("sample.krona", "other.krona"),
            "sample identity": BASE_HTML.replace('node name="S1"', 'node name="S2"', 1),
            "unclassified allocation": BASE_HTML.replace(
                'node name="Unclassified"', 'node name="Unclassified reads"', 1
            ),
        }
        for name, html in mutations.items():
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.verify(html)


if __name__ == "__main__":
    unittest.main()
