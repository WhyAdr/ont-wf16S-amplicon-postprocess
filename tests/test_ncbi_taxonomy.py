import hashlib
import importlib.util
import json
import os
import pathlib
import csv
import tempfile
import unittest
from unittest import mock
from urllib import request

ROOT = pathlib.Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "analysis" / "utils" / "ncbi_taxonomy.py"
SPEC = importlib.util.spec_from_file_location("ncbi_taxonomy", MODULE_PATH)
taxonomy = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(taxonomy)


class TaxonomyResolverTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.work = pathlib.Path(self.tempdir.name)
        self.lineage = "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus subtilis"
        self.abundance = self.work / "abundance.tsv"
        self.abundance.write_text(
            "tax\tS1\ttotal\n"
            "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown\t1\t1\n"
            f"{self.lineage}\t2\t2\n",
            encoding="utf-8",
        )
        self.assignment = self.work / "assignment.tsv"
        self.assignment.write_text(
            "C\tread1\t1423\t0|1500\tBacteria|Bacillota|Bacilli|Bacillales|Bacillaceae|Bacillus|Bacillus subtilis\n"
            "C\tread2\t1423\t1501\tBacteria|Bacillota|Bacilli|Bacillales|Bacillaceae|Bacillus|Bacillus subtilis\n"
            "U\tread3\t0\t1490\tUnclassified\n",
            encoding="utf-8",
        )
        parts = self.lineage.split(";")
        self.cache = self.work / "cache.json"
        cache_payload = {";".join(parts[:depth]): depth for depth in range(1, 8)}
        cache_payload[self.lineage] = 0
        self.cache.write_text(json.dumps(cache_payload), encoding="utf-8")

    def tearDown(self):
        self.tempdir.cleanup()

    def args(self, mode="cache_only"):
        return [
            "ncbi_taxonomy.py", "--abundance", str(self.abundance),
            "--assignments", str(self.assignment), "--cache", str(self.cache),
            "--resolved-cache", str(self.work / "resolved.json"), "--mode", mode,
            "--unresolved-policy", "warn", "--unresolved-tsv", str(self.work / "unresolved.tsv"),
            "--conflicts-tsv", str(self.work / "conflicts.tsv"),
            "--resolution-sources-tsv", str(self.work / "resolution_sources.tsv"),
            "--provenance", str(self.work / "provenance.json"),
        ]

    def test_cache_only_uses_assignments_without_mutating_source_cache(self):
        before = hashlib.sha256(self.cache.read_bytes()).hexdigest()
        with mock.patch.object(request, "urlopen") as urlopen, \
             mock.patch("sys.argv", self.args()):
            self.assertEqual(taxonomy.main(), 0)
        urlopen.assert_not_called()
        after = hashlib.sha256(self.cache.read_bytes()).hexdigest()
        self.assertEqual(before, after)
        resolved = json.loads((self.work / "resolved.json").read_text(encoding="utf-8"))
        self.assertEqual(resolved[self.lineage], "1423")
        self.assertEqual((self.work / "unresolved.tsv").read_text(encoding="utf-8").count("\n"), 1)
        sources = (self.work / "resolution_sources.tsv").read_text(encoding="utf-8")
        self.assertIn(f"{self.lineage}\t1423\tassignment", sources)
        provenance = json.loads((self.work / "provenance.json").read_text(encoding="utf-8"))
        self.assertEqual(provenance["resolution_source_counts"]["assignment"], 1)
        self.assertEqual(provenance["resolution_source_counts"]["unresolved"], 0)
        self.assertEqual(provenance["source_cache_sha256_before"], before)
        self.assertEqual(provenance["source_cache_sha256_committed"], before)
        self.assertRegex(provenance["source_cache_sha256_candidate"], r"^[0-9a-f]{64}$")

    def test_refresh_failure_preserves_source_cache_and_returns_nonzero(self):
        before = hashlib.sha256(self.cache.read_bytes()).hexdigest()
        with mock.patch.dict(os.environ, {"NCBI_EMAIL": "test@example.org"}), \
             mock.patch.object(taxonomy, "read_assignment_taxids", return_value=({}, [])), \
             mock.patch.object(taxonomy, "query_exact_scientific_name", return_value=(0, "simulated failure", [])), \
             mock.patch("sys.argv", self.args(mode="refresh")):
            self.assertEqual(taxonomy.main(), 1)
        after = hashlib.sha256(self.cache.read_bytes()).hexdigest()
        self.assertEqual(before, after)

    def test_assignment_tie_break_is_deterministic(self):
        self.assignment.write_text(
            "C\tread1\t200\t1500\tBacteria|Example\n"
            "C\tread2\t100\t1500\tBacteria|Example\n",
            encoding="utf-8",
        )
        resolved, conflicts = taxonomy.read_assignment_taxids([self.assignment])
        self.assertEqual(resolved["Bacteria|Example"], 100)
        self.assertEqual(conflicts[0]["winner_taxid"], 100)

    def test_assignment_read_ids_reset_between_files(self):
        second = self.work / "assignment-second.tsv"
        second.write_text(
            "C\tread1\t1423\t1500\tBacteria|Bacillota\n", encoding="utf-8"
        )
        resolved, conflicts = taxonomy.read_assignment_taxids([self.assignment, second])
        self.assertEqual(resolved["Bacteria|Bacillota|Bacilli|Bacillales|Bacillaceae|Bacillus|Bacillus subtilis"], 1423)
        self.assertEqual(conflicts, [])

    def test_assignment_rejects_zero_read_length(self):
        self.assignment.write_text(
            "C\tread1\t1423\t0\tBacteria|Bacillota\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(ValueError, r"0 < length"):
            taxonomy.read_assignment_taxids([self.assignment])

    def test_expected_input_hash_is_checked_before_read(self):
        with self.assertRaisesRegex(ValueError, "input changed"):
            taxonomy.read_abundance_paths(self.abundance, "tax", "0" * 64)

    def test_ambiguous_exact_name_query_is_not_silently_selected(self):
        payload = json.dumps({"esearchresult": {"idlist": ["22", "11", "22"]}}).encode()
        response = mock.MagicMock()
        response.__enter__.return_value.read.return_value = payload
        response.__exit__.return_value = False
        with mock.patch.object(request, "urlopen", return_value=response):
            taxid, error, ambiguous = taxonomy.query_exact_scientific_name(
                "Example", "test@example.org", None, attempts=1
            )
        self.assertEqual(taxid, 0)
        self.assertIsNone(error)
        self.assertEqual(ambiguous, [11, 22])

    def test_exact_name_query_validates_rank_and_ancestry(self):
        esearch = json.dumps({"esearchresult": {"idlist": ["11"]}}).encode()
        efetch = b"""<TaxaSet><Taxon><TaxId>11</TaxId><ScientificName>Bacillus</ScientificName><Rank>genus</Rank><LineageEx><Taxon><ScientificName>Bacteria</ScientificName></Taxon></LineageEx></Taxon></TaxaSet>"""

        def response(payload):
            item = mock.MagicMock()
            item.__enter__.return_value.read.return_value = payload
            item.__exit__.return_value = False
            return item

        with mock.patch.object(request, "urlopen", side_effect=[response(esearch), response(efetch)]):
            taxid, error, ambiguous = taxonomy.query_exact_scientific_name(
                "Bacillus", "test@example.org", None, attempts=1,
                expected_rank="genus", ancestor_names=["Bacteria"]
            )
        self.assertEqual(taxid, 11)
        self.assertIsNone(error)
        self.assertEqual(ambiguous, [])

    def test_refresh_records_all_resolution_source_labels(self):
        parts = self.lineage.split(";")
        cache_payload = {";".join(parts[:depth]): depth for depth in range(1, 7)}
        cache_payload[";".join(parts[:7])] = 0
        cache_payload[self.lineage] = 0
        self.cache.write_text(json.dumps(cache_payload), encoding="utf-8")
        with mock.patch.dict(os.environ, {"NCBI_EMAIL": "test@example.org"}), \
             mock.patch.object(
                 taxonomy,
                 "query_exact_scientific_name",
                 side_effect=[(777, None, [])],
             ), \
             mock.patch("sys.argv", self.args(mode="refresh")):
            self.assertEqual(taxonomy.main(), 0)
        with (self.work / "resolution_sources.tsv").open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        by_path = {row["TaxonPath"]: row["ResolutionSource"] for row in rows}
        self.assertIn("source_cache", by_path.values())
        self.assertEqual(by_path[";".join(parts[:7])], "ncbi_refresh")
        self.assertEqual(by_path[self.lineage], "assignment")

    def test_unresolved_source_is_reported(self):
        self.assignment.write_text("U\tread1\t0\t1500\tUnclassified\n", encoding="utf-8")
        with mock.patch("sys.argv", self.args()):
            self.assertEqual(taxonomy.main(), 0)
        sources = (self.work / "resolution_sources.tsv").read_text(encoding="utf-8")
        self.assertIn(f"{self.lineage}\t0\tunresolved", sources)

    def test_ordered_ancestry_validation(self):
        record = {
            "scientific_name": "Bacillus subtilis",
            "rank": "species",
            "lineage": ["Bacteria", "Bacillota", "Bacilli", "Bacillales", "Bacillaceae", "Bacillus"],
        }
        # Positive case
        err = taxonomy.validate_taxonomy_context(
            record, "Bacillus subtilis", expected_rank="species",
            ancestor_names=["Bacteria", "Bacillota", "Bacillus"]
        )
        self.assertIsNone(err)

        # Reversed order
        err = taxonomy.validate_taxonomy_context(
            record, "Bacillus subtilis", expected_rank="species",
            ancestor_names=["Bacillus", "Bacteria"]
        )
        self.assertIsNotNone(err)
        self.assertIn("out of order or reversed", err)

        # Missing ancestor
        err = taxonomy.validate_taxonomy_context(
            record, "Bacillus subtilis", expected_rank="species",
            ancestor_names=["Archaea"]
        )
        self.assertIsNotNone(err)
        self.assertIn("missing ancestor context", err)

        # Wrong rank
        err = taxonomy.validate_taxonomy_context(
            record, "Bacillus subtilis", expected_rank="genus",
            ancestor_names=["Bacteria"]
        )
        self.assertIsNotNone(err)
        self.assertIn("rank mismatch", err)

        # Wrong name
        err = taxonomy.validate_taxonomy_context(
            record, "Escherichia coli", expected_rank="species",
            ancestor_names=["Bacteria"]
        )
        self.assertIsNotNone(err)
        self.assertIn("does not exactly match", err)

        # Duplicate ancestor ambiguity
        record_dup = {
            "scientific_name": "Bacillus subtilis",
            "rank": "species",
            "lineage": ["Bacteria", "Bacillus", "Bacillaceae", "Bacillus"],
        }
        err = taxonomy.validate_taxonomy_context(
            record_dup, "Bacillus subtilis", expected_rank="species",
            ancestor_names=["Bacteria", "Bacillus"]
        )
        self.assertIsNotNone(err)
        self.assertIn("ambiguity", err)

    def test_concurrent_cache_mutation_aborts_refresh(self):
        parts = self.lineage.split(";")
        cache_payload = {";".join(parts[:depth]): depth for depth in range(1, 7)}
        cache_payload[";".join(parts[:7])] = 0
        cache_payload[self.lineage] = 0
        self.cache.write_text(json.dumps(cache_payload), encoding="utf-8")

        def mutate_cache(*args, **kwargs):
            self.cache.write_text(json.dumps({"external_mutation": 9999}), encoding="utf-8")
            return 777, None, []

        with mock.patch.dict(os.environ, {"NCBI_EMAIL": "test@example.org"}), \
             mock.patch.object(taxonomy, "query_exact_scientific_name", side_effect=mutate_cache), \
             mock.patch("sys.argv", self.args(mode="refresh")):
            exit_code = taxonomy.main()
            self.assertEqual(exit_code, 1)

        # External mutation must be preserved, not overwritten by candidate
        cache_content = json.loads(self.cache.read_text(encoding="utf-8"))
        self.assertIn("external_mutation", cache_content)
        self.assertEqual(cache_content["external_mutation"], 9999)

    def test_conflicting_expected_input_rejected(self):
        args = self.args() + [
            "--expected-input", f"{self.abundance}\t{'0'*64}",
            "--expected-input", f"{self.abundance}\t{'1'*64}",
        ]
        with mock.patch("sys.argv", args):
            self.assertEqual(taxonomy.main(), 1)

    def test_cache_lock_timeout_returns_busy(self):
        parts = self.lineage.split(";")
        cache_payload = {";".join(parts[:depth]): depth for depth in range(1, 7)}
        cache_payload[";".join(parts[:7])] = 0
        self.cache.write_text(json.dumps(cache_payload), encoding="utf-8")

        # Hold the lock externally
        with taxonomy.acquire_cache_lock(str(self.cache)):
            # Nested attempt should fail with E_TAXONOMY_CACHE_BUSY
            with self.assertRaises(SystemExit) as cm:
                with taxonomy.acquire_cache_lock(str(self.cache), timeout=0.2, poll_interval=0.05):
                    pass
            self.assertIn("E_TAXONOMY_CACHE_BUSY", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
