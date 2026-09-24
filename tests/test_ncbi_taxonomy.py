import hashlib
import importlib.util
import json
import os
import pathlib
import csv
import tempfile
import threading
import unittest
from http import client
from unittest import mock
from urllib import error, request

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
             mock.patch.object(taxonomy, "query_exact_scientific_name", return_value=taxonomy.LookupOutcome(
                 "request_failed", code="E_NCBI_REQUEST", message="simulated failure")), \
             mock.patch("sys.argv", self.args(mode="refresh")):
            self.assertEqual(taxonomy.main(), 1)
        after = hashlib.sha256(self.cache.read_bytes()).hexdigest()
        self.assertEqual(before, after)

    def test_deferred_refresh_writes_candidate_without_mutating_source_cache(self):
        before = hashlib.sha256(self.cache.read_bytes()).hexdigest()
        args = self.args(mode="refresh") + [
            "--defer-cache-commit", "--cache-lock-held",
            "--cache-lock-owner-pid", str(os.getpid()),
        ]
        with mock.patch.dict(os.environ, {"NCBI_EMAIL": "test@example.org"}), \
             mock.patch("sys.argv", args):
            self.assertEqual(taxonomy.main(), 0)
        self.assertEqual(hashlib.sha256(self.cache.read_bytes()).hexdigest(), before)
        resolved = json.loads((self.work / "resolved.json").read_text(encoding="utf-8"))
        self.assertEqual(resolved[self.lineage], "1423")
        provenance = json.loads((self.work / "provenance.json").read_text(encoding="utf-8"))
        self.assertTrue(provenance["source_cache_commit_deferred"])
        self.assertFalse(provenance["source_cache_updated"])
        self.assertEqual(provenance["source_cache_sha256_committed"], before)

    def test_assignment_tie_break_is_deterministic(self):
        self.assignment.write_text(
            "C\tread1\t200\t1500\tBacteria|Example\n"
            "C\tread2\t100\t1500\tBacteria|Example\n",
            encoding="utf-8",
        )
        resolved, conflicts = taxonomy.read_assignment_taxids([self.assignment])
        self.assertEqual(resolved["Bacteria|Example"], 100)
        self.assertEqual(conflicts[0]["winner_taxid"], 100)

    def test_assignment_conflict_is_explicit_in_resolution_provenance(self):
        assignment_lineage = "|".join(
            [self.lineage.split(";")[0]] + self.lineage.split(";")[2:]
        )
        self.assignment.write_text(
            f"C\tread1\t200\t1500\t{assignment_lineage}\n"
            f"C\tread2\t100\t1500\t{assignment_lineage}\n",
            encoding="utf-8",
        )
        with mock.patch("sys.argv", self.args()):
            self.assertEqual(taxonomy.main(), 0)
        sources = (self.work / "resolution_sources.tsv").read_text(encoding="utf-8")
        self.assertIn(f"{self.lineage}\t100\tassignment_conflict", sources)
        provenance = json.loads((self.work / "provenance.json").read_text(encoding="utf-8"))
        self.assertEqual(provenance["resolution_source_counts"]["assignment_conflict"], 1)

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
        payload = json.dumps({"esearchresult": {"count": "2", "idlist": ["22", "11"]}}).encode()
        response = mock.MagicMock()
        response.__enter__.return_value.read.return_value = payload
        response.__exit__.return_value = False
        with mock.patch.object(request, "urlopen", return_value=response):
            outcome = taxonomy.query_exact_scientific_name(
                "Example", "test@example.org", None, attempts=1
            )
        self.assertEqual(outcome.status, "ambiguous")
        self.assertEqual(outcome.candidates, (11, 22))

    def test_exact_name_query_validates_rank_and_ancestry(self):
        esearch = json.dumps({"esearchresult": {"count": "1", "idlist": ["11"]}}).encode()
        efetch = b"""<TaxaSet><Taxon><TaxId>11</TaxId><ScientificName>Bacillus</ScientificName><Rank>genus</Rank><LineageEx><Taxon><ScientificName>Bacteria</ScientificName></Taxon></LineageEx></Taxon></TaxaSet>"""

        def response(payload):
            item = mock.MagicMock()
            item.__enter__.return_value.read.return_value = payload
            item.__exit__.return_value = False
            return item

        with mock.patch.object(request, "urlopen", side_effect=[response(esearch), response(efetch)]):
            outcome = taxonomy.query_exact_scientific_name(
                "Bacillus", "test@example.org", None, attempts=1,
                expected_rank="genus", ancestor_names=["Bacteria"]
            )
        self.assertEqual(outcome.status, "resolved")
        self.assertEqual(outcome.taxid, 11)
        self.assertEqual(outcome.rank_rule, "exact")

    def test_current_cellular_domain_rank_is_narrowly_accepted(self):
        for name in ("Bacteria", "Archaea", "Eukaryota"):
            record = {"scientific_name": name, "rank": "domain", "lineage": []}
            code, message, rule = taxonomy.validate_taxonomy_context(
                record, name, expected_rank="superkingdom", ancestor_names=[]
            )
            self.assertIsNone(code)
            self.assertIsNone(message)
            self.assertEqual(rule, "cellular_domain_legacy_alias")

        rejected = [
            ({"scientific_name": "Viruses", "rank": "domain", "lineage": []},
             "Viruses", "superkingdom"),
            ({"scientific_name": "Bacteria", "rank": "realm", "lineage": []},
             "Bacteria", "superkingdom"),
            ({"scientific_name": "Bacillus", "rank": "domain", "lineage": ["Bacteria"]},
             "Bacillus", "genus"),
        ]
        for record, name, rank in rejected:
            code, message, rule = taxonomy.validate_taxonomy_context(
                record, name, expected_rank=rank, ancestor_names=[]
            )
            self.assertEqual(code, "E_TAXONOMY_CONTEXT")
            self.assertIn("rank mismatch", message)
            self.assertIsNone(rule)

    def test_missing_search_fields_are_fatal_response_invalid(self):
        response = mock.MagicMock()
        response.__enter__.return_value.read.return_value = b'{"esearchresult": {}}'
        response.__exit__.return_value = False
        with mock.patch.object(request, "urlopen", return_value=response):
            outcome = taxonomy.query_exact_scientific_name(
                "Example", "test@example.org", None, attempts=1
            )
        self.assertEqual(outcome.status, "response_invalid")
        self.assertEqual(outcome.code, "E_NCBI_RESPONSE")

    def test_interrupted_http_body_retries_then_writes_failure_summary_through_main(self):
        before = hashlib.sha256(self.cache.read_bytes()).hexdigest()
        diagnostics = self.work / "interrupted diagnostics"
        diagnostics.mkdir()
        transaction_id = "tx-" + "b" * 64
        args = self.args(mode="refresh") + [
            "--diagnostics-dir", str(diagnostics), "--transaction-id", transaction_id
        ]

        def interrupted_response():
            response = mock.MagicMock()
            response.__enter__.return_value = response
            response.__exit__.return_value = False
            response.read.side_effect = client.IncompleteRead(b"partial", 99)
            return response

        with mock.patch.dict(os.environ, {"NCBI_EMAIL": "test@example.org"}), \
             mock.patch.object(taxonomy, "read_assignment_taxids", return_value=({}, [])), \
             mock.patch.object(request, "urlopen", side_effect=[
                 interrupted_response(), interrupted_response(), interrupted_response()
             ]) as urlopen, \
             mock.patch.object(taxonomy.time, "sleep"), \
             mock.patch("sys.argv", args):
            self.assertEqual(taxonomy.main(), 1)

        self.assertEqual(urlopen.call_count, taxonomy.MAX_REQUEST_ATTEMPTS)
        self.assertEqual(hashlib.sha256(self.cache.read_bytes()).hexdigest(), before)
        summary = json.loads((diagnostics / "taxonomy_failure.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["code"], "E_NCBI_REQUEST")
        self.assertEqual(summary["outcome"], "request_failed")
        self.assertIn("IncompleteRead", summary["message"])
        events = (diagnostics / "taxonomy_events.jsonl").read_text(encoding="utf-8")
        self.assertEqual(events.count('"event": "request_retry"'), 2)
        self.assertIn('"event": "lookup_failure"', events)

    def test_unsupported_xml_encoding_is_response_invalid_through_main(self):
        before = hashlib.sha256(self.cache.read_bytes()).hexdigest()
        diagnostics = self.work / "encoding diagnostics"
        diagnostics.mkdir()
        transaction_id = "tx-" + "c" * 64
        args = self.args(mode="refresh") + [
            "--diagnostics-dir", str(diagnostics), "--transaction-id", transaction_id
        ]
        esearch = json.dumps({
            "esearchresult": {"count": "1", "idlist": ["11"]}
        }).encode()
        efetch = b'<?xml version="1.0" encoding="x-not-installed"?><TaxaSet/>'

        def response(payload):
            item = mock.MagicMock()
            item.__enter__.return_value = item
            item.__exit__.return_value = False
            item.read.return_value = payload
            return item

        with mock.patch.dict(os.environ, {"NCBI_EMAIL": "test@example.org"}), \
             mock.patch.object(taxonomy, "read_assignment_taxids", return_value=({}, [])), \
             mock.patch.object(request, "urlopen", side_effect=[response(esearch), response(efetch)]) as urlopen, \
             mock.patch("sys.argv", args):
            self.assertEqual(taxonomy.main(), 1)

        self.assertEqual(urlopen.call_count, 2)
        self.assertEqual(hashlib.sha256(self.cache.read_bytes()).hexdigest(), before)
        summary = json.loads((diagnostics / "taxonomy_failure.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["code"], "E_NCBI_RESPONSE")
        self.assertEqual(summary["outcome"], "response_invalid")
        self.assertEqual(summary["message"], "efetch returned malformed XML")

    def test_search_count_mismatch_is_not_treated_as_unique(self):
        payload = json.dumps({"esearchresult": {"count": "2", "idlist": ["11"]}}).encode()
        response = mock.MagicMock()
        response.__enter__.return_value.read.return_value = payload
        response.__exit__.return_value = False
        with mock.patch.object(request, "urlopen", return_value=response):
            outcome = taxonomy.query_exact_scientific_name(
                "Example", "test@example.org", None, attempts=1
            )
        self.assertEqual(outcome.status, "response_invalid")
        self.assertIn("count", outcome.message)

    def test_efetch_retry_does_not_repeat_successful_esearch(self):
        esearch = json.dumps({"esearchresult": {"count": "1", "idlist": ["11"]}}).encode()
        efetch = b"""<TaxaSet><Taxon><TaxId>11</TaxId><ScientificName>Bacillus</ScientificName><Rank>genus</Rank><LineageEx><Taxon><ScientificName>Bacteria</ScientificName></Taxon></LineageEx></Taxon></TaxaSet>"""

        def response(payload):
            item = mock.MagicMock()
            item.__enter__.return_value.read.return_value = payload
            item.__exit__.return_value = False
            return item

        calls = [response(esearch), error.URLError("temporary"), response(efetch)]
        with mock.patch.object(request, "urlopen", side_effect=calls) as urlopen:
            outcome = taxonomy.query_exact_scientific_name(
                "Bacillus", "test@example.org", None, attempts=2,
                expected_rank="genus", ancestor_names=["Bacteria"], sleep=mock.Mock()
            )
        self.assertEqual(outcome.status, "resolved")
        self.assertEqual(urlopen.call_count, 3)
        endpoints = [call.args[0].full_url for call in urlopen.call_args_list]
        self.assertEqual(sum("esearch.fcgi" in url for url in endpoints), 1)
        self.assertEqual(sum("efetch.fcgi" in url for url in endpoints), 2)

    def test_top_level_api_errors_are_invalid_even_with_zero_or_valid_results(self):
        payloads = [
            {"ERROR": "simulated service failure",
             "esearchresult": {"count": "0", "idlist": []}},
            {"Error": "simulated service failure",
             "esearchresult": {"count": "1", "idlist": ["11"]}},
        ]
        for payload in payloads:
            response = mock.MagicMock()
            response.__enter__.return_value = response
            response.__exit__.return_value = False
            response.read.return_value = json.dumps(payload).encode()
            with mock.patch.object(request, "urlopen", return_value=response):
                outcome = taxonomy.query_exact_scientific_name(
                    "Example", "test@example.org", None, attempts=1
                )
            self.assertEqual(outcome.status, "response_invalid")
            self.assertEqual(outcome.code, "E_NCBI_RESPONSE")

    def test_permanent_http_failure_is_not_retried_or_leaked(self):
        failure = error.HTTPError(
            "https://example.invalid/?api_key=SECRET&email=user@example.org",
            401, "unauthorized SECRET", {}, None
        )
        with mock.patch.object(request, "urlopen", side_effect=failure) as urlopen:
            outcome = taxonomy.query_exact_scientific_name(
                "Example", "test@example.org", "SECRET", attempts=3
            )
        self.assertEqual(outcome.status, "request_failed")
        self.assertEqual(urlopen.call_count, 1)
        self.assertNotIn("SECRET", outcome.message)
        self.assertNotIn("user@example.org", outcome.message)

    def test_retry_limiter_runs_before_every_request_attempt(self):
        limiter = mock.Mock()
        limiter.wait_for_request_slot = mock.Mock()
        events = mock.Mock()
        failure = error.HTTPError("https://example.invalid", 503, "down", {}, None)
        with self.assertRaises(taxonomy.SafeRequestFailure):
            taxonomy.request_payload(
                "esearch", {}, limiter, events, attempts=3,
                opener=mock.Mock(side_effect=failure), sleep=mock.Mock()
            )
        self.assertEqual(limiter.wait_for_request_slot.call_count, 3)

    def test_local_refresh_validation_reports_pending_without_network_or_writes(self):
        args = self.args(mode="refresh") + ["--validate-only"]
        before = {path.name for path in self.work.iterdir()}
        with mock.patch.dict(os.environ, {"NCBI_EMAIL": "test@example.org"}), \
             mock.patch.object(taxonomy, "read_assignment_taxids", return_value=({}, [])), \
             mock.patch.object(request, "urlopen") as urlopen, \
             mock.patch("sys.argv", args), \
             mock.patch("builtins.print") as printer:
            self.assertEqual(taxonomy.main(), 0)
        urlopen.assert_not_called()
        self.assertEqual(before, {path.name for path in self.work.iterdir()})
        self.assertTrue(any("pending_online" in str(call) for call in printer.call_args_list))

    def test_online_preflight_fatal_request_fails_under_both_policies(self):
        for policy in ("warn", "error"):
            args = self.args(mode="refresh") + ["--validate-only", "--online-preflight"]
            policy_index = args.index("--unresolved-policy") + 1
            args[policy_index] = policy
            with mock.patch.dict(os.environ, {"NCBI_EMAIL": "test@example.org"}), \
                 mock.patch.object(taxonomy, "read_assignment_taxids", return_value=({}, [])), \
                 mock.patch.object(
                     taxonomy, "query_exact_scientific_name",
                     return_value=taxonomy.LookupOutcome(
                         "request_failed", code="E_NCBI_REQUEST", message="esearch timeout"
                     )
                 ), \
                 mock.patch("sys.argv", args):
                self.assertEqual(taxonomy.main(), 1)

    def test_online_preflight_flag_requires_validate_only(self):
        args = self.args(mode="refresh") + ["--online-preflight"]
        with mock.patch("sys.argv", args), mock.patch.object(request, "urlopen") as urlopen:
            self.assertEqual(taxonomy.main(), 1)
        urlopen.assert_not_called()

    def test_semantic_mismatch_warns_and_continues_but_error_policy_fails(self):
        parts = self.lineage.split(";")
        self.cache.write_text(json.dumps({
            ";".join(parts[:depth]): depth for depth in range(1, 7)
        }), encoding="utf-8")
        outcomes = [
            taxonomy.LookupOutcome(
                "context_mismatch", code="E_TAXONOMY_CONTEXT",
                message="expected genus, returned family",
                expected_rank="genus", returned_rank="family"
            ),
            taxonomy.LookupOutcome(
                "resolved", taxid=777, expected_rank="species",
                returned_rank="species", rank_rule="exact"
            ),
        ]
        with mock.patch.dict(os.environ, {"NCBI_EMAIL": "test@example.org"}), \
             mock.patch.object(taxonomy, "read_assignment_taxids", return_value=({}, [])), \
             mock.patch.object(taxonomy, "query_exact_scientific_name", side_effect=outcomes), \
             mock.patch("sys.argv", self.args(mode="refresh")):
            self.assertEqual(taxonomy.main(), 0)
        provenance = json.loads((self.work / "provenance.json").read_text(encoding="utf-8"))
        self.assertEqual(len(provenance["taxonomy_mismatches"]), 1)
        self.assertEqual(provenance["query_failures"], [])
        resolved = json.loads((self.work / "resolved.json").read_text(encoding="utf-8"))
        self.assertEqual(resolved[self.lineage], "777")
        self.assertNotIn(";".join(parts[:7]), resolved)

        self.cache.write_text(json.dumps({
            ";".join(parts[:depth]): depth for depth in range(1, 7)
        }), encoding="utf-8")
        error_args = self.args(mode="refresh")
        error_args[error_args.index("--unresolved-policy") + 1] = "error"
        with mock.patch.dict(os.environ, {"NCBI_EMAIL": "test@example.org"}), \
             mock.patch.object(taxonomy, "read_assignment_taxids", return_value=({}, [])), \
             mock.patch.object(taxonomy, "query_exact_scientific_name", side_effect=outcomes), \
             mock.patch("sys.argv", error_args):
            self.assertEqual(taxonomy.main(), 1)

    def test_fatal_lookup_stops_before_later_nodes(self):
        self.cache.write_text("{}", encoding="utf-8")
        query = mock.Mock(return_value=taxonomy.LookupOutcome(
            "response_invalid", code="E_NCBI_RESPONSE", message="malformed payload"
        ))
        with mock.patch.dict(os.environ, {"NCBI_EMAIL": "test@example.org"}), \
             mock.patch.object(taxonomy, "read_assignment_taxids", return_value=({}, [])), \
             mock.patch.object(taxonomy, "query_exact_scientific_name", query), \
             mock.patch("sys.argv", self.args(mode="refresh")):
            self.assertEqual(taxonomy.main(), 1)
        self.assertEqual(query.call_count, 1)

    def test_execution_failure_diagnostics_survive_and_redact(self):
        diagnostics = self.work / "diagnostics"
        diagnostics.mkdir()
        transaction_id = "tx-" + "a" * 64
        args = self.args(mode="refresh") + [
            "--diagnostics-dir", str(diagnostics), "--transaction-id", transaction_id
        ]
        with mock.patch.dict(os.environ, {
                 "NCBI_EMAIL": "private@example.org", "NCBI_API_KEY": "SECRETKEY"
             }), \
             mock.patch.object(taxonomy, "read_assignment_taxids", return_value=({}, [])), \
             mock.patch.object(
                 taxonomy, "query_exact_scientific_name",
                 return_value=taxonomy.LookupOutcome(
                     "request_failed", code="E_NCBI_REQUEST",
                     message="https://example.invalid/?api_key=SECRETKEY&email=private@example.org"
                 )
             ), \
             mock.patch("sys.argv", args):
            self.assertEqual(taxonomy.main(), 1)
        combined = "\n".join(path.read_text(encoding="utf-8") for path in diagnostics.iterdir())
        self.assertIn("E_NCBI_REQUEST", combined)
        self.assertNotIn("SECRETKEY", combined)
        self.assertNotIn("private@example.org", combined)

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
                 side_effect=[taxonomy.LookupOutcome("resolved", taxid=777,
                                                     expected_rank="genus",
                                                     returned_rank="genus",
                                                     rank_rule="exact")],
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
        code, err, rule = taxonomy.validate_taxonomy_context(
            record, "Bacillus subtilis", expected_rank="species",
            ancestor_names=["Bacteria", "Bacillota", "Bacillus"]
        )
        self.assertIsNone(err)
        self.assertEqual(rule, "exact")

        # Reversed order
        _, err, _ = taxonomy.validate_taxonomy_context(
            record, "Bacillus subtilis", expected_rank="species",
            ancestor_names=["Bacillus", "Bacteria"]
        )
        self.assertIsNotNone(err)
        self.assertIn("out of order or reversed", err)

        # Missing ancestor
        _, err, _ = taxonomy.validate_taxonomy_context(
            record, "Bacillus subtilis", expected_rank="species",
            ancestor_names=["Archaea"]
        )
        self.assertIsNotNone(err)
        self.assertIn("missing ancestor context", err)

        # Wrong rank
        _, err, _ = taxonomy.validate_taxonomy_context(
            record, "Bacillus subtilis", expected_rank="genus",
            ancestor_names=["Bacteria"]
        )
        self.assertIsNotNone(err)
        self.assertIn("rank mismatch", err)

        # Wrong name
        _, err, _ = taxonomy.validate_taxonomy_context(
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
        _, err, _ = taxonomy.validate_taxonomy_context(
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
            return taxonomy.LookupOutcome("resolved", taxid=777,
                                          expected_rank="genus", returned_rank="genus",
                                          rank_rule="exact")

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

        lock_path = pathlib.Path(f"{self.cache}.lock")
        self.assertTrue(lock_path.is_file())
        self.assertEqual(lock_path.read_bytes(), b"")
        with taxonomy.acquire_cache_lock(str(self.cache), timeout=0.2):
            self.assertTrue(lock_path.is_file())

    def test_interrupted_lock_acquisition_cleans_process_guard(self):
        lock_identity = os.path.normcase(os.path.realpath(f"{self.cache}.lock"))
        with mock.patch.object(taxonomy, "_try_lock_fd", side_effect=OSError("busy")), \
             mock.patch.object(taxonomy.time, "sleep", side_effect=KeyboardInterrupt):
            with self.assertRaises(KeyboardInterrupt):
                with taxonomy.acquire_cache_lock(str(self.cache), timeout=1.0, poll_interval=0.01):
                    pass

        self.assertNotIn(lock_identity, taxonomy._ACTIVE_CACHE_LOCKS)
        with taxonomy.acquire_cache_lock(str(self.cache), timeout=0.2):
            self.assertIn(lock_identity, taxonomy._ACTIVE_CACHE_LOCKS)

    def test_non_os_error_during_lock_poll_cleans_process_guard(self):
        lock_identity = os.path.normcase(os.path.realpath(f"{self.cache}.lock"))
        with mock.patch.object(taxonomy, "_try_lock_fd", side_effect=OSError("busy")), \
             mock.patch.object(taxonomy.time, "sleep", side_effect=RuntimeError("injected poll failure")):
            with self.assertRaisesRegex(RuntimeError, "injected poll failure"):
                with taxonomy.acquire_cache_lock(str(self.cache), timeout=1.0, poll_interval=0.01):
                    pass

        self.assertNotIn(lock_identity, taxonomy._ACTIVE_CACHE_LOCKS)
        with taxonomy.acquire_cache_lock(str(self.cache), timeout=0.2):
            pass

    def test_relative_and_absolute_cache_aliases_share_lock_identity(self):
        previous = os.getcwd()
        os.chdir(self.work.parent)
        try:
            relative = os.path.relpath(self.cache, self.work.parent)
            self.assertEqual(
                os.path.normcase(os.path.realpath(os.path.abspath(relative) + ".lock")),
                os.path.normcase(os.path.realpath(os.path.abspath(str(self.cache)) + ".lock")),
            )
            with taxonomy.acquire_cache_lock(relative, timeout=0.2):
                with self.assertRaises(SystemExit):
                    with taxonomy.acquire_cache_lock(str(self.cache), timeout=0.2):
                        pass
        finally:
            os.chdir(previous)

    def test_final_cache_symlink_shares_lock_identity(self):
        alias = self.work / "cache-alias.json"
        try:
            os.symlink(self.cache, alias)
        except (OSError, NotImplementedError) as exc:
            self.skipTest(f"cache symlink unavailable: {exc}")
        expected_lock = pathlib.Path(
            os.path.realpath(os.path.abspath(str(alias))) + ".lock"
        )
        target_lock = pathlib.Path(
            os.path.realpath(os.path.abspath(str(self.cache))) + ".lock"
        )
        self.assertEqual(
            os.path.normcase(str(expected_lock)),
            os.path.normcase(str(target_lock)),
        )
        with taxonomy.acquire_cache_lock(str(alias), timeout=0.2):
            self.assertTrue(target_lock.is_file())
            self.assertFalse(pathlib.Path(str(alias) + ".lock").exists())
            with self.assertRaises(SystemExit):
                with taxonomy.acquire_cache_lock(str(self.cache), timeout=0.2):
                    pass

    def test_parent_directory_symlink_shares_lock_identity(self):
        real_dir = self.work / "real-cache-dir"
        alias_dir = self.work / "alias-cache-dir"
        real_dir.mkdir()
        real_cache = real_dir / "cache.json"
        real_cache.write_text("{}", encoding="utf-8")
        try:
            os.symlink(real_dir, alias_dir, target_is_directory=True)
        except (OSError, NotImplementedError) as exc:
            self.skipTest(f"directory symlink unavailable: {exc}")
        alias_cache = alias_dir / "cache.json"
        self.assertEqual(
            os.path.normcase(os.path.realpath(str(alias_cache) + ".lock")),
            os.path.normcase(os.path.realpath(str(real_cache) + ".lock")),
        )

    def test_same_process_threads_do_not_poison_active_guard(self):
        entered = threading.Event()
        release = threading.Event()
        first_errors = []

        def holder():
            try:
                with taxonomy.acquire_cache_lock(str(self.cache), timeout=1.0):
                    entered.set()
                    release.wait(timeout=2.0)
            except BaseException as exc:  # pragma: no cover - diagnostic path
                first_errors.append(exc)

        thread = threading.Thread(target=holder)
        thread.start()
        self.assertTrue(entered.wait(timeout=2.0))
        with self.assertRaisesRegex(SystemExit, "E_TAXONOMY_CACHE_BUSY"):
            with taxonomy.acquire_cache_lock(str(self.cache), timeout=0.2):
                pass
        release.set()
        thread.join(timeout=2.0)
        self.assertFalse(thread.is_alive())
        self.assertEqual(first_errors, [])
        with taxonomy.acquire_cache_lock(str(self.cache), timeout=0.2):
            pass


if __name__ == "__main__":
    unittest.main()
