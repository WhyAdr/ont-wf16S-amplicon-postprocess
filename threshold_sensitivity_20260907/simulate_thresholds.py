#!/usr/bin/env python3
"""Simulate relaxed wf-16s minimap2 thresholds from saved read statistics.

This is a sensitivity reconstruction, not a replacement wf-16s run. It leaves
all native workflow outputs unchanged and assigns mapped TaxID-0 reads through
their recorded best reference only when they pass the simulated thresholds.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import statistics
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent
RUNS = {
    "BAER": ROOT / "wf-16s_BAER-NCBI",
    "BANAE": ROOT / "wf-16s_BANAE-NCBI",
    "BLEA": ROOT / "wf-16s_BLEA-NCBI",
    "BGRN": ROOT / "wf-16s_BGRN-NCBI",
}
RANKS = (
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)


def one(path: Path, pattern: str) -> Path:
    matches = sorted(path.glob(pattern))
    if len(matches) != 1:
        raise ValueError(f"Expected one {pattern!r} under {path}, found {matches}")
    return matches[0]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite_float(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def fmt(value: object) -> object:
    if isinstance(value, float):
        return f"{value:.4f}"
    return value


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key, "")) for key in fields})


def median(rows: list[dict[str, object]], key: str) -> float | None:
    values = [row[key] for row in rows if isinstance(row.get(key), float)]
    return statistics.median(values) if values else None


def load_run(label: str, run_dir: Path) -> tuple[list[dict[str, object]], dict[str, object]]:
    params_path = run_dir / "output" / "params.json"
    abundance_path = run_dir / "output" / "abundance_table_species.tsv"
    assignment_path = one(run_dir / "output" / "reads_assignments", "*.assignments.tsv")
    bamstats_path = one(
        run_dir / "output" / "bams", "*.bamstats_results/bamstats.readstats.tsv.gz"
    )
    alignment_path = one(run_dir / "output" / "alignment_tables", "*-alignment-stats.tsv")

    params = json.loads(params_path.read_text(encoding="utf-8"))
    if params.get("classifier") != "minimap2":
        raise ValueError(f"{label}: expected minimap2, found {params.get('classifier')}")
    baseline_id = float(params["min_percent_identity"])
    baseline_cov = float(params["min_ref_coverage"])

    with abundance_path.open(encoding="utf-8", newline="") as handle:
        abundance_rows = list(csv.DictReader(handle, delimiter="\t"))
    if not abundance_rows:
        raise ValueError(f"{label}: empty abundance table")
    sample_columns = [
        field for field in abundance_rows[0] if field not in ("tax", "total")
    ]
    if len(sample_columns) != 1:
        raise ValueError(f"{label}: expected one abundance sample column, found {sample_columns}")
    sample_column = sample_columns[0]
    abundance_total = sum(float(row[sample_column]) for row in abundance_rows)
    abundance_unclassified = sum(
        float(row[sample_column])
        for row in abundance_rows
        if row["tax"].startswith("Unclassified;")
    )

    assignments: dict[str, tuple[str, int, str]] = {}
    with assignment_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) != 5:
                raise ValueError(f"{label}: assignment row does not have five fields")
            status, read_id, taxid, _mapping, lineage = row
            if read_id in assignments:
                raise ValueError(f"{label}: duplicate assignment read ID {read_id}")
            assignments[read_id] = (status, int(taxid), lineage)

    references: dict[str, dict[str, object]] = {}
    with alignment_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row["reference"]
            record: dict[str, object] = {"inferred_taxid": int(row["taxid"])}
            record.update({rank: row.get(rank, "") for rank in RANKS})
            previous = references.get(accession)
            if previous is not None and previous != record:
                raise ValueError(f"{label}: conflicting taxonomy for reference {accession}")
            references[accession] = record

    readstats: dict[str, dict[str, object]] = {}
    with gzip.open(bamstats_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            read_id = row["name"]
            if read_id in readstats:
                raise ValueError(f"{label}: duplicate bamstats read ID {read_id}")
            readstats[read_id] = {
                "best_reference": row["ref"],
                "identity": finite_float(row["iden"]),
                "reference_coverage": finite_float(row["ref_coverage"]),
                "query_coverage": finite_float(row["coverage"]),
                "read_length": finite_float(row["read_length"]),
                "mean_quality": finite_float(row["mean_quality"]),
            }

    if assignments.keys() != readstats.keys():
        missing_stats = assignments.keys() - readstats.keys()
        missing_assignments = readstats.keys() - assignments.keys()
        raise ValueError(
            f"{label}: read-ID mismatch; no stats={len(missing_stats)}, "
            f"no assignment={len(missing_assignments)}"
        )

    rows: list[dict[str, object]] = []
    classified_mismatches = 0
    for read_id, (status, taxid, original_lineage) in assignments.items():
        stats = readstats[read_id]
        ref = str(stats["best_reference"])
        mapped = (
            status == "C"
            and ref != "*"
            and isinstance(stats["identity"], float)
            and isinstance(stats["reference_coverage"], float)
        )
        taxonomy = references.get(ref, {}) if mapped else {}
        if mapped and not taxonomy:
            raise ValueError(f"{label}: mapped reference absent from alignment table: {ref}")
        inferred_taxid = int(taxonomy.get("inferred_taxid", 0))
        if taxid > 0 and inferred_taxid != taxid:
            classified_mismatches += 1

        identity = stats["identity"]
        ref_cov = stats["reference_coverage"]
        if status == "U":
            failure = "raw_unmapped"
        elif taxid > 0:
            failure = "classified"
        else:
            fail_id = isinstance(identity, float) and identity < baseline_id
            fail_cov = isinstance(ref_cov, float) and ref_cov < baseline_cov
            if fail_id and fail_cov:
                failure = "identity_and_reference_coverage"
            elif fail_id:
                failure = "identity_only"
            elif fail_cov:
                failure = "reference_coverage_only"
            else:
                failure = "unexplained_taxid0"

        record: dict[str, object] = {
            "dataset": label,
            "read_id": read_id,
            "raw_status": status,
            "baseline_taxid": taxid,
            "original_lineage": original_lineage,
            **stats,
            "baseline_failure_mode": failure,
            "inferred_taxid": inferred_taxid,
        }
        record.update({rank: taxonomy.get(rank, "") for rank in RANKS})
        rows.append(record)

    if classified_mismatches:
        raise ValueError(f"{label}: {classified_mismatches} classified TaxID mismatches")

    classified = [row for row in rows if int(row["baseline_taxid"]) > 0]
    c0 = [row for row in rows if row["baseline_failure_mode"] not in ("classified", "raw_unmapped")]
    if any(
        float(row["identity"]) < baseline_id
        or float(row["reference_coverage"]) < baseline_cov
        for row in classified
    ):
        raise ValueError(f"{label}: classified read below recorded baseline threshold")
    if any(row["baseline_failure_mode"] == "unexplained_taxid0" for row in c0):
        raise ValueError(f"{label}: mapped TaxID-0 read passes both baseline thresholds")
    if abundance_total != len(rows):
        raise ValueError(
            f"{label}: abundance total {abundance_total} != assignment total {len(rows)}"
        )
    if abundance_unclassified != len(rows) - len(classified):
        raise ValueError(
            f"{label}: abundance unclassified {abundance_unclassified} != "
            f"assignment unclassified {len(rows) - len(classified)}"
        )

    provenance = {
        "dataset": label,
        "run_directory": str(run_dir.relative_to(ROOT)),
        "classifier": params.get("classifier"),
        "database_set": params.get("database_set"),
        "baseline_min_percent_identity": baseline_id,
        "baseline_min_ref_coverage": baseline_cov,
        "min_len": params.get("min_len"),
        "max_len": params.get("max_len"),
        "min_read_qual": params.get("min_read_qual"),
        "input_files": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in (
                params_path,
                abundance_path,
                assignment_path,
                bamstats_path,
                alignment_path,
            )
        },
    }
    return rows, provenance


def main() -> None:
    all_rows: list[dict[str, object]] = []
    provenance: list[dict[str, object]] = []
    breakdown_rows: list[dict[str, object]] = []
    grid_rows: list[dict[str, object]] = []
    selected_rows: list[dict[str, object]] = []
    taxa_rows: list[dict[str, object]] = []
    diagnostics: list[dict[str, object]] = []

    selected = {
        "baseline_90_90": (90, 90),
        "identity_89": (89, 90),
        "identity_88": (88, 90),
        "identity_87": (87, 90),
        "identity_85": (85, 90),
        "refcov_89": (90, 89),
        "refcov_88": (90, 88),
        "refcov_87": (90, 87),
        "refcov_85": (90, 85),
        "balanced_89": (89, 89),
        "balanced_88": (88, 88),
        "balanced_87": (87, 87),
        "balanced_85": (85, 85),
    }

    report_sections: list[str] = []
    for label, run_dir in RUNS.items():
        rows, run_provenance = load_run(label, run_dir)
        all_rows.extend(rows)
        provenance.append(run_provenance)
        total = len(rows)
        baseline_classified = sum(int(row["baseline_taxid"]) > 0 for row in rows)
        baseline_unclassified = total - baseline_classified
        modes = Counter(str(row["baseline_failure_mode"]) for row in rows)
        mapped_c0 = [
            row
            for row in rows
            if int(row["baseline_taxid"]) == 0 and row["raw_status"] == "C"
        ]
        diagnostics.extend(row for row in rows if int(row["baseline_taxid"]) == 0)

        breakdown = {
            "dataset": label,
            "total_reads": total,
            "baseline_classified": baseline_classified,
            "baseline_classified_pct": 100 * baseline_classified / total,
            "baseline_unclassified": baseline_unclassified,
            "baseline_unclassified_pct": 100 * baseline_unclassified / total,
            "raw_unmapped": modes["raw_unmapped"],
            "mapped_taxid0": len(mapped_c0),
            "mapped_taxid0_pct_of_unclassified": 100 * len(mapped_c0) / baseline_unclassified,
            "identity_only": modes["identity_only"],
            "reference_coverage_only": modes["reference_coverage_only"],
            "identity_and_reference_coverage": modes["identity_and_reference_coverage"],
            "identity_implicated": modes["identity_only"] + modes["identity_and_reference_coverage"],
            "reference_coverage_implicated": modes["reference_coverage_only"]
            + modes["identity_and_reference_coverage"],
            "median_classified_length": median(
                [row for row in rows if int(row["baseline_taxid"]) > 0], "read_length"
            ),
            "median_mapped_taxid0_length": median(mapped_c0, "read_length"),
            "median_classified_quality": median(
                [row for row in rows if int(row["baseline_taxid"]) > 0], "mean_quality"
            ),
            "median_mapped_taxid0_quality": median(mapped_c0, "mean_quality"),
        }
        breakdown_rows.append(breakdown)

        scenario_lookup: dict[tuple[int, int], dict[str, object]] = {}
        for min_identity in range(90, 79, -1):
            for min_ref_cov in range(90, 79, -1):
                passing = [
                    row
                    for row in rows
                    if row["raw_status"] == "C"
                    and isinstance(row["identity"], float)
                    and isinstance(row["reference_coverage"], float)
                    and float(row["identity"]) >= min_identity
                    and float(row["reference_coverage"]) >= min_ref_cov
                ]
                classified_count = len(passing)
                rescued = classified_count - baseline_classified
                scenario = {
                    "dataset": label,
                    "min_identity": min_identity,
                    "min_ref_coverage": min_ref_cov,
                    "classified_reads": classified_count,
                    "newly_rescued_reads": rescued,
                    "remaining_unclassified": total - classified_count,
                    "classified_pct": 100 * classified_count / total,
                    "rescued_pct_of_baseline_unclassified": 100
                    * rescued
                    / baseline_unclassified,
                }
                grid_rows.append(scenario)
                scenario_lookup[(min_identity, min_ref_cov)] = scenario

        for scenario_name, cutoffs in selected.items():
            scenario = dict(scenario_lookup[cutoffs])
            scenario["scenario"] = scenario_name
            selected_rows.append(scenario)
            if cutoffs == (90, 90):
                continue
            rescued_rows = [
                row
                for row in mapped_c0
                if float(row["identity"]) >= cutoffs[0]
                and float(row["reference_coverage"]) >= cutoffs[1]
            ]
            taxa = Counter(
                (
                    int(row["inferred_taxid"]),
                    str(row["superkingdom"]),
                    str(row["kingdom"]),
                    str(row["phylum"]),
                    str(row["class"]),
                    str(row["order"]),
                    str(row["family"]),
                    str(row["genus"]),
                    str(row["species"]),
                )
                for row in rescued_rows
            )
            for taxonomy, count in sorted(taxa.items(), key=lambda item: (-item[1], item[0])):
                taxa_rows.append(
                    {
                        "dataset": label,
                        "scenario": scenario_name,
                        "min_identity": cutoffs[0],
                        "min_ref_coverage": cutoffs[1],
                        "inferred_taxid": taxonomy[0],
                        **dict(zip(RANKS, taxonomy[1:])),
                        "newly_rescued_reads": count,
                    }
                )

        report_sections.append(
            f"| {label} | {total:,} | {baseline_classified:,} "
            f"({100 * baseline_classified / total:.2f}%) | {baseline_unclassified:,} | "
            f"{modes['raw_unmapped']:,} | {len(mapped_c0):,} | "
            f"{modes['identity_only']:,} | {modes['reference_coverage_only']:,} | "
            f"{modes['identity_and_reference_coverage']:,} |"
        )

    breakdown_fields = list(breakdown_rows[0])
    grid_fields = list(grid_rows[0])
    selected_fields = ["dataset", "scenario"] + [
        field for field in grid_fields if field != "dataset"
    ]
    diagnostic_fields = [
        "dataset",
        "read_id",
        "raw_status",
        "baseline_taxid",
        "best_reference",
        "identity",
        "reference_coverage",
        "query_coverage",
        "read_length",
        "mean_quality",
        "baseline_failure_mode",
        "inferred_taxid",
        *RANKS,
        "original_lineage",
    ]
    taxa_fields = [
        "dataset",
        "scenario",
        "min_identity",
        "min_ref_coverage",
        "inferred_taxid",
        *RANKS,
        "newly_rescued_reads",
    ]

    write_tsv(OUT / "baseline_failure_breakdown.tsv", breakdown_rows, breakdown_fields)
    write_tsv(OUT / "threshold_grid_80_to_90.tsv", grid_rows, grid_fields)
    write_tsv(OUT / "selected_threshold_scenarios.tsv", selected_rows, selected_fields)
    write_tsv(OUT / "selected_scenario_rescued_taxa.tsv", taxa_rows, taxa_fields)
    with gzip.open(OUT / "unclassified_read_diagnostics.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=diagnostic_fields)
        writer.writeheader()
        for row in diagnostics:
            writer.writerow({key: fmt(row.get(key, "")) for key in diagnostic_fields})

    (OUT / "provenance.json").write_text(
        json.dumps(
            {
                "method": "Read-ID join of wf-16s assignments, bamstats, and alignment tables",
                "interpretation": (
                    "Sensitivity reconstruction using the recorded best alignment; "
                    "not a native wf-16s rerun"
                ),
                "threshold_grid": {
                    "min_percent_identity": "integer values 80 through 90 inclusive",
                    "min_ref_coverage": "integer values 80 through 90 inclusive",
                    "comparison": "greater than or equal to cutoff",
                },
                "runs": provenance,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    selected_for_report = {
        "identity_85",
        "refcov_85",
        "balanced_89",
        "balanced_88",
        "balanced_85",
    }
    scenario_report_rows = []
    for row in selected_rows:
        if row["scenario"] not in selected_for_report:
            continue
        scenario_report_rows.append(
            f"| {row['dataset']} | {row['scenario']} | "
            f"{row['min_identity']} | {row['min_ref_coverage']} | "
            f"{int(row['newly_rescued_reads']):,} | "
            f"{float(row['classified_pct']):.2f}% |"
        )

    cause_report_rows = []
    for row in breakdown_rows:
        dominant = (
            "identity"
            if int(row["identity_implicated"]) > int(row["reference_coverage_implicated"])
            else "reference coverage"
        )
        cause_report_rows.append(
            f"- **{row['dataset']}**: {float(row['mapped_taxid0_pct_of_unclassified']):.2f}% "
            f"of unclassified reads were mapped but threshold-rejected; {dominant} was "
            f"the more frequent implicated threshold."
        )

    report = """# BAER, BANAE, BLEA, and BGRN threshold sensitivity

This is a read-level sensitivity reconstruction from the saved native wf-16s
outputs. It does not alter or replace those outputs and does not rerun minimap2.
The coverage threshold is **reference coverage**, matching `min_ref_coverage`.

## Baseline 90/90 read accounting

| Dataset | Total | Classified | Unclassified | Raw unmapped | Mapped TaxID 0 | Identity only | Ref. coverage only | Both |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
""" + "\n".join(report_sections) + """

Every baseline-classified read passed both recorded thresholds. Every mapped
TaxID-0 read failed identity, reference coverage, or both. All existing
classified TaxIDs were reproduced exactly by joining the bamstats best reference
to the alignment table.

## Main cause

""" + "\n".join(cause_report_rows) + """

Read length alone does not explain the rejected reads; median classified and
mapped-TaxID-0 lengths are reported in `baseline_failure_breakdown.tsv`.

## Selected simulations

The scenario name and both numeric columns make the cutoff order explicit:
identity first, reference coverage second.

| Dataset | Scenario | Min identity | Min ref. coverage | Newly rescued | Classified total |
|---|---|---:|---:|---:|---:|
""" + "\n".join(scenario_report_rows) + """

## Outputs

- `baseline_failure_breakdown.tsv`: baseline cause counts and QC medians.
- `threshold_grid_80_to_90.tsv`: complete 121-scenario grid per dataset.
- `selected_threshold_scenarios.tsv`: compact set of one-axis and balanced scenarios.
- `unclassified_read_diagnostics.tsv.gz`: every baseline-unclassified read with
  alignment metrics, failure mode, and inferred best-reference taxonomy.
- `selected_scenario_rescued_taxa.tsv`: inferred taxa introduced by selected scenarios.
- `provenance.json`: parameters, input paths, and SHA-256 identities.

The inferred TaxIDs are suitable for cutoff selection and sensitivity analysis.
Use a native wf-16s rerun for canonical workflow reports and final abundance
artifacts after selecting a threshold.
"""
    (OUT / "README.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
