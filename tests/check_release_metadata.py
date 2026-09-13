#!/usr/bin/env python3
"""Independently validate release metadata and optional annotated tag identity."""

from __future__ import annotations

import argparse
import datetime as dt
import re
import subprocess
from pathlib import Path


VERSION_RE = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+$")
DATE_RE = re.compile(r"^[0-9]{4}-[0-9]{2}-[0-9]{2}$")


def fail(message: str) -> None:
    raise SystemExit(f"release metadata check failed: {message}")


def git(repo_root: Path, *args: str) -> str:
    command = ["git", "-c", f"safe.directory={repo_root}", "-C", str(repo_root), *args]
    try:
        return subprocess.check_output(command, text=True, stderr=subprocess.STDOUT).strip()
    except subprocess.CalledProcessError as exc:
        fail(f"git {' '.join(args)} failed: {exc.output.strip()}")
    raise AssertionError("unreachable")


def read_strict_version(path: Path) -> str:
    try:
        raw = path.read_bytes()
    except OSError as exc:
        fail(f"could not read VERSION: {exc}")
    if not raw.endswith(b"\n") or raw.count(b"\n") != 1:
        fail("VERSION must contain exactly one newline-terminated value")
    try:
        value = raw[:-1].decode("utf-8")
    except UnicodeDecodeError as exc:
        fail(f"VERSION is not UTF-8: {exc}")
    if not VERSION_RE.fullmatch(value):
        fail("VERSION must contain one SemVer value")
    return value


def read_cff_fields(path: Path) -> tuple[str, str]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        fail(f"could not read CITATION.cff: {exc}")
    version = next((re.sub(r"^version:\s*", "", line).strip() for line in lines
                    if re.match(r"^version:\s*", line)), None)
    date_value = next((re.sub(r"^date-released:\s*", "", line).strip().strip('"\'')
                       for line in lines if re.match(r"^date-released:\s*", line)), None)
    if version is None or date_value is None or not VERSION_RE.fullmatch(version) or not DATE_RE.fullmatch(date_value):
        fail("CITATION.cff has invalid top-level version/date-released fields")
    return version, date_value


def read_first_release_heading(path: Path) -> tuple[str, str]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        fail(f"could not read CHANGELOG.md: {exc}")
    heading = re.compile(r"^## \[([^]]+)\] - ([0-9]{4}-[0-9]{2}-[0-9]{2})$")
    for line in lines:
        match = heading.fullmatch(line)
        if match and match.group(1) != "Unreleased":
            version, date_value = match.groups()
            if not VERSION_RE.fullmatch(version) or not DATE_RE.fullmatch(date_value):
                fail("first non-Unreleased CHANGELOG heading is invalid")
            return version, date_value
    fail("CHANGELOG has no dated non-Unreleased release heading")
    raise AssertionError("unreachable")


def assert_readme_anchor(path: Path, version: str) -> None:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        fail(f"could not read README.md: {exc}")
    expected = f"Current pipeline version: **{version}**."
    if sum(line == expected for line in lines) != 1:
        fail("README does not contain exactly one current-version anchor")


def parse_tagger_date(value: str) -> dt.date:
    try:
        parsed = dt.datetime.fromisoformat(value)
    except ValueError as exc:
        fail(f"annotated tagger timestamp is invalid: {exc}")
    if parsed.tzinfo is None:
        fail("annotated tagger timestamp lacks a timezone")
    return parsed.astimezone(dt.timezone.utc).date()


def check(repo_root: Path, required_tag: str | None = None) -> str:
    version = read_strict_version(repo_root / "VERSION")
    cff_version, release_date = read_cff_fields(repo_root / "CITATION.cff")
    changelog_version, changelog_date = read_first_release_heading(repo_root / "CHANGELOG.md")
    if not (version == cff_version == changelog_version):
        fail(f"version mismatch: VERSION={version}, CFF={cff_version}, CHANGELOG={changelog_version}")
    if release_date != changelog_date:
        fail(f"release date mismatch: CFF={release_date}, CHANGELOG={changelog_date}")
    assert_readme_anchor(repo_root / "README.md", version)

    if required_tag is not None:
        expected_tag = f"v{version}"
        if required_tag != expected_tag:
            fail(f"--require-tag must equal {expected_tag}")
        tag_ref = f"refs/tags/{expected_tag}"
        if git(repo_root, "cat-file", "-t", tag_ref) != "tag":
            fail("required tag is not an annotated tag object")
        tagged_commit = git(repo_root, "rev-parse", f"{tag_ref}^{{commit}}")
        head = git(repo_root, "rev-parse", "HEAD")
        if tagged_commit != head:
            fail(f"tag target {tagged_commit} does not equal HEAD {head}")
        if git(repo_root, "status", "--porcelain", "--untracked-files=all"):
            fail("exact-tag mode requires a clean checkout")
        tagger_value = git(repo_root, "for-each-ref", "--format=%(taggerdate:iso-strict)", tag_ref)
        if parse_tagger_date(tagger_value).isoformat() != release_date:
            fail("annotated tagger UTC date does not equal release metadata date")
    return version


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-tag", default=None)
    args = parser.parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    version = check(repo_root, args.require_tag)
    print(f"Release metadata check passed for v{version}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
