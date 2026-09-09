#!/usr/bin/env python3
"""Check text files changed by HEAD without platform-specific Git heuristics."""

from __future__ import annotations

import os
import pathlib
import subprocess
import sys


GIT = ["git", "-c", f"safe.directory={pathlib.Path.cwd().resolve().as_posix()}"]


def valid_commit(value: str | None) -> bool:
    if not value or value == "0" * 40:
        return False
    probe = subprocess.run(
        [*GIT, "rev-parse", "--verify", value],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    return probe.returncode == 0


def empty_tree() -> str:
    # The well-known Git empty-tree object is stable across repositories.
    return "4b825dc642cb6eb9a060e54bf8b69288fbee4904"


def choose_base() -> str:
    diff_base = os.environ.get("WF16S_DIFF_BASE")
    if diff_base and valid_commit(diff_base):
        return diff_base

    event = os.environ.get("GITHUB_EVENT_NAME", "")
    before = os.environ.get("GITHUB_EVENT_BEFORE")
    if event == "push" and valid_commit(before):
        return before  # The push range is exactly before..HEAD.

    if event == "pull_request":
        base_ref = os.environ.get("GITHUB_BASE_REF")
        if base_ref:
            remote_base = f"origin/{base_ref}"
            if valid_commit(remote_base):
                merge_base = subprocess.run(
                    [*GIT, "merge-base", "HEAD", remote_base],
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.DEVNULL,
                ).stdout.decode("ascii", errors="ignore").strip()
                if valid_commit(merge_base):
                    return merge_base

    for candidate in ("origin/main", "HEAD^"):
        if valid_commit(candidate):
            resolved = subprocess.run(
                [*GIT, "rev-parse", "--verify", candidate],
                check=True,
                stdout=subprocess.PIPE,
            ).stdout.decode("ascii").strip()
            if resolved != subprocess.run([*GIT, "rev-parse", "HEAD"], check=True,
                                           stdout=subprocess.PIPE).stdout.decode("ascii").strip():
                return resolved
    return empty_tree()


def changed_paths() -> list[pathlib.PurePosixPath]:
    base = choose_base()
    command = [*GIT, "diff", "--name-only", "-z", "--diff-filter=AM", base, "HEAD"]
    result = subprocess.run(command, check=True, stdout=subprocess.PIPE)
    return [pathlib.PurePosixPath(name.decode("utf-8")) for name in result.stdout.split(b"\0") if name]


def check_file(path: pathlib.PurePosixPath) -> list[str]:
    # Vendored Krona assets are pinned by byte hash in SOURCE.json.  Their
    # upstream JavaScript intentionally contains whitespace, so do not rewrite
    # it to satisfy this repository's text-style check.
    if path.as_posix().startswith("analysis/vendor/krona-2.8.1/"):
        return []
    data = subprocess.run(
        [*GIT, "show", f"HEAD:{path.as_posix()}"],
        check=True,
        stdout=subprocess.PIPE,
    ).stdout
    if b"\0" in data:
        return []

    errors: list[str] = []
    for line_number, line in enumerate(data.splitlines(keepends=True), start=1):
        body = line.rstrip(b"\r\n")
        if body.endswith((b" ", b"\t")):
            errors.append(f"{path}:{line_number}: trailing whitespace")
        indent = body[: len(body) - len(body.lstrip(b" \t"))]
        if b" " in indent and b"\t" in indent:
            errors.append(f"{path}:{line_number}: mixed spaces and tabs in indentation")

    normalized = data.replace(b"\r\n", b"\n")
    if normalized.endswith(b"\n\n"):
        errors.append(f"{path}: blank line at end of file")
    return errors


def main() -> int:
    errors = [error for path in changed_paths() for error in check_file(path)]
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("Committed whitespace check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
