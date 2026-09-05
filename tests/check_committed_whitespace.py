#!/usr/bin/env python3
"""Check text files changed by HEAD without platform-specific Git heuristics."""

from __future__ import annotations

import pathlib
import subprocess
import sys


GIT = ["git", "-c", f"safe.directory={pathlib.Path.cwd().resolve().as_posix()}"]


def changed_paths() -> list[pathlib.PurePosixPath]:
    parent = subprocess.run(
        [*GIT, "rev-parse", "--verify", "HEAD^"],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    shallow = subprocess.run(
        [*GIT, "rev-parse", "--is-shallow-repository"],
        check=True,
        stdout=subprocess.PIPE,
    )
    if parent.returncode != 0 and shallow.stdout.strip() == b"true":
        raise RuntimeError(
            "HEAD parent is unavailable in a shallow checkout; "
            "configure checkout fetch-depth >= 2."
        )

    command = [
        *GIT,
        "diff-tree",
        "--no-commit-id",
        "--name-only",
        "--diff-filter=AM",
        "-z",
        "-r",
    ]
    command.extend(["HEAD^", "HEAD"] if parent.returncode == 0 else ["--root", "HEAD"])
    result = subprocess.run(command, check=True, stdout=subprocess.PIPE)
    return [pathlib.PurePosixPath(name.decode("utf-8")) for name in result.stdout.split(b"\0") if name]


def check_file(path: pathlib.PurePosixPath) -> list[str]:
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
