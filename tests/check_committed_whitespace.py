#!/usr/bin/env python3
"""Check text files changed by HEAD without platform-specific Git heuristics."""

from __future__ import annotations

import pathlib
import subprocess
import sys


def changed_paths() -> list[pathlib.Path]:
    git = ["git", "-c", f"safe.directory={pathlib.Path.cwd().resolve().as_posix()}"]
    parent = subprocess.run(
        [*git, "rev-parse", "--verify", "HEAD^"],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    command = [*git, "diff-tree", "--no-commit-id", "--name-only", "-z", "-r"]
    command.extend(["HEAD^", "HEAD"] if parent.returncode == 0 else ["--root", "HEAD"])
    result = subprocess.run(command, check=True, stdout=subprocess.PIPE)
    return [pathlib.Path(name.decode("utf-8")) for name in result.stdout.split(b"\0") if name]


def check_file(path: pathlib.Path) -> list[str]:
    if not path.is_file():
        return []
    data = path.read_bytes()
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
