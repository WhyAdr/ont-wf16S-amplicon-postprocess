"""Deterministic Git-range tests for ``check_committed_whitespace``."""

from __future__ import annotations

import contextlib
import importlib.util
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock


SCRIPT = Path(__file__).with_name("check_committed_whitespace.py")
SPEC = importlib.util.spec_from_file_location("check_committed_whitespace", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
CHECKER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECKER)


class ChooseBaseTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.repo = Path(self.tempdir.name)
        self.git("init", "-q")
        self.git("config", "user.name", "Whitespace Test")
        self.git("config", "user.email", "whitespace-test@example.invalid")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def git(self, *args: str) -> str:
        result = subprocess.run(
            ["git", *args], cwd=self.repo, check=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True,
        )
        return result.stdout.strip()

    def commit(self, message: str) -> str:
        index = len(list(self.repo.glob("commit-*.txt"))) + 1
        (self.repo / f"commit-{index}.txt").write_text(f"{message}\n", encoding="utf-8")
        self.git("add", ".")
        self.git("commit", "-qm", message)
        return self.git("rev-parse", "HEAD")

    @contextlib.contextmanager
    def checker_context(self, **env: str):
        values = {
            "WF16S_DIFF_BASE": "",
            "GITHUB_EVENT_NAME": "",
            "GITHUB_EVENT_BEFORE": "",
            "GITHUB_BASE_REF": "",
            **env,
        }
        git = ["git", "-c", f"safe.directory={self.repo.as_posix()}"]
        with (
            contextlib.chdir(self.repo),
            mock.patch.object(CHECKER, "GIT", git),
            mock.patch.dict(os.environ, values, clear=False),
        ):
            yield

    def test_root_commit_falls_back_to_the_empty_tree(self) -> None:
        self.commit("root")

        with self.checker_context():
            self.assertEqual(CHECKER.choose_base(), CHECKER.empty_tree())

    def test_push_prefers_explicit_base_then_uses_event_before(self) -> None:
        first = self.commit("first")
        second = self.commit("second")
        self.commit("third")

        with self.checker_context(
            WF16S_DIFF_BASE=first,
            GITHUB_EVENT_NAME="push",
            GITHUB_EVENT_BEFORE=second,
        ):
            self.assertEqual(CHECKER.choose_base(), first)

        with self.checker_context(GITHUB_EVENT_NAME="push", GITHUB_EVENT_BEFORE=second):
            self.assertEqual(CHECKER.choose_base(), second)

    def test_all_zero_new_branch_before_sha_falls_back_to_parent(self) -> None:
        parent = self.commit("parent")
        self.commit("new branch tip")

        with self.checker_context(
            GITHUB_EVENT_NAME="push", GITHUB_EVENT_BEFORE="0" * 40,
        ):
            self.assertEqual(CHECKER.choose_base(), parent)

    def test_pull_request_uses_merge_base_of_origin_base_branch(self) -> None:
        base = self.commit("main base")
        self.git("branch", "-M", "main")
        self.git("update-ref", "refs/remotes/origin/main", base)
        self.git("switch", "-qc", "feature")
        self.commit("feature change")

        with self.checker_context(GITHUB_EVENT_NAME="pull_request", GITHUB_BASE_REF="main"):
            self.assertEqual(CHECKER.choose_base(), base)


if __name__ == "__main__":
    unittest.main()
