"""Smoke driver: parse the ATAC-sm Snakefile under both alignment backends.

This is the workflow-level parse/lint smoke test. It uses the public
Snakemake 9.23.1 API with a ``conda_inject`` stub (the project workflow has
no conda installed locally) and asserts that:

- ``backend=legacy`` produces the legacy DAG and excludes
  ``align_tachyon_upstream``.
- ``backend=tachyon_upstream`` produces the Tachyon producer rule and excludes
  ``fastp`` and ``prealign_reads``.
- ``backend=tachyon_upstream`` + ``prealign.enabled=true`` aborts at parse
  time (1-shot downstream guard in ``processing.smk``).

The parse verifies producer selection without submitting jobs.
"""
from __future__ import annotations

import os
import sys
import types
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SNAKEFILE = REPO / "workflow" / "Snakefile"
TEST_CFG = REPO / "config" / "config_test.yaml"


def _build_api_config(backend: str, prealign: bool) -> dict:
    return {
        "alignment": {
            "backend": backend,
            "prealign": {
                "enabled": prealign,
                "indices": [{"name": "chrM", "path": "/tmp/chrM"}],
            },
        }
    }


def _parse(backend: str, prealign: bool):
    """Return ``(exit_code, rules_set, message)`` parsing the workflow."""
    stub = types.ModuleType("conda_inject")
    stub.inject_env = stub.inject_env_file = lambda *a, **kw: None
    stub._get_envs = lambda *a, **kw: {}

    class _N:
        def __init__(self, *a, **kw):
            pass

        def __call__(self, *a, **kw):
            return self

        def __getattr__(self, _):
            return _N()

    stub.Environment = stub.PackageManager = _N
    sys.modules["conda_inject"] = stub

    from snakemake.api import SnakemakeApi
    from snakemake.settings.types import (
        ConfigSettings,
        DeploymentSettings,
        ResourceSettings,
        StorageSettings,
        WorkflowSettings,
    )

    with SnakemakeApi() as api:
        wfapi = api.workflow(
            snakefile=str(SNAKEFILE),
            workdir=REPO,
            config_settings=ConfigSettings(
                configfiles=[TEST_CFG],
                config=_build_api_config(backend, prealign),
            ),
            resource_settings=ResourceSettings(),
            storage_settings=StorageSettings(),
            workflow_settings=WorkflowSettings(),
            deployment_settings=DeploymentSettings(deployment_method=frozenset()),
        )
        try:
            wfapi.lint()
        except SystemExit as e:
            return e.code or 1, set(), f"SystemExit({e.code})"
        except Exception as e:
            return 2, set(), f"{type(e).__name__}: {e}"
        rules = {r.name for r in wfapi._workflow_store.rules}
    return 0, rules, ""


class SmokeParseTest(unittest.TestCase):
    def test_legacy_backend(self):
        code, rules, msg = _parse("legacy", prealign=True)
        self.assertEqual(code, 0, msg=msg)
        self.assertIn("align_bwa_mem", rules)
        self.assertIn("prealign_reads", rules)
        self.assertIn("fastp", rules)
        self.assertNotIn("align_tachyon_upstream", rules)

    def test_tachyon_backend(self):
        code, rules, msg = _parse("tachyon_upstream", prealign=False)
        self.assertEqual(code, 0, msg=msg)
        self.assertIn("align_tachyon_upstream", rules)
        self.assertIn("bwa_mem2_index", rules)
        self.assertNotIn("fastp", rules)
        self.assertNotIn("prealign_reads", rules)

    def test_tachyon_plus_prealign_aborts(self):
        # Parse-time guard inside processing.smk forces SystemExit(1) when both
        # backend=tachyon_upstream and prealign.enabled=true are active.
        code, _rules, _msg = _parse("tachyon_upstream", prealign=True)
        self.assertEqual(code, 1)


if __name__ == "__main__":
    unittest.main()
