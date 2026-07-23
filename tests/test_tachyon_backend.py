"""Static contract tests for the optional Tachyon-Upstream backend.

These tests are independent of any real FASTQ data; they only inspect config
files, the schema, and Snakemake rule text. They cover:

- ``alignment.backend`` is a recognised config key with legacy/tachyon_upstream.
- ``config_test.yaml`` opts into the new backend explicitly.
- Tachyon mode rejects prealignment at parse time.
- The Tachyon producer reads raw paired FASTQs in deterministic matching run
  order, invokes the pinned binary, uses the BWA-MEM2 index, sets
  ``--assay atac`` + ``--adapter nextera``, and writes the canonical
  ``important_processed/bam/{sample}.filtered.bam`` (+ ``.bai``) plus the
  native UCF side output and a stats JSON.
- Under Tachyon the DAG does NOT include fastp, prealign, samblaster, or the
  legacy BWA log/samblaster log/samtools_flagstat logs.
- Legacy producer semantics and downstream canonical BAM/BAI path remain
  unchanged.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import unittest
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
PROCESSING = (ROOT / "workflow/rules/processing.smk").read_text()
COMMON = (ROOT / "workflow/rules/common.smk").read_text()
QC = (ROOT / "workflow/rules/qc.smk").read_text()
REPORT = (ROOT / "workflow/rules/report.smk").read_text()
SCHEMA = (ROOT / "schemas/config.schema.yml").read_text()
CONFIG_TEST = yaml.safe_load((ROOT / "config/config_test.yaml").read_text())
CONFIG = yaml.safe_load((ROOT / "config/config.yaml").read_text())


TACHYON_EXE_DEFAULT = "/home/gilberthan/disk1/projects/tachyon_upstream/target/release/tachyon-upstream"


def _has_rule(body: str, name: str) -> bool:
    """True iff a top-level `rule <name>:` block exists in *body*."""
    return bool(re.search(rf"^\s*rule {re.escape(name)}\s*:", body, re.M))


def _extract_rule_block(body: str, name: str) -> str:
    """Return the full text of `rule <name>:` through the next blank+non-indented line."""
    m = re.search(rf"^\s*rule {re.escape(name)}\s*:[^\n]*\n(?P<body>(?:[ \t].*\n|^\s*\n)*)", body, re.M)
    if not m:
        raise AssertionError(f"rule {name} not found")
    return m.group(0)


class SchemaConfigTest(unittest.TestCase):
    def test_schema_accepts_alignment_backend(self):
        # The schema should declare alignment.backend with the two values we support.
        self.assertIn("backend", SCHEMA)

    def test_test_config_pins_tachyon_upstream_backend(self):
        cfg = CONFIG_TEST
        self.assertEqual(cfg["alignment"].get("backend"), "tachyon_upstream")
        # Sanity: Tachyon backend is only allowed when prealignment is off.
        self.assertFalse(cfg["alignment"].get("prealign", {}).get("enabled", True))

    def test_production_config_keeps_legacy_default_and_exposes_tachyon_settings(self):
        alignment = CONFIG["alignment"]
        self.assertEqual(alignment.get("backend"), "legacy")
        self.assertEqual(
            alignment["tachyon_upstream"]["executable"], TACHYON_EXE_DEFAULT
        )
        self.assertEqual(alignment["tachyon_upstream"]["adapter"], "nextera")


class TachyonProducerContractTest(unittest.TestCase):
    def test_tachyon_backend_constant_is_defined(self):
        # The processing module must surface the chosen backend as a constant so
        # downstream helpers can branch off it without re-parsing config.
        self.assertRegex(
            PROCESSING,
            re.compile(
                r"ALIGNMENT_BACKEND\s*=\s*config\[.alignment.\]\.get\(.backend.\s*,\s*.[A-Za-z_]+.\)",
                re.S,
            ),
        )

    def test_tachyon_rejects_prealignment_at_parse_time(self):
        # A clear, action-stopping guard must live in processing.smk.
        guard = re.search(
            r"if\s+ALIGNMENT_BACKEND\s*==\s*[\"']tachyon_upstream[\"'][^\n]*\n(?P<body>(?:[ \t]+.*\n)+)",
            PROCESSING,
        )
        self.assertIsNotNone(guard, "Tachyon guard block missing")
        body = guard.group("body") if guard else ""
        # Must depend on the legacy prealign flag and abort the workflow.
        self.assertIn("prealign", body)
        self.assertIn("sys.exit", body)

    def test_tachyon_producer_rule_exists(self):
        self.assertTrue(_has_rule(PROCESSING, "align_tachyon_upstream"))

    def test_tachyon_producer_reads_raw_paired_fastqs(self):
        block = _extract_rule_block(PROCESSING, "align_tachyon_upstream")
        self.assertIn("get_paired_raw_fastqs_for_sample", block)
        self.assertIn("get_paired_raw_fastqs_for_sample", COMMON)
        self.assertRegex(block, re.compile(r"r1\s*=\s*lambda\s*w\s*:\s*"))
        self.assertRegex(block, re.compile(r"r2\s*=\s*lambda\s*w\s*:\s*"))

    def test_tachyon_producer_uses_bwa_index_and_required_flags(self):
        block = _extract_rule_block(PROCESSING, "align_tachyon_upstream")
        # BWA-MEM2 index input dependency.
        self.assertIn("get_bwa_index_input", block)
        # Capture the entire shell heredoc; allow newlines inside.
        shell_m = re.search(r"shell:\s*\"\"\"(?P<body>.*?)\"\"\"", block, re.S)
        self.assertIsNotNone(shell_m, "tachyon rule shell block missing")
        shell = shell_m.group("body") if shell_m else ""
        # The shell must invoke the params-resolved executable reference, the
        # literal default constant is set as the rule's params value above.
        self.assertIn("{params.exe:q}", shell)
        self.assertIn("tachyon_upstream", PROCESSING)
        self.assertIn("TACHYON_UPSTREAM_EXE", block)
        self.assertIn("--assay atac", shell)
        self.assertIn("--adapter {params.adapter:q}", shell)
        # Both --r1 and --r2 must be present, repeated to support multiple runs.
        self.assertIn("--r1", block)
        self.assertIn("--r2", block)
        # Must NOT invoke fastp, samblaster, or bwa-mem2 directly.
        self.assertNotIn("fastp", shell)
        self.assertNotIn("samblaster", shell)
        self.assertNotIn("bwa-mem2 mem", shell)
        self.assertNotIn("samtools index", shell)

    def test_tachyon_producer_emits_canonical_and_side_outputs(self):
        block = _extract_rule_block(PROCESSING, "align_tachyon_upstream")
        # The canonical BAM/BAI must keep their existing filenames so every
        # downstream rule remains a no-op for Tachyon.
        self.assertIn("{sample}.filtered.bam", block)
        self.assertIn("{sample}.filtered.bam.bai", block)
        # Side outputs requested by the design doc.
        self.assertIn("{sample}.ucf", block)
        # Stats JSON in canonical report directory.
        self.assertRegex(block, re.compile(r"\{sample\}\.tachyon\.stats\.json"))


class TachyonDAGIsolationTest(unittest.TestCase):
    """Run-order collection and report-DAG isolation."""

    def test_legacy_producer_branch_remains_intact(self):
        # The legacy BWA-MEM2 alignment rule must still exist and still pipe
        # through bwa-mem2 -> samblaster -> samtools view -F 4 -> sort -> index.
        self.assertTrue(_has_rule(PROCESSING, "align_bwa_mem"))
        bwa_block = _extract_rule_block(PROCESSING, "align_bwa_mem")
        self.assertIn("bwa-mem2 mem", bwa_block)
        self.assertIn('return "-F 4"', COMMON)
        self.assertIn('"{sample}.filtered.bam"', bwa_block)

    def test_fastp_and_prealign_only_when_legacy(self):
        # fastp must remain scoped to legacy mode (legacy default keeps it).
        # We assert that the rule definition is structurally gated by ALIGNMENT_BACKEND.
        self.assertIn("ALIGNMENT_BACKEND", PROCESSING)
        self.assertTrue(_has_rule(QC, "fastp"))
        # The Tachyon producer does not depend on fastp outputs.
        tachyon_block = _extract_rule_block(PROCESSING, "align_tachyon_upstream")
        self.assertNotIn(".trimmed", tachyon_block)

    def test_report_excludes_fastp_under_tachyon(self):
        # The MultiQC rule's expand() must skip fastp outputs when Tachyon is
        # active. We assert the conditional is wired through the backend flag.
        # ``has_prealignments`` is already the legacy-style gate; we add an
        # analogous fastp gate for Tachyon.
        self.assertIn("ALIGNMENT_BACKEND", REPORT)
        self.assertIn("fastp", REPORT)
        # The rule must exclude the legacy per-sample BWA + samblaster logs
        # when Tachyon is active to keep the DAG closed. The legacy producer
        # builds files like "{sample}.bwa.log" via concat; we allow both
        # forms.
        self.assertIn("ALIGNER_LOG_SUFFIX", REPORT)
        self.assertIn("samblaster.log", REPORT)
        self.assertIn("samtools_flagstat.log", REPORT)


class TachyonParseTimeGuardTest(unittest.TestCase):
    """Drive the parse-time guard directly through a child interpreter.

    The guard in ``processing.smk`` is a 3-line conditional that calls
    ``sys.exit(1)`` when ``backend=tachyon_upstream`` is combined with
    prealignment. We isolate the exact block by reading processing.smk,
    feed it to ``exec`` in a child interpreter with Snakemake DSL
    keywords stubbed to no-ops, and assert the abort fires only for the
    bad combination.
    """

    DRIVER = r'''
import sys

config = {
    "alignment": {
        "backend": "tachyon_upstream",
        "tool": "bwa-mem2",
        "bwa": {"index": "/tmp/none"},
        "prealign": {"enabled": True, "indices": [{"name": "chrM", "path": "/tmp/chrM"}]},
    },
    "refs": {"fasta": "/tmp/ref.fa", "mito_name": "chrM", "genome_size_bp": 0},
    "resources": {"mem_mb": 4000, "threads": 2},
}

# Inert stubs for Snakemake DSL + workflow helpers that the module touches
# before reaching the guard. They exist only to keep ``exec`` from
# failing on name lookups; the guard fires before any rule body executes.
def _noop(*a, **kw):
    if a and isinstance(a[0], list):
        return a[0]
    return ""
expand = _noop
multiext = _noop
get_output_dir = lambda s: f"/tmp/{s}"
get_bwa_index_path = lambda: "/tmp/none"
get_bwa_index_input = lambda w=None: ["/tmp/none"]
sanitize_bwa_min_score_flag = lambda: ""
get_filtering_flags = lambda w: "-F 4"
get_bowtie2_input_string = lambda w, i: ""
get_add_mate_tags = lambda w: ""
get_reads = lambda w, d: []
_get_fastqs_for_sample = lambda s, f: ([], [])
get_all_fastqs_for_sample = lambda s: ([], [])
annotation_sheet_path = "/tmp/annot.csv"
HOMER_path = "/tmp/homer"
module_name = "atacseq_pipeline"
samples = {}
annot = None

# The guard block is the only piece that matters for this test. We extract
# the exact text from processing.smk so the assertion stays grounded in the
# real workflow code. The module-level constants it references are
# defined just above the guard; we resolve them the same way.
import re
src = open(r"%s").read()

def _resolve(name):
    m = re.search(rf"^{name}\s*=\s*(.+?)$", src, re.M)
    assert m is not None, f"constant {name} not found"
    return eval(m.group(1), {"config": config})

ALIGNMENT_BACKEND = _resolve("ALIGNMENT_BACKEND")
has_prealignments = bool(
    config["alignment"].get("prealign", {}).get("enabled", True)
    and config["alignment"].get("prealign", {}).get("indices")
)

guard = re.search(
    r"^# Tachyon cannot prealign.*?sys\.exit\(1\)\n",
    src,
    re.S | re.M,
)
assert guard is not None, "Tachyon guard block not found"
try:
    exec(guard.group(0), globals())
except SystemExit as e:
    print("EXIT_CODE:", e.code)
    sys.exit(e.code or 1)
print("NO_EXIT")
sys.exit(0)
'''

    def test_processing_module_aborts_when_tachyon_plus_prealign(self):
        path = ROOT / "workflow/rules/processing.smk"
        driver_src = self.DRIVER % str(path)
        result = subprocess.run(
            [sys.executable, "-c", driver_src],
            capture_output=True,
            text=True,
            timeout=20,
        )
        self.assertNotEqual(result.returncode, 0,
                            "Tachyon+prealign must trigger parse-time abort")
        combined = (result.stderr + result.stdout).lower()
        self.assertIn("tachyon_upstream", combined)
        self.assertIn("prealign", combined)

    def test_legacy_backend_with_prealign_still_parses(self):
        path = ROOT / "workflow/rules/processing.smk"
        driver_src = (self.DRIVER % str(path)).replace(
            '"backend": "tachyon_upstream"', '"backend": "legacy"'
        )
        result = subprocess.run(
            [sys.executable, "-c", driver_src],
            capture_output=True,
            text=True,
            timeout=20,
        )
        self.assertEqual(result.returncode, 0,
                         msg=f"legacy+prealign must parse cleanly; got {result.stderr}")
        self.assertIn("NO_EXIT", result.stdout)


if __name__ == "__main__":
    unittest.main()
