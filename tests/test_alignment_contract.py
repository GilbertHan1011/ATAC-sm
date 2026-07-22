import re
import unittest
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
INDEX = "{database_dir}/hg38/indices_for_BWA-MEM2/hg38"


class AlignmentContractTest(unittest.TestCase):
    def test_configs_use_isolated_bwa_mem2_index(self):
        for path in (ROOT / "config/config.yaml", ROOT / "config/config_test.yaml"):
            config = yaml.safe_load(path.read_text())
            self.assertEqual(config["alignment"]["bwa"]["index"], INDEX)
            self.assertEqual(config["refs"]["fasta"], "{database_dir}/hg38/hg38.fa")

    def test_workflow_retains_all_mapped_records(self):
        common = (ROOT / "workflow/rules/common.smk").read_text()
        body = re.search(
            r"def get_filtering_flags\(wildcards\):\n(?P<body>(?:    .*\n)+)", common
        ).group("body")
        self.assertIn('return "-F 4"', body)
        for strict_flag in ("-q 30", "2316", "-L", "-f 2"):
            self.assertNotIn(strict_flag, body)

    def test_index_rule_uses_expected_prefix_and_files(self):
        processing = (ROOT / "workflow/rules/processing.smk").read_text()
        self.assertIn('bwa-mem2 index -p {params.index_prefix} {input.fasta}', processing)
        self.assertRegex(
            processing,
            re.compile(r"rule bwa_mem2_index:.*?mem_mb\s*=\s*120000", re.S),
        )
        for suffix in (".amb", ".ann", ".pac", ".bwt.2bit.64", ".0123"):
            self.assertIn(suffix, processing)
        self.assertIn('"{sample}.filtered.bam"', processing)

    def test_tss_coverage_accepts_extra_reference_contigs(self):
        qc = (ROOT / "workflow/rules/qc.smk").read_text()
        self.assertIn("samtools view -bh {input.bam}", qc)
        self.assertIn("bedtools coverage -a \"$TSS_SLOPPED\" -b stdin -d -sorted", qc)


if __name__ == "__main__":
    unittest.main()
