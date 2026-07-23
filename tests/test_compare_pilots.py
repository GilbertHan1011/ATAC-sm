from pathlib import Path


def test_fragment_histogram_uses_proper_first_mates():
    source = (Path(__file__).parents[1] / "tests/hpc_validation/compare_pilots.py").read_text()
    assert '[SAMTOOLS, "view", "-f", "66", str(bam)]' in source


def test_fragment_histogram_limits_qc_range_to_one_kilobase():
    source = (Path(__file__).parents[1] / "tests/hpc_validation/compare_pilots.py").read_text()
    assert "if 0 < length <= 1000:" in source


def test_bam_metrics_checks_bam_and_index():
    source = (Path(__file__).parents[1] / "tests/hpc_validation/compare_pilots.py").read_text()
    assert '"quickcheck", "-v", str(bam)' in source
    assert 'str(bam), "chr1:1-100000"' in source
