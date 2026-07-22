# B56-1_ATAC19 real-data validation

Validation output: `/storage/zhangkaiLab/hanlitian/macrophage/process/20260722_ATAC_reference_validation`.

## Verified

- Built BWA-MEM2 2.2.1 index from the CUT&Tag FASTA at
  `hg38/indices_for_BWA-MEM2/hg38` (16 GB); Bowtie2 index was untouched.
- Canonical BAM contains 99,908,555 mapped records and retains 98,405
  secondary alignments and 1,070,670 non-proper pairs.
- CD99 `chrX:2690000-2750000` has 9,714 retained records. The legacy strict
  BAM has 0 records in this interval.
- New BAM, BAI, tagAlign, MACS peak calls, bigWig, and TSS histogram were
  created successfully. Peak count is 300,486 and FRiP is 0.56863366.
- The TSS histogram has 4,002 positions and finite enrichment values
  (maximum 4.55125).
- Strict filtering applied after retention produced 80,882,757 records.
  A direct strict remap with the same BWA thread count, FASTQ streaming,
  samblaster, and filters produced exactly the same record count and SAM MD5:
  `a5cf785e017d59c1ee734404b0f3a81e`.

## Interpretation

The retained BAM preserves information needed for later strict filtering; the
strict derivative is bit-for-bit equivalent at SAM-record level to direct
strict filtering from the same new reference alignment. The legacy strict BAM
has 80,777,235 records; its small count difference from the new strict BAM is
expected because the reference sequence/index changed.

## Scope note

ATAQV was not run in this HPC validation because the compute nodes lack the
required cached environment and cannot reach conda channels. This does not
affect the BAM, tagAlign, MACS, bigWig, or TSS validation above.
