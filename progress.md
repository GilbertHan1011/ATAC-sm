# Progress

## 2026-07-22

- Diagnosed CD99 hard-zero signal as the interaction of a different ATAC BWA
  reference and MAPQ/flag/whitelist filtering.
- Confirmed BWA-MEM2 supports an independent `-p` index prefix.
- Wrote the implementation and real-data validation plan.
- Added `tests/test_alignment_contract.py` and observed three expected failures.
- Updated both configs to use the isolated BWA-MEM2 prefix.
- Relaxed canonical BAM filtering to `-F 4` and removed the whitelist input dependency.
- Updated the BWA-MEM2 index rule to build the isolated prefix with `-p`.
- Regression check passes: 3 tests.
- `git diff --check` passes.
- Snakemake 9.23.1 parses the workflow and lists all rules, including
  `bwa_mem2_index` and `align_bwa_mem`.
- Full local dry-run reaches DAG construction but stops at expected HPC-only
  FASTQ paths; full data-aware dry-run is deferred to HPC.
- Deployed an isolated workflow copy at
  `/storage/zhangkaiLab/hanlitian/software/ATAC-sm-validation-20260722` and added
  a one-sample `B56-1_ATAC19` validation configuration.
- HPC dry-run resolved the expected 12-rule DAG. The first real job (2845615)
  stopped before computation because compute nodes cannot access conda channels;
  validation was adjusted to reuse installed HPC tools without online installs.
- Job 2846463 completed fastp but BWA-MEM2 index construction exceeded 64 GB;
  the index rule now explicitly requests 120 GB and the validation job 128 GB.
- The isolated 16 GB BWA-MEM2 index built successfully. Job 2851854 completed
  BAM, bigWig, and tagAlign but exposed a TSS `-sorted` assumption incompatible
  with retained random-contig reads. The unsorted fallback reached 128 GB, so
  TSS coverage now streams only TSS contigs and retains the sorted sweep.
- Job 2852843 completed successfully after the TSS fix: BAM, bigWig, tagAlign,
  MACS peaks, and TSS histogram exist. Retained BAM contains 99,908,555 reads
  and CD99 has 9,714 reads. A strict-filter control rerun remains pending with
  exactly the pipeline's BWA thread/input settings.
- Job 2852861 confirmed strict-filter deferral exactly: the retained-BAM
  strict derivative and direct strict remap each have 80,882,757 records and
  the same SAM-record MD5. Results are in `tests/hpc_validation/RESULTS.md`.
- Formal HPC `rule all` dry-run parsed the 22-row archived sample annotation
  but stopped at the missing project HOMER executable. An explicit formal BAM
  dry-run resolved the changed alignment path for all inputs: 22 fastp, 18
  prealign, 18 align, and one BWA index jobs.
- Corrected the HPC path: the real checkout is
  `/storage/zhangkaiLab/hanlitian/macrophage/script/ATAC-sm`, which has HOMER.
  Its complete `rule all --dry-run` succeeds (12 pending downstream jobs for
  its current ATAC26 configuration); no computation was submitted.

## 2026-07-23: Tachyon-Upstream Integration

- Recorded an opt-in backend plan. The integration boundary is raw paired FASTQ
  to the existing canonical filtered BAM/BAI; all downstream rules remain shared.
- Local implementation and dry-run/fixture checks precede HPC execution. Full
  real-data validation remains blocked until a frozen raw FASTQ manifest and
  pinned HPC executable/reference hashes are available.
- Schema gained `alignment.backend ∈ {legacy, tachyon_upstream}` with default
  `legacy`. `config/config_test.yaml` opts into Tachyon mode for the unit
  tests; default `config/config.yaml` keeps legacy.
- `processing.smk` added a parse-time guard that aborts on
  `backend=tachyon_upstream` + `prealign.enabled=true`, plus a new
  `rule align_tachyon_upstream` that reads paired, PEP-expanded raw FASTQs,
  shells out to the configured Tachyon executable with `--assay atac`, the
  configured adapter, repeated `--r1`/`--r2`, and BWA-MEM2 `--index`. It
  writes canonical `important_processed/bam/{sample}.filtered.bam` (+`.bai`),
  `middle_files/ucf/{sample}.ucf`, and `report/align/{sample}.tachyon.stats.json`.
  The legacy `align_bwa_mem`, `align_bowtie2`, and `prealign_reads` block is
  gated by `ALIGNMENT_BACKEND == "legacy"`.
- `qc.smk` `rule fastp:` is gated to legacy. `report.smk` `rule multiqc:`
  stops expanding fastp_html/fastp_json/aligner_log/samblaster_log/flagstat_log
  when the backend is Tachyon, so MultiQC keeps a clean input list.
- Regression suite is green: 21 tests cover legacy, Tachyon, and invalid
  Tachyon-plus-prealign workflow parsing. The apparent `bam_to_bed` parser
  error exposed a malformed new legacy conditional; it was fixed while keeping
  the rule's existing command semantics.
- Real local CLI dry-runs parse both configurations: Tachyon reaches its BWA
  index prerequisite and stops only because the local hg38 FASTA is HPC-only;
  legacy stops only because its archived `batch_18-19.csv` is absent locally.
- Rebuilt Tachyon at source commit `d30dc24`; `tachyon-upstream 0.1.0` SHA256
  is `19a2011d0ecd1a38ccbb172692eece50c99ad85222cfbb4ca2386581fa8d3046`.
  The ATAC fixture command generated an indexed mapped-only BAM (4,124 records,
  zero unmapped), UCF (4,112 rows), and stats exactly matching the prior native
   fixture baseline.
- Isolated HPC pilot completed for `B56-1_ATAC19`: legacy without chrM
  prealignment and Tachyon both produced valid indexed canonical BAMs and all
  requested downstream targets. BAM records differ by 0.43%, fragment Spearman
  is 0.99974, peak Jaccard 0.9532, BigWig Pearson 0.99989, and TSS Pearson
  0.9999997. Duplicate rates differ by 4.34 pp, so the backend remains opt-in
  pending expanded-cohort confirmation. Details: `tests/hpc_validation/TACHYON_RESULTS.md`.
