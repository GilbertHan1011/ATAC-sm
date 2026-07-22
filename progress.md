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
