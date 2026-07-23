# ATAC reference and BAM retention plan

## Goal

Use the same `hg38.fa` sequence as CUT&Tag with an isolated BWA-MEM2 index,
retain all mapped BAM records, keep downstream file paths and rules unchanged,
and validate the change on real data.

## Acceptance criteria

- BWA-MEM2 index lives under `hg38/indices_for_BWA-MEM2/` and does not touch
  the Bowtie2 index directory.
- Alignment uses the CUT&Tag FASTA at `{database_dir}/hg38/hg38.fa`.
- The retained BAM applies only `samtools view -F 4` after duplicate marking.
- Existing `{sample}.filtered.bam` paths and downstream rules remain unchanged.
- Local regression checks and a Snakemake dry-run pass.
- `B56-1_ATAC19` completes through BAM, tagAlign, peaks, TSS, ATAQV, peak QC,
  and bigWig on HPC.
- CD99 (`chrX:2690000-2750000`) is no longer a hard zero.
- A strict BAM derived from the retained BAM reproduces direct strict filtering
  on the same new alignment stream.

## Phases

1. **Plan and inspect** — complete
2. **Regression check (red)** — complete
   - Assert isolated index prefix and BWA-MEM2 index suffixes.
   - Assert filtering is exactly `-F 4`.
   - Assert downstream BAM paths are unchanged.
3. **Minimal implementation (green)** — complete
   - Update production and test configs.
   - Update BWA index dependency/path handling.
   - Update index rule to use `bwa-mem2 index -p`.
   - Relax BAM filtering to `-F 4`.
4. **Local verification** — complete
   - Run regression check.
   - Validate configuration/schema and Snakemake DAG.
5. **HPC real-data validation** — complete
   - Build index with Slurm.
   - Run `B56-1_ATAC19` to selected downstream targets in a new output root.
   - Check BAM integrity, flags, CD99 reads, and downstream artifacts.
6. **Controlled metric comparison** — complete
   - Compare old versus retained BAM descriptively.
   - Compare direct strict filtering versus strict filtering derived from the
     retained BAM on the same new alignment.
   - Record results and remaining differences.

## Files expected to change

- `config/config.yaml`
- `config/config_test.yaml`
- `workflow/rules/common.smk`
- `workflow/rules/processing.smk`
- one minimal regression check under `tests/`

## Errors encountered

| Error | Resolution |
|---|---|
| Root-level `Snakefile` not found | Workflow entry point is `workflow/Snakefile`. |
| `python -m unittest tests/test_alignment_contract.py` could not import a non-package path | Run the test file directly. |
| Default workflow profile requested an unavailable local Slurm executor plugin | Disabled the workflow profile for local parsing. |
| Full local dry-run could not see HPC-only FASTQs | Used local rule/DAG parsing; full dry-run will run on HPC. |
| First HPC dry-run could not resolve the optional project-specific HOMER installation needed by `peak_annotation` | Excluded HOMER annotation from this focused validation and will calculate peak count/FRIP directly from the unchanged BAM and MACS2 outputs. |
| Slurm rejected `amd-ep2` from login05 | Switched the validation job to available `amd-ep5,intel-sc3` partitions. |
| Job 2845615 failed before computation because compute nodes cannot reach conda channels | Reused the existing HPC `py311` tool environment, supplied the pinned samblaster binary locally, aliased installed MACS3 as `macs2`, and disabled per-rule conda creation for the validation run. ATAQV is excluded because it is not installed in the shared environment. |
| Locally compiled samblaster required a newer glibc and some moved-environment Python entry points had stale shebangs | Compile samblaster 0.1.24 on the HPC login node and use small wrappers that invoke MACS3/deepTools with the current environment Python. |
| BWA-MEM2 index construction was killed at 64 GB (exit 137) | Added a 120 GB memory resource to the index rule and raised the validation allocation to 128 GB. |
| TSS coverage aborted because the retained BAM has reads on random contigs absent from the TSS BED | Stream only contigs present in the slopped TSS BED from the indexed BAM, retaining memory-efficient `bedtools -sorted`. |
| First direct strict-remap control used 12 BWA threads while pipeline alignment uses 6 | Rerun the control with the exact pipeline thread count and process-substitution FASTQ input before comparing records. |
| Standalone validation copy lacked its project HOMER installation | The actual HPC checkout at `/storage/zhangkaiLab/hanlitian/macrophage/script/ATAC-sm` includes HOMER and its complete `rule all --dry-run` succeeds. The all-sample formal BAM dry-run in the validation copy also resolves 59 jobs: 22 fastp, 18 prealign, 18 align, and one BWA index. |

---

# Optional Tachyon-Upstream Backend Plan

## Goal

Add an opt-in `alignment.backend: tachyon_upstream` that replaces only raw
paired FASTQ preprocessing and alignment with the pinned
`/home/gilberthan/disk1/projects/tachyon_upstream` executable. Keep `legacy`
as the default and preserve the canonical filtered BAM/BAI contract for every
existing downstream rule.

## Implementation phases

1. Add static contract tests for backend selection, canonical outputs, ordered
   paired raw inputs, and explicit rejection of Tachyon with chrM prealignment. — complete
2. Add the smallest producer branch: raw FASTQs -> Tachyon -> canonical BAM,
   BAI, native UCF side output, and stats JSON. Do not change downstream QC. — complete
3. Run local legacy/tachyon dry-run, rulegraph, summary, lint, and fixture
   producer smoke checks. — local parse and fixture smoke complete; full local
   DAG remains blocked by HPC-only reference/sample metadata.
4. Completed pilot on HPC: compare isolated legacy-no-prealign and Tachyon
   output roots. Expanded validation remains required before promotion.

## Acceptance gates

- Tachyon DAG schedules no fastp, prealign, legacy alignment, or samblaster.
- Both branches produce a mapped-only, coordinate-sorted, indexed canonical BAM.
- Existing peak, track, TSS, FRiP, ATAQV, quantification, and reproducibility
  rules consume unchanged canonical paths.
- Real-data validation compares ATAC against legacy with chrM prealignment
  disabled; chrM prealignment remains a separate descriptive policy arm.
- Promotion stays opt-in until the HPC pilot and expanded cohort meet the
  documented mapping, fragment, peak, track, and QC agreement gates.
