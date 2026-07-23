# Findings

- CUT&Tag uses `{database_dir}/hg38/hg38.fa`; its Bowtie2 index is isolated at
  `{database_dir}/hg38/indices_for_Bowtie2/hg38`.
- The current ATAC config points BWA-MEM2 at a different reference under
  `{shared_genome_dir}/bwa-mem2/GRCh38`.
- BWA-MEM2 supports `bwa-mem2 index -p <prefix> <in.fasta>`, so the same FASTA
  can feed an isolated index without copying or modifying it.
- BWA-MEM2 2.3 index outputs used by the existing rule are `.amb`, `.ann`,
  `.pac`, `.bwt.2bit.64`, and `.0123`.
- Current filtering is `-q 30 -F 2316 -L whitelist` plus `-f 2` for paired
  samples. The requested retained BAM filter is only `-F 4`.
- Downstream rules consistently consume
  `important_processed/bam/{sample}.filtered.bam`; retaining this filename
  avoids broad unrelated edits.
- Existing downstream metrics include mitochondrial fraction, peak count,
  FRiP, regulatory fraction, TSS coverage, and ATAQV metrics.
- Exact equality with the old production run is not expected because both the
  reference sequence and filter change. Exact reproducibility will instead be
  tested between two strict-filter paths using the same new alignment.

## Tachyon-Upstream Backend

- The source CLI at `/home/gilberthan/disk1/projects/tachyon_upstream` accepts
  paired raw `--r1`/`--r2`, `--assay atac`, BWA-MEM2 index, BAM, required UCF,
  and optional stats JSON. It writes coordinate-sorted, indexed mapped-only BAM.
- It cannot emit fastp reports or apply the existing chrM prealignment policy;
  `alignment.backend: tachyon_upstream` therefore requires prealignment off.
- The existing 2,500-pair fixture recovered every legacy no-prealign ATAC
  fragment, with 92.67% within 5 bp and chromosome Pearson 0.99961. BAM bytes
  and aligner-specific tags are intentionally not parity criteria.

## Verification status

- `pytest -q tests`: 21 passed. The suite covers backend config, raw paired
  FASTQ collection, producer selection, legacy preservation, MultiQC inputs,
  and real Snakemake parsing for legacy, Tachyon, and rejected prealignment.
- Local CLI dry-run parses both branches. Tachyon reaches `bwa_mem2_index` and
  stops only because `/storage/.../hg38/hg38.fa` is unavailable locally; legacy
  stops only because archived `config/batch_18-19.csv` is absent.
- Rebuilt the source-pinned Tachyon executable at `d30dc24`:
  `tachyon-upstream 0.1.0`, SHA256
  `19a2011d0ecd1a38ccbb172692eece50c99ad85222cfbb4ca2386581fa8d3046`.
  Its ATAC fixture run yields the exact prior native stats: 2,500 input pairs,
  4,124 mapped BAM records, zero unmapped records in BAM, 12 duplicate pairs,
  and 4,112 UCF rows.
