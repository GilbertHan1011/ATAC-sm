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
