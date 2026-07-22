#!/bin/bash
#SBATCH --job-name=atac-ref-val
#SBATCH --partition=amd-ep5,intel-sc3
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/storage/zhangkaiLab/hanlitian/software/ATAC-sm-validation-20260722/tests/hpc_validation/slurm-%j.out
#SBATCH --error=/storage/zhangkaiLab/hanlitian/software/ATAC-sm-validation-20260722/tests/hpc_validation/slurm-%j.err

set -euo pipefail

repo=/storage/zhangkaiLab/hanlitian/software/ATAC-sm-validation-20260722
out=/storage/zhangkaiLab/hanlitian/macrophage/process/20260722_ATAC_reference_validation
sample=B56-1_ATAC19
snake=/storage/zhangkaiLab/hanlitian/software/micromamba/micromamba/envs/py311/bin/snakemake
tools=/storage/zhangkaiLab/hanlitian/software/micromamba/micromamba/envs/py311/bin
export PATH="$repo/tests/hpc_validation/bin:$tools:$PATH"

targets=(
  "$out/important_processed/bam/$sample.filtered.bam"
  "$out/middle_files/prealigned/${sample}_prealigned_1.fq.gz"
  "$out/middle_files/prealigned/${sample}_prealigned_2.fq.gz"
  "$out/middle_files/bed/$sample.tagAlign.gz"
  "$out/important_processed/peaks/${sample}_peaks.narrowPeak"
  "$out/report/tss_coverage/$sample.tss_histogram.csv"
  "$out/important_processed/tracks/$sample.bw"
)

"$snake" \
  --directory "$repo" \
  --snakefile "$repo/workflow/Snakefile" \
  --workflow-profile none \
  --executor local \
  --cores "$SLURM_CPUS_PER_TASK" \
  --resources mem_mb=120000 \
  --notemp \
  --rerun-incomplete \
  "${targets[@]}"

bam="$out/important_processed/bam/$sample.filtered.bam"
old_bam=/storage/zhangkaiLab/hanlitian/macrophage/data/ATAC/bam/$sample.filtered.bam
peaks="$out/important_processed/peaks/${sample}_peaks.narrowPeak"
whitelist=/storage/zhangkaiLab/hanlitian/database/algin_ref/atac_pipeline/hg38/hg38-ENCODE_whitelist.sorted.v3.bed
index=/storage/zhangkaiLab/hanlitian/database/algin_ref/atac_pipeline/hg38/indices_for_BWA-MEM2/hg38
fq1="$out/middle_files/prealigned/${sample}_prealigned_1.fq.gz"
fq2="$out/middle_files/prealigned/${sample}_prealigned_2.fq.gz"
validation="$out/validation"
mkdir -p "$validation"
old_bam_local="$validation/old_strict.bam"
ln -sfn "$old_bam" "$old_bam_local"
"$tools/samtools" index "$old_bam_local"

"$tools/samtools" quickcheck -v "$bam"
"$tools/samtools" flagstat "$bam" > "$validation/retained.flagstat.txt"
"$tools/samtools" stats "$bam" > "$validation/retained.stats.txt"

"$tools/samtools" view -q 30 -F 2316 -f 2 -L "$whitelist" -bh "$bam" | \
  "$tools/samtools" sort -@ 4 -o "$validation/derived_strict.bam" -
"$tools/samtools" index "$validation/derived_strict.bam"

samblaster="$repo/tests/hpc_validation/bin/samblaster"
RG="@RG\tID:$sample\tSM:$sample\tPL:illumina"
"$tools/bwa-mem2" mem -M -R "$RG" -t 6 "$index" <(cat "$fq1") <(cat "$fq2") 2> "$validation/direct_strict.bwa.log" | \
  "$samblaster" --addMateTags 2> "$validation/direct_strict.samblaster.log" | \
  "$tools/samtools" view -q 30 -F 2316 -f 2 -L "$whitelist" -bhS - | \
  "$tools/samtools" sort -@ 4 -o "$validation/direct_strict.bam" -
"$tools/samtools" index "$validation/direct_strict.bam"

{
  printf 'metric\told_strict\tnew_retained\tderived_strict\tdirect_strict\n'
  printf 'records\t%s\t%s\t%s\t%s\n' \
    "$("$tools/samtools" view -c "$old_bam_local")" \
    "$("$tools/samtools" view -c "$bam")" \
    "$("$tools/samtools" view -c "$validation/derived_strict.bam")" \
    "$("$tools/samtools" view -c "$validation/direct_strict.bam")"
  printf 'cd99_records\t%s\t%s\t%s\t%s\n' \
    "$("$tools/samtools" view -c "$old_bam_local" chrX:2690000-2750000)" \
    "$("$tools/samtools" view -c "$bam" chrX:2690000-2750000)" \
    "$("$tools/samtools" view -c "$validation/derived_strict.bam" chrX:2690000-2750000)" \
    "$("$tools/samtools" view -c "$validation/direct_strict.bam" chrX:2690000-2750000)"
  printf 'peak_count\tNA\t%s\tNA\tNA\n' "$(wc -l < "$peaks")"
  printf 'frip\tNA\t%s\tNA\tNA\n' "$(perl -e '$n=shift; $t=shift; printf "%.8f", $t ? $n/$t : 0' "$("$tools/samtools" view -c -L "$peaks" "$bam")" "$("$tools/samtools" view -c "$bam")")"
} > "$validation/metrics.tsv"

"$tools/samtools" view "$validation/derived_strict.bam" | md5sum > "$validation/derived_strict.records.md5"
"$tools/samtools" view "$validation/direct_strict.bam" | md5sum > "$validation/direct_strict.records.md5"
diff -u "$validation/derived_strict.records.md5" "$validation/direct_strict.records.md5" > "$validation/strict_records_md5.diff" || true
