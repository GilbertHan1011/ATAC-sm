#!/bin/bash
#SBATCH --job-name=atac-tachyon-pilot
#SBATCH --partition=intel-sc3
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/storage/zhangkaiLab/hanlitian/macrophage/validation/20260723_tachyon_backend/logs/atac-%j.out
#SBATCH --error=/storage/zhangkaiLab/hanlitian/macrophage/validation/20260723_tachyon_backend/logs/atac-%j.err

set -euo pipefail

root=/storage/zhangkaiLab/hanlitian/macrophage/validation/20260723_tachyon_backend
repo=$root/code/ATAC-sm
snake=/storage/zhangkaiLab/hanlitian/software/micromamba/micromamba/envs/py311/bin/snakemake
sample=B56-1_ATAC19
export CONDA_PKGS_DIRS=$root/conda-pkgs
export TMPDIR=$root/tmp
export XDG_CACHE_HOME=$root/xdg-cache

run_arm() {
  local arm=$1 config=$2
  local out=$root/results/atac/$arm
  "$snake" \
    "$out/important_processed/bam/$sample.filtered.bam" \
    "$out/middle_files/bed/$sample.tagAlign.gz" \
    "$out/important_processed/peaks/${sample}_peaks.narrowPeak" \
    "$out/report/tss_coverage/$sample.tss_histogram.csv" \
    "$out/important_processed/tracks/$sample.bw" \
    "$out/report/ataqv/$sample.ataqv.json" \
    --directory "$repo" --snakefile "$repo/workflow/Snakefile" \
    --workflow-profile none --executor local --cores "$SLURM_CPUS_PER_TASK" \
    --resources mem_mb=120000 --use-conda --conda-prefix "$root/conda" \
    --rerun-incomplete --notemp --configfile "$repo/tests/hpc_validation/$config"
}

run_arm legacy_no_prealign config_legacy_no_prealign.yaml
run_arm tachyon config_tachyon.yaml
