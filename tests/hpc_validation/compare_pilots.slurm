#!/bin/bash
#SBATCH --job-name=tachyon-compare
#SBATCH --partition=intel-sc3
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/storage/zhangkaiLab/hanlitian/macrophage/validation/20260723_tachyon_backend/logs/compare-%j.out
#SBATCH --error=/storage/zhangkaiLab/hanlitian/macrophage/validation/20260723_tachyon_backend/logs/compare-%j.err

set -euo pipefail
root=/storage/zhangkaiLab/hanlitian/macrophage/validation/20260723_tachyon_backend
export XDG_CACHE_HOME="$root/xdg-cache"
export TMPDIR="$root/tmp"
py=/storage/zhangkaiLab/hanlitian/software/micromamba/micromamba/envs/py311/bin/python
sam=/storage/zhangkaiLab/hanlitian/software/micromamba/micromamba/envs/py311/bin/samtools
"$py" "$root/code/ATAC-sm/tests/hpc_validation/compare_pilots.py" "$root" "$sam"
