#!/bin/bash
#SBATCH --job-name=hsc-icl-pixext
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/pixel-extinction-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/pixel-extinction-%j.err

set -euo pipefail
module purge
module load gnu13
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

REPO=/home/kjhan/BACKUP/Astrafek/OGFinder
SAMPLE=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13
PYTHON=/home/kjhan/miniconda3/bin/python3
WORK=${WORK:-$SAMPLE/icl_pixel_extinction_2026-10-19_01}
ARGS=()
if [[ -n "${LIMIT:-}" ]]; then
  ARGS+=(--limit "$LIMIT")
fi
if [[ -n "${FILTER:-}" ]]; then
  ARGS+=(--filter "$FILTER")
fi
if [[ -e "$WORK" ]]; then
  echo "Refusing to overwrite $WORK" >&2
  exit 2
fi
cd "$REPO"
exec flock /tmp/astrafex_heavy.lock "$PYTHON" \
  analysis/icl_wgl_membership/hsc_debug/measure_pixel_extinction_profiles.py \
  --data "$SAMPLE" --work "$WORK" --replicas 5000 "${ARGS[@]}"
