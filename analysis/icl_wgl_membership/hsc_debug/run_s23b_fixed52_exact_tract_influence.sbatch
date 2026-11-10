#!/usr/bin/env bash
#SBATCH --job-name=s23b-tract-exact
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl-s23b-exact-tract-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl-s23b-exact-tract-%j.err
set -euo pipefail

module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1

cd /home/kjhan/BACKUP/Astrafek/OGFinder
GEOMETRY_DIR="${HSC_PDR3_GEOMETRY_DIR:-/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/control_patches_2026-10-14_02/public_geometry}"
OUTPUT_DIR="${HSC_S23B_EXACT_TRACT_OUTPUT_DIR:-}"
test -d "$GEOMETRY_DIR" || {
  echo "Missing existing HSC PDR3 geometry directory: $GEOMETRY_DIR" >&2
  exit 2
}
test -n "$OUTPUT_DIR" || {
  echo "Set HSC_S23B_EXACT_TRACT_OUTPUT_DIR to a fresh output path" >&2
  exit 2
}

exec /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/analyze_s23b_fixed52_tract_region_influence.py \
  --geometry-dir "$GEOMETRY_DIR" \
  --output "$OUTPUT_DIR"
