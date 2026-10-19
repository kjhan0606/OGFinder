#!/usr/bin/env bash
# CPU-only stack/covariance calculation; reserve an explicit shared-node slice.
# Supply paths through sbatch --export=ALL,... after the random row export is verified.
# This job performs no HSC/CAS query and stores no credentials.
#SBATCH --job-name=wl-random-stack
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_stack_%j.out

set -euo pipefail

ASTRAFEK_REPO_DIR="/home/kjhan/BACKUP/Astrafek/OGFinder"
ASTRAFEK_PYTHON_BIN="/home/kjhan/miniconda3/bin/python3"
module purge
module load python/3.13
cd "$ASTRAFEK_REPO_DIR"

: "${RANDOM_ROWS_DIR:?set RANDOM_ROWS_DIR to the completed random shape-row export}"
: "${RANDOM_PLAN_DIR:?set RANDOM_PLAN_DIR to the immutable random row-plan directory}"
: "${LENS_ROWS_DIR:?set LENS_ROWS_DIR to the verified cluster shape-row export}"
: "${LENS_PREFLIGHT_DIR:?set LENS_PREFLIGHT_DIR to the verified cluster count/preflight directory}"
: "${LENS_PLAN_DIR:?set LENS_PLAN_DIR to the verified cluster row-plan directory}"
: "${LENS_STACK_JSON:?set LENS_STACK_JSON to the verified cluster-only stack summary.json}"
: "${NZ_SACC:?set NZ_SACC to the local HSC Y3 SACC file}"
: "${RANDOM_STACK_OUTPUT:?set RANDOM_STACK_OUTPUT to a new output JSON path}"

if [[ ! -x "$ASTRAFEK_PYTHON_BIN" ]]; then
  echo "Astrafek Python is not executable: $ASTRAFEK_PYTHON_BIN" >&2
  exit 2
fi
for ASTRAFEK_INPUT in \
  "$RANDOM_ROWS_DIR/row_status.json" \
  "$RANDOM_ROWS_DIR/row_summary.json" \
  "$RANDOM_PLAN_DIR/random_shape_row_plan.json" \
  "$LENS_ROWS_DIR/row_status.json" \
  "$LENS_PREFLIGHT_DIR/count_status.json" \
  "$LENS_PLAN_DIR/row_export_plan.json" \
  "$LENS_STACK_JSON" "$NZ_SACC"; do
  if [[ ! -f "$ASTRAFEK_INPUT" ]]; then
    echo "Required verified input is missing: $ASTRAFEK_INPUT" >&2
    exit 2
  fi
done
if [[ -e "$RANDOM_STACK_OUTPUT" ]]; then
  echo "Refusing to overwrite existing output: $RANDOM_STACK_OUTPUT" >&2
  exit 2
fi

read -r -a ASTRAFEK_RADIAL_EDGES <<< "${RADIAL_EDGES_MPC:-0.1 0.2 0.4 0.7 1.0 1.5 2.0 3.0}"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

"$ASTRAFEK_PYTHON_BIN" analysis/icl_wgl_membership/hsc_debug/run_wl_random_shape_stack.py \
  --rows "$RANDOM_ROWS_DIR" \
  --plan-dir "$RANDOM_PLAN_DIR" \
  --nz-sacc "$NZ_SACC" \
  --b-mode-masked-value "${B_MODE_MASKED_VALUE:-0}" \
  --radial-edges-mpc "${ASTRAFEK_RADIAL_EDGES[@]}" \
  --output "$RANDOM_STACK_OUTPUT" \
  --lens-stack "$LENS_STACK_JSON" \
  --lens-rows "$LENS_ROWS_DIR" \
  --lens-preflight "$LENS_PREFLIGHT_DIR" \
  --lens-plan "$LENS_PLAN_DIR" \
  --block-width-deg "${BLOCK_WIDTH_DEG:-5.0}" \
  --h0 "${H0_KM_S_MPC:-70.0}" \
  --omega-m "${OMEGA_M:-0.3}"
