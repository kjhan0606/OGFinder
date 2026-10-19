#!/bin/bash
#SBATCH --job-name=icl-independent-support
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-independent-support-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-independent-support-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLBACKEND=Agg

REPO=/home/kjhan/BACKUP/Astrafek/OGFinder
SAMPLE=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13
PYTHON=/home/kjhan/miniconda3/bin/python3
WORK25=$SAMPLE/icl_2d_independent_support_k2p5_2026-10-19_01
WORK35=$SAMPLE/icl_2d_independent_support_k3p5_2026-10-19_01
WORK50=$SAMPLE/icl_2d_independent_support_k5p0_2026-10-19_01
COMMON25=$SAMPLE/icl_2d_matched_stack_2026-10-19_02
COMMON35=$SAMPLE/icl_2d_matched_stack_2026-10-19_k3p5
COMMON50=$SAMPLE/icl_2d_matched_stack_2026-10-19_k5p0

for work in "$WORK25" "$WORK35" "$WORK50"; do
  if [[ -e "$work" ]]; then
    echo "Refusing to overwrite $work" >&2
    exit 2
  fi
done

cd "$REPO"
exec flock /tmp/astrafex_heavy.lock bash -c '
set -euo pipefail
python=$1; repo=$2; sample=$3
shift 3
while (( $# )); do
  factor=$1; common=$2; work=$3
  shift 3
  echo "=== independent support K=$factor ===" >&2
  "$python" "$repo/analysis/icl_wgl_membership/hsc_debug/stack_hsc_icl_independent_support.py" \
    --sample "$sample" --mask-factor "$factor" --common-work "$common" --work "$work"
  "$python" "$repo/analysis/icl_wgl_membership/hsc_debug/verify_icl_independent_support.py" \
    --work "$work"
done
' bash "$PYTHON" "$REPO" "$SAMPLE" \
  2.5 "$COMMON25" "$WORK25" \
  3.5 "$COMMON35" "$WORK35" \
  5.0 "$COMMON50" "$WORK50"
