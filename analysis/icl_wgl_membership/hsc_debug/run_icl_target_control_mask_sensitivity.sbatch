#!/bin/bash
#SBATCH --job-name=hsc-icl-mask-sens
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-mask-sensitivity-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-mask-sensitivity-%j.err

set -euo pipefail
module purge
module load gnu13
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

REPO=/home/kjhan/BACKUP/Astrafek/OGFinder
SAMPLE=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13
WORK35=$SAMPLE/icl_2d_matched_stack_2026-10-19_k3p5
WORK50=$SAMPLE/icl_2d_matched_stack_2026-10-19_k5p0
PYTHON=/home/kjhan/miniconda3/bin/python3

for work in "$WORK35" "$WORK50"; do
  if [[ -e "$work" ]]; then
    echo "Refusing to overwrite $work" >&2
    exit 2
  fi
done

cd "$REPO"
exec flock /tmp/astrafex_heavy.lock bash -c '
set -euo pipefail
python=$1
repo=$2
sample=$3
work35=$4
work50=$5
cd "$repo"
for spec in "3.5:$work35" "5.0:$work50"; do
  factor=${spec%%:*}
  work=${spec#*:}
  "$python" analysis/icl_wgl_membership/hsc_debug/stack_hsc_icl_target_control_2d.py \
    --sample "$sample" --work "$work" --mask-factor "$factor" --kpc-max 500 --n-side 201
  "$python" analysis/icl_wgl_membership/hsc_debug/verify_icl_target_control_2d_stack.py \
    --sample "$sample" --work "$work"
  "$python" analysis/icl_wgl_membership/hsc_debug/evaluate_icl_target_control_2d_stack.py \
    --work "$work"
done
' bash "$PYTHON" "$REPO" "$SAMPLE" "$WORK35" "$WORK50"
