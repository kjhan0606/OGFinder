#!/bin/bash
#SBATCH --job-name=hsc-icl-2d-stack
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-2d02-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-2d02-%j.err

set -euo pipefail
module purge
module load gnu13
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

REPO=/home/kjhan/BACKUP/Astrafek/OGFinder
SAMPLE=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13
WORK=$SAMPLE/icl_2d_matched_stack_2026-10-19_k2p5
PYTHON=/home/kjhan/miniconda3/bin/python3

if [[ -e "$WORK" ]]; then
  echo "Refusing to overwrite completed output in $WORK" >&2
  exit 2
fi
cd "$REPO"
exec flock /tmp/astrafex_heavy.lock "$PYTHON" \
  analysis/icl_wgl_membership/hsc_debug/stack_hsc_icl_target_control_2d.py \
  --sample "$SAMPLE" --work "$WORK" --mask-factor 2.5 --kpc-max 500 --n-side 201
