#!/bin/bash
#SBATCH --job-name=hsc-icl-pilot
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/logs/hsc-icl-pilot-%j.out
set -euo pipefail
module purge
module load gnu13
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
export PYTHONDONTWRITEBYTECODE=1
PILOT_REPO=${ASTRAFEX_OGF_ROOT:-/home/kjhan/BACKUP/Astrafek/OGFinder}
PILOT_PYTHON=${ASTRAFEK_PYTHON:-/home/kjhan/miniconda3/bin/python3}
cd "$PILOT_REPO"
if [ "${1:-}" = "--test" ]; then
    exec flock /tmp/astrafex_heavy.lock "$PILOT_PYTHON" -m pytest -p no:cacheprovider \
      icl/tests/test_pilot.py icl/tests/test_wgl.py plugins/icl/tests/test_icl_stack.py -q
fi
if [ "$#" != 2 ]; then
    printf 'Usage: sbatch pilot_job.sh MANIFEST NEW_OUTPUT_DIRECTORY\n' >&2
    exit 2
fi
exec flock /tmp/astrafex_heavy.lock "$PILOT_PYTHON" \
    analysis/icl_wgl_membership/hsc_pilot.py measure --manifest "$1" --work "$2"
