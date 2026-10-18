#!/usr/bin/env bash
#SBATCH --job-name=icl-plane-inj
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-plane-inj-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-plane-inj-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLBACKEND=Agg
cd /home/kjhan/BACKUP/Astrafek/OGFinder
exec flock /tmp/astrafex_icl_plane_transfer.lock /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/audit_icl_plane_injection.py \
  --sample /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13 \
  --input /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_pixel_extinction_2026-10-19_03/pixel_extinction_satellite_association.tsv \
  --work /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_plane_injection_2026-10-19_02 \
  --replicas 5000 --seed 20261019 --amplitude 1.0
