#!/usr/bin/env bash
#SBATCH --job-name=icl-psf-i2-proxy
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-psf-i2-proxy-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-psf-i2-proxy-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLBACKEND=Agg
cd /home/kjhan/BACKUP/Astrafek/OGFinder
data=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13
exec flock /tmp/astrafex_icl_redetection_injection.lock /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/extended_psf_injection.py \
  --reference-dir "$data/extended_psf_reference_2026-10-14_01" \
  --reference-psf "$data/extended_psf_reference_2026-10-14_01/HSC_PSF_i.fits" \
  --output-dir "$data/icl_extended_psf_i2_proxy_pilot_2026-10-21_01" \
  --cluster-id 10 \
  --filter i2
