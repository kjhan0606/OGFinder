#!/usr/bin/env bash
#SBATCH --job-name=icl-wide-psf-sens
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-wide-psf-sens-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-wide-psf-sens-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLBACKEND=Agg
cd /home/kjhan/BACKUP/Astrafek/OGFinder
data=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13
runner=analysis/icl_wgl_membership/hsc_debug/extended_psf_injection.py
reference_dir="$data/extended_psf_reference_2026-10-14_01"
exec flock /tmp/astrafex_icl_redetection_injection.lock bash -c '
  set -euo pipefail
  /home/kjhan/miniconda3/bin/python3 "$1" --reference-psf "$2/HSC_PSF_i.fits" \
    --reference-dir "$2" --output-dir "$3/icl_extended_psf_symmetric_2026-10-19_01"
  /home/kjhan/miniconda3/bin/python3 "$1" --reference-psf "$2/HSC_PSF_i_asymmetric.fits" \
    --reference-dir "$2" --output-dir "$3/icl_extended_psf_asymmetric_2026-10-19_01"
' bash "$runner" "$reference_dir" "$data"
