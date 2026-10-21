#!/usr/bin/env bash
#SBATCH --job-name=icl-corewing-pilot
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-11%4
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/corewing-pilot-%A_%a.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/corewing-pilot-%A_%a.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLBACKEND=Agg

repo=/home/kjhan/BACKUP/Astrafek/OGFinder
data=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13
reference_dir="$data/extended_psf_reference_2026-10-14_01"
runner="$repo/analysis/icl_wgl_membership/hsc_debug/extended_psf_injection.py"
validator="$repo/analysis/icl_wgl_membership/hsc_debug/validate_extended_psf_injection_batch.py"
python=/home/kjhan/miniconda3/bin/python3

# Predeclared 3 x 4 grid: splice radii 1.8/2.4/3.0 arcsec;
# wing amplitudes 0.50/0.65/0.80/1.00; fixed smooth blend width 0.6 arcsec.
radii=(1.8 2.4 3.0)
amplitudes=(0.50 0.65 0.80 1.00)
task=${SLURM_ARRAY_TASK_ID:?SLURM array task ID is required}
(( task >= 0 && task < 12 )) || { echo "invalid array index: $task" >&2; exit 2; }
radius=${radii[$((task / 4))]}
amplitude=${amplitudes[$((task % 4))]}
radius_tag=${radius/./p}
amplitude_tag=${amplitude/./p}
output="$data/icl_hybrid_psf_i_pilot_2026-10-21_r${radius_tag}_a${amplitude_tag}"

cd "$repo"
"$python" "$runner" \
  --reference-psf "$reference_dir/HSC_PSF_i.fits" \
  --reference-dir "$reference_dir" \
  --output-dir "$output" \
  --filter i --cluster-id 1056 \
  --hybrid-splice-radius-arcsec "$radius" \
  --hybrid-blend-width-arcsec 0.6 \
  --hybrid-wing-amplitude "$amplitude"
"$python" "$validator" \
  --batch-dir "$output" --filter i --expected-frames 1 \
  --output "$output/validation.json"
