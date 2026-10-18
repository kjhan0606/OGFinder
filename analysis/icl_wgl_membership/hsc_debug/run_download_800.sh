#!/usr/bin/env bash
set -euo pipefail
case "$(hostname | tr '[:upper:]' '[:lower:]')" in
  lageunha*) ;;
  *) echo 'Run this CPU-only download on lageunha.' >&2; exit 1 ;;
esac
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec python3 /home/kjhan/BACKUP/Astrafek/OGFinder/analysis/icl_wgl_membership/hsc_debug/download_pilot.py \
  --user kjhan0606 --resume \
  --manifest /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/sample_manifest.json \
  --work /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/download_01
