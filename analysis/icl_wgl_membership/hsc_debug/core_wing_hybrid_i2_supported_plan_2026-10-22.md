# Supported i2-frame hybrid PSF sensitivity plan — 2026-10-22

## Purpose

Run the same frozen local-core/public-HSC-I-wing injection sensitivity grid
already applied to 13 PSF-supported `i` frames on the 44 supported `i2` frames.
The `i2` sample contains the exploratory ICL–satellite endpoint, so this tests
how that endpoint's measurement operator responds to the same deliberately
varied PSF splice models. It does not validate the public wing as an `i2` PSF,
estimate an astrophysical association, or provide a correction factor.

## Frozen scope

- Inputs: existing HSC cutouts, masks, position-specific PSF cores, and the
  public HSC-I extended reference already present under `/gpfs`.
- Frames: all 44 `i2` cluster/filter frames with the approved local-core map.
- Grid: splice radii 1.8, 2.4, and 3.0 arcsec; 0.6-arcsec blend; public-wing
  amplitudes 0.50, 0.65, 0.80, and 1.00 (12 configurations).
- Injections: diffuse ICL, satellite envelopes, and separate noiseless unit
  point-source profiles through the saved masks, following the existing
  `extended_psf_injection.py` estimator paths.
- No HSC service, CAS query, credential, or new download is involved.
- Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_core_wing_hybrid_i2_supported44_2026-10-22_01/`.
- Driver: [`run_core_wing_hybrid_grid.py`](run_core_wing_hybrid_grid.py);
  Slurm entry point:
  [`run_core_wing_hybrid_i2_supported.sbatch`](run_core_wing_hybrid_i2_supported.sbatch).

## Interpretation gate

Every configuration must contain the identical 44 cluster IDs, complete
without frame errors, and pass the existing response arithmetic validator.
Summaries are across this PSF-supported subset and are model-sensitivity
descriptions, not confidence intervals. The generic HSC-I wing remains
position-, epoch-, and seeing-unmatched to these `i2` coadds. Preserve
`science_measurement=false`, `psf_approved=false`,
`observational_validation_passed=false`, and
`flux_covariance_available=false` regardless of whether the grid is stable.

The result can identify whether the ICL and satellite transfer operator is
sensitive to splice choices on the primary association filter. It cannot close
the sky-transfer, source-photon-noise, mask/completeness, or matched-control
null gates.
