# Supported i-frame hybrid PSF grid — 2026-10-22

## Scope

Extend the frozen local-core/public-HSC-I-wing grid from cluster 1056 to all
13 existing K2.5/40%-coverage `i` frames with approved position-PSF mappings.
The 12 configurations remain splice radii 1.8/2.4/3.0 arcsec, 0.6-arcsec
blend, and public-wing amplitudes 0.50/0.65/0.80/1.00. Each configuration
uses the same 13 cluster/filter pairs. The cluster-center hybrid is used for
diffuse ICL; each satellite envelope uses that satellite's mapped local PSF
core joined to the public wing.

The v3 injection runner also computes a separate unit point-source profile
through each existing saved ICL mask. This noiseless profile is not added to
the redetection scenarios and is not a source-detection recovery test. The
redetection and fixed-mask branches for diffuse/satellite light remain
separate outcomes.

No HSC data request is included. Observed coadds and masks are read-only.
The grid estimates frame-to-frame sensitivity within this PSF-supported
subset; it is not the 800-cluster population estimate and cannot validate the
public wing, sky transfer, coadd covariance, or source completeness.

## Execution

Driver: [`run_core_wing_hybrid_grid.py`](run_core_wing_hybrid_grid.py), using
`--all-supported --filter i`. The model summaries and independent
frame/response/point-source validation are assembled by
[`summarize_core_wing_hybrid_grid.py`](summarize_core_wing_hybrid_grid.py).
Slurm wrapper: [`run_core_wing_hybrid_i_supported.sbatch`](run_core_wing_hybrid_i_supported.sbatch).

Job `432803` completed on `syntax`, `a10` partition, in 2:07:02 using one
CPU and 32 GB requested RAM, with no GPU request. All 12 configurations
completed across the same 13 frames; independent stored-response arithmetic
validation reported maximum error 0. Output root:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_core_wing_hybrid_i_supported13_2026-10-22_01/`.
The final assessment and v2 frame-set audit are in
[`core_wing_hybrid_i_supported_assessment_2026-10-22.md`](core_wing_hybrid_i_supported_assessment_2026-10-22.md).

Scientific state remains `science_measurement=false`, `psf_approved=false`,
`observational_validation_passed=false`, and `flux_covariance_available=false`.
