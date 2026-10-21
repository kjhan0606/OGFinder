# Local-core/public-wing hybrid PSF injection pilot — 2026-10-21

## Predeclared scope

This pilot applies the existing HSC-I extended reference to the median-support
`i` frame, cluster `1056` (10 injected satellite centers among the 13-frame
support-selected `i` set). It uses 12 predeclared combinations: splice radii
1.8, 2.4, and 3.0 arcsec; a fixed 0.6-arcsec smooth blend; and public-wing
amplitudes 0.50, 0.65, 0.80, and 1.00. The amplitude range brackets the
previously measured overlap scales. The code constructs one hybrid kernel
from the cluster-center local core for diffuse ICL injection and a separate
hybrid kernel from each satellite's own local core for its envelope injection.

The existing image-level runner measures diffuse-to-ICL response in
200–300 kpc and satellite-envelope response in 50–75 kpc, with both
redetection and fixed-baseline-mask branches. It does not inject isolated point
sources in this pilot; point-source wing recovery remains a separate required
validation before any PSF approval. The public model is symmetric HSC-I and
is not matched to the exact visit/epoch. This is a sensitivity pilot, not an
observational PSF calibration or an ICL science result.

No HSC data request is part of this run. Original coadds and masks are
read-only. Each grid cell has a unique output directory and the driver writes
an atomic checkpoint after each validated configuration, allowing safe resume
without overwriting results.

## Execution and verification state

Astrafek driver: [`run_core_wing_hybrid_grid.py`](run_core_wing_hybrid_grid.py).
Slurm wrapper: [`run_core_wing_hybrid_i_pilot.sbatch`](run_core_wing_hybrid_i_pilot.sbatch).
The independent batch validator also checks every satellite's local-core
filename, checksum, hybrid parameters, signed normalization, and false
approval flags. Focused tests for the PSF constructor, injection paths,
validator, and grid definition passed (22 tests).

Slurm job `432792` was submitted on `syntax` to partition `a10` with one CPU,
32 GB RAM, no GPU request, and a four-hour limit. The initial status check
reported `RUNNING` on `syn04`. Output root:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_core_wing_hybrid_i_pilot_cluster1056_2026-10-21_01/`.
The final report will append per-configuration completion, arithmetic
validation, response ranges, and file hashes after the job reaches a terminal
state.

`science_measurement=false`, `psf_approved=false`,
`observational_validation_passed=false`, and
`flux_covariance_available=false` remain required invariants.
