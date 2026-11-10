# S23B fixed-52 tract influence audit

This CPU analysis reuses only the previously approved 20-file HSC S19A shape-row export and pinned HSC Y3 SACC n(z). It makes no CAS request, additional download, or image read. `science_measurement=false`.

Official HSC PDR3 Wide patch polygons assign the center and unique source endpoints to their parent tract IDs. Same-field patch overlaps use the nearest containing patch center; active endpoints that are unassigned or cross-field ambiguous stop the run. Each lens/control pair is removed if either endpoint maps to the deleted tract. The calibrated pooled-control and equal-control-center profiles are recomputed under the same four-control floor.

| Region | Lens centers | Control centers | Affected lens sets | Common supported sets after deletion | Minimum positive controls |
|---|---:|---:|---:|---:|---:|
| AUTUMN:tract-10180 | 0 | 1 | 1 | 51 | 4 |
| AUTUMN:tract-10183 | 1 | 7 | 2 | 50 | 4 |
| AUTUMN:tract-10419 | 1 | 1 | 1 | 51 | 4 |
| AUTUMN:tract-10420 | 0 | 4 | 1 | 51 | 4 |
| AUTUMN:tract-10422 | 1 | 4 | 1 | 51 | 4 |
| AUTUMN:tract-10423 | 0 | 0 | 1 | 52 | 4 |
| AUTUMN:tract-8280 | 1 | 1 | 1 | 51 | 4 |
| AUTUMN:tract-8281 | 0 | 0 | 1 | 52 | 4 |
| AUTUMN:tract-8284 | 0 | 0 | 1 | 52 | 4 |
| AUTUMN:tract-8285 | 0 | 0 | 1 | 52 | 4 |
| AUTUMN:tract-8520 | 0 | 3 | 1 | 51 | 4 |
| AUTUMN:tract-8521 | 0 | 4 | 1 | 51 | 4 |
| AUTUMN:tract-8522 | 0 | 0 | 1 | 52 | 4 |
| AUTUMN:tract-8525 | 0 | 4 | 1 | 51 | 4 |
| AUTUMN:tract-8526 | 2 | 6 | 2 | 50 | 4 |
| AUTUMN:tract-8762 | 1 | 2 | 1 | 51 | 4 |
| AUTUMN:tract-8763 | 0 | 0 | 1 | 52 | 4 |
| AUTUMN:tract-9226 | 0 | 2 | 1 | 51 | 4 |
| AUTUMN:tract-9456 | 0 | 2 | 1 | 51 | 4 |
| AUTUMN:tract-9457 | 0 | 1 | 1 | 52 | 4 |
| AUTUMN:tract-9466 | 0 | 0 | 1 | 52 | 4 |
| AUTUMN:tract-9469 | 1 | 2 | 1 | 51 | 4 |
| AUTUMN:tract-9470 | 0 | 1 | 1 | 52 | 4 |
| AUTUMN:tract-9471 | 1 | 4 | 1 | 51 | 4 |
| AUTUMN:tract-9692 | 0 | 3 | 1 | 51 | 4 |
| AUTUMN:tract-9693 | 0 | 1 | 1 | 52 | 4 |
| AUTUMN:tract-9698 | 1 | 0 | 1 | 51 | 4 |
| AUTUMN:tract-9699 | 2 | 11 | 4 | 48 | 4 |
| AUTUMN:tract-9700 | 2 | 7 | 3 | 49 | 4 |
| AUTUMN:tract-9707 | 1 | 4 | 3 | 51 | 4 |
| AUTUMN:tract-9708 | 1 | 8 | 3 | 50 | 4 |
| AUTUMN:tract-9709 | 1 | 5 | 1 | 51 | 4 |
| AUTUMN:tract-9710 | 1 | 1 | 1 | 51 | 4 |
| AUTUMN:tract-9711 | 0 | 5 | 2 | 50 | 4 |
| AUTUMN:tract-9714 | 0 | 0 | 1 | 52 | 4 |
| AUTUMN:tract-9934 | 0 | 1 | 1 | 52 | 4 |
| AUTUMN:tract-9935 | 1 | 0 | 1 | 51 | 4 |
| AUTUMN:tract-9941 | 1 | 7 | 3 | 50 | 4 |
| AUTUMN:tract-9942 | 0 | 0 | 2 | 52 | 4 |
| AUTUMN:tract-9949 | 0 | 2 | 2 | 51 | 4 |
| AUTUMN:tract-9950 | 1 | 1 | 1 | 51 | 4 |
| HECTOMAP:tract-15820 | 0 | 2 | 1 | 51 | 4 |
| HECTOMAP:tract-15821 | 0 | 0 | 1 | 52 | 4 |
| HECTOMAP:tract-15832 | 0 | 2 | 1 | 51 | 4 |
| HECTOMAP:tract-15998 | 0 | 0 | 1 | 52 | 4 |
| HECTOMAP:tract-15999 | 1 | 3 | 1 | 51 | 4 |
| HECTOMAP:tract-16009 | 1 | 0 | 1 | 51 | 4 |
| HECTOMAP:tract-16010 | 1 | 8 | 2 | 50 | 4 |
| SPRING:tract-9090 | 0 | 0 | 1 | 52 | 4 |
| SPRING:tract-9091 | 0 | 1 | 1 | 52 | 4 |
| SPRING:tract-9107 | 0 | 1 | 1 | 52 | 4 |
| SPRING:tract-9112 | 0 | 0 | 1 | 52 | 4 |
| SPRING:tract-9113 | 0 | 0 | 2 | 52 | 4 |
| SPRING:tract-9117 | 0 | 1 | 1 | 52 | 4 |
| SPRING:tract-9118 | 0 | 0 | 1 | 52 | 4 |
| SPRING:tract-9127 | 0 | 0 | 1 | 52 | 4 |
| SPRING:tract-9128 | 0 | 0 | 1 | 52 | 4 |
| SPRING:tract-9133 | 0 | 0 | 1 | 52 | 4 |
| SPRING:tract-9317 | 1 | 4 | 1 | 51 | 4 |
| SPRING:tract-9333 | 1 | 4 | 1 | 51 | 4 |
| SPRING:tract-9334 | 0 | 0 | 1 | 52 | 4 |
| SPRING:tract-9339 | 1 | 5 | 1 | 51 | 4 |
| SPRING:tract-9341 | 1 | 5 | 1 | 51 | 4 |
| SPRING:tract-9343 | 1 | 0 | 1 | 51 | 4 |
| SPRING:tract-9344 | 0 | 0 | 1 | 52 | 4 |
| SPRING:tract-9348 | 0 | 1 | 1 | 52 | 4 |
| SPRING:tract-9349 | 2 | 8 | 2 | 50 | 4 |
| SPRING:tract-9350 | 0 | 0 | 1 | 52 | 4 |
| SPRING:tract-9351 | 0 | 0 | 1 | 52 | 4 |
| SPRING:tract-9352 | 0 | 0 | 1 | 52 | 4 |
| SPRING:tract-9354 | 0 | 3 | 1 | 51 | 4 |
| SPRING:tract-9355 | 1 | 5 | 2 | 50 | 4 |
| SPRING:tract-9356 | 1 | 2 | 1 | 51 | 4 |
| SPRING:tract-9358 | 0 | 2 | 1 | 51 | 4 |
| SPRING:tract-9359 | 0 | 2 | 1 | 51 | 4 |
| SPRING:tract-9360 | 1 | 4 | 1 | 51 | 4 |
| SPRING:tract-9370 | 1 | 4 | 1 | 51 | 4 |
| SPRING:tract-9371 | 1 | 2 | 2 | 51 | 4 |
| SPRING:tract-9372 | 0 | 3 | 1 | 51 | 4 |
| SPRING:tract-9375 | 1 | 5 | 2 | 51 | 4 |
| SPRING:tract-9376 | 1 | 5 | 1 | 51 | 4 |
| SPRING:tract-9556 | 0 | 4 | 1 | 51 | 4 |
| SPRING:tract-9560 | 0 | 1 | 1 | 51 | 4 |
| SPRING:tract-9583 | 1 | 0 | 1 | 51 | 4 |
| SPRING:tract-9584 | 0 | 5 | 1 | 51 | 4 |
| SPRING:tract-9586 | 0 | 4 | 2 | 51 | 4 |
| SPRING:tract-9587 | 1 | 6 | 3 | 50 | 4 |
| SPRING:tract-9588 | 1 | 5 | 1 | 51 | 4 |
| SPRING:tract-9589 | 1 | 5 | 3 | 51 | 4 |
| SPRING:tract-9590 | 2 | 9 | 2 | 50 | 4 |
| SPRING:tract-9594 | 1 | 4 | 1 | 51 | 4 |
| SPRING:tract-9595 | 0 | 1 | 1 | 52 | 4 |
| SPRING:tract-9600 | 1 | 4 | 1 | 51 | 4 |
| SPRING:tract-9601 | 1 | 4 | 2 | 51 | 4 |
| SPRING:tract-9602 | 1 | 3 | 1 | 51 | 4 |
| SPRING:tract-9609 | 1 | 2 | 1 | 51 | 4 |
| SPRING:tract-9610 | 0 | 3 | 1 | 51 | 4 |
| SPRING:tract-9612 | 2 | 10 | 2 | 50 | 4 |
| SPRING:tract-9613 | 0 | 0 | 1 | 52 | 4 |
| SPRING:tract-9799 | 1 | 1 | 1 | 51 | 4 |
| SPRING:tract-9831 | 0 | 0 | 1 | 52 | 4 |
| SPRING:tract-9832 | 0 | 1 | 1 | 52 | 4 |

The reported regional shifts are influence diagnostics, not errors or significance. Individual adjacent HSC tracts are not independent blocks, so this runner does not estimate covariance. `region_profile_shifts.tsv` gives the per-bin pooled, equal-center, and paired-contrast shifts.

The fixed sample is selected for radial support, and the retained common lens sets can change after deletion. Other unresolved limits include n(z) uncertainty, reduced shear, magnification, member dilution/boost, shear-selection systematics, and PSF uncertainty. Keep `science_measurement=false`; do not infer a mass or detection.

- Runner SHA-256: `cead8c69dc852fb5c67a419cd78a4266c957732cbc1d482600b9cf6b3e0ecd80`
- Approved row-plan SHA-256: `7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`
- Reference estimator/covariance summary SHA-256: `00c8427ef6810b78f0900da1d515e3fca2ed77233bd012ebbeddbf5bc5e5c033`
- Official geometry directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/control_patches_2026-10-14_02/public_geometry`
- Official geometry parser SHA-256: `d15a79ec1ee085ec27b5e4b0c06c3bbd638a54841af2c6d4ecf6345e4516ae05`
- Result directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed52_exact_tract_influence_2026-11-11_02`
