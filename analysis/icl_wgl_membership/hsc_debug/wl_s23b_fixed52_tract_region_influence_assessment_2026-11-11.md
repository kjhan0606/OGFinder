# S23B fixed-52 approximate tract-region influence (2026-11-11)

## Scope and method

This CPU analysis reused the approved 20-file HSC S19A shape-row export, the
pinned HSC Y3 SACC `n(z)`, and the saved HSC patch-QA table. No CAS query,
download, or image read was made. The fixed sample contains 52 lens sets and
260 selected controls. After the source-bin and calibration-row filters, 259
control centers have retained source rows.

The 840 patch-QA tracts were grouped by proximity of their spherical-mean
centers at 2 degrees. Lenses, controls, and source endpoints were assigned to
their nearest tract center. A pair was removed when either endpoint belonged
to the component being deleted. Both pooled-control and equal-control-center
profiles were recalibrated with the same minimum of four positive controls in
each radial bin. The leave-one-region-out common sample was then selected
again, so changes combine endpoint deletion and sample-composition change.

## Results

The tract-center graph has four components in the full patch table; this sample
occupies three. Its lens centers are distributed 28, 21, and 3 across them.
After filtering, the corresponding row-backed control counts are 139, 105,
and 15. Every deletion retains at least four positive controls per supported
lens set, but only 24, 31, and 49 of the original 52 sets remain.

Tangential shifts below are relative to the full 52-set estimator. The scale is
`10^13 M☉ Mpc⁻²`.

| Deleted component | Lens sets / row-backed controls | Common sets after deletion | Pooled largest fractional shift | Equal-center largest fractional shift |
|---:|---:|---:|---:|---:|
| 0 | 28 / 139 | 24 | −117.6% at 1.5–2.0 Mpc (−0.7938) | −101.1% at 1.5–2.0 Mpc (−0.7813) |
| 1 | 21 / 105 | 31 | +116.1% at 1.5–2.0 Mpc (+0.7838) | +99.6% at 1.5–2.0 Mpc (+0.7696) |
| 2 | 3 / 15 | 49 | −16.7% at 0.4–0.7 Mpc (−0.6700) | −17.3% at 0.4–0.7 Mpc (−0.7055) |

The two larger components materially change the outer profile. Removing either
one shifts the 1.5–2.0 Mpc value by about its full-sample amplitude, and the
sign changes between the two deletions. The smaller component has less
influence. Pooled and equal-center estimators show the same regional pattern.

## Assessment

Three occupied components do not support a stable full survey-region
covariance. These shifts are influence diagnostics, not uncertainties or
significance. Because the common lens-set sample changes from 52 to 24 or 31
after the dominant-region deletions, the shifts also include sample-selection
effects.

The source-row export contains coordinates but no tract or patch IDs. Nearest
tract-center assignment is an approximation; the patch-center graph is not an
official field map or polygon-overlap calculation. The analysis therefore
does not resolve the spatial-covariance gate. Keep `science_measurement=false`
and do not fit a halo mass or claim a detection.

## Reproducibility

- Runner: [`analyze_s23b_fixed52_tract_region_influence.py`](analyze_s23b_fixed52_tract_region_influence.py)
- Runner SHA-256: `d21057f9a5176f22fdd01d1bf18dbb6302ca26effd6892b6528cc6b3df015911`
- Approved row-plan SHA-256: `7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`
- Reference estimator/covariance summary SHA-256: `00c8427ef6810b78f0900da1d515e3fca2ed77233bd012ebbeddbf5bc5e5c033`
- Patch-QA SHA-256: `5c982c954b1a4026e5fe4ab75e978d4734b67d41e04aeddf8725eb5351ea3c64`
- Output directory (outside Git): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed52_tract_region_influence_2026-11-11_01/`
- Summary SHA-256: `afc151a6b189af182f975b7a70a2e8bc176b85fa62ab09bff7a2b8718cb41326`
- Region-deletion table SHA-256: `5724b6fadbbec70d1ecc806dacba390147c57b014ebd8b9134ae1aaa6ffa1833`
- Per-bin profile-shift table SHA-256: `b2f1e5ac95c42354cca8c3e593c8207527fcc23dae92e136c03f6593f5640d69`
- Generated assessment SHA-256: `d09b391ea3cbab14355039a0b8ffd52bf1cf33e2a2330dec0ae95ab00d341dd7`
