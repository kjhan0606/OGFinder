# DESI HSC PDR3 `n(z)` sensitivity in the cluster weak-lensing stack — 2026-10-20

## Data and method

The user approved only the four author-published `_v_1p1` NPZ/JSON products at
commit `abc70af29910e49fce17a63ef3083160d515fb6c` (549,818 bytes total). The
files are stored outside Git under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/desi_y3_hsc_nz_abc70af2_2026-10-20/`.
The manifest pins each raw URL, size, and SHA-256. No image, new shape row,
individual photo-z, or other repository file was fetched. The NPZ schema has
four source bins on a shared 1,000-point `z=0..3` grid, with mean/median and
pointwise lower/upper curves plus magnification-corrected variants. It has no
joint posterior samples.

The stack used the existing, verified 800-lens S19A shape export unchanged:
3,069,544 lens–source rows, of which 2,254,646 have assigned HSC Y3 source bins
1–4. Each model used the same 5-degree, 52-cell spatial jackknife and the same
shape-selection, responsivity, and selection-bias calibration as the stored
HSC Y3 SACC baseline. The baseline was reused only after exact comparison of
all 32 shape-row batch hashes, plan hash, masked value, SACC hash, radial bins,
cosmology, and occupied-cell count. The lower/upper curves were evaluated as
profile sensitivities without a covariance; they are not independent posterior
realizations. This run has no random-lens subtraction and remains
`science_measurement=false`.

## Results

Relative to SACC, the central DESI curves raise the median `DeltaSigma` profile
in source bin 1 by 12.5–14.8% across the seven radial bins, while lowering bins
2–4 by roughly 1–3.5%. The sign and size pattern is stable across both scale
cuts and the raw/magnification-corrected central curves.

| DESI scale cut | Curve | Bin 1 | Bin 2 | Bin 3 | Bin 4 |
|---|---|---:|---:|---:|---:|
| 0.3–3 Mpc/h | mean | +12.9% | −2.1% | −3.1% | −2.9% |
| 0.3–3 Mpc/h | mean, magnification corrected | +14.8% | −1.1% | −2.7% | −2.8% |
| 1–5 Mpc/h | mean | +12.5% | −3.5% | −2.4% | −3.1% |
| 1–5 Mpc/h | mean, magnification corrected | +13.9% | −2.9% | −2.0% | −3.0% |

The DESI-vs-SACC median source-redshift shifts also match the direction in the
paper: bin 1 is lower by 0.023–0.026, bin 2 higher by 0.005–0.024, bin 3 higher
by 0.038–0.043, and bin 4 higher by 0.050–0.058, depending on scale cut and
magnification correction.

The 4-bin tomography-consistency diagnostic improves numerically under the
DESI central curves but still shows disagreement at 0.1–0.7 Mpc:

| Calibration | 0.1–0.2 Mpc | 0.2–0.4 Mpc | 0.4–0.7 Mpc |
|---|---:|---:|---:|
| HSC Y3 SACC baseline | 1.05×10⁻⁵ | 8.60×10⁻¹⁰ | 9.48×10⁻⁵ |
| DESI 0.3–3, mean | 1.03×10⁻⁴ | 6.55×10⁻⁸ | 7.93×10⁻⁴ |
| DESI 0.3–3, mean μ-corrected | 1.29×10⁻⁴ | 9.57×10⁻⁸ | 1.17×10⁻³ |
| DESI 1–5, mean | 1.05×10⁻⁴ | 4.36×10⁻⁸ | 5.69×10⁻⁴ |
| DESI 1–5, mean μ-corrected | 1.24×10⁻⁴ | 5.59×10⁻⁸ | 7.27×10⁻⁴ |

These are per-radius raw 4-bin chi-square diagnostic p-values under the
jackknife covariance Gaussian approximation. Every central curve retains a
small p-value at all three radii; in particular, the 0.2–0.4 Mpc discrepancy
remains. The DESI calibration changes the stack materially in bin 1 but does
not reconcile the source-bin profiles. The p-values compare each calibration
to internal tomography consistency; they are not a formal test that the DESI
profiles differ significantly from SACC, since a paired covariance of the
model differences was not computed.

Each of the four central models retained 52 occupied cells, a 56×56 covariance
of rank 51, and identical per-bin shape-row counts. Relative covariance
asymmetry was below `4.1×10⁻¹⁷`. An independent validator checked the pinned
manifest and all four artifact hashes, baseline/shape provenance, 12 scenario
profiles, stored comparison arithmetic, covariance rank/symmetry, and profile
counts. Its status is `validated_cluster_only_sensitivity`.

## Interpretation and limits

The DESI redshift calibration is important for the amplitude comparison across
source bins, particularly because source bin 1 moves the opposite direction
from bins 2–4. It reduces the source-bin inconsistency statistic relative to
SACC at these radii, but the internal discrepancy remains strong. This cannot
yet be called a lensing detection or a correction to the cluster masses:
random-lens subtraction is absent, the published lower/upper curves are only
marginal pointwise envelopes, the equal-area cells can join disconnected HSC
fields, and remaining shear/redshift calibration effects are not all
propagated.

The author repository reports no SPDX license; the cited paper is CC BY 4.0.
The imported results retain their source and pinned revision. This note does
not make a legal determination about repository artifacts.

## Reproducibility

- Runner: [`audit_wl_desi_nz_distribution_sensitivity.py`](audit_wl_desi_nz_distribution_sensitivity.py)
- Validator: [`validate_wl_desi_nz_distribution_sensitivity.py`](validate_wl_desi_nz_distribution_sensitivity.py)
- Tests: [`test_audit_wl_desi_nz_distribution_sensitivity.py`](test_audit_wl_desi_nz_distribution_sensitivity.py)
- Slurm job: `430481`, completed successfully in 37m 41s.
- Result directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/desi_y3_hsc_nz_abc70af2_2026-10-20/spatial_jackknife_2026-10-20_01/`.
- `summary.json` SHA-256: `f0da334af3393ac1ffe27a5948dc0be385c856094a842c5fe11e11e54f979113`.
- `validation.json` SHA-256: `dcc1a49fc2c8dda85b0ccf307a52569131fc6a2766ccd3fc4e1a50a71606e4a3`.
