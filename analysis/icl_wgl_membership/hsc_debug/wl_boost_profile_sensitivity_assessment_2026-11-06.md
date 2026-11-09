# Weak-lensing boost-profile sensitivity (2026-11-06)

## Scope and method

This diagnostic uses only the already approved and locally saved HSC S19A
shape rows for the 435 common-support cluster lenses and 2,175 planned random
centers (five per lens). It makes no HSC request, downloads no data, reads no
images, and leaves the saved weak-lensing profile unchanged.

For each radial bin, the candidate profile is the exact S19A estimator-weight
ratio, `B_W(R)`, multiplied by the saved random-subtracted tangential and
cross `DeltaSigma` profiles. `B_W` is calculated from the response-adjusted
S19A estimator weight, including lens-distance-dependent `Sigma_crit^-2`;
the raw source-count ratio is not used. The nonlinear product is recomputed
inside each paired delete-one spatial-block replicate, preserving the
covariance between boost and lensing. Results are shown for the existing 2.5°
and 5° block choices.

This is a sensitivity calculation only. Source excess can arise from member
dilution, magnification, depth, masking, PSF, sky transfer, or source-selection
effects. The calculation does not identify or calibrate those mechanisms.

## Results

The exact estimator-weight ratio ranges from 1.125 to 1.163. Its paired
delete-one-block jackknife standard error is about 0.03–0.05, depending on
radius and block width. The candidate tangential profile is 12.5–16.3% higher
than the saved random-subtracted profile at every radius. The 2.5° and 5°
candidate central values agree to stored precision.

| Proper radius (Mpc) | `B_W` | Saved `DeltaSigma` (`10^13 M_sun Mpc^-2`) | Candidate `B_W DeltaSigma` (`10^13 M_sun Mpc^-2`) | Candidate jackknife error (`10^13 M_sun Mpc^-2`, 2.5° / 5°) |
|---:|---:|---:|---:|---:|
| 0.1–0.2 | 1.139 | 16.01 | 18.23 | 1.78 / 1.96 |
| 0.2–0.4 | 1.142 | 11.17 | 12.76 | 1.03 / 0.97 |
| 0.4–0.7 | 1.163 | 7.33 | 8.52 | 0.69 / 0.77 |
| 0.7–1.0 | 1.132 | 5.17 | 5.85 | 0.53 / 0.54 |
| 1.0–1.5 | 1.125 | 3.94 | 4.43 | 0.33 / 0.34 |
| 1.5–2.0 | 1.134 | 2.68 | 3.03 | 0.31 / 0.36 |
| 2.0–3.0 | 1.126 | 1.68 | 1.89 | 0.20 / 0.19 |

The candidate cross profile remains small compared with its paired spatial
jackknife errors. No significance claim is made from the conditional
chi-square diagnostic in the machine summary: it tests this constructed
sensitivity profile against zero while holding the unvalidated boost model
fixed, so it is not a detection test for a physically calibrated signal.

## Interpretation and remaining gates

The result quantifies the numerical effect of applying the measured
estimator-weight ratio; it does not justify applying that ratio as a physical
member-dilution correction. Full-aperture matched depth, mask, PSF, and sky
response are still unverified. Magnification and source-selection effects,
uncertainty in the source-redshift distribution, reduced shear, source-Poisson
noise, and regional survey covariance are not fully propagated. The current
HSC region partition has only a few occupied connected components.

Do not update the saved science profile. Keep
`boost_correction_applied=false` and `science_measurement=false`. Revisit the
correction only after the physical dilution interpretation and coverage,
selection, and covariance gates are validated.

## Reproduction and provenance

Runner: [`analyze_wl_source_density_boost_blocks.py`](analyze_wl_source_density_boost_blocks.py).
The run uses existing local shape rows and stores the full paired nonlinear
jackknife covariance and seven-bin table under:

`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_density_boost_spatial_2026-11-05_06/`

The summary records `status=complete_diagnostic`,
`hsc_network_used=false`, `images_read=false`,
`boost_correction_applied=false`, and `science_measurement=false`.

SHA-256:

- Runner: `779037a82edab9d94f80b3623875c8e88a458e3bd4bc67fca546b83936ab4360`
- Summary: `c98ad74a4e6201894d8f2f62e229cfd94a5e5a0f47f73f72c83b2f5345c0cbe3`
- Seven-bin table: `5796a119087e1f19f0bc738402f25b2565de756749d15943c95afe0a23cdd231`
