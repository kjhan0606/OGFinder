# HSC S19A raw tangential-shape profile diagnostic (2026-10-19)

## Scope

The approved 800-lens source-row export was used to calculate a descriptive,
shape-weighted mean raw tangential ellipticity `e_t` and cross component `e_x`.
The calculation applies the HSC Y3 source cuts `|e|<2`,
`0<=sigma_e<=0.4`, `resolution>=0.3`, `mag_A<=25.5`, and source tomographic
bins 1–4. It uses proper projected radii in a flat `H0=70 km/s/Mpc`,
`Omega_m=0.3` cosmology, with radial edges 0.1, 0.2, 0.4, 0.7, 1.0, 1.5,
2.0, and 3.0 Mpc.

This is **not** a calibrated shear or `DeltaSigma` measurement. It applies no
responsivity, multiplicative/additive, or selection-bias correction; it has no
source `n(z)` lensing efficiency or foreground-dilution correction, no random-
lens subtraction, and no covariance for repeated sources or shared sky blocks.
The stored `b_mode_mask` polarity remains undecoded, so both possible exclusions
are reported alongside the unmasked diagnostic. No data were downloaded for
this calculation.

## Raw profile

For source position angle `phi` measured east of the lens toward north,
`e_t=-e1*cos(2 phi)-e2*sin(2 phi)` and
`e_x=e1*sin(2 phi)-e2*cos(2 phi)`. The table gives the shape-weighted means.
Counts are lens–source pairs; one source can contribute to multiple lenses.

| Proper R (Mpc) | Pairs, no B-mode exclusion | Lenses | mean raw `e_t` | mean raw `e_x` | mean raw `e_t`, excluding raw 0 | mean raw `e_t`, excluding raw 1 |
|---|---:|---:|---:|---:|---:|---:|
| 0.1–0.2 | 11,681 | 561 | 0.034340 | −0.000757 | 0.035105 | 0.014120 |
| 0.2–0.4 | 39,333 | 586 | 0.029459 | −0.000545 | 0.029010 | 0.041421 |
| 0.4–0.7 | 96,297 | 600 | 0.022836 | −0.001068 | 0.022837 | 0.022792 |
| 0.7–1.0 | 136,858 | 609 | 0.019620 | −0.000429 | 0.019522 | 0.022234 |
| 1.0–1.5 | 315,199 | 617 | 0.014335 | −0.000226 | 0.014397 | 0.012616 |
| 1.5–2.0 | 428,222 | 622 | 0.010671 | +0.000370 | 0.010665 | 0.010831 |
| 2.0–3.0 | 1,202,278 | 632 | 0.005649 | +0.000075 | 0.005646 | 0.005734 |

The shape cuts and tomographic assignment retain 2,234,510 of 3,069,544
exported pairs before radial and B-mode selection. Counts by source bin 1–4 are
604,274, 712,085, 590,222, and 327,929. In the no-mask diagnostic, the raw
`e_t` mean decreases with projected radius and `e_x` stays near zero. This is
qualitatively compatible with a tangential alignment pattern, but it is not
evidence of a calibrated lensing detection: source-bin foreground dilution,
shape response, additive terms, random sky structure, shared sources, and field
covariance have not been removed or propagated. The bin-1 inner mean is close
to zero while higher source bins have larger inner means, which also shows why
the calibrated redshift distributions and dilution correction are necessary.

## B-mode polarity sensitivity

The raw catalog bit counts are 105,692 zeros and 2,963,852 ones. Within the
0.1–3.0 Mpc radial range, the no-exclusion profile contains 2,229,868 selected
pairs; excluding raw zero retains 2,153,004, while excluding raw one retains
76,864. Excluding raw zero changes the all-bin `e_t` profile only slightly.
Excluding raw one leaves only 20–26 lenses per radial bin and the cross
component becomes appreciable in several bins.
These are sensitivity scenarios, not evidence for which bit value is masked;
the frequency imbalance does not decode the bit. A release schema or validated
catalog convention is still needed before a masked science stack.

The calculation verifies that the pair counts from the two opposite exclusions
sum to the no-mask count in every radial and source-bin group. The final product
contains 105 profile rows (3 B-mode scenarios × 5 source-bin groups × 7 radial
bins), and its counts and status were checked independently.

## Reproducibility and status

The analysis is `audit_wl_raw_shape_profiles.py`; the relevant adapter change
allows the export-integrity loader to validate the binary mask without assigning
its polarity. Twenty-three focused tests pass, including HSC selection,
mask-polarity, angle-sign, and radial-edge checks.

The hash-recorded run is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_raw_shape_profiles_2026-10-19_02/`.
It records the 32 source-batch hashes, row plan/count provenance, code hashes,
profile hash, and `science_measurement=false`. `profiles.tsv` SHA-256 is
`8251fdb3ee19ada63d3b4ec5e4c2b9f20dfad88886ae86933cc47572e099e779`;
`summary.json` SHA-256 is
`572d8e02209c278164b9c918e18963b280c40d4e12577a5ae1fb072c90f396ec`.

The next calibrated weak-lensing step still requires the official HSC Y3 source
`n(z)`, verified mask polarity, random-lens pairs, per-source response and
selection corrections, and covariance that accounts for source overlap and
spatial blocks. The calibrated stack and cluster-mass inference remain pending.
