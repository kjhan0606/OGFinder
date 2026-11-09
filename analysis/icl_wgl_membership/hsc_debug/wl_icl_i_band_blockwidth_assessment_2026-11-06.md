# i-band ICL split weak-lensing assessment — 2026-11-06

## Purpose and method

The first ICL-to-weak-lensing split used only `i2`. This follow-up applies the
same Astrafek analysis to the saved `i`-band common-support sample, then changes
only the approximate equal-area delete-one block width (2.5°, 5°, 10°). It
uses HSC source tomography bins 2–4, the 200–300 kpc target-minus-control ICL
median split, one-to-one redshift/richness matching (`|Δz|≤0.03`,
`|Δln richness|≤0.25`), five matched random apertures per lens, and the same
random-subtracted shape rows at all three block widths. No new HSC request,
download, or image-pixel read was made.

The verified ICL/WL overlap contains 37 usable `i` lenses. The median split
has 19 low-ICL and 18 high-ICL lenses; matching retains only 10 pairs. The
matching loses nine low-ICL and eight high-ICL objects. The resulting groups
each contain 10 lenses. The random stack has 100 requested centers, with rows
for 97. The split threshold is 2.119 nJy arcsec⁻²; median 200–300 kpc ICL
contrasts are −1.053 and +8.035 nJy arcsec⁻² in the low/high groups.

## Results

The high-minus-low random-subtracted `DeltaSigma` profile is radially mixed.
Values and delete-one errors are in `10^12 M☉ Mpc⁻²`:

| Proper radius (Mpc) | High − low | Error, 2.5° | Error, 5° | Error, 10° |
|---:|---:|---:|---:|---:|
| 0.1–0.2 | −235.1 | 175.6 | 106.9 | 102.9 |
| 0.2–0.4 | −130.5 | 81.0 | 86.3 | 103.6 |
| 0.4–0.7 | −68.0 | 82.9 | 97.7 | 91.9 |
| 0.7–1.0 | −41.6 | 28.2 | 30.8 | 26.9 |
| 1.0–1.5 | +57.4 | 17.1 | 14.0 | 18.5 |
| 1.5–2.0 | +10.6 | 18.5 | 13.3 | 16.9 |
| 2.0–3.0 | −0.4 | 14.3 | 10.1 | 10.1 |

The point profile, exact matched-pair list, and target source-pair counts are
identical across widths (matched-pair SHA-256
`025840dba31db0e671916f752efcb735ba4064b2f1b3aa20ab605fa7622c562b`). The
covariance changes materially: the conditional seven-bin tangential Wald
diagnostic is χ²=15.79, rank 7, nominal p=0.027 at 2.5°; χ²=161.92, rank 7,
nominal p=1.27×10⁻³¹ at 5°; and χ²=101.28, rank 7, nominal p=5.87×10⁻¹⁹ at
10°. The corresponding cross-component null diagnostics are χ²=24.70
(nominal p=8.58×10⁻⁴), 90.97 (7.82×10⁻¹⁷), and 105.50 (7.87×10⁻²⁰).

These nominal values are not calibrated detection probabilities. There are
only 9–24 approximate blocks for 10 matched cluster pairs; at the 5° and 10°
widths, the number of blocks is close to the seven-dimensional profile rank.
The large cross-null diagnostic at every width flags strong cross-component
structure or unmodeled systematics. The steep change in tangential χ² with
block width shows that the apparent contrast is dominated by the assumed
spatial covariance. The radially mixed profile and this covariance behavior
do not support an ICL-strength–halo-mass interpretation.

## Assessment

Adding `i` does not strengthen the ICL–weak-lensing evidence. Its matched
sample is too small for stable spatial covariance, and its cross-component
null fails strongly under all three approximate block assignments. Together
with the existing `i2` result (29 matched pairs; tangential nominal p changes
from 2.42×10⁻⁴ at 5° to 0.099 at 2.5°, while cross-null is consistent with
zero under those approximate covariances), the two filters do not produce a
robust, mutually supported association.

This is an exploratory sample-size and systematics diagnostic, not a
population measurement. The ICL split remains conditional on imaging support;
source-Poisson correction, validated contiguous HSC covariance, extended-PSF
transfer, sky transfer, photo-z/shear calibration, and random-aperture
depth/mask/PSF matching remain open. Keep `science_measurement=false`.

## Approximate footprint balance check

The saved parent-sample tract-component labels place the 38-lens `i` common
sample in two of three occupied components (18 and 20 lenses). The 10 matched
pairs have low/high group counts 4/3 in component 0 and 6/7 in component 1;
five pairs are within one component and five cross components. For comparison,
the 70-lens `i2` common sample has 48/18/4 lenses across components 0/1/2; its
29 matched pairs have low/high counts 22/21, 4/7, and 3/1, with 17 same-region
and 12 cross-region pairs. Thus the `i` high/low groups are not grossly
imbalanced between the two components that contain them, while half the
matched pairs cross components. These approximate nearest-tract components
are not official HSC fields. This lens-only check does not map source or
random-center coordinates and cannot explain the `i` cross null.

The Astrafek audit is
[`audit_wl_icl_split_tract_region_balance.py`](audit_wl_icl_split_tract_region_balance.py)
(SHA-256 recorded in both output summaries). Its outputs are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_icl_split_tract_region_balance_2026-11-06_02/`
and
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_icl_split_tract_region_balance_i2_2026-11-06_02/`.

## Source- and random-position cross-component decomposition

The next Astrafek audit reapplied the exact HSC Y3 source cuts to the approved
shape-row export and stacked the low/high ICL groups and their five random
centers separately within each lens/center tract component. It reproduces
53,399 used target pairs and 241,258 used random pairs from the matched 10-pair
sample. For every used target and random pair, the source and its lens/center
are assigned to the same approximate tract component (fraction 1.0). Thus
cross-component source leakage does not explain the cross-null.

The random-subtracted high-minus-low cross profiles (in `10^12 M☉ Mpc⁻²`) are:

| Proper radius (Mpc) | Component 0 | Component 1 |
|---:|---:|---:|
| 0.1–0.2 | −101.3 | −190.2 |
| 0.2–0.4 | +101.1 | +179.7 |
| 0.4–0.7 | +44.4 | −16.9 |
| 0.7–1.0 | +19.7 | +94.8 |
| 1.0–1.5 | +36.4 | −72.0 |
| 1.5–2.0 | −53.2 | −2.6 |
| 2.0–3.0 | +2.6 | −13.0 |

The oscillatory cross contrast appears in both components; its first two bins
have the same signs and large amplitudes. This argues against a single
component or boundary-crossing source pairs being the sole cause. The regional
profiles have no covariance, and the groups contain only 3–7 lenses per
component, so they cannot identify a physical or instrumental origin.

Runner:
[`audit_wl_icl_split_source_region_cross_null.py`](audit_wl_icl_split_source_region_cross_null.py)
(SHA-256 `961874ffdeb80efa9abe09c519adac680ea35adf7bdc83adceff5043fbfd4fdb`).
Summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_icl_split_source_region_cross_null_2026-11-06_01/summary.json`
(SHA-256 `fb90c982da784dff60b270efdfe986d0fbe5b8f3a656e15dd5bc14629810862c`).

## Target-versus-random cross decomposition and radial cuts

The endpoint arithmetic shows that the first two radial bins of the
high-minus-low cross profile are primarily present in the target stacks. In
component 0, target high-minus-low cross values are −142.1 and +120.2, while
the corresponding random-center differences are −40.9 and +19.1. In component
1 the target values are −114.4 and +159.5; random differences are +75.8 and
−20.2, which reinforce the random-subtracted contrast. Values are in
`10^12 M☉ Mpc⁻²`. Random stacks therefore modify the result, but they do not
remove the target-side cross pattern.

A radial-cut audit then subselects the same saved profile and covariance.
After the post-hoc cut `R≥1 Mpc`, the three-bin tangential conditional nulls
have χ²/p = 11.50/0.0093 (2.5°), 21.97/6.63×10⁻⁵ (5°), and
22.45/5.27×10⁻⁵ (10°); the cross nulls are 3.61/0.307, 3.64/0.303, and
6.09/0.107. The inner cross inconsistency is thus concentrated below 1 Mpc,
while an outer tangential high/low contrast remains in these conditional
fits. This is an exploratory lead, not a detection: the cut was chosen after
seeing the cross profile, only 10 matched lens pairs are used, and the
covariance comes from 9–24 approximate blocks, not validated HSC fields. It
does not reproduce in the larger `i2` split: for `R≥1 Mpc`, the tangential
null has nominal p=0.172 at 2.5° (54 blocks) and p=0.114 at 5° (26 blocks),
while the cross null remains consistent with zero (p=0.639 and 0.550). The
outer tangential candidate is therefore filter/sample-specific in this
exploratory comparison.

Astrafek endpoint arithmetic:
[`audit_wl_icl_split_target_random_cross_decomposition.py`](audit_wl_icl_split_target_random_cross_decomposition.py),
runner SHA-256 `051bc67f3be1f23f0134d7509287fa201b2da7b2689d86302fbbf1b3e2b4f42c`,
summary SHA-256 `ffc453ac4b78e7f25a354d19469b31411856cc6d4e0ab41274f6940f331fe5af`.
Radial-cut audit:
[`audit_wl_icl_split_radial_cut_nulls.py`](audit_wl_icl_split_radial_cut_nulls.py),
runner SHA-256 `00a2f24b6a4c7e63b13f84c393d4778224db03b07a4ce5a6316f3b3f70717cf2`,
summary SHA-256 `583426fef8b60c011ea392dd80c3a7783d1572fa06006d72b729bebea1eade45`.

## Reproducibility

- Astrafek runner:
  `analysis/icl_wgl_membership/hsc_debug/analyze_wl_random_subtracted_icl_split.py`
  (SHA-256 `705d6389c70f1a7d6dce6ff2c07bc1517ede631b0a507757b81ef39417bb4a3f`).
- Exact matched sample: 10 pairs; matched-pair-list SHA-256
  `025840dba31db0e671916f752efcb735ba4064b2f1b3aa20ab605fa7622c562b`.
- 2.5° summary SHA-256 `5c62835fcae76a8e9251f05218c03c73629579bf15ef8460257aabfc1a1d5baf`.
- 5° summary SHA-256 `61bbf1d605ed15b7666aee3490ed91051fa22f7452956e839b3d4b41535733ff`.
- 10° summary SHA-256 `ad38c8b5d7747d264199eaf5a420e48206be8d5c1054b82971fafc98d822007f`.
- All summaries are under
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_icl_split_i_2026-11-06_{01,02,03}/`.
- Input provenance hashes (common lens set, random-subtracted rows, random
  plan, target stack, lens manifest, and HSC Y3 `n(z)`) are embedded in each
  summary.
