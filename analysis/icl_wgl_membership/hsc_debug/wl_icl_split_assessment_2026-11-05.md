# ICL-split, random-subtracted weak-lensing assessment (2026-11-05)

## Scope and method

Astrafek split the existing `i2` ICL/WL common-support sample at the median
target-minus-control ICL surface-brightness contrast in the 200–300 kpc ring.
It used only source tomography bins 2–4, excluding bin 1 because earlier
weak-lensing diagnostics showed strong source-density/boost sensitivity there.
High- and low-ICL lenses were paired one-to-one with `|Δz| ≤ 0.03` and
`|Δln(richness)| ≤ 0.25`, without using the shear profiles. Each group used the
same saved lens shape rows and five previously matched random apertures per
lens. The reported contrast is high-ICL minus low-ICL after subtracting the
corresponding random stack from each group.

The exact common sample has 70 `i2` lenses (35 per median-split group). Matching
retained 29 pairs (58 lenses); six from each group had no partner within both
calipers. The mean and maximum matched absolute redshift differences are
0.00958 and 0.02583. Median CAMIRA richness is 18.97 in the low-ICL group and
19.47 in the high-ICL group. Median ICL contrasts are −0.886 and
6.360 nJy arcsec⁻², respectively. The five random apertures per lens yield 290
requested centers, 267 of which contain usable rows in source bins 2–4.

## Results

The high-minus-low tangential contrast is radially mixed. Values and delete-one
block errors are in `10^12 Msun Mpc^-2`:

| Proper radius (Mpc) | High minus low ΔΣ | Error, 5° blocks | Error, 2.5° blocks |
|---:|---:|---:|---:|
| 0.1–0.2 | −39.9 | 97.8 | 87.1 |
| 0.2–0.4 | +5.1 | 44.4 | 56.0 |
| 0.4–0.7 | +25.6 | 31.8 | 29.7 |
| 0.7–1.0 | −35.4 | 20.0 | 22.8 |
| 1.0–1.5 | −9.0 | 15.2 | 15.5 |
| 1.5–2.0 | −9.5 | 13.5 | 12.0 |
| 2.0–3.0 | −16.6 | 7.0 | 7.6 |

The same matched lenses, source rows, and point profile were recomputed with
5° and 2.5° spatial blocks. The tangential contrast changes from `χ²=27.77`
(`p=2.42×10⁻⁴`, 26 blocks) to `χ²=12.05` (`p=0.099`, 54 blocks). Cross-null
values are `χ²=4.54`, `p=0.716` and `χ²=5.22`, `p=0.633`, respectively. The
profile points are identical by construction; the tangential error bars shift
by factors 0.89–1.26. The apparent 5° joint contrast therefore does not pass
this block-width sensitivity. Neither nominal probability is a calibrated
detection significance, and the radially mixed profile does not support a
simple monotonic high-ICL halo-mass interpretation.

![Exploratory ICL-split random-subtracted weak-lensing profile with 5-degree blocks](</gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_icl_split_2026-11-05_02/icl_split_profile_5deg.png>)

![The same profile with 2.5-degree blocks](</gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_icl_split_2026-11-05_03/icl_split_profile_2p5deg.png>)

## Limitations and next gates

This is a selected 70-lens overlap, reduced to 29 matched pairs. It does not
represent the 800-cluster image parent or the full 224-pair ICL-support set.
The ICL contrast uncertainty is not propagated into the split. The image-side
measurement still lacks a physical injected-source Poisson-noise correction,
validated extended-PSF wing transfer, sky-transfer and completeness
corrections. On the lensing side, equal-area blocks are not validated
contiguous HSC regions, and redshift/photo-z calibration, shear calibration,
miscentering, magnification, mass-distribution effects, and matched-random
depth/mask/PSF/sky differences remain. The `i` filter has only 38 ICL/WL common
lenses and was not included in this first calculation.

Treat this as a useful Astrafek analysis path and a sample-size diagnostic,
not as evidence that ICL strength causes or securely traces weak-lensing mass.
Every output retains `science_measurement=false`. No HSC query, download, or
image-pixel read occurred in this analysis.

## Reproducibility

- Runner: `analyze_wl_random_subtracted_icl_split.py`, SHA-256
  `705d6389c70f1a7d6dce6ff2c07bc1517ede631b0a507757b81ef39417bb4a3f`.
- Block-width audit: `compare_wl_icl_split_blockwidth.py`, SHA-256
  `e9b669dd980d83362f4a08875467e2a24a9d899c397a67f8783735af1119ab4b`.
- Plotter: `plot_wl_random_subtracted_icl_split.py`, SHA-256
  `e008ee7150a7da21a67f433b3eaf27ca236a6f56834dd65165a8c55f908f90bd`.
- Result: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_icl_split_2026-11-05_02/summary.json`, SHA-256
  `1da415abdc45720294056e7379ceff049ecb2baff66f67921a69ec8ddc1ed495`.
- The preliminary `_01` run is superseded by `_02`, which records the cross-component null statistic and matches the final runner hash above.
- 2.5° result: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_icl_split_2026-11-05_03/summary.json`, SHA-256
  `809dc2362adc8a3b345fe469abdb50339c9aeb7baf1b8b45e3859a8f43464722`.
- Block-width comparison: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_icl_split_blockwidth_2026-11-05_01/summary.json`, SHA-256
  `6d9be0508c1e89aecb4a2b7f87e69764af21a9534d4449eb9cc028ab8bf0d442`.
- Figures: `icl_split_profile_5deg.png` SHA-256
  `1566f601e065eb70fa3ee93d63a09f67c5a003686831c41980722ac3cb7cf39d` and
  `icl_split_profile_2p5deg.png` SHA-256
  `64aed42cef186cd7499aac16749be92d4332c75c97023f874362e12a89a19d3e`,
  in their corresponding result directories.
- ICL support table SHA-256: `86e6d453df847ba4c9316c4fdc188565b1bc292bb41b0c5ba8ad310864f02aa5`.
- ICL/WL overlap summary SHA-256: `98d94834b895f1275337f494a801aab7e021048ef502d4f7074635bb26b25e12`.
