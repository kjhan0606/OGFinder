# ICL–satellite sensitivity to strict HSC warp support

## Question and scope

Does requiring at least three valid HSC visit warps and shared mask-plane support in the 200–300 proper-kpc ICL aperture materially change the measured ICL–satellite association?

This offline test uses the already approved/downloaded 94-pair, 188-center HSC-I warp plan (plan SHA-256 `bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a`). No new HSC query or download was made. It keeps the original saved K=2.5 source masks and the same 40% target/control annulus coverage requirement. Only the ICL aperture's valid-pixel support changes; the satellite 50–75 kpc values remain the saved measurements.

## Exact baseline reproduction

Before applying warp support, the runner remeasures each target and matched-control 200–300 kpc aperture from the coadd image and K=2.5 masks. Every tested target/control surface brightness and coverage reproduces the saved annulus-profile table within 2×10⁻⁶ nJy arcsec⁻² / coverage fraction. This pins the comparison to the existing Astrafek measurement.

## Results

- Four approved center requests were repeat HTTP 404, so 90 pairs had both warp cutouts: 24 `i`, 66 `i2`.
- Applying the strict warp-valid mask and then retaining the original 40% coverage floor leaves 81 pairs: 23 `i`, 58 `i2`. The additional nine-pair loss is one `i` and eight `i2` pairs. Median strict target/control annulus coverage is 56.0%/57.3% (`i`) and 58.2%/58.9% (`i2`).
- On the same retained pairs, redshift/richness-adjusted partial correlation changes from 0.0779 to 0.0799 in `i` (paired change +0.00203; connected-block 95% interval −0.00168 to +0.01487).
- In `i2`, it changes from 0.3200 to 0.3153 (paired change −0.00478; interval −0.01293 to +0.00362).
- Both paired intervals include zero. The apparent difference between `i2` baseline r=0.320 on this retained subset and r=0.293 on the original 68 pairs is sample attrition; it is not a mask-only effect. The mask-only comparison uses the same 58 pairs on both sides.

The bootstrap resamples connected target/control 5-degree spatial blocks recomputed on each strict retained filter sample (12 blocks for `i`, 32 for `i2`), with 5,000 replicates and fixed per-filter seeds.

## Method and limitations

For each returned visit, the runner checks the archive hash and required visit inventory, fits a local affine warp-to-coadd WCS map on a 7×7 grid, rejects maps with maximum residual above 0.05 pixel, and nearest-neighbor samples mask/validity onto coadd pixels. A strict ICL pixel must be valid in the coadd, have at least three finite positive-variance visits, and pass the shared-mask criteria. The original K=2.5 source masks are then applied.

Nearest-neighbor mask sampling is not the exact HSC resampling kernel. The 4 HTTP 404 center gaps and 9 further coverage losses make the retained set selected; this sensitivity cannot establish the full 94-pair or 800-parent relation. It is not a source-Poisson, PSF-wing, sky-transfer, or full covariance correction. Keep `science_measurement=false`.

## Reproducibility

- Runner: `analysis/icl_wgl_membership/hsc_debug/analyze_icl_supported94_warp_mask_association.py`
- Runner SHA-256: `b40eff47b5d2bf2449f32b0141d4cf570c97278698f8ac616ce72809c8d50519`
- Pair table: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported94_warp_mask_association_2026-11-06_02/strict_mask_association_pairs.tsv`
- Summary: same directory, `summary.json`
- Summary SHA-256: `f7d62ae88864157636d9d7be64ad09d4907060a2cefd7d9ce9b57faf09ef69ad`
