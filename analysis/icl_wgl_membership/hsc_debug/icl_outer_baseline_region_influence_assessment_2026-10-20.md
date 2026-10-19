# ICL outer-baseline regional influence assessment — 2026-10-20

## Result

The positive `i`-band 400–500 kpc target-minus-control baseline is present in
all three approximate HSC regions occupied by this support-selected sample.
Deleting any one region leaves the baseline positive, between 1.08 and 2.20
nJy/arcsec² depending on inner annulus. The `i2` full-sample baseline is near
zero, while leave-one-region-out values straddle zero. This rules out a single
region as the sole source of the positive `i` offset under this grouping. It
does not identify whether the offset comes from residual sky, genuine extended
cluster light, or another spatially varying effect.

| Filter | Inner annulus (kpc) | Pairs | Full outer baseline | Leave-one-region-out outer baseline | Full outer-subtracted profile | Leave-one-region-out outer-subtracted profile |
|---|---:|---:|---:|---:|---:|---:|
| i | 50–100 | 157 | 1.658 | 1.221–1.801 | 38.544 | 35.368–40.320 |
| i | 100–200 | 177 | 1.668 | 1.120–1.956 | 10.133 | 9.843–10.412 |
| i | 200–300 | 189 | 1.783 | 1.079–2.196 | 2.530 | 1.862–3.176 |
| i2 | 50–100 | 279 | −0.008 | −0.595–0.759 | 46.098 | 44.081–47.208 |
| i2 | 100–200 | 338 | −0.163 | −0.464–0.378 | 11.251 | 9.629–11.744 |
| i2 | 200–300 | 368 | −0.007 | −0.221–0.638 | 3.168 | 2.715–3.805 |

All surface-brightness values are nJy/arcsec². Means are equal-weighted over
the listed target-control pairs. The outer-subtracted quantity is the inner
pair difference minus that same pair's 400–500 kpc difference.

## Method and limits

[`audit_icl_outer_baseline_regions.py`](audit_icl_outer_baseline_regions.py)
re-reads the saved outer-annulus result and its source stack, verifies their
SHA-256 checksums and target/control mapping, then assigns each target to the
nearest component built from HSC patch-QA tract centers with 2° adjacency. Of
840 patch-QA tracts, this approximate graph has four components; three contain
the analyzed targets. Labels use target coordinates only. The matched-control
coordinates are not used to assign regions.

The same per-annulus rows selected by the 5% common-support floor are held fixed
for the full mean and each leave-one-region-out calculation. This is an
influence diagnostic over three occupied approximate regions, not a stable
regional covariance or confidence interval. Tract-center components are not
official HSC footprint polygons. The selection differs between annuli, and
400–500 kpc can contain real ICL, so subtracting it may remove astrophysical
signal. This analysis estimates neither pixel covariance nor a science
detection; `science_measurement=false`.

The verified JSON is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_outer_baseline_region_influence_2026-10-20_01.json`
(SHA-256 `ec941eb765827de518d8428b4d26256cac30d5756cdaa445b0b73110cb789aac`).
Its input paths and hashes are recorded in the JSON. No HSC query or download
was made.
