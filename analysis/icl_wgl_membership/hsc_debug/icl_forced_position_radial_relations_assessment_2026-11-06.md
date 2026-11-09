# Forced unmatched-position radial sensitivity (2026-11-06)

## Question and scope

This exploratory analysis asks whether the existing 200–300 proper-kpc ICL
target-minus-control residual covaries with light at unmatched CAMIRA positions
as a function of projected cluster-centric radius. It reuses saved HSC cutout
profiles only: no network access, new download, or image-pixel read occurred in
the analysis. The profiles use the post-hoc 20/20% target/control annular
coverage floor and at least three rotated controls.

The positions do not have unique SEP detections and are not confirmed isolated
satellite galaxies. The response is therefore called **unmatched-position
light**, not satellite-galaxy luminosity or satellite outer-envelope light.
`science_measurement=false` remains in force.

## Estimator and inference

Within each cluster, filter, and radial bin, position-level target-minus-control
surface-brightness measurements are averaged with equal position weights. The
association model gives each supported cluster/filter row equal weight and
regresses that mean light against the continuous ICL residual, standardized
within filter and radial bin, adjusting for standardized redshift and log
CAMIRA richness. Standard errors use CR1 cluster-robust covariance over
connected paired target/control 5-degree sky blocks.

Only bins with at least 12 cluster rows and 6 connected blocks enter inference.
The nine eligible filter/radius tests use 10,000 restricted HC2 Webb six-point
wild-cluster bootstrap-t replicates, with shared block weights across models,
and report both Holm and single-step max-T familywise adjustment. The
familywise scope is the nine supported tests in this run; it does not include
other outcomes, radial definitions, or prior exploratory choices.

## Support and results

| Filter | Radius (kpc) | Clusters | Positions | 5° blocks | Slope per SD ICL (nJy arcsec⁻²) | CR1 95% CI | Webb p | Holm p | max-T p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| i | 20–30 | 16 | 16 | 9 | −3.893 | [−12.079, 4.293] | 0.280 | 1.000 | 0.987 |
| i | 30–40 | 18 | 20 | 9 | −0.231 | [−4.149, 3.686] | 0.862 | 1.000 | 1.000 |
| i | 40–50 | 25 | 29 | 12 | +0.065 | [−4.632, 4.762] | 0.980 | 1.000 | 1.000 |
| i | 50–75 | 34 | 46 | 12 | −3.287 | [−5.452, −1.123] | 0.030 | 0.274 | 0.059 |
| i2 | 15–20 | 14 | 14 | 12 | +9.884 | [1.555, 18.213] | 0.045 | 0.358 | 0.223 |
| i2 | 20–30 | 27 | 29 | 19 | +0.689 | [−7.974, 9.352] | 0.862 | 1.000 | 1.000 |
| i2 | 30–40 | 43 | 47 | 29 | +1.433 | [−4.324, 7.190] | 0.651 | 1.000 | 1.000 |
| i2 | 40–50 | 53 | 71 | 29 | +1.699 | [−2.940, 6.337] | 0.468 | 1.000 | 0.999 |
| i2 | 50–75 | 78 | 121 | 35 | +0.943 | [−1.201, 3.087] | 0.395 | 1.000 | 0.997 |

No max-T adjusted p-value is below 0.05. The i 50–75 kpc slope has a nominal
negative association, but max-T p=0.059 and Holm p=0.274. The i2 15–20 kpc
nominal result is also unresolved after familywise adjustment. The intervals
are broad in several bins, especially at small radius.

## Interpretation and limitations

These data do not resolve a radial ICL–satellite association. That is not
evidence that the physical association is absent: support is small in inner
bins, block counts are modest, and the forced-position sample is selected.
The earlier support audit found that retained 20/20 pairs have higher median
ICL residual than excluded pairs; no transport correction is available.

The analysis does not propagate position/profile measurement uncertainty,
source Poisson noise, exact coadd covariance, PSF wings, sky transfer, or
injection response. It does not make a correction or science promotion. The
next substantive gate is to establish position identity and unbiased radial
support, then propagate validated image-estimator uncertainty and selection
before testing this relation on a broader, spatially representative sample.

## Reproducibility

Runner: `analyze_icl_forced_position_radial_relations.py`.

Machine output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_forced_position_radial_relations_2026-11-06_01/`.

The summary records hashes for the five source artifacts and three result
tables. The report builder validates those hashes and the support/test counts
before including this analysis in the integrated evidence report.
