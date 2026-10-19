# ICL matched-stack tract-region influence check — 2026-10-20

## Scope and validation

This analysis tests whether the matched target-minus-control ICL profile is
stable to removing one approximate connected HSC survey region. It reuses the
existing 783 matched target/control pairs (302 `i`, 481 `i2`), the saved
per-pair annular profiles, the CAMIRA manifest, and HSC patch-QA tract centers.
No image or catalogue request was made. The result remains exploratory and
`science_measurement=false`.

The 840 HSC tracts form four connected components from spherical-mean patch
centers at a 2-degree adjacency threshold. The target sample occupies three;
the partition and region counts remain stable over 1.6–3.0 degrees. The 800
manifest targets fall in components of 499, 246, and 55 clusters. Of the 783
targets with a same-filter matched control, 487, 244, and 52 fall in those
components. Each saved broad-annulus pair value is re-averaged with equal pair
weight, after deleting one region. All twelve full-sample filter/annulus means
and pair counts reproduce `broad_radial_profiles.tsv`.

## Results

The table reports the deletion of the largest region (487 matched targets).
It shows absolute changes because relative shifts become misleading when the
outer profile is near zero.

| Filter | Annulus (kpc) | Valid pairs | Full excess (nJy/arcsec²) | Median common-pixel support | Excess after largest-region deletion | Shift |
|---|---:|---:|---:|---:|---:|---:|
| `i` | 50–100 | 235 | 41.39 | 15.0% | 40.66 | −0.74 (−1.8%) |
| `i2` | 50–100 | 378 | 44.59 | 19.0% | 40.49 | −4.10 (−9.2%) |
| `i` | 100–200 | 277 | 12.92 | 14.5% | 14.24 | +1.32 (+10.2%) |
| `i2` | 100–200 | 447 | 10.75 | 17.7% | 9.69 | −1.06 (−9.8%) |
| `i` | 200–300 | 283 | 2.08 | 11.7% | 5.46 | +3.38; middle-region deletion changes sign |
| `i2` | 200–300 | 467 | 3.18 | 15.4% | 3.79 | +0.61 |

For the `i` 200–300 kpc annulus, deleting the middle region instead changes
the mean from `+2.08` to `−0.12 nJy/arcsec²`. In `i2`, the 300–400 and
400–500 kpc means are close to zero and change sign under regional deletion.
Their large percentage shifts are not meaningful as relative error measures;
the sign changes show that the outer excess is not stable to the sampled sky
regions.

## Interpretation and limits

The 50–100 kpc mean is less sensitive to the largest region than the 100–200
kpc and outer profiles. Even in the better-supported annuli, only 15–19% of
geometric pixels are common per pair at the median. The 100–200 kpc excess
changes by about 10% when the largest region is removed, and the low-level
outer profile can change sign. Thus this spatial check does not strengthen an
ICL–satellite detection; it reinforces the existing warnings about limited
pixel support and regional dependence.

The tract-center graph is approximate rather than a tract-polygon overlap or
official HSC field map. Only three regions are occupied, so the deletions are
influence diagnostics, not a covariance estimate, confidence interval, or
significance test. Existing mask, sky, extended-PSF, and coadd-noise limits
remain. Filters stay separate, and no correction is applied to the measured
profile.

## Reproduction

- Astrafek audit: [`audit_icl_tract_region_influence.py`](audit_icl_tract_region_influence.py)
- Unit tests: [`test_audit_icl_tract_region_influence.py`](test_audit_icl_tract_region_influence.py), 2 passed
- Output JSON: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_tract_region_influence_2026-10-20_01.json`
- Output SHA-256: `4beaaec655567629a113f21c497d431328b276c1238eb96900025bd33be65ace`
- Input stack summary SHA-256: `9c804cfb12ff2b86a5bdc8fdbfea94aa1dc473711e6f7ac6f4d233cc6b853173`
- Paired profile SHA-256: `ddd19d9b633b9def4b47fb4021a5124482a8a49a18a5293d50a9b15c3dcd9187`
- Independent input-stack verification: `independent_verification_2026-10-20.json`, status `verified`, SHA-256 `0da3238fa4105118667d5e67392e09c07dc1b59db68adc4344aada464447d740`
