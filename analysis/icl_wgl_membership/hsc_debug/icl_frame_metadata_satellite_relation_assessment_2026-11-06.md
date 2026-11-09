# Frame metadata sensitivity of the ICL--satellite relation (2026-11-06)

## Question and scope

This offline analysis asks whether the observed ICL--satellite partial
correlation changes when coarse target/control CCD-frame metadata offsets are
added to the adjustment set. It joins the already audited 67,071-row HSC frame
metadata export to the fixed 94-pair support table (26 `i`, 68 `i2`). The
baseline estimator matches the existing analysis: partial Pearson
correlation between the ICL 200--300 kpc target-minus-control contrast and
satellite 50--75 kpc excess, controlling linearly for cluster redshift and
CAMIRA richness.

Five metadata offsets are added one at a time: log target/control unique-frame
count, unique-visit count, a median sky-RMS/FLUXMAG0 proxy, median FLUXMAG0,
and geometric-mean amplifier gain. A shared set of 5,000 connected 5-degree
spatial-block bootstrap resamples estimates each paired change from baseline.
A single-step max-T adjustment covers all five changes in both filters.

## Results

| Filter | Baseline partial r | Metadata adjustment | Adjusted r | Change from baseline (pointwise 95% interval) | max-T adjusted p |
|---|---:|---|---:|---:|---:|
| `i` | 0.079 | frame count | 0.151 | +0.073 [−0.025, +0.265] | 0.895 |
| `i` | 0.079 | visit count | 0.086 | +0.007 [−0.070, +0.070] | 1.000 |
| `i` | 0.079 | sky-RMS proxy | 0.085 | +0.006 [−0.128, +0.204] | 1.000 |
| `i` | 0.079 | FLUXMAG0 | 0.086 | +0.007 [−0.155, +0.280] | 1.000 |
| `i` | 0.079 | gain geometric mean | 0.046 | −0.033 [−0.199, +0.035] | 0.990 |
| `i2` | 0.293 | frame count | 0.286 | −0.007 [−0.059, +0.033] | 1.000 |
| `i2` | 0.293 | visit count | 0.309 | +0.017 [−0.015, +0.092] | 0.975 |
| `i2` | 0.293 | sky-RMS proxy | 0.290 | −0.003 [−0.051, +0.006] | 1.000 |
| `i2` | 0.293 | FLUXMAG0 | 0.289 | −0.004 [−0.048, +0.036] | 1.000 |
| `i2` | 0.293 | gain geometric mean | 0.300 | +0.007 [−0.046, +0.039] | 1.000 |

All familywise adjusted p-values are at least 0.895, and every simultaneous
interval includes zero. The `i2` association is nearly unchanged by these
coarse frame-summary adjustments. The `i` estimate moves more under frame
count and gain adjustments, but its paired intervals are broad. This analysis
does not identify a resolved metadata explanation for the association and
does not establish that the `i2` relation is real.

## Limits

The metadata aggregate pools contributing frames at the patch level and does
not identify which exposure or amplifier contributed to each measured pixel.
The converted sky-RMS is only a proxy; none of these predictors is a calibrated
coadd-noise estimate. The 26/68 pair sample is selected and small, especially
in `i`. The 5-degree endpoint blocks are approximate spatial regions. The
analysis does not modify the ICL profiles or correct their uncertainties.
Source Poisson, pixel coadd weights, sky transfer, extended PSF wings, and
selection uncertainty remain unresolved; `science_measurement=false` remains.

## Reproducibility

- Astrafek runner: `analyze_icl_frame_metadata_satellite_relation.py`, SHA-256
  `cf2aa6cfc75bc8f0b9e6db5ac694d52315bfbbe3b4e97d6e402d3eefca0d07d3`.
- Result directory:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_satellite_relation_2026-11-06_03/`.
- Summary SHA-256:
  `cd2b17e5b9e3579ae580b0e7531b54f5ee6e9508477f77d1b3b78bb7c9a02643`.
- Joined-pair table SHA-256:
  `368137cae49ff728f97fc97562721dd3cbcaaa5765aa83d32c65d664ccec5950`.
- Inputs are the fixed supported-pair table and the independently re-audited
  frame-metadata aggregates from plan SHA-256
  `49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3`.
- An initial output directory ending `_01` used log richness in the baseline
  and is superseded. `_02` has the correct estimator but its summary omits the
  runner digest; use `_03` as the final reproducible output. Earlier
  directories are retained only as provenance records.
- The run used no HSC network request and read no image pixels.
