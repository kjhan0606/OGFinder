# Extended-scale paired-field ICL plane-injection audit — 2026-10-20

## Purpose and method

This extends Astrafek's paired target/control test of the downstream local-plane
estimator to a 500 kpc exponential stress profile. It reuses the fixed 74-pair
support sample (i: 21; i2: 53), with 100, 250, and 500 kpc exponential scales
and axis ratios 1.0 and 0.65. The injected surface brightness is 1 nJy/arcsec²
at elliptical radius 100 kpc; the position angle is fixed at 37 degrees.

For each field, the calculation measures the 200–300 kpc annulus before and
after rerunning Astrafek's 600–900 kpc tile-plane estimate. The same K2.5 mask
and annular pixels are used. Target and matched-control responses are paired
through the connected 5-degree sky blocks; intervals use 5,000 identical block
bootstrap draws. Control masks are derived from the unmodified control image
and then held fixed. The test is downstream of HSC's delivered global-sky
coadd/bg product: it does not rerun HSC's upstream sky subtraction.

No PSF convolution, source Poisson noise, mask re-detection, or change to the
delivered diagonal variance is included. The 500 kpc profile is a stress case,
not a claim that it describes the true ICL profile.

## Results

Entries are mean transfer fractions after local-plane subtraction. The paired
difference is target minus matched control; parentheses show its 95% connected
block-bootstrap interval. The raw-annulus response is unity to numerical
precision for all models.

| Filter | Exponential scale | q | Pairs | Target recovery | Control recovery | Paired difference (95% interval) |
|---|---:|---:|---:|---:|---:|---:|
| i | 100 kpc | 1.00 | 21 | 0.9911 | 0.9912 | −0.0001 (−0.00065, 0.00068) |
| i | 100 kpc | 0.65 | 21 | 0.9951 | 0.9948 | +0.0004 (−0.00056, 0.00149) |
| i | 250 kpc | 1.00 | 21 | 0.8596 | 0.8600 | −0.0004 (−0.00456, 0.00467) |
| i | 250 kpc | 0.65 | 21 | 0.9058 | 0.9020 | +0.0038 (−0.00476, 0.01383) |
| i | 500 kpc | 1.00 | 21 | 0.6303 | 0.6309 | −0.0006 (−0.00599, 0.00639) |
| i | 500 kpc | 0.65 | 21 | 0.7084 | 0.7023 | +0.0061 (−0.00765, 0.02239) |
| i2 | 100 kpc | 1.00 | 53 | 0.9911 | 0.9911 | +0.0001 (−0.00019, 0.00039) |
| i2 | 100 kpc | 0.65 | 53 | 0.9954 | 0.9950 | +0.0004 (−0.00012, 0.00100) |
| i2 | 250 kpc | 1.00 | 53 | 0.8598 | 0.8593 | +0.0005 (−0.00123, 0.00269) |
| i2 | 250 kpc | 0.65 | 53 | 0.9077 | 0.9055 | +0.0022 (−0.00228, 0.00771) |
| i2 | 500 kpc | 1.00 | 53 | 0.6305 | 0.6299 | +0.0006 (−0.00164, 0.00355) |
| i2 | 500 kpc | 0.65 | 53 | 0.7110 | 0.7084 | +0.0026 (−0.00518, 0.01143) |

For the 500 kpc profiles, the plane estimator retains about 63% of the
circular exponential signal and 70–71% of the flatter profile in both target
and control fields. Thus the estimator removes roughly 29–37% of this broad
injected light. Every paired target-minus-control interval includes zero; the
test does not show a measurable differential transfer between these matched
fields for the assumed profiles. It also cannot rule out differential bias at
the precision available from 13 i and 28 connected bootstrap blocks, or 28
i2 blocks.

This scale dependence is a methodological result: a local-plane estimator can
strongly attenuate broad diffuse components even when matched target/control
transfer is similar. It does not supply a universal correction factor for
observed ICL, because the true spatial profile, HSC upstream processing,
mask/redetection response, PSF wings, and full noise covariance are not
modeled. In particular, this calculation cannot be used to correct the
observed ICL–satellite relation.

## Reproducibility and status

The Astrafek implementation is `audit_icl_plane_injection.py`, with the
extended run wrapper `run_icl_plane_injection_extended.sbatch`. The
shape-grouping regression tests in `test_audit_icl_plane_injection.py` passed
before execution. Slurm job 428956 completed with exit 0 in 6 min 23 sec; peak
RSS was 608,268 KiB. The output contains 74 pair-level summaries and 444
filter × shape × pair transfer rows. The result is in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_plane_injection_extended_2026-10-20_01/`.

- `summary.json` SHA-256: `387f1236346a85d15e43bc296fedf4d4bbc8a880756093576fce98755f367dca`
- `transfer_rows.tsv` SHA-256: `e605fec6de9f687d27b6e5a71842e54ba0a3fb875c284da1568bea8143abdf25`
- `audit_icl_plane_injection.py` SHA-256: `b1c9426ff7d7b5667dbb5615d5d8010abd6e5d4c5c5f751d26072c12c72bdd28`

An independent row-level recomputation reproduced the stored mean responses
and paired block-bootstrap intervals from the frozen input rows. The batch
remains exploratory: `science_measurement=false` and
`observational_validation_passed=false`.

## Next ICL validation gate

The next Astrafek run should carry these injected profiles through the full
source-detection and mask-regeneration path, with target and control fields
treated symmetrically. Add source-noise and correlated-coadd-noise treatment
only after the image metadata and empirical covariance model support it. Keep
filter-specific results, failures, and injection recovery visible, and do not
turn the measured plane-transfer fraction into an observed-data correction.
