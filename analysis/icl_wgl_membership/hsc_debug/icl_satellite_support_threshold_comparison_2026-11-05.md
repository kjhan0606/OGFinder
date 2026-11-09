# ICL–satellite support threshold comparison

## Result

Compared the completed 40% satellite-annulus-coverage run with the earlier 60% run on the same 57 PSF-supported frames (13 `i`, 44 `i2`; 57 cluster IDs). Per-frame image hashes, all position-dependent PSF hashes, and the re-detection outputs match. The diffuse profile, injected amplitudes, radial bins, and redetection configuration also match. The changed selection is the satellite support floor and the control-profile requirement.

| Filter | Median satellites/frame, 60% → 40% | Median controls/frame, 60% → 40% | Satellite unit-injection response, median 60% → 40% | Paired median change (40%−60%), cluster bootstrap 95% interval |
|---|---:|---:|---:|---:|
| i | 2 → 5 | 7 → 26 | 0.338 → −0.144 | −0.033 [−0.768, 0.487] |
| i2 | 2 → 6 | 9 → 34.5 | 0.409 → 0.572 | +0.098 [−0.264, 0.526] |
| Combined | 2 → 6 | 8 → 33 | 0.373 → 0.464 | +0.009 [−0.250, 0.416] |

The response intervals include zero in both filters. At the 50–75 kpc profile measurement, the `both` injection scenario has a combined median satellite-minus-control change of 1.05 nJy/arcsec² at 60% and 1.58 nJy/arcsec² at 40%; the paired median change is +0.11 nJy/arcsec² (95% cluster-bootstrap interval −0.73 to +0.51). This remains broad and does not establish a recovered satellite signal. The diffuse-unit response in the ICL 200–300 kpc annulus is exactly unchanged between runs (median 0.215 in `i`, 0.225 in `i2`), as expected because that estimator is unaffected by the satellite profile support threshold.

## Interpretation and limits

Lowering the coverage floor increases the median number of usable satellite profiles by roughly three per frame and provides substantially more azimuthal controls. Detection and masking were re-run and are identical between paired runs, so this yield increase comes from profile support eligibility rather than a changed source catalog. The satellite injection-transfer response remains noisy and compatible with no difference between thresholds. The small `i` subset (13 frames) is especially uncertain.

This is a support/estimator sensitivity check on the existing 57-frame PSF subset. It does not estimate population completeness, correct radial selection, calibrate coadd covariance, or constitute an astrophysical ICL–satellite detection. Injection metadata mark `science_measurement=false`; the artificial profiles do not alter the native FITS or variance maps.

## Reproducibility

- Runner: [`compare_icl_satellite_support_thresholds.py`](compare_icl_satellite_support_thresholds.py), SHA-256 `55b8a89b40fdc9191412f30602d48218015d63c73157431aada0c114d84f66d6`.
- Paired result directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_support_threshold_comparison_2026-11-05_02`.
- `paired_frames.tsv` SHA-256 `c86368c4dfb8aa4811f1d32bf9f133b341c2d14a3189a31fdc0baab084ca0add`.
- `summary.json` SHA-256 `c2722311dc4c5dbfaf6592751a6ae58358152baa0376e27a3cb88f8e2529ffe2`.
- `assessment.md` SHA-256 `bae8b17c56961b608c6658a88ff8b1ad32194519b92f73ce2326d2cf39310697`.

The summary uses 10,000 paired cluster-ID bootstrap resamples with seed `20261105`. The intervals quantify paired sample uncertainty for this diagnostic, not uncertainty in a population correction.
