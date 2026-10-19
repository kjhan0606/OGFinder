# HSC published n(z) mean-shift grid assessment — 2026-10-20

## Scope and input handling

`audit_wl_published_nz_shift_sensitivity.py` ran five cluster-only sensitivity
scenarios against the already downloaded HSC Y3 SACC `n(z)`: the central
published mean shifts and four corners of the quoted marginal 68% intervals
for source bins 3 and 4. Bin 3 was translated by 0.021–0.063 and bin 4 by
0.038–0.062 in redshift. The same approved 2,153,004 lens–source rows, four
source bins, seven radial bins, and 52 sky-cell jackknife blocks were used in
all scenarios. No HSC query or data download occurred.

The immutable run summary is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_published_nz_shift_sensitivity_2026-10-20_02/summary.json`
(SHA-256 `cd16cb81c766877781424818533ab4114d3aca0a41d6ef84818c7d59736d8fdb`).
Astrafek's [`validate_wl_published_nz_shift_grid.py`](validate_wl_published_nz_shift_grid.py)
independently checked the run code and SACC/baseline hashes, all five requested
and realized shifts, each 56×56 covariance's symmetry and rank, the seven-bin
test tables, and the stored envelope. It passed with status
`validated_approximate_sensitivity_grid`; the validation JSON is next to the
summary at `validation_2026-10-20.json` (SHA-256
`1b13d25ed6b90651ac2f7eab6215bde20608bce69b95e121637dcb43cb9ec515`).

## Result

The four-bin profile inconsistency remains in all five scenarios at 0.1–0.7
Mpc. Across the central point and four interval corners, Holm-adjusted
p-values over the seven radial tests range as follows:

| Radius (Mpc) | Holm-adjusted p-value range |
|---|---:|
| 0.1–0.2 | 7.59×10⁻⁵–8.61×10⁻⁵ |
| 0.2–0.4 | 9.22×10⁻⁹–1.71×10⁻⁸ |
| 0.4–0.7 | 7.41×10⁻⁴–1.04×10⁻³ |
| 0.7–1.0 | 1.0 in all scenarios |
| 1.0–1.5 | 1.0 in all scenarios |
| 1.5–2.0 | 0.354–0.507 |
| 2.0–3.0 | 1.0 in all scenarios |

The shift grid changes the three inner-bin consistency statistics only modestly
and does not remove their discrepancy. This makes the exact published DESI
calibrated distributions worth testing against the SACC baseline; the
approximate translations do not establish that the new calibration explains
the mismatch.

## Limits and next gate

The scenarios translate and renormalize the SACC curves; they do not propagate
the DESI spline `n(z)` posterior, its covariance, or the calibration's
selection uncertainty. The analysis remains cluster-only and lacks
random-lens subtraction and covariance aligned to the connected HSC footprint.
The Holm-adjusted tests diagnose inconsistency among source-bin profiles; they
are not lens-detection significances. The grid remains
`science_measurement=false`.

The next calibration step is to obtain the specifically identified published
DESI distribution/metadata files only after explicit download approval, verify
their license, format, provenance, and checksums, then rerun the same paired
jackknife alongside SACC. A true random-subtracted profile and connected-field
covariance remain independent requirements.
