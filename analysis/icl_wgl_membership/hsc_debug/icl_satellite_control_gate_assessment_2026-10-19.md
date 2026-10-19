# ICL–satellite control-gate attrition audit — 2026-10-19

## Scope and reproducibility

This audit explains why 398 satellite positions that passed the target-side
coverage gate were rejected by the complete seven-annulus estimator for having
fewer than three valid azimuthal controls. It uses only the checksum-verified
candidate diagnostics from the existing full-profile and per-annulus runs. No
images, pixels, queries, or downloads were added. The audit extends
[`audit_icl_satellite_sample_retention.py`](audit_icl_satellite_sample_retention.py)
and preserves the original complete-profile selection as the reference.

The machine-readable output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_sample_retention_2026-10-19_02/summary.json`
(SHA-256 `81e36f8d0f1b23aac4cd4f6ad614bca5f1910b94475354349dcf66ed0d5df052`).
Both source diagnostic tables matched the hashes recorded in their run
summaries. The 50–75 kpc candidate set was independently cross-matched against
the saved per-annulus satellite profiles.

## Why the complete profiles fail

The 398 failures comprise 118 i-band and 280 i2-band positions. Candidate
generation supplied eight azimuthal positions for 396 of them; one had seven
and one had six. The loss therefore comes mainly after candidates are generated:
each control must have at least 60% usable area in **all seven** 5–75 kpc
annuli to count as one complete profile.

| Valid controls across all seven annuli | Failed positions |
|---:|---:|
| 0 | 186 |
| 1 | 128 |
| 2 | 84 |

The number of failures with at least three controls in each individual annulus
falls from 263 in 5–10 kpc to 249, 226, 176, 153, 134, and 102 in successive
annuli through 50–75 kpc. This points to control-profile coverage and masking
across the full radial range, rather than a general inability to generate
azimuthal positions. The diagnostics do not retain a per-position causal
breakdown for every rejected control (edge, mask, or source clearance), so the
relative contribution of those mechanisms is not identified here.

## Primary outer-annulus recovery

The 50–75 kpc bin has at least three valid controls for 25 i-band and 77 i2-band
members among the 398 complete-profile failures (102 total). All 102 appear in
the saved independent-annulus profiles, with exact member/filter/detection-key
agreement. This reproduces the intended behavior of the radial-bin estimator:
a satellite can contribute to the prespecified outer bin even when it lacks a
complete seven-bin profile. It does not add profiles by relaxing the coverage
or control threshold, and it does not make the full seven-bin profile complete.

The separate full-sample per-annulus stack contains 70 i and 308 i2 satellites
in 50–75 kpc; those totals also include candidates that failed target coverage
in one or more other annuli and therefore were never in the 398 complete-profile
control-failure set. The independent-bin estimator remains the appropriate
reference for the outer-satellite association, with its bin-specific sample
selection stated explicitly.

This is a sample-support audit, not a new flux measurement or a science
detection. `science_measurement=false` remains in force. Tests are in
[`test_audit_icl_satellite_sample_retention.py`](test_audit_icl_satellite_sample_retention.py).
