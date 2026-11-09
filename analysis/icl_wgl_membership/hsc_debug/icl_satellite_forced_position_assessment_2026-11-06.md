# Forced-position satellite profile assessment — 2026-11-06

## Purpose and input scope

This offline analysis asked whether the existing 800 HSC-I target cutouts could
extend the CAMIRA satellite-position stack to member positions missed by the
saved SEP detection catalog. It did not query or download HSC data. The
catalog-only inventory found 3,197 positions at 50–500 proper kpc with no
unique SEP match within 1 arcsec and no SEP detection within 3 arcsec. The
measurement then used those existing coadds and preparation masks.

Runner: `hsc_debug/measure_forced_satellite_profiles.py`.
Support audit: `hsc_debug/audit_satellite_forced_position_support.py`.
Machine output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_forced_outer_profile_2026-11-06_01/`.
The support inventory is in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_forced_position_support_2026-11-06_03/`.

## Measurement and support result

The pipeline forced the center at each CAMIRA coordinate, masked saved source
detections, HSC quality/bright-object bits, and other unmatched CAMIRA
positions, then compared the profile with up to eight rotated controls at the
same cluster-centric radius. A radial bin was retained only when target
coverage was at least 60% and at least three controls independently met the
same coverage floor.

Only 14/3,197 candidates (0.44%) passed in any radial bin, and only 9/3,197
(0.28%) contributed to 50–75 kpc. There were 20 retained candidate-bin
measurements across 721 processed target frames. 3,183 candidates passed no
annulus. At 50–75 kpc, target coverage alone reached 60% for only 8/1,292 `i`
positions and 10/1,905 `i2` positions; only four and five, respectively,
passed the additional control-support rule. Most candidates had eight valid
rotated control profiles, so the dominant bottleneck was target-annulus
coverage, not the number of available control placements. The forced-position
stack therefore has too little support to repair the detected-member
selection problem.

The 50–75 kpc forced-only ICL-association sample contains 2 pairs in `i` and
4 in `i2`, both below the analysis minimum; no forced-only association estimate
is reported. Combining the few retained forced profiles with the saved
detection-selected profiles leaves 26 `i` and 68 `i2` target/control pairs.
The redshift/log-richness-adjusted, connected-5-degree-block exploratory
partial correlations are `r=0.131` (95% interval −0.211 to 0.556) for `i` and
`r=0.287` (0.037 to 0.490) for `i2`. The latter remains a selected-sample
exploratory interval; it is not a science detection, does not correct
completeness, and has no familywise adjustment across the broader analysis.
The change from the earlier detected-only `i2` estimate (`r=0.293`, interval
0.043 to 0.498) is small compared with its uncertainty.

## Interpretation and limitations

These are forced measurements around catalog positions, not verified
individual satellite galaxies. The 3,197-position parent set is still
CAMIRA-selected; unmatched positions can include blends, contaminated members,
or positions with no recoverable galaxy signal. Strongly position-dependent
coverage makes the 14 retained profiles a highly selected subset. This
pipeline does not estimate the missed-position completeness function and does
not justify inverse-probability weighting.

The existing masks, BCG handling, PSF wings, residual-sky transfer,
source-Poisson noise, and total coadd covariance are not science-certified.
The connected 5-degree block resampling is approximate. All outputs retain
`science_measurement=false`; the ICL quantity remains BCG plus diffuse light,
not a separated pure ICL component.

## Reproducibility

Analysis summary reports 3,197 candidates, 721 processed frames, 20 retained
bins, 14 candidates with any retained bin, and 9 candidates at 50–75 kpc. Key
output SHA-256 values:

- `forced_satellite_profiles.tsv`: `d211d27bfafbbf15a1d0bde06c3813ebd43d50d097af91be3723d31a05abe3aa`
- `candidate_diagnostics.tsv`: `73eded921f44966f7860a810aa4ea68c906458004ba8557c5336d23404f49992`
- `cluster_profiles.tsv`: `a770d8affeba5e63eda68036ab02dceade24fc2fddb99e5182d8d192cb6c7e8b`
- `icl_satellite_association.tsv`: `c28c30881e1781e7546a0fd5a0f46d9a2d2d46c590b7fffe5cb1ea790dde30d5`

The machine-readable summary records all input hashes and the runner hash.
