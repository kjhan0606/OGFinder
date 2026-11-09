# Host-centric satellite-position mix sensitivity for the ICL relation

## Sample and method

The analysis uses the fixed 94-pair common-support sample, with 26 `i` pairs
and 68 `i2` pairs. Each stack contains detected members with a measured
satellite-centric 50–75 proper-kpc profile. The analysis recomputed the
contributing-member count and the median host-centric proper radius from the
individual member rows, then verified the member count against the saved
cluster-level profile table.

The baseline partial Pearson correlation controls lens redshift and linear
CAMIRA richness. The sensitivity model also controls the log of the median
host-centric radius of the contributing members. Uncertainty uses a joint
10,000-draw resampling of connected 5-degree target/control sky blocks across
both filters. A total of 9,989 draws produced finite estimates in both filters.
The intervals are pointwise and exploratory.

## Results

Across clusters, the distribution of per-cluster median host-centric radii has
median 311.7 kpc in `i` (16th–84th percentiles 199.4–405.3 kpc) and 326.8 kpc
in `i2` (264.7–394.5 kpc).

For `i2`, the baseline ICL–satellite partial correlation is 0.293. Adding the
position-mix covariate raises the estimate to 0.322. The paired change is
0.030 with a 95% block-bootstrap interval of −0.019 to 0.101. The partial
correlation of log median radius with ICL excess is 0.205 (−0.063, 0.432),
and with satellite excess is −0.103 (−0.348, 0.078).

For `i`, the baseline partial correlation is 0.079. Adding the position-mix
covariate lowers the estimate to 0.029. The paired change is −0.050 with a
95% interval of −0.330 to 0.044. The partial correlations of log median
radius with ICL excess and satellite excess are 0.153 (−0.214, 0.659) and
0.338 (−0.057, 0.547), respectively. All position-mix association and
model-change intervals include zero.

The measured host-centric position mix does not resolve the ICL–satellite
association in either filter. The `i2` point estimate rises modestly, while
the paired interval permits either sign. The `i` estimate shifts downward
with a broad interval. The position summary is noisy for stacks with few
members and does not recover undetected satellites. The analysis remains a
selected-sample sensitivity, not a population or causal correction. Source
Poisson noise, PSF wings, sky transfer, total image covariance, and selection
completeness remain unresolved. Keep `science_measurement=false`.

## Reproducibility

- Astrafek runner: [`analyze_icl_satellite_clustercentric_position_adjustment.py`](analyze_icl_satellite_clustercentric_position_adjustment.py).
- Runner SHA-256: `65e8557516e6df8633a85e75ad519bf343f997f52a5e223424f7e5e4f9e55f86`.
- Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_clustercentric_position_adjustment_2026-11-06_01/`.
- Summary SHA-256: `78317905c3d0b2b0ada48c8b39a27ae564caaeb0bd744f365c12fedc0a7d8971`.
- Joined-pair table SHA-256: `9966cbb4eb978e0944f6aa164f52c747bd5148e153e96427c88700056f679195`.
- Filter summary SHA-256: `bdcdff169313f10af2cb05418982b26ddc3b05e19b0860b90e75f9bb3eba501e`.
- Existing local profile tables and coordinate manifests were used. No HSC request, download, or image-pixel read occurred.
