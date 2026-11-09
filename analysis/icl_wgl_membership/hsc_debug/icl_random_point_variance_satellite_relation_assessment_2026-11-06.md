# Random-point local conditions and the ICL–satellite relation (2026-11-06)

## Question and data

Does the observed association between the target-minus-control ICL
200–300 proper-kpc contrast and satellite 50–75 kpc excess change after
accounting for local conditions measured from HSC random points in the same
annulus? The fixed supported sample contains 26 `i` and 68 `i2` pairs. Every
pair has positive target and control random-point medians for `i_pix_variance`,
`i_sky_std`, `i_inputcount_value`, and the local PSF moment trace.

The analysis uses the already approved and audited random-point export. It
reads no HSC endpoint and no image pixels. For each pair it forms target/control
log-ratio and pair-mean log values for those four local metrics. Each is added
one at a time to the baseline partial-correlation model controlling redshift
and CAMIRA richness. A shared 5,000-draw connected 5-degree block bootstrap
uses a single max-T family across all eight adjustments in both filters (16
comparisons).

## Results

The baseline partial correlations reproduce the existing association:
`r=0.079` in `i` and `r=0.293` in `i2`. In `i2`, adding local variance
log-ratio changes r by −0.011 (max-T 95% interval −0.176 to +0.155); adding
mean log variance changes r by −0.033 (−0.225 to +0.160). Sky-RMS, input-count,
and PSF-trace adjustments likewise leave r between 0.260 and 0.295. In `i`,
adjusted estimates range from −0.021 to +0.130 with broad intervals. Across all
16 tests, the minimum max-T adjusted p is 0.951, and every simultaneous interval
for the change in r includes zero.

The local random-point metrics do not produce a resolved change in this
selected-sample ICL–satellite association. This is not evidence that local noise
or image quality has no effect: the `i` sample is small, random-point medians
vary in precision with 11–39 points per center (q16–q84), and the 94-pair
support sample is not representative of the 800-cluster parent. The test does
not propagate random-point median uncertainty and does not measure pixel
covariance, source Poisson noise, PSF wings, or sky transfer. No flux correction
or science promotion follows; `science_measurement=false`.

## Reproducibility

- Runner: `analyze_icl_random_point_variance_satellite_relation.py`, SHA-256
  `738703bc9cadd561595440eaa9e395a783d1b67ffcabfb4c68e0515085a4b826`.
- Output directory:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_point_variance_satellite_relation_2026-11-06_01/`.
- Summary SHA-256:
  `0c75520ad6b04dd500ba66ae10982e30d6fa9758379615a4372e7aa129569f0f`.
- Joined 94-pair table SHA-256:
  `e9d84e1ee6e78d483685893ade99849850019e1fd66713907b0a627d064207cf`.
- Pinned ICL–satellite parent SHA-256:
  `79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da`.
- Random-point audit summary SHA-256:
  `96245d565375635a994c93ec3294fc581d007161c4479232d5091bbd3cf52389`.
- HSC network used: false; image pixels read: false; data requested: false.
