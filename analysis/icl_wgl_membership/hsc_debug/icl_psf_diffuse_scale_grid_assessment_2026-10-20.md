# ICL diffuse-profile transfer across scale length — 2026-10-20

## Result

Astrafek extended the existing 100/500 kpc injection comparison with 200 and
300 kpc exponential scale lengths. Each run used the same 57 approved
cluster/filter frames (13 `i`, 44 `i2`), the local HSC-I model-PSF core at the
injection positions, and fresh SEP detection/K2.5 source masks for every
scenario. The 200 and 300 kpc batches each completed all 57 frames with zero
errors; the independent PSF-to-position audits verified all frame mappings.

The table reports the median diffuse ICL response in the 200–300 kpc annulus
before and after Astrafek's downstream local-plane estimator, along with the
per-frame corrected/raw retention ratio. Brackets are the 16th–84th percentiles
across frames, not confidence intervals.

| Filter | Exponential scale length (kpc) | Raw response | Plane-corrected response | Retention ratio |
|---|---:|---:|---:|---:|
| i | 100 | 0.215 [0.170, 0.232] | 0.206 [0.163, 0.239] | 1.002 [0.954, 1.066] |
| i | 200 | 0.456 [0.406, 0.485] | 0.420 [0.367, 0.443] | 0.913 [0.888, 0.946] |
| i | 300 | 0.599 [0.524, 0.618] | 0.478 [0.412, 0.515] | 0.816 [0.777, 0.840] |
| i | 500 | 0.746 [0.682, 0.758] | 0.476 [0.411, 0.498] | 0.644 [0.611, 0.679] |
| i2 | 100 | 0.224 [0.183, 0.252] | 0.227 [0.174, 0.256] | 0.997 [0.957, 1.055] |
| i2 | 200 | 0.481 [0.436, 0.509] | 0.437 [0.390, 0.477] | 0.918 [0.896, 0.942] |
| i2 | 300 | 0.604 [0.571, 0.644] | 0.490 [0.445, 0.536] | 0.812 [0.783, 0.829] |
| i2 | 500 | 0.739 [0.710, 0.770] | 0.465 [0.423, 0.510] | 0.633 [0.588, 0.658] |

Across these fixed profiles, the median retention declines smoothly with
diffuse scale length in both filters. The paired median retention change from
the 100 kpc profile was −0.064 and −0.075 at 200 kpc, −0.170 and −0.182 at
300 kpc, and −0.341 and −0.362 at 500 kpc for `i` and `i2`, respectively.
These paired changes are descriptive; their q16–q84 ranges describe frame
scatter and are not inferential intervals.

## Interpretation and limits

The result shows that the downstream plane estimator attenuates broader
diffuse profiles even after local model-PSF convolution and mask redetection.
For the injected 200 and 300 kpc scale lengths, median retention is about
0.81–0.92; the 500 kpc stress profile retains about 0.63–0.64. These values
are conditional on the selected 57 frames, the synthetic exponential profiles,
and the current masks. They are not universal transfer functions and must not
be applied as multiplicative corrections to the observed ICL stack.

The test does not model HSC upstream sky processing, extended PSF wings,
injected-source Poisson noise, or correlated coadd noise. The 300 and 500 kpc
`i2` runs each lack the diffuse-to-satellite cross-response for one frame whose
satellite aperture lost support after mask redetection; the diffuse ICL
response itself remains valid for all 44 `i2` frames. No HSC query or download
was made. `science_measurement=false` remains in force.

The next ICL gate remains a validated PSF-wing and sky-transfer model applied
symmetrically to target and matched-control data, followed by joint covariance
estimation with spatial blocks. Until that gate is met, the observed ICL and
satellite association remains exploratory.

## Reproduction and provenance

The new CPU-only runner is
[`run_icl_psf_scale_bracketing_cpu.sh`](run_icl_psf_scale_bracketing_cpu.sh);
it runs on `lageunha` and refuses to overwrite existing outputs. The paired
scale-grid summarizer is
[`summarize_icl_psf_diffuse_scale_grid.py`](summarize_icl_psf_diffuse_scale_grid.py),
with three focused tests in
[`test_summarize_icl_psf_diffuse_scale_grid.py`](test_summarize_icl_psf_diffuse_scale_grid.py).

The two new outputs and their independent PSF/redetection audits are outside
Git:

| Scale | Batch summary SHA-256 | Independent audit SHA-256 |
|---:|---|---|
| 200 kpc | `1407fbb95ecc89114a948501cb597fcef0b63933a76b09879c3f560a741595d0` | `bfe4b2afee000c59fb0165fe01c8f4af168c179eaf8f8476c2f95ee314bc6483` |
| 300 kpc | `c060ed71d9542815a44984728cc1589e92cd776a8bb2b51c3355b3c7f52da84f` | `d0747064e8d55e35940cde150088ac89e25a6bbd4a2a0fdf5b1a355fbb72cd97` |

The combined 100/200/300/500 kpc summary is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_diffuse_scale_grid_2026-10-20_01/summary.json`
(SHA-256 `5cdde4204365e696c92e6c83bd2463911bdd3d71578c9cc75200c03e6653ccb6`).
It verifies identical frame support, runtime code hashes, and scientific inputs
across all four scale runs. The original 100/500 kpc analysis remains in
[`icl_psf_redetection_scale_comparison_assessment_2026-10-20.md`](icl_psf_redetection_scale_comparison_assessment_2026-10-20.md).
