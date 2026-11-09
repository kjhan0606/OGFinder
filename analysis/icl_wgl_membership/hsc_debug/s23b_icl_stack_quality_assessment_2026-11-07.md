# S23B ICL target/control stack quality assessment — 2026-11-07

## Status

The approved v2 matched-control download, independent FITS/WCS audit, control
mask preparation, mask-only paired-support audit, and exploratory flux stack
are complete. The stack remains `science_measurement=false` and
`images_ready_for_science=false`.

## Provenance

- Approved control manifest SHA-256:
  `ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df`.
- Controls downloaded and verified: 1,041/1,041 FITS files; 58,271,158,080
  bytes; no failed files or partial downloads. The independent header/WCS/file
  audit found 248 `i` and 793 `i2` controls.
- One valid MASK header interchanges `MP_NOT_DEBLENDED` and `MP_UNMASKEDNAN`
  between bits 11 and 12. The download record preserves the per-file map and
  all mask logic uses the actual header map. HSC's [pipeline workshop
  documentation](https://hsc.mtk.nao.ac.jp/pipedoc/DataWorkshop2022/_downloads/2682edf0a494c8b005acf4efa51952ae/hsc_kaiseki_20221128.pdf)
  says to read the MASK header for bit/event correspondence.
- Control source-mask preparation: 1,041/1,041, 0 errors. These are automatic
  source masks and support diagnostics; no surface brightness was measured at
  this step.
- The v2 mask-only target/control audit is
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_pair_common_support_2026-11-07_02/summary.json`
  (SHA-256 `376d9c076a54166bf8bc0f30f1cdb26ce65be9c2c43743bf5b60d0367867f91c`).
  It reads headers and prepared masks, not science-image pixels.
- The exploratory stack summary is
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_target_control_stack_2026-11-07_01/summary.json`
  (SHA-256 `c8bac8dbdd0a886604c216b2410725eb97883b0cdc1a91fd47e7a7cfcaaaa20b`).
  It pins 4,105 input hashes and records `science_measurement=false`.

## Pair and filter accounting

The 1,041 control positions yield 1,022 filter-matched target/control pairs:
241 in `i` and 781 in `i2`. Seventeen pairs are excluded for filter mismatch
and two because the paired target has zero usable pixels. Two other targets
(`S23B_10647`, `S23B_8226`) have no matched control. No pair has zero common
mask support across the full 0–1,000 proper-kpc region. The stack keeps `i` and
`i2` separate.

## Usable-area assessment

Coverage is the intersection of target and control valid pixels after warping
to the 401×401 proper-kpc grid. “Supported pairs” means a pair has at least one
common pixel in the shell; this is not a minimum-area criterion.

| Shell (proper kpc) | `i` supported / 241 | `i` median common area | `i2` supported / 781 | `i2` median common area |
|---|---:|---:|---:|---:|
| 0–50 | 128 (53%) | 0.7% | 445 (57%) | 1.3% |
| 50–100 | 189 (78%) | 13.1% | 633 (81%) | 18.5% |
| 100–300 | 232 (96%) | 14.1% | 756 (97%) | 19.7% |
| 300–800 | 241 (100%) | 11.1% | 779 (99.7%) | 16.5% |
| 800–1,000 | 241 (100%) | 10.2% | 781 (100%) | 15.8% |

Only 23% of `i` pairs and 25% of `i2` pairs have any common support in the
innermost 0–25 kpc bin. The median common area there is zero. At 25–50 kpc the
median common area is 0.8% (`i`) and 1.7% (`i2`). The share of pairs with any
annular support reaches 90% only at about 175 kpc (`i`) or 150 kpc (`i2`). The
2-D maps have median per-pixel pair counts of 30 (`i`) and 134 (`i2`), with maxima of
59 and 214, respectively. The central region is therefore particularly weak,
and even outer bins use only a fraction of the nominal annular area.

## Radial stack (descriptive pair SEM)

The values below are equal-pair means of each pair's common-support annular
target-minus-control profile. The uncertainty shown is the between-pair SEM;
it is not the full error budget.

| Shell (proper kpc) | `i`: N, ΔSB ± pair SEM | `i2`: N, ΔSB ± pair SEM |
|---|---:|---:|
| 0–50 | 128; 584 ± 237 | 445; 541 ± 124 |
| 50–100 | 189; 30.7 ± 3.1 | 633; 37.3 ± 2.0 |
| 100–300 | 232; 4.47 ± 1.14 | 756; 4.77 ± 0.95 |
| 300–800 | 241; −0.04 ± 0.42 | 779; 0.15 ± 0.63 |
| 800–1,000 | 241; −0.56 ± 0.46 | 781; 0.84 ± 0.42 |

ΔSB is in nJy arcsec⁻². A target excess appears in 50–300 kpc and declines
toward zero at larger radii in both filters. This is a candidate diffuse-light
signal, not an ICL detection: the very large central value is dominated by
compact central light and has sparse common coverage. The outer profile is a
diagnostic of target/control balance, not a validated residual-sky correction.

## Qualitative and scientific limits

The 2-D maps show a bright, centrally concentrated component plus patchy
pixel-scale structure. The `i` map also has a few localized asymmetric
residuals. The radial signal is stable in sign and broadly similar between the
two filters at 50–300 kpc, while the outskirts fluctuate around zero. The
map-level pair counts are much smaller than the nominal sample size because
the stack deliberately uses only the strict target/control common mask.

The output does not yet account for Galactic extinction, PSF-wing matching,
residual-sky transfer, coadd-pixel covariance, field/region covariance, or
manual BCG identity review. Automatic masks are proposals, and the bright
central galaxies remain part of the measured central light. The pair SEM omits
those shared and instrumental uncertainties. No coverage threshold was chosen
after inspecting the result; no frames were removed based on measured flux.

**Use:** suitable for pipeline validation and for locating a candidate
50–300-kpc excess. The center and the present uncorrected profile are not ready
for a paper-level ICL measurement. Before science use, correct Galactic
extinction, match PSF wings, estimate and transfer residual sky, quantify
coadd covariance and spatial jackknife uncertainty, inspect BCG/source masks,
and publish sensitivity to mask size and common-versus-independent support.

## Outputs

- FITS stacks:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_target_control_stack_2026-11-07_01/s23b_icl_target_control_stack_i.fits`
  and `..._i2.fits`.
- Pair radial profiles:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_target_control_stack_2026-11-07_01/paired_radial_profiles.tsv`.
- Stacked radial profiles:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_target_control_stack_2026-11-07_01/stacked_radial_profiles.tsv`.
- Diagnostic figure:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_stack_quality_assessment_2026-11-07_01/stack_quality_diagnostic.png`.
