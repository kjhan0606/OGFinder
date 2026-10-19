# HSC target–control 2-D ICL stack assessment — 2026-10-19

## Result and audit

Astrafek's `icl.wgl.warp_metric` was used to stack the existing HSC PDR3 global-
sky coadd/bg cutouts on a common ±500 proper-kpc grid with 5-kpc pixels. Target
and matched-control images were kept in separate `i` and `i2` stacks, converted
from `FLUXMAG0` to nJy/arcsec², and differenced only on their shared valid
interpolation support. The same pair weights enter the target, control, and
difference maps. The run hashed and rechecked every included input image and target mask.

Of 800 targets, 798 had downloaded matched controls. Fifteen of those pairs had
different primary `FILTER` values (`i` versus `i2`) and were excluded. The
remaining 783 pairs are 302 `i` and 481 `i2`; two targets have no matched
control. The output consists of two 201×201 FITS cubes, per-pair and stacked
radial profiles, and diagnostic plots. The independent verifier confirmed the
pair selection, output checksums, all 40 radial-summary rows, map support/count
agreement, and `DELTA = TARGET − CONTROL` to float32 precision.

The summary initially serialized `science_measurement=false` as JSON integer
`0`. The FITS products already carried `SCIVAL=false`, and no data arrays were
changed. The exact script used for the run is preserved in the data output and
hashed there; `repair_icl_stack_summary_serialization.py` corrected only that
JSON type and recorded both pre- and post-correction checksums. The current
stack source now serializes booleans correctly and has a regression test.

## Radial result

The table reports the equal-weight mean of each pair's annular mean. The quoted
uncertainty is the between-pair standard error only.

| Filter | Radius (kpc) | Pairs with support | Median common-area coverage | ΔSB (nJy/arcsec²) | Pair SEM |
|---|---:|---:|---:|---:|---:|
| i | 0–50 | 141/302 | 0.062 | 521.6 | 89.9 |
| i | 50–100 | 235/302 | 0.150 | 41.4 | 2.90 |
| i | 100–200 | 277/302 | 0.145 | 12.9 | 1.49 |
| i | 200–300 | 283/302 | 0.117 | 2.08 | 1.56 |
| i | 300–400 | 295/302 | 0.103 | 4.36 | 1.48 |
| i | 400–500 | 299/302 | 0.101 | 4.34 | 2.65 |
| i2 | 0–50 | 240/481 | 0.067 | 412.0 | 43.2 |
| i2 | 50–100 | 378/481 | 0.190 | 44.6 | 2.32 |
| i2 | 100–200 | 447/481 | 0.177 | 10.7 | 1.01 |
| i2 | 200–300 | 467/481 | 0.154 | 3.18 | 0.78 |
| i2 | 300–400 | 476/481 | 0.130 | 0.14 | 0.73 |
| i2 | 400–500 | 480/481 | 0.127 | −0.30 | 0.76 |

The central pixel is dominated by the BCG/core (52,514 nJy/arcsec² in `i`,
43,991 in `i2`), so it is not an ICL-only measurement. The broader profile
shows positive target-minus-control light at 50–200 kpc in both filter subsets,
then values near zero by 300–500 kpc in `i2`. The `i` profile instead has a
roughly 4.3 nJy/arcsec² positive outer offset. That filter-dependent outer
baseline is a warning that residual sky or calibration differences can mimic
low-surface-brightness light. Do not combine `i` and `i2` or interpret their
unsubtracted offsets as a physical color or ICL detection.

## Image-level review and scientific status

The target maps show the central BCG and a declining diffuse-looking halo; the
control map is noise-like. The difference maps retain positive structure near
the center and mottled/asymmetric low-level residuals farther out. The residual
shape is qualitatively consistent with a central excess, but the pairwise
common-area support is sparse: among pairs contributing to each annulus, the
median covered fraction is only about 6–19%; roughly half the pairs have no
common support at all inside 50 kpc. The median coverage and the set of
contributing pairs change with radius.

This is an exploratory stacking product (`science_measurement=false`), not a
confirmed ICL detection or a measurement of ICL–satellite correlation. The
between-pair SEM does not include the spatial covariance of coadd noise,
residual-sky transfer, PSF mismatch, extinction/calibration, or astrophysical
covariance. The available position-dependent PSFs cover only the core and do
not constrain extended wings. No background-plane correction was imposed in
this run. The large descriptive mean/SEM ratios in some bins are therefore not
to be read as detection significance.

The next ICL gate is to resolve the sparse and radius-changing mask support,
measure a matched-filter residual-sky baseline and its uncertainty, and test
the stack against PSF-wing and mask-radius choices. Only after those checks
should the cluster-centered ICL profile be joined to the satellite-centered
profiles and their shared-pixel covariance. This result does not change the
separate weak-lensing pipeline or its random-point subtraction requirement.

## Reproducible artifacts

- Astrafek core: `icl/wgl.py` (`warp_metric`).
- Stack: `stack_hsc_icl_target_control_2d.py`; Slurm entry point:
  `run_icl_target_control_2d_stack.sbatch`.
- Independent audit: `verify_icl_target_control_2d_stack.py`.
- Profile and figure evaluation: `evaluate_icl_target_control_2d_stack.py`.
- Run data: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_2d_matched_stack_2026-10-19_02/`.
- The run summary, input/code hashes, filter exclusions, correction record,
  independent verification, FITS stacks, radial tables, and plots are retained
  in that directory. No new HSC request or download was made.
