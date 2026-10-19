# HSC ICL target–control 2-D stack

`stack_hsc_icl_target_control_2d.py` constructs an exploratory, matched-
control-subtracted 2-D stack from the downloaded HSC PDR3 coadd/bg cutouts. It
uses the Astrafek `icl.wgl.warp_metric` resampler and keeps the native `i` and
`i2` products in separate outputs.

## Inputs and selection

- 800 target cutouts and their saved `MASK_K2P5` products from
  `prestack_cluster_masks_2026-10-16_01`.
- Up to 798 matched-control cutouts from `control_download_01`, with the target
  redshift assigned to each control only for angular-distance scaling.
- A pair is stacked only if target and control have the same primary `FILTER`
  (`i` or `i2`). Filter evidence is checked against both the download status
  and FITS headers. Mismatches are listed as exclusions in `summary.json`.
- Control K=2.5 masks are regenerated from the existing science, variance, and
  HSC mask planes with the same SEP/Kron mask utility used by control QA.

## Measurement

Each map is sampled onto a 201 × 201 grid spanning ±500 proper kpc (5 kpc per
pixel). HSC `FLUXMAG0` converts the input to nJy/pixel; Astrafek converts that
to nJy/arcsec² while resampling. Target and control maps contribute only where
both interpolated samples have full valid-pixel support. The stack gives each
matched pair equal weight at each pixel and records the number of contributing
pairs. The radial table first measures one annular mean per pair, then reports
the equal-weight pair mean, between-pair SD/SEM, common-area coverage, and
radius-dependent pair count.

## Run and products

Submit `run_icl_target_control_2d_stack.sbatch` from the repository root. It
refuses to reuse an existing output directory. The default data output is
`icl_2d_matched_stack_2026-10-19_k2p5/` under the HSC sample directory. Select
the alternative target/control Kron masks with `--mask-factor 3.5` or `5.0` and
write each run to a new `--work` directory. The paired K=3.5/K=5.0 sensitivity
job is `run_icl_target_control_mask_sensitivity.sbatch`.

Each filter FITS contains target, matched-control, target-minus-control,
pixelwise pair count, between-pair SEM maps, and a radius grid. The TSV outputs
are `paired_radial_profiles.tsv` and `stacked_radial_profiles.tsv`; `summary.json`
contains input and code SHA-256 values, filter exclusions, processing
definitions, and limitations. After independent verification,
`evaluate_icl_target_control_2d_stack.py` also writes broad-annulus summaries
and per-filter map/radial diagnostic figures. The FITS primary header records
the mask factor, and all outputs remain separate by filter and mask choice.

`bootstrap_icl_2d_spatial_blocks.py` adds pair and connected target/control
sky-block bootstrap intervals for the 50–500 kpc profile. It checks four cell
widths (2.5°, 5°, 10°, 15°), verifies the source stack headers and manifests,
and retains the per-annulus covariance matrices and diagnostic figures outside
Git. See `icl_2d_spatial_block_bootstrap_assessment_2026-10-19.md` for the
interpretation; these intervals remain exploratory.

`stack_hsc_icl_independent_support.py` measures each pair's target and control
annular means on their own valid pixels, then compares that result with the
strict common-support estimator on the same pairs and the exact same 25 kpc
annuli. Its Slurm sensitivity run
and output verifier are `run_icl_independent_support_sensitivity.sbatch` and
`verify_icl_independent_support.py`. The paired comparison distinguishes the
effect of relaxing pixel intersection from the effect of admitting additional
pairs; the alternate estimator remains a sensitivity because target/control
masks may remove different source environments.

## Interpretation limits

The output is marked `science_measurement=false`. It does not apply extinction,
residual sky-transfer correction, or PSF-core homogenization; the available
position-dependent PSF stamps do not constrain extended wings. The pixelwise
between-pair SEM is not the full covariance-aware uncertainty and does not
include sky, PSF, calibration, or astrophysical covariance. Target and control
noise independence is an assumption. Inspect pair-count and coverage maps
before interpreting the outskirts. This ICL stack is independent of the
background-source weak-lensing stack.
