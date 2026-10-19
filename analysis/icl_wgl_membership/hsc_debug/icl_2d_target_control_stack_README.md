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
`icl_2d_matched_stack_2026-10-19_02/` under the HSC sample directory.

Each filter FITS contains target, matched-control, target-minus-control,
pixelwise pair count, between-pair SEM maps, and a radius grid. The TSV outputs
are `paired_radial_profiles.tsv` and `stacked_radial_profiles.tsv`; `summary.json`
contains input and code SHA-256 values, filter exclusions, processing
definitions, and limitations. After independent verification,
`evaluate_icl_target_control_2d_stack.py` also writes broad-annulus summaries
and per-filter map/radial diagnostic figures.

## Interpretation limits

The output is marked `science_measurement=false`. It does not apply extinction,
residual sky-transfer correction, or PSF-core homogenization; the available
position-dependent PSF stamps do not constrain extended wings. The pixelwise
between-pair SEM is not the full covariance-aware uncertainty and does not
include sky, PSF, calibration, or astrophysical covariance. Target and control
noise independence is an assumption. Inspect pair-count and coverage maps
before interpreting the outskirts. This ICL stack is independent of the
background-source weak-lensing stack.
