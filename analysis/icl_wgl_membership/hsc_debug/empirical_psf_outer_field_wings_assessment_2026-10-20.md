# Matched outer-field empirical PSF-wing pilot — 2026-10-20

## Purpose and input scope

The existing 333 Gaia-supported HSC star stamps stop at 96 pixels, corresponding
to a median 88.4 kpc, while the ICL analysis extends to 200–300 kpc. This pilot
uses full-size target coadds already present in the approved 800-cluster sample
and their 798 already-downloaded matched control coadds. It downloads no data.

The selector found 18 usable star/target/control pairs among the 333 complete
stamps. Each star is more than 800 proper kpc from its cluster center, each
target and matched control has geometric coverage through 300 kpc, the filter
matches, the HSC detection has zero flags, and the K2.5 central 12-pixel aperture
is at least 99% open. The pairs span 7 `i` and 11 `i2` frames. Their Gaia DR3
variability status is `NOT_AVAILABLE`; this is not evidence that they are
non-variable, and epoch matching remains unresolved.

For each pair, Astrafek measures 25-kpc radial bins to 300 kpc, applies the same
Kron-2.5 source masking and 120–150 or 160–190 pixel local-sky annulus to target
and control, normalizes both by the target star's annulus-specific central flux,
and forms target-minus-control profiles. A bin is retained only when at least
50% of its total area is unmasked and each azimuth sector has at least 25%
coverage. Bootstrap resampling is by matched target/control frame pair, not by pixel.

## Result and use decision

The computation completed as Slurm job `430734`. The output contains 864 unique
profile rows (18 candidates × 2 roles × 2 sky annuli × 12 radial bins) and 24
paired-stack bins. An independent validator recomputed every paired mean from
the profile table and verified output hashes.

The paired sample is too small for a science correction. Across either sky
choice, the number of complete target/control pairs in a radial bin ranges from
zero to four. No pair survives at 250–300 kpc. At 200–225 kpc there is one pair
per sky choice. This falls far short of a stable outer-ICL PSF-wing estimate.
The measured paired differences change sign between bins and sky choices; they
are residual diagnostics, not a detection of the stellar PSF wing. Therefore
the output is not suitable for correcting the ICL profile or evaluating an
ICL–satellite correlation.

The principal bottleneck is background and radial support after source and HSC
quality masks: only 6–8 of 18 target profiles pass the local-sky coverage gate,
and only 0–4 matched controls pass any given radial bin. Relaxing those gates
would risk leaving azimuthally biased source contamination in the wing profile.
The next PSF step needs a larger star sample with robust target/control support,
or an independently validated effective coadd PSF; merely extending the current
stamp is insufficient.

An earlier run in `_01` duplicated each control profile once per sky annulus and
is superseded. `_02` has the corrected one-row-per-pair/annulus/radius structure.
Both runs remain explicitly `science_measurement=false` and `psf_approved=false`.

## Reproduction and files

Astrafek implementation:

- [`stack_empirical_psf_outer_field_wings.py`](stack_empirical_psf_outer_field_wings.py)
- [`run_empirical_psf_outer_field_wings.sbatch`](run_empirical_psf_outer_field_wings.sbatch)
- [`validate_empirical_psf_outer_field_wings.py`](validate_empirical_psf_outer_field_wings.py)
- [`test_stack_empirical_psf_outer_field_wings.py`](test_stack_empirical_psf_outer_field_wings.py)

Reproduce on the cluster with a new, absent output directory:

```bash
OUTPUT_DIR=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_outer_field_wings_2026-10-20_03 \
  sbatch analysis/icl_wgl_membership/hsc_debug/run_empirical_psf_outer_field_wings.sbatch
```

Authoritative output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_outer_field_wings_2026-10-20_02/`.
The summary SHA-256 is `c8085c8a7459a34ab03e551de6d2d37cccf637a80156bff4f70ba61b0dc46369`;
the independent `validation.json` SHA-256 is
`c70a144aae21922044d14a230404c8633a7d9be755580b976d7752a0fad860f6`.
The reproduction script SHA-256 is
`2041a7da07968734158509d647b55c4217ef9bdc9dd14e64047755e60f3dacda`; the
pipeline rechecked the source-image SHA-256 values before measuring them.
No target or control image was downloaded in this step.
