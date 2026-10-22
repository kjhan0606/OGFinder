# S23B candidate overlaps in matched controls: ICL sensitivity

## Question and method

Seven centers from the pinned 1,333-object S23B CAMIRA expansion table fell
inside existing matched-control coadd WCS footprints. This test asks whether
removing the corresponding control pairs changes the HSC-I target-minus-control
profile. It is a metadata-triggered sensitivity test: center containment does
not establish that a candidate is detected in the coadd or that it contributes
diffuse light.

The original K=2.5 stack used 783 matched pairs (302 `i`, 481 `i2`). One of the
seven control overlaps, target 1746 / candidate 3509, was already excluded by
the baseline target/control filter mismatch. It was not counted twice. The
alternate stack therefore removes the other six pairs and retains 777 pairs
(300 `i`, 477 `i2`). Both stacks use the same 201 × 201 grid over ±500 proper
kpc, 5 kpc pixels, masks, input images, and separate filter products.

The Astrafek verifier checked the alternate pair list, all six additional
exclusions, both filter maps, radial products, summary provenance, and source
hashes. Each map has 40,401 covered pixels; the maximum `DELTA - (TARGET -
CONTROL)` residual is 0.00391 nJy/arcsec² for `i` and 0.00098 nJy/arcsec² for
`i2`. The output is still marked `science_measurement=false`.

## Profile changes

The table gives annulus-area-weighted means from the 25 kpc radial bins. Units
are nJy/arcsec².

| Filter | Radius (kpc) | Baseline | Exclude six controls | Change | Relative change |
|---|---:|---:|---:|---:|---:|
| i | 50–200 | 18.05296 | 18.11743 | +0.06447 | +0.357% |
| i | 200–300 | 3.95750 | 4.00977 | +0.05227 | +1.321% |
| i2 | 50–200 | 18.43218 | 18.42854 | −0.00363 | −0.020% |
| i2 | 200–300 | 3.35652 | 3.29386 | −0.06266 | −1.867% |

The largest per-annulus changes were 0.08947 nJy/arcsec² (0.877%) for `i`
150–175 kpc; 0.08874 (2.633%) for `i` 250–275 kpc; 0.24161 (1.452%) for `i2`
75–100 kpc; and 0.20313 (4.611%) for `i2` 200–225 kpc. Radial pair counts
fall by one in each listed `i` bin, by two to four in `i2` 50–200 kpc, and by
four in `i2` 200–300 kpc, reflecting the existing pairwise valid-pixel support.

The broad means move by less than 1.9%, while a single outer `i2` annulus moves
by 4.6%. The sign is not the same in the two filters. This does not support a
large, coherent control-contamination offset in these broad averages, but it
also does not prove the controls are clean or that ICL is detected. The center
overlap could be incidental, could have no detectable flux effect, or could
interact with the coadd background subtraction. No uncertainty or significance
for the *difference* is reported: the alternate stack is a subset of the
baseline, and their covariance has not been estimated.

## Reproduction and limits

The WCS audit read FITS headers only. The stack itself read existing FITS data
and wrote to a new directory; it requested no data from HSC and left the
baseline untouched. Reproduce the stack with
`run_icl_s23b_control_exclusion_sensitivity.sbatch`, verify it with
`verify_icl_target_control_2d_stack.py`, then use
`compare_icl_s23b_control_exclusion.py` to create the comparison JSON.

The run directories and provenance are:

- Baseline stack: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_2d_matched_stack_2026-10-19_02/`; summary SHA-256 `9c804cfb12ff2b86a5bdc8fdbfea94aa1dc473711e6f7ac6f4d233cc6b853173`.
- Alternate stack: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_2d_s23b_control_exclusion_2026-10-22_01/`; summary SHA-256 `45353a85ff864d5ced9e8334157207d6892160ccae9fa88d281b672bec2cdf21`.
- Alternate independent verification: `independent_verification.json`, SHA-256 `f0fb9a4ba6e9048366e481507c73074cebb88cccd6bb87cfa842963bb69535ca`.
- Comparison: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_control_exclusion_comparison_2026-10-22_01/comparison.json`, SHA-256 `6efdfa9b9df334f5e038e98a677ed4ce55f40e4c15178528223145c826bb1f7b`.
- Input WCS audit: `icl_s23b_existing_cutout_support_2026-10-22_01/summary.json`, SHA-256 `071377813fcfbf15f298a5d48f3e7a3f0b74f60899d88c46eafb7d5bb58ffebb`.

This sensitivity does not resolve residual sky transfer, PSF wings, extinction,
coadd-noise covariance, or the effect of satellite masking on ICL flux. Those
gates remain before an ICL–satellite association can be treated as a science
measurement. It is separate from the background-source weak-lensing stack.
