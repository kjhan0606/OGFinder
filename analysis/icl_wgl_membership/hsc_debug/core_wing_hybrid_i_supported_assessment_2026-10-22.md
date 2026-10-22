# Supported i-frame hybrid PSF sensitivity — 2026-10-22

## Execution and independent checks

The frozen 12-model grid was run on all 13 existing K2.5 / 40%-ICL-coverage
`i` frames with approved position-PSF mappings. Each configuration uses splice
radii 1.8/2.4/3.0 arcsec, a 0.6-arcsec blend, and public-wing amplitudes
0.50/0.65/0.80/1.00. Diffuse light uses the center-position local core; each
satellite envelope uses its own position-specific core. Each frame also
stores a separate noiseless unit point-source profile convolved with the same
hybrid and evaluated through the saved ICL exclusion mask.

Slurm job 432803 completed on `syntax`/`syn04` in 2:07:02, requested one CPU
and 32 GB, used no GPU, and reported 4,363,184 KiB maximum RSS. All 12
configurations completed across the same 13 frames (156 frame records total),
with no run errors. Independent `audit_batch` validation recomputed each
stored injection response: maximum arithmetic error was 0 for every frame.
The grid assembler also checked the checkpoint/summary/frame SHA-256 values,
the point-profile schema, and all approval-state invariants. It did not
independently rebuild point-source profiles from the FITS inputs.

## Main response results

The table gives the median across the 13 selected frames for each estimator;
the model interval is the minimum–maximum of those 13-frame medians across the
12 PSF alternatives. These are descriptive sensitivity summaries, not
confidence intervals. The 13 frames are the position-PSF-supported subset,
not a random sample of all 800 clusters.

| Injection → measurement | Fixed-mask median; PSF-grid envelope | Redetection median; PSF-grid envelope | Median of paired redetection − fixed |
|---|---:|---:|---:|
| Diffuse → ICL 200–300 kpc | 0.2250; 0.2250–0.2252 | 0.2141; 0.2141–0.2150 | −0.0128; grid −0.0128 to −0.0115 |
| Satellite envelopes → ICL 200–300 kpc | 0.4903; 0.4896–0.4917 | 0.3701; 0.3610–0.3852 | −0.1022; grid −0.1045 to −0.0825 |
| Diffuse → satellite 50–75 kpc | 0.00073; 0.00069–0.00074 | −0.0394; −0.0395 to −0.0392 | −0.05135; grid −0.05135 to −0.05135 |
| Satellite envelopes → satellite 50–75 kpc | 0.2632; 0.2597–0.2743 | 0.2874; 0.2852–0.3093 | −0.0536; grid −0.0977 to −0.0333 |

Frame-to-frame spread is much wider than the shift among model medians. For
example, with the central `r2.4_a0.65` model, fixed-mask diffuse→ICL has
16th–84th percentiles 0.2219–0.2272, but its redetection branch reaches
−2.7339 on one frame; the paired difference spans −2.9554 to 0.0464. In the
satellite-envelope→satellite measurement, the paired difference spans
−3.0894 to +2.1677 across these frames (clusters 430 and 4982 are the extreme
negative and positive cases). Its 16th–84th percentile interval is
−0.3158–+0.2065. The median paired shift therefore does not support a
single global transfer correction. Rerunning detection and masks changes
which source pixels and satellites enter the estimator, sometimes in either
direction.

The measured ICL–satellite association is still not established by these
injections. They only test the measurement pipeline's response to known
synthetic light; they do not measure the astrophysical correlation or its
significance.

## Point-source wing profile

For the central `r2.4_a0.65` model, the median fraction of unit point-source
flux falling in valid (unexcluded) pixels was 0.00259 in the 50–100 kpc
annulus and 0.000957 in 200–300 kpc. Across the 12 PSF alternatives, the
13-frame median ranges were 0.00200–0.00391 and 0.000741–0.001446,
respectively. The median valid-area mask coverage was 50.1% and 67.2% in the
inner and outer annuli; these coverage values are set by the masks and do not
vary with PSF model.

This calculation is a noiseless profile of a unit point source through the
saved mask, not source extraction or a noise-injected recovery experiment. It
does not include source Poisson noise, coadd covariance, sky subtraction,
observed-source subtraction, or spatial/epoch/seeing matching of the public
HSC-I reference. The profile is a useful wing-leakage diagnostic, not an ICL
correction.

## Reproduction and state

Astrafek sources:
[`run_core_wing_hybrid_grid.py`](run_core_wing_hybrid_grid.py),
[`run_core_wing_hybrid_i_supported.sbatch`](run_core_wing_hybrid_i_supported.sbatch),
[`summarize_core_wing_hybrid_grid.py`](summarize_core_wing_hybrid_grid.py),
and [`validate_extended_psf_injection_batch.py`](validate_extended_psf_injection_batch.py).

Output root: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_core_wing_hybrid_i_supported13_2026-10-22_01/`.
The first `grid_audit.json` summary is retained. The superseding v2 audit also
checks that all 12 configurations contain the same 13 unique cluster IDs;
its response and point-source envelopes exactly match the first summary.
`grid_audit_v2.json` SHA-256:
`c0545e9f36c1b0d47589744abfcf21264461e8bf467f6f63fd379436ee35ba2c`.
`grid_status.json` SHA-256:
`91f5aed2dbff64fd5ae07290c30b8f77b69fde15535c881d1822352ba25747b6`.
Audit script SHA-256:
`f8fd10adbf5bdc7ea8f5e11d972cd31804a2ea58072c5d76ce7be39f784c6974`.

No additional HSC data was requested. The original science images were not
modified. `science_measurement=false`, `psf_approved=false`,
`observational_validation_passed=false`, and `flux_covariance_available=false`
remain in force. The next science gate is to test the pixel/noise/sky response
on matched controls and to obtain a PSF-wing treatment validated over the
relevant physical radii; do not apply these medians as corrections.
