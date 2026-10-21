# Full-sample ICL sensitivity to target/control PSF-core mismatch — 2026-10-21

## Scope and result

Astrafek joined the newly received HSC-I model-PSF cores at the 800 target and
798 matched-control centers, then tested whether the target-minus-control
major-axis core-width difference covaries with the existing ICL
200–300 kpc target-minus-control profile. This is a PSF-core sensitivity
diagnostic, not a corrected ICL measurement or a detection claim.

All 1,598 center positions now have a checksum-verified PSF core: 1,542 from
the approved target/control-center batch and 56 from the earlier approved
injection PSF batch. Of the 798 matched pairs, the median target-minus-control
moment FWHM difference is −0.0070 arcsec, with a position-to-position q16–q84
range of −0.0572 to +0.0392 arcsec and median absolute difference 0.0319
arcsec. The target and control major-axis FWHM medians are 1.0670 and 1.0734
arcsec, respectively. These are descriptive distributions across locations,
not confidence intervals. The model stamps are only 41–43 pixels across.

The ICL join contains 137 common-area-qualified profile pairs from the
800-cluster parent: 44 `i` and 93 `i2`. The rows use the existing K2.5,
50–300 kpc selection with at least 40% common area in each radial bin. The
outcome is the saved 200–300 kpc surface-brightness difference; the PSF
predictor is the target-minus-control moment-based major-axis core FWHM.
Connected 5-degree target/control sky blocks were resampled 5,000 times. The
percentile ranges below are descriptive block-bootstrap sensitivity ranges;
no p-values were calculated.

| Stratum | N / blocks | Model | Partial correlation | Block-bootstrap 2.5–97.5 percentile | Leave-one-block-out range |
|---|---:|---|---:|---:|---:|
| i | 44 / 15 | Unadjusted | 0.106 | [−0.175, 0.419] | [0.026, 0.222] |
| i | 44 / 15 | z, richness, common coverage | 0.154 | [−0.129, 0.453] | [0.062, 0.240] |
| i | 44 / 15 | plus central-mask flag | 0.016 | [−0.349, 0.388] | [−0.155, 0.117] |
| i2 | 93 / 37 | Unadjusted | 0.066 | [−0.121, 0.270] | [0.036, 0.168] |
| i2 | 93 / 37 | z, richness, common coverage | 0.038 | [−0.173, 0.250] | [−0.006, 0.146] |
| i2 | 93 / 37 | plus central-mask flag | 0.036 | [−0.164, 0.270] | [−0.004, 0.155] |

Every block-bootstrap interval includes zero. In `i`, the point estimate changes
substantially after accounting for the central-mask flag; in `i2`, adjusted
associations stay close to zero. The present sample therefore gives no stable
evidence that model-core mismatch explains the measured ICL excess. This does
not rule out extended-PSF-wing or sky-transfer systematics: the model core is a
seeing proxy and does not measure the PSF at 200–300 kpc.

## Reproducibility and limits

The center-level metrics and pair table are produced by
[`analyze_hsc_target_control_psf_cores.py`](analyze_hsc_target_control_psf_cores.py).
The ICL sensitivity is reproduced by
[`analyze_full_sample_icl_psf_core_dependence.py`](analyze_full_sample_icl_psf_core_dependence.py)
and independently checked by
[`validate_full_sample_icl_psf_core_dependence.py`](validate_full_sample_icl_psf_core_dependence.py).
The focused audit, join, and regression tests pass (9 tests total).

The analysis output is outside Git at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_core_dependence_full_sample_2026-10-21_01/`.
Its `summary.json` SHA-256 is
`f67019b8cb6f3c096b45e8640d489022db19f2b1d29f7b4db1ab89441936c34e`; the joined
137-row table SHA-256 is
`d6e0f0b669e3eb2c0f7a57c37669260da4ad17216daaf66fda2c60af7223fcf5`.
The independent validator reproduced all six point correlations, the six
5,000-replica block-bootstrap percentile triplets (30,000 resampled
correlations), and all leave-one-block-out ranges.

The core-center diagnostic is saved at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/target_control_psf_core_diagnostic_2026-10-21_01/`.
Its `summary.json` SHA-256 is
`cb9244c39e5608932434eec8eb0e7aca0a50ebac585cd769a19f996937843bb0`; its
1,598-row center table SHA-256 is
`7c61821b2c73e7a7c3a43bc8c4e58be53e5c7e53ad564cc1d076734a8cb83a27`; its
798-row paired-core table SHA-256 is
`957f425d68b2ed1ccdbcc0686b338cb7f02a0c752bab06fcfa7ec6909c786c90`.

The analysis did not modify image pixels, subtract a PSF correction, request
additional HSC data, or declare a scientific result. It leaves
`psf_approved=false` and `science_measurement=false`. The extended PSF wings,
local sky-subtraction transfer, extinction, photometric calibration, and
mask/noise systematics still need separate validation before interpreting the
ICL–satellite relation scientifically.
