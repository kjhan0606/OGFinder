# Position-specific PSF core plus public-wing sensitivity — implementation note

## Current state

Astrafek now has `position_psf.build_core_wing_hybrid_psf`, which takes a
signed, unit-normalized local HSC coadd PSF core and a signed, unit-normalized
public extended reference on the same pixel grid. It zero-pads the local core,
uses a smooth radial transition around a caller-selected splice radius, scales
the public contribution by a caller-selected amplitude, and normalizes the
resulting signed kernel to unit sum. It reports the pre-normalization sum and
keeps `psf_approved=false` and `science_measurement=false` in its metadata.

This is only the kernel-construction primitive. It has not been applied to
science images and does not establish that the public HSC-I wing is matched to
any target/control coadd. In particular, HSC-I remains a proxy for `i2`.

## Validation and one-frame pilot

The focused hybrid/injection/validator suites pass 23 tests. `py_compile`,
`bash -n` for the CPU-only Slurm array script, and `git diff --check` pass.
The pilot used existing files only; it made no HSC request and changed no
science coadd.

One i-band frame (cluster 1056; 10 satellite positions) completed with splice
radius 2.4 arcsec, blend width 0.6 arcsec, and public-wing amplitude 0.65.
The local core at the center is 43×43; satellite cores include 41×43 and
43×43 files. All 10 per-position hybrids and the center hybrid have signed
sum 1 within floating-point precision, preserve negative PSF pixels, and carry
`psf_approved=false`/`science_measurement=false`. The independent saved-JSON
validator passed with maximum response arithmetic error 0.

With fixed masks, unit diffuse injection recovered 0.230718 in the ICL
aperture; the sum of satellite envelopes contributed 0.506852 there. A unit
satellite-envelope injection recovered 0.314188 in the satellite aperture;
diffuse-to-satellite leakage was −0.012496. Relative to the existing public
PSF fixed-mask result for this same cluster, the changes are approximately
+0.000105, −0.000554, −0.000093, and −0.012785 for those four responses. The
single-frame values are sensitivity diagnostics, not population estimates.

A separate noiseless unit point-source diagnostic was measured through the
existing saved masks. The valid masked area covers 46.3% of the 50–100 kpc
annulus and 55.4% of the 200–300 kpc annulus; the masked annular flux fractions
are 0.001996 and 0.000834 per unit input flux, respectively. This is not a
noise-injected recovery test and does not model source Poisson noise, coadd
covariance, sky transfer, or observed-source subtraction.

The result is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_hybrid_psf_i_pilot_2026-10-21_r2p4_a0p65_perpos_point_01/`.
`summary.json`, `cluster_1056.json`, and `validation.json` SHA-256 values are
`8a2fa26dd197aae9f29fafaaaac7e91afeab0bcc9c1d2299b2a73704a367d16e`,
`ddb650883c7e664fb9c026d28a83768c195d7b5ccaa867b6b5c58d8048666711`, and
`d547fba7bfbae4869f48c6946a316190c860afe4dd711b301c2205e55d03c568`.

## Next execution gate

Run the predeclared 12-point splice-radius/amplitude grid on the same pilot
frame, with point-source diagnostics kept separate from diffuse/envelope
injections. Compare both fixed-mask and redetection outputs, and independently
validate every frame before broadening to supported i-band frames. Keep i2 as
a separately labeled HSC-I proxy sensitivity. Report sky transfer and coadd
noise limitations independently; neither this one-frame run nor the grid can
approve a PSF or turn on science measurements.
