# Public HSC PDR3 extended-PSF reference assessment — 2026-10-21

## Decision

The local bright-star stack is not a usable ICL-wing correction: its strict
filter-separated support is effectively zero over 200–300 kpc. A directly
relevant public reference exists: the published HSC-SSP PDR3 extended PSFs
were built for the intermediate `coadd/bg` product, reach 5.6 arcmin, and are
provided as 4001×4001-pixel FITS models. That angular range exceeds the
200–300 kpc science annulus for the current cluster redshift range.

This is a promising external reference and a major improvement over extrapolating
the current short star stamps. It is **not yet an approved correction** for our
stack: the published model is a survey-level empirical PSF, not a local PSF for
each target coadd/frame. HSC's PDR3 processing documentation says the similar,
not identical `i` and `i2` exposures are combined during coaddition and an
object-specific effective-response afterburner converts photometry. Thus this
public `i` model is a plausible coadd-level wing reference, but its weighting
and response still need to be matched to the Astrafek cutouts. Keep the
analysis' `i`/`i2` star-support diagnostics separate; do not interpret those
labels as proof that there are two independent coadd images. Use the public PSF
as a radial-wing/forward-scattered-light sensitivity template and retain
`psf_approved=false` until local core/wing transfer, photometric normalization,
and sky/mask effects are tested. Do not describe a template subtraction as a
measured ICL profile.

## Exact public artifacts proposed for the next gate

The authors' GitHub README points to Zenodo record
[`10.5281/zenodo.14915578`](https://zenodo.org/records/14915578), version
`v4.0.0`, whose metadata declares CC BY 4.0. The record lists two `i`-band
4001×4001 FITS variants:

| Artifact | Bytes | Zenodo MD5 |
|---|---:|---|
| `HSC_PSF_i.fits` (rotationally symmetric) | 9,809,280 | `5bf3317bcc41a83c81c25acab2faab09` |
| `HSC_PSF_i_asymmetric.fits` | 9,809,280 | `665e7b459cda203d6610133c3626194d` |
| **Total** | **19,618,560** | |

The user approved and Astrafek downloaded only these two files to
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/extended_psf_public_reference_2026-10-21_01/`.
Both byte counts and Zenodo MD5 values match. The manifest SHA-256 is
`91ce5b5495bc9c306497d6e7bd9c01472d3e82b102b6db86af9a6dfc9778e8ea`;
the symmetric and asymmetric FITS SHA-256 values are respectively
`8897af89767e5fbb8d5f9c70965126bc89409f6c65ecd8d117c7995f34be2d85` and
`b2ab473c6af4da11c1df74850bdd0a8cc0c1da4d109c118d1f4033b4235027fd`.
No HSC image, star catalogue, shape row, source-redshift file, reconstruction
script, or other Zenodo artifact was fetched. The asymmetric model preserves
azimuthal structure; the symmetric model is the radial-profile reference.
Their difference can contribute a model-form sensitivity, but not a complete
uncertainty budget.

## Compatibility and scientific use

The paper says the final aggressive local-sky subtraction is absent from the
`coadd/bg` inputs used to construct these PSFs. That matches the project's
HSC intermediate cutouts and is more appropriate for diffuse-light work than
the final detection-optimized coadd. Its published construction uses saturated
bright stars for the outer wings and fainter stars for inner regions, combining
the four regions into a band PSF. The model reaches roughly 2000 pixels, or
5.6 arcmin at the HSC pixel scale.

Planned local validation gates:

1. Verify exact byte counts and Zenodo MD5; compute and record SHA-256 locally.
2. Inspect FITS HDUs, image dimensions, pixel scale, center, finite values,
   radial coverage, and encircled-flux normalization without changing the
   source files.
3. Reproduce radial and azimuthal profiles from both variants and compare
   them with existing `i` and `i2` empirical star/control profiles over their
   actual overlapping support; do not pool the two subfilters.
4. Convolve representative satellite/BCG light models with the two PSFs and
   quantify scattered-light leakage into 50–100 and 200–300 kpc annuli. Keep
   background/sky transfer as an independent term; a PSF model cannot repair
   oversubtracted or spatially varying sky.
5. Accept an ICL correction only after local-frame agreement, mask/coverage
   treatment, and injected-profile recovery pass predeclared limits. Otherwise
   publish the two PSF variants only as sensitivity bounds and preserve
   `science_measurement=false`.

## First FITS characterization

Astrafek's local-only implementation is
[`audit_public_extended_psf_reference.py`](audit_public_extended_psf_reference.py).
It checks FITS extension 1, shape, finite pixels, peak offset, total flux,
annular profile, physical aperture leakage, and eight-sector flux shares, and
writes checksummed CSV/JSON without setting either approval flag. Six synthetic
tests pass. The actual run used Python/Astropy on lageunha CPU, without a GPU or
another HSC request. The output directory is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/extended_psf_public_reference_audit_2026-10-21_03/`.
The summary SHA-256 is
`3e94417b453eb7cf2ebb98f51ba044478b114ff9c7b5a95604fdec1cc728c8d8`.
The independent [`validate_public_extended_psf_reference_audit.py`](validate_public_extended_psf_reference_audit.py)
reopened both raw FITS files, verified Zenodo MD5 and SHA-256, and independently
recomputed physical aperture fractions. Its status is
`validated_reference_characterization_only`; `validation.json` SHA-256 is
`258d2716144e50ef1803792ea1bfcb00b13c6e846fe48e2fb4fed45444ef70e4`.

Both FITS images are 4001×4001, peak exactly at pixel (2000,2000), and have
raw sums within 1.2×10⁻⁸ of unity. At 0.168 arcsec/pixel, the radial profile
reaches 5.6 arcmin and encloses 99.9908% of total signed flux. Under the
project's flat `H0=70`, `Omega_m=0.3` cosmology, the 200–300 kpc annulus is
44.9–67.4 arcsec at `z=0.30` and 31.2–46.8 arcsec at `z=0.55`. A point source
places 0.2221% (`z=0.30`) to 0.2415% (`z=0.55`) of its model flux in that
annulus; the 50–100 kpc annulus contains 0.7826–0.9721%. The symmetric and
asymmetric models' radial annular fractions differ by less than 2.2×10⁻⁷
fractionally in those ranges; eight-sector shares are within about 0.1% RMS of
equal. This characterizes the template's wing scale but does not establish
contamination in cluster images, which also depends on source light profiles,
local seeing, coadd weights, masks, and sky transfer.

The first `_01` audit attempt is preserved as a failed partial output after a
header-serialization bug; `_02` lacks the physical-aperture diagnostics. Use
`_03` as the complete characterization. `psf_approved=false` and
`science_measurement=false` remain in force.

## Transfer check against the approved local PSF cores

The central 43×43 pixels of each public model were renormalized only for a
shape comparison with the previously approved, checksum-verified 744-position
PDR3 PSF Picker batch. The crop retains 95.7803% of the extended model's signed
flux. The public core's major-axis Gaussian-equivalent FWHM is 1.3168 arcsec,
21.5% above the local-core median of 1.0838 arcsec. Its `R50`, `R80`, and `R90`
are 0.385, 0.725, and 1.087 arcsec, compared with local medians 0.363, 0.653,
and 0.913 arcsec. The public model's core therefore should not replace the
position-specific core directly.

The repeatable comparison is
[`compare_public_extended_psf_to_picker_cores.py`](compare_public_extended_psf_to_picker_cores.py),
which checks both public-file hashes and the approved local-core summary hash.
Its output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/extended_psf_public_vs_picker_core_2026-10-21_01/summary.json`
with SHA-256
`03e9623c56eec6a2b14cbb7218463a166d0f02f9f1bb60ed53e7c5e333c9f787`.
Five synthetic tests pass. This core mismatch motivates a **hybrid sensitivity
test**, not a correction: preserve the signed local PSF core, append the
published extended wing only beyond a validated overlap radius, and measure
how the inferred tail scale varies over multiple matching annuli and local
frames. The current local 43×43 stamps and sparse outer-star support may not
constrain that transfer tightly enough. Do not construct or apply a hybrid
kernel until the radial overlap and injected-profile recovery pass. Keep
`psf_approved=false` and `science_measurement=false`.

The next Astrafek check measures the radial overlap from 1.5 to 3.3 arcsec
between each of the 744 local models and the public symmetric PSF. It fits a
single local/public amplitude in six annuli while keeping signed model pixels.
The median fitted scale is 0.739 (5th–95th percentiles 0.467–1.045), with a
median weighted shape residual of 12.6% (5th–95th percentiles 8.2–22.8%). The
frame-filter groups remain separate: `i` has 158 positions and median scale
0.685; `i2` has 586 and median scale 0.748. This gives a data-driven bracket
for a possible local-core/public-wing splice, but does not validate the outer
wing amplitude past 3.3 arcsec. The matching factor varies by position and the
profiles are not identical, so use the scale distribution and fit residuals as
sensitivity terms rather than applying one global factor.

The run is
[`audit_public_psf_picker_core_overlap.py`](audit_public_psf_picker_core_overlap.py)
and its synthetic tests pass. The 744-position summary is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/extended_psf_local_overlap_2026-10-21_02/summary.json`
with SHA-256
`32e24b98d32eddbf026786c558f48f63dcdc435ec70b289136b2e5e22f885dd1`.
Next, bracket the splice radius and amplitude with injected point sources and
diffuse profiles in the real target/control masks. Keep the observed-star
comparison and sky-transfer term independent; do not call this a calibrated
PSF or an ICL detection.

## Coverage of target and matched-control centers

The approved 744-position PSF Picker batch covers the 56 sampled target-frame
centers in its original 57-frame injection set (the 57th center is the
separately received test PSF) and 688 satellite-injection positions. A local
coordinate-only cross-match against the full manifests finds 800 unique target
centers and 798 unique matched-control centers. Of their 1,598 unique combined
centers, 56 target centers already have an exact approved PSF coordinate;
no matched-control center does. Thus 1,542 exact center positions remain
uncovered. The 688 satellite-injection PSFs are useful for the existing
injection sample but do not provide local center PSFs for the rest of the
800-pair sample.

An offline request plan is prepared at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_target_control_centers_plan_2026-10-21_01/`.
It lists only those 1,542 missing coordinates for centered `pdr3_wide`
coadd `HSC-I` model PSFs. At the observed batch median of 20,160 bytes per
stamp, the estimated total is 31,086,720 bytes (29.6 MiB); the observed mean
implies 28,908,023 bytes (27.6 MiB). The coordinate-list SHA-256 is
`78bb8b426ad3aa99df68d757b8981e0893517c9e246667aa15d0c891d071d99c`, and the
position-mapping SHA-256 is
`4c2f060f350b65fc93dbb0a8b7c12421517b9be54ddba09c866b08a401dc0212`.
The plan is reproducible with
[`prepare_hsc_psf_picker_target_control_plan.py`](prepare_hsc_psf_picker_target_control_plan.py).
The script has no network path; it checks manifest counts, unique IDs, valid
coordinates, and the received PSF file count/bytes against the saved audit
before writing a plan. Its two focused tests pass, and a run against the real
manifests reproduced the coordinate list and mapping byte-for-byte. This plan
has not been submitted: its scope is separate from the previous 744-position
approval. It requests no images, shape rows, or photo-z values and makes no
science measurement. Even if received, the core models would support local
core matching only; they would not establish extended wings or approve an ICL
correction.

## Exploratory HSC-I proxy sensitivity on `i2` frames

The injection runner now keeps `i` and `i2` selections separate. The HSC-I
public extended model is allowed on `i2` only as a labeled proxy sensitivity;
the output explicitly says this is not an `i2` PSF calibration. This expands
the model-dependence check to the 44 `i2` frames with approved position-core
mapping while retaining the original 13-frame `i` path.

A single-frame pilot on cluster 10 completed in 4.08 seconds using the
symmetric public reference and eight satellite injection centers. Recomputed
responses are 0.255 for diffuse injection into the ICL annulus, 0.757 for all
injected satellite envelopes into that annulus, 1.096 for diffuse injection
into the satellite estimator, and 0.997 for satellite injection into the
satellite estimator. This is only one frame, so it measures neither a stable
filter response nor population completeness. The output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_extended_psf_i2_proxy_pilot_2026-10-21_01/`;
its summary SHA-256 is
`acd7fada5bc1caf1e3beec8bc8da112a18e7b68ac434043057837076b6431e21`, and the
independent response-arithmetic validation SHA-256 is
`6546b76c35e592e41585ad9b6834bc6861d3ab348ba933691e3055827ba6a4bb`.

The 44-frame symmetric and asymmetric proxy runs are defined by
[`run_extended_psf_i2_proxy_sensitivity.sbatch`](run_extended_psf_i2_proxy_sensitivity.sbatch)
and completed sequentially as Slurm job 430903 in 9 min 17 sec; peak RSS was
1,459,504 KiB. Both batches have 44/44 completed frames and zero errors. The
independent validator recomputed all stored transfer responses with maximum
arithmetic error 0 for both variants. The symmetric-model descriptive
responses are:

| Injected component → measured estimator | Median | 16th–84th percentile | Frame range |
|---|---:|---:|---:|
| Smooth diffuse → ICL 200–300 kpc | 0.224 | 0.169–0.255 | 0.004–0.318 |
| Satellite envelopes → ICL 200–300 kpc | 0.485 | 0.198–0.855 | −0.047–1.424 |
| Smooth diffuse → satellite 50–75 kpc | 0.017 | −0.184–0.174 | −2.443–1.179 |
| Satellite envelopes → satellite 50–75 kpc | 0.401 | −0.546–1.111 | −5.933–10.728 |

The broad frame-to-frame ranges, including negative responses in some
measurements, show strong mask/sample dependence and prohibit applying these
numbers as scalar corrections. The output's summary SHA-256 is
`5f2e86d6c3499148b1dbb3151e261daf827de5c434b2b32da7c9af1be846c59d`; the
independent 44-frame arithmetic validation SHA-256 is
`31b8443bbc46ff8c57c711eeb1ac886f21ca63e0056631f535c002a05dd8e5f5`.

The independent
[`validate_extended_psf_injection_batch.py`](validate_extended_psf_injection_batch.py)
checks complete frame support, checksums, filter/proxy labels, science flags,
and recomputes the saved transfer-response arithmetic without rereading the
large source FITS files. Its focused tests pass. The paired
[`compare_extended_psf_injection_variants.py`](compare_extended_psf_injection_variants.py)
comparison confirms the same 44 frames, positions, and injection amplitudes.
The asymmetric-minus-symmetric median differences are below `4×10⁻⁹` for the
ICL responses and zero at printed precision for the satellite responses; the
largest absolute difference is `9.54×10⁻⁷` for satellite injection into the
satellite estimator. Thus the two public variants are indistinguishable for
these radial metrics at the current numerical precision. This does not
validate either model against local `i2` wings. The asymmetric summary
SHA-256 is
`03bdf8f12fde0f9e84ab28ae6c39819f4129564aab0dd4e0c3bb3c1e0340194e`, its
independent validation SHA-256 is
`2e7c32c6cb0e68a6aa69f25d272d42c4f1d683d1a204ee48b15ed1e5e9feea2b`, and
the paired comparison SHA-256 is
`17ab4d0db5ac73ce4b772a1ef0dc176f2796d0bdecf6f4f4b114188cec70526d`.
All results remain
`science_measurement=false`; the proxy cannot establish local `i2` wings, sky
transfer, or an ICL correction.

## Sources and provenance

- Garate-Nuñez et al., *The Hyper Suprime-Cam extended Point Spread Functions
  and applications*, [MNRAS 531, 2517](https://academic.oup.com/mnras/article/531/2/2517/7676869),
  [arXiv:2309.16244](https://arxiv.org/abs/2309.16244).
- Authors' [HSC_PSFs README](https://github.com/luciagarate/HSC_PSFs) describes
  the FITS dimensions, symmetric/asymmetric variants, and `coadd/bg` suitability.
- Zenodo record [10.5281/zenodo.14915578](https://zenodo.org/records/14915578),
  version `v4.0.0`; record metadata reports CC BY 4.0 and the two file checksums
  above.

The HSC PDR3 processing notes that `i` and `i2` transmission curves differ but
are combined at coaddition, with object-specific effective-response corrections
available for photometry ([official processing documentation](https://hsc-release.mtk.nao.ac.jp/doc/index.php/processing__pdr3/)).
The project currently retains `psf_approved=false` and
`science_measurement=false`.
