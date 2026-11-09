# ICL-supported PSF redetection injection assessment — 2026-11-06

## Scope and validation

The explicitly approved HSC PDR3 Wide centered HSC-I coadd PSF request returned
all 2,005 requested files (38,088,000 FITS bytes). The download auditor matched
all coordinates and validated finite, near-unit-sum PSF cores. The combined
bank audit then verified 224 PSF-complete cluster/filter frames, 2,693 satellite
injection roles, and 2,917 unique PSF files, with request-approval hash
`1d7c0189721600fc31e8524a4258dc70a5c0f0d436df45b57b9a0c3dd0058b65`.

The runner selected 137 of those 224 pairs after applying the existing
200–300 kpc ICL common-coverage threshold of 40%: 44 `i` and 93 `i2`. The other
87 pairs were below that measured profile-support threshold. For each selected
pair the run generated baseline, ICL-only, satellite-only, and combined
scenarios, convolved injections with the position-specific HSC model-PSF core,
and reran SEP detection and masks under the conservative `mask_all` unresolved
BCG policy. The run completed all 137 pairs with zero errors in 671.5 s. The
original FITS images and their variance planes were not modified.

An independent audit checked each frame JSON checksum, selected-pair coverage,
scenario support bookkeeping, code provenance, and the exact PSF filenames and
hashes against the combined approved bank. It verified 137/137 frame results.
The auditor reports that at least one satellite-profile support gap occurs in
79/137 frames. Gaps by scenario are baseline 74, ICL-only 73, satellite-only 76,
and both 75. These gaps are retained as missing support rather than converted
to zero flux.

## Response summary

Responses below are per-unit injected amplitude. For the diffuse component, one
unit is 1 nJy arcsec⁻² at 100 kpc with a 100 kpc exponential scale. Satellite
envelopes use one unit per galaxy at 50 kpc with a 15 kpc scale. Values are
frame-level injected-minus-baseline measurements summarized by median and the
16th/84th percentiles; they are not calibration factors.

| Measured response | Filter | Valid frames | Median | 16th–84th percentile |
|---|---|---:|---:|---:|
| ICL 200–300 kpc ← diffuse injection | `i` | 44 | 0.226 | 0.178–0.249 |
| ICL 200–300 kpc ← diffuse injection | `i2` | 93 | 0.226 | 0.194–0.257 |
| Plane-corrected ICL ← diffuse injection | `i` | 44 | 0.223 | 0.164–0.255 |
| Plane-corrected ICL ← diffuse injection | `i2` | 93 | 0.227 | 0.182–0.257 |
| ICL 200–300 kpc ← all satellite envelopes | `i` | 44 | 0.342 | 0.146–0.685 |
| ICL 200–300 kpc ← all satellite envelopes | `i2` | 93 | 0.449 | 0.182–0.838 |
| Satellite profile ← smooth diffuse injection | `i` | 17 | −0.031 | −0.096–0.067 |
| Satellite profile ← smooth diffuse injection | `i2` | 45 | 0.019 | −0.181–0.236 |
| Satellite profile ← satellite injection | `i` | 14 | 0.349 | −0.306–0.955 |
| Satellite profile ← satellite injection | `i2` | 44 | 0.409 | −0.384–1.263 |

## Assessment

Under this particular detection/masking and compact-PSF setup, the measured
200–300 kpc ICL annulus recovers about 23% of the nominal smooth injection in
both filters. This is a useful warning that the image-side response is
substantial and similar medians across filters do not remove it. Satellite
envelopes also contribute to the measured ICL annulus, with median responses of
0.34 (`i`) and 0.45 (`i2`) per unit input amplitude.

The satellite-profile response is not stable enough to calibrate: support gaps
reduce the valid counts to 14/44 in `i` and 44/93 in `i2`, and both percentile
ranges span negative values and values near or above unity. The near-zero
diffuse-to-satellite medians also have broad intervals crossing zero. These
are sensitivity diagnostics for the selected image sample, not evidence for
or against an astrophysical ICL–satellite relation.

This result is not ready to correct the science stack. Injected-source Poisson
noise was not added; coadd variance was held fixed; upstream HSC sky transfer
and extended PSF wings were not modeled; and `mask_all` can remove real BCG
light when the central detection is ambiguous. The full 224-pair bank and 137
pair injection subset also remain selected by measured coverage. Keep
`science_measurement=false`, `observational_validation_passed=false`, and
`flux_covariance_available=false`.

## Reproduction and artifacts

- Runner: [`run_joint_icl_satellite_batch.py`](run_joint_icl_satellite_batch.py).
- Combined-bank verifier: [`verify_supported_pairs_position_psf_bank.py`](verify_supported_pairs_position_psf_bank.py).
- Batch auditor: [`audit_icl_supported_pairs_psf_redetection_batch.py`](audit_icl_supported_pairs_psf_redetection_batch.py), SHA-256 `0868c6f189e4c188af8e755676e6d0a836ccd7127e1a704640bcc07cf3ab8300`.
- Download audit: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_satellite_support224_download_2026-11-06_01/psf_picker_audit.json`, SHA-256 `79f559c0c85720a7281cf743083517dffe8d15d78af38ace538d72431d843f54`.
- Combined PSF bank audit: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported_pairs_position_psf_bank_audit_2026-11-06_01/summary.json`, SHA-256 `9c0002344fc7642a5a224ef108b209bb65ca468948ed0309420e8e12e4fac5ac`.
- Injection batch: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/joint_icl_satellite_psf_redetection_icl_supported224_2026-11-06_01/summary.json`, SHA-256 `be3c2e8cb22185c06768d9022ae3e6192072f180fa54d55e21b85920e5572750`.
- Independent audit: the same directory's `psf_redetection_audit.json`, SHA-256 `b91783ee4e1582c24676d9a9bd85aacec2cebafa46d76980e909fcbfecbc0bd1`.

The existing generic 744-position PSF/redetection auditor intentionally only
accepts the original 57-frame role set; this batch has its own auditor so its
expanded PSF manifest is checked without weakening that earlier gate.
