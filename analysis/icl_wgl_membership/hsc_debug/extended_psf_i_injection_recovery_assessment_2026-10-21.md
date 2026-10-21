# HSC `i`-band ICL/satellite injection recovery — 2026-10-21

## Scope

The existing Astrafek injection runner was applied to the 13 `i`-filter
cluster/control frames that pass the K2.5 and 40% ICL 200–300 kpc support cuts
and have approved satellite-position PSF mappings. It injected a smooth
circular exponential profile (100 kpc scale, normalized to 1 nJy/arcsec² at
100 kpc) and detected-member satellite envelopes (15 kpc scale, normalized to
1 nJy/arcsec² at 50 kpc per member). Source detection, the ICL mask, and the
satellite masks were rerun on each synthetic image. The measured estimators
were ICL 200–300 kpc and satellite outer light 50–75 kpc.

Synthetic light was convolved with the public symmetric or asymmetric 4001²
HSC-I reference. These are full-reference model sensitivities in real image
masks, not position-specific PSF calibration and not the proposed local-core /
public-wing splice. The measured science images were not changed. No injected
source shot noise or correlated coadd noise was added.

## Results

Both model variants completed all 13 frames with zero errors. For the symmetric
reference, median transfer responses and frame q16–q84 ranges were:

| Injected component → measured estimator | Median | q16–q84 | Full frame range |
|---|---:|---:|---:|
| Smooth diffuse → ICL 200–300 kpc | 0.215 | 0.161–0.231 | −2.734–0.271 |
| Satellite envelopes → ICL 200–300 kpc | 0.353 | 0.150–0.611 | 0.024–1.269 |
| Smooth diffuse → satellite 50–75 kpc | −0.039 | −0.100–0.134 | −0.915–4.324 |
| Satellite envelopes → satellite 50–75 kpc | 0.294 | −0.234–0.726 | −3.024–1.895 |

The asymmetric-minus-symmetric paired medians were effectively zero. Across the
13 common frames, the largest absolute difference was `4.77×10⁻⁷` for the
satellite-envelope response in the satellite estimator; the largest ICL
estimator difference was `2.27×10⁻⁸`. This numerical agreement only says these
two public templates yield the same aperture responses in this setup. They
share the same source-star sample and neither is matched to each local frame.

The wide frame ranges, negative recovery values, and unstable cross-talk show
that a single transfer coefficient is not supported. The 13-frame subset is
small and support-selected; its quantiles are not population confidence
intervals. These diagnostics do not measure completeness, spatial covariance,
or observational PSF/sky transfer. No ICL–satellite scientific measurement is
approved (`science_measurement=false`, `psf_approved=false`).

### Fixed-baseline-mask sensitivity

A follow-up run remeasured all four injection scenarios with the baseline
detection, ICL mask, and satellite targets held fixed, alongside the original
per-scenario redetection branch. The independent validator recomputed the
saved arithmetic for all 13 frames with maximum error 0. With fixed masks,
diffuse-to-ICL response is tightly clustered (median `0.22495`, q16–q84
`0.22175–0.22717`), while diffuse-to-satellite response is near zero (median
`0.00062`, q16–q84 `−0.01161–0.01260`). Satellite-envelope responses remain
variable: satellite-to-ICL median `0.494` (q16–q84 `0.165–0.713`) and
satellite-to-satellite median `0.278` (`−0.174–0.509`).

The redetection-minus-fixed-mask difference isolates the effect of rerunning
detection and masks in this synthetic setup. Its diffuse-to-ICL median is
`−0.0115`, but the frame range reaches `−2.955`; diffuse-to-satellite has
median `−0.0514` and range `−0.900–4.326`. The extreme recovery values are
therefore strongly affected by reselection in a few frames, rather than
providing a stable PSF transfer coefficient. This remains a model sensitivity
test: it does not add source noise, reproduce coadd covariance, validate sky
transfer, or establish completeness. All scientific approval flags remain
false.

The frame-level influence audit identifies cluster `4982` as the largest
reselection outlier. Its matched satellite support falls from 4 objects / 16
controls at baseline to 2 / 8 after diffuse-only and joint injections, while
ICL 200–300 kpc coverage falls from `0.5395` to `0.4982`. This frame has
redetection-minus-fixed-mask responses `−2.955` for diffuse→ICL and `+4.326`
for diffuse→satellite. Removing it changes the 13-frame diffuse→ICL median
from `0.2150` to `0.2181` (fixed-mask: `0.22495` to `0.22522`), so that
median is not driven by this single extreme. The satellite-envelope response
median changes from `0.294` to `0.218` under the same leave-one-frame-out
check, illustrating substantial sample influence. These leave-one-out shifts
are descriptive for the selected 13 frames, not population confidence
intervals. The reusable audit is
[`audit_extended_psf_injection_reselection_influence.py`](audit_extended_psf_injection_reselection_influence.py);
its saved record also verifies per-frame hashes and response arithmetic.

## Reproduction and provenance

The run used [`extended_psf_injection.py`](extended_psf_injection.py). The
filter-aware paired comparator is
[`compare_extended_psf_injection_variants.py`](compare_extended_psf_injection_variants.py);
it now accepts both `i` and `i2` while preserving the existing `i2` default.
Eleven focused tests covering the comparator, validator, and injection model
passed; syntax and whitespace checks passed.

Outputs are under `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/`:

- Symmetric run: `icl_extended_psf_i_sensitivity_2026-10-21_01/`, 13/13 frames. Summary SHA-256 `e9853c6626cace7b630d6db17e3d9a8587a26007d726b884801608f2719be228`; validation SHA-256 `aa1c1b3fd80f044d76f66a25ca2da19517b90c2e511ef48d3d2c3033b43df2d1`.
- Asymmetric run: `icl_extended_psf_i_asymmetric_sensitivity_2026-10-21_01/`, 13/13 frames. Summary SHA-256 `dd3f40c56c55840a6a268eb9192cfaf7acede0183a138951df978bf63f882527`; validation SHA-256 `b62255417a19e069bbb53879b002bde3932a40d10a98cf2f3f8cba80db382aa8`.
- Paired comparison: symmetric run's `variant_comparison_filter_aware.json`, SHA-256 `9281dca86f14aaa4ce7e63537b021e2ce861e1c3663f5b4a17b7ff86f5fba575`. It records the filter as `i` and uses the filter-neutral paired-sensitivity status; the earlier comparison artifact is retained but superseded.
- Fixed-mask run: `icl_extended_psf_i_fixed_masks_2026-10-21_01/`, 13/13 frames, zero errors. Summary SHA-256 `a732fbc5fdaf41932f5685a329143f8d794f543885b5a5c0917ea321c3c264c5`; validation SHA-256 `2bcd49fa584c18e2415be56ed4fa75fa4894f2762b839ec9cdf689b2147c0b45`.
- Reselection influence audit: `icl_extended_psf_i_reselection_influence_2026-10-21_01/reselection_influence_audit.json`, SHA-256 `09dc009995b8dc8bbbc29cff0b392812f47b5351e399c57c943a26d8a56cee29`; the audit script SHA-256 is `6b9f2a1e9f9f97dcbc1c27a21e1383524a42aed8a023c5a8b5fe086c8129ccad`.

Both independent validators recomputed all saved transfer-response arithmetic
with maximum error 0. The public reference SHA-256 values are
`8897af89767e5fbb8d5f9c70965126bc89409f6c65ecd8d117c7995f34be2d85`
(symmetric) and
`b2ab473c6af4da11c1df74850bdd0a8cc0c1da4d109c118d1f4033b4235027fd`
(asymmetric).
