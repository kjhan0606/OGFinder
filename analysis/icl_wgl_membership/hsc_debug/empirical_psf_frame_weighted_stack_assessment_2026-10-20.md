# Frame-weighted stack of existing empirical-star profiles — 2026-10-20

## Purpose and method

The existing Gaia-supported diagnostic contains 20 complete 401×401 HSC-I star
stamps in 15 target frames. The earlier per-star diagnostics had mixed signs in
the outer profile. This Astrafek run tests whether an equal-frame radial stack
changes that result without requesting new data.

[`stack_empirical_psf_star_profiles.py`](stack_empirical_psf_star_profiles.py)
checks the ledger and every stamp SHA-256, then measures the two existing source
masks (Kron 2.5 and 3.5) and two existing equal-azimuth local-sky annuli
(120–150 and 160–190 pixels). It subtracts the measured background, normalizes
each radial-bin mean by the measured flux inside 12 pixels, averages stars
within each HSC frame, and gives frames equal weight. The uncertainty band is a
deterministic 10,000-draw bootstrap of frames; it does not treat stars in one
frame as independent. Native radial bins, sectors, coverage rules, and the
background estimator match the original diagnostic.

Reproduce on `lageunha` with an unused output path:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/stack_empirical_psf_star_profiles.py \
  --input /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_2026-10-14_01 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_frame_weighted_stack_NEW
```

The script has two focused tests for frame weighting/reproducibility and signed
profile retention. Both passed.

## Results

The run read 20 hash-matched stamps from 15 frames and produced 880 star-level
method/radial-bin rows plus 44 frame-weighted stack rows. The normalized
profile is a signed mean surface brightness divided by the inner-core flux, in
inverse-pixel units; it is not a normalized PSF encircled-energy curve.

| Mask / sky annulus | 48–64 px: frames; mean [16th, 84th bootstrap] | 64–96 px: frames; mean [16th, 84th bootstrap] |
|---|---:|---:|
| Kron 2.5 / 120–150 px | 8; −2.48×10⁻⁶ [−3.35, −1.61]×10⁻⁶ | 9; −2.02×10⁻⁶ [−2.47, −1.57]×10⁻⁶ |
| Kron 2.5 / 160–190 px | 5; −1.80×10⁻⁶ [−2.54, −1.16]×10⁻⁶ | 6; −1.62×10⁻⁶ [−2.03, −1.21]×10⁻⁶ |
| Kron 3.5 / 120–150 px | 2; −2.08×10⁻⁶ [−2.21, −1.95]×10⁻⁶ | 5; −1.40×10⁻⁶ [−2.10, −0.69]×10⁻⁶ |
| Kron 3.5 / 160–190 px | 2; −1.84×10⁻⁶ [−1.94, −1.73]×10⁻⁶ | 4; −1.51×10⁻⁶ [−2.11, −0.90]×10⁻⁶ |

The consistently negative net profile beyond roughly 48 pixels (about 8 arcsec
at 0.168 arcsec/pixel) shows that these local-background-subtracted stamps do
not provide a positive extended-wing measurement. The apparent agreement
between sky-annulus choices does not distinguish PSF wings from a sky offset:
both estimates subtract a local level that may contain real wings, cluster
light, undetected sources, or residual background structure. The small frame
counts for Kron 3.5 also make its bootstrap bands descriptive only.

## Interpretation and limits

This stack does not approve an empirical extended PSF, a wing correction, or an
ICL surface-brightness correction. The 15 frames are a sparse, cluster-field
subset rather than a representative sample of the 800 coadds; star variability,
Gaia match identity at the exact coadd epoch, HSC coadd covariance, and the effective
background transfer remain unresolved. The output explicitly retains
`science_measurement=false` and `psf_approved=false`.

Output directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_frame_weighted_stack_2026-10-20_02/`.
Summary SHA-256: `274aed73330cd64b9787662ac61af071949bb4841be6345dc0f529f0dd6406c5`.
The profile table, stack table, and plot SHA-256 values are in `summary.json`.
The implementation SHA-256 is `558af1f3f37f9e446e16aeb5734430db5057eb6437d42c38436a8bf2f441b07e`.

## Expanded Gaia-supported sample update (2026-10-20)

The expanded Gaia-supported ledger was subsequently used to test whether the
small initial result persisted. It contains 480 candidates, of which 333 have
complete hash-verified stamps from 227 distinct frames. The parent 800-frame
download manifest confirms all 800 science cutouts are HSC-I. This larger
diagnostic stack contains 14,652 radial-profile rows and 44 frame-weighted
radial-bin combinations. A deterministic 10,000-draw bootstrap resamples frames.

| Mask / sky annulus | 48–64 px: frames; mean [16th, 84th bootstrap] | 64–96 px: frames; mean [16th, 84th bootstrap] |
|---|---:|---:|
| Kron 2.5 / 120–150 px | 83; −5.20×10⁻⁷ [−9.75, −0.47]×10⁻⁷ | 90; −0.77×10⁻⁷ [−3.02, +1.48]×10⁻⁷ |
| Kron 2.5 / 160–190 px | 66; −9.93×10⁻⁷ [−17.28, −2.94]×10⁻⁷ | 70; −6.96×10⁻⁷ [−9.17, −4.79]×10⁻⁷ |
| Kron 3.5 / 120–150 px | 32; −9.90×10⁻⁷ [−20.53, +0.88]×10⁻⁷ | 48; +2.99×10⁻⁷ [−0.53, +6.55]×10⁻⁷ |
| Kron 3.5 / 160–190 px | 26; −18.70×10⁻⁷ [−36.30, −0.84]×10⁻⁷ | 30; −6.02×10⁻⁷ [−10.22, −1.83]×10⁻⁷ |

The larger sample removes the earlier impression that all outer profiles are
strictly negative: two of the four 64–96 pixel estimates are compatible with
zero, and one has a positive mean. The other two remain negative. The sign
depends on the selected local-sky annulus and source-mask radius, so these data
still do not isolate extended PSF light from local sky, ICL, or contaminating
sources. This is not a positive-wing detection and does not support an ICL
correction. Exact Gaia-to-HSC cross-match timing and stellar variability checks
remain open. HSC `X_IMAGE/Y_IMAGE` already centers each extracted stamp, so Gaia
propagation is not required for image-space centering. `science_measurement=false`
and `psf_approved=false`.

The authoritative full-sample diagnostic is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_frame_weighted_full_gaia_2026-10-20_02/`.
Its summary SHA-256 is
`0be06f1a827823bcd04ce23db4113dca1bd368e1a2b7205d0bab0ed76911a7ba`; profile,
stack, and plot checksums are in that summary. It supersedes the 20-stamp result
for interpretation while preserving that earlier run unchanged.

후속 header-only 점검에서 Gaia 우표를 포함하는 227개 HSC target coadd 어디에도 표준 DATE/MJD epoch 카드가 없음을 확인했다. 따라서 ±2/5/10년 분석은 실제 coadd epoch 보정이 아니다. 기여 exposure/visit 시각이 확보되기 전까지 경험 프로필을 ICL PSF 보정으로 사용하지 않는다. 재현기·행별 헤더 digest·해시는 [HSC coadd epoch 감사](empirical_psf_hsc_coadd_epoch_audit_2026-10-20.md)에 기록했다.
