# S23B ICL Galactic-extinction sensitivity — 2026-11-07

## Scope and method

This analysis applies separate nominal HSC `i` and `i2` SFD corrections to
the saved target and matched-control annular means, then re-forms their
target-minus-control profiles on the existing 1,022-pair common support. The
conversion is `f = 10**(0.4 k E(B−V))`, using the hscPipe coefficients
`k_i=1.641` and `k_i2=1.617` listed in the [hscPipe 7 manual](https://hsc.mtk.nao.ac.jp/pipedoc/pipedoc_7_e/colorterms.html).
These are treated as nominal filter coefficients for a sensitivity, not as a
new validation of the PDR3 coadd calibration.

Within each pair, 25-kpc bins are combined using their common valid-pixel
counts. Pair outcomes are then averaged equally. Uncertainty uses 5,000 draws
of the connected target/control sky blocks in six 2.5°, 5°, and 10° layouts
with zero/half-cell phases. The RA origin and block counts reproduce the
pinned spatial-block audit. The uncorrected relative contrasts reproduce the
previous S23B contrast table exactly in all 180 checked values (36 rows × 5
summary fields). No image pixels or network data were used.

## Primary result

Entries are equal-pair means in nJy arcsec⁻², with 95% connected-block
bootstrap intervals. The contrast is the inner annulus minus 300–800 kpc.

| Filter | Inner radius (kpc) | Pairs / blocks | Raw | SFD-corrected | Correction |
|---|---:|---:|---:|---:|---:|
| `i` | 50–300 | 232 / 26 | 6.561 [3.844, 8.667] | 6.930 [4.043, 9.168] | +0.369 [0.190, 0.516] |
| `i` | 100–300 | 232 / 26 | 4.528 [1.901, 6.537] | 4.774 [1.984, 6.918] | +0.246 [0.078, 0.384] |
| `i` | 200–300 | 231 / 26 | 3.253 [0.498, 5.441] | 3.452 [0.526, 5.777] | +0.199 [0.017, 0.345] |
| `i2` | 50–300 | 756 / 39 | 7.611 [5.722, 8.877] | 8.055 [6.117, 9.385] | +0.444 [0.337, 0.569] |
| `i2` | 100–300 | 756 / 39 | 4.686 [3.239, 5.707] | 4.965 [3.460, 6.058] | +0.279 [0.197, 0.400] |
| `i2` | 200–300 | 752 / 39 | 1.908 [0.401, 3.035] | 2.024 [0.450, 3.217] | +0.116 [0.029, 0.215] |

For the 50–300 kpc contrast, the nominal extinction correction raises the
primary estimate by 5.6% in `i` and 5.8% in `i2`. Its paired block interval
stays positive in all six spatial layouts. The 100–300 kpc correction also
stays positive in all six layouts. For 200–300 kpc, at least one layout's
correction interval crosses zero in each filter. Median absolute target/control
SFD mismatch is 0.00071/0.00073 mag in `i`/`i2` (95th percentile
0.00264/0.00317 mag).

## Interpretation and limits

The matched-control stack's relative 50–300 kpc excess changes modestly but
measurably under this nominal correction. This does not make the stack a
paper-ready ICL measurement. Multiplication also rescales non-astrophysical
residual sky in the delivered coadds; the 300–800 kpc reference can contain
cluster light; and this sensitivity does not address PSF wings, upstream sky
transfer, source-Poisson noise, coadd covariance, mask/center review, or member
completeness. Coefficient and dust-map systematics are not in the block
intervals. Keep `science_measurement=false` and `images_ready_for_science=false`.

The separate S23B CAMIRA member catalog has not been downloaded; direct
member-centered ICL/satellite analysis still awaits its separate explicit
download approval. The existing weak-lensing estimate also remains
exploratory because its spatial support occupies only three connected regions.

## Reproducibility

Runner: `analyze_s23b_icl_extinction_sensitivity.py` (SHA-256
`746228c31f8e80c8b183ef7b50e91b356afec4463001126c12bf547176cc7b2c`).
It pins the target/control manifests, stack summary, pair-profile table, and
spatial-block summary by hash, then asserts the RA origin, block counts, and
raw-contrast reproduction before writing output.

Output directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_extinction_sensitivity_2026-11-07_04/`.

- Summary SHA-256: `5dfb730f26043f2073790d21326156fd7e3d79730f457e7c0bc738745b32e8fe`.
- Block-bootstrap table SHA-256: `003d3c574867e60ec3d056bbe41fc255ba455bf7ccb2253315b747723097b0b2`.
- Pair table SHA-256: `caf73e45bd03ead10e07afb4b658c759c9d468bba4d7f5c8e3293cf10190c232`.
- Figure SHA-256: `6a2bbec11aa01b6a5b0ca19c6abab856c2db50d7c44c1fd8de23fdf9c6e259c4`.

No HSC request, download, image-pixel read, commit, or push was performed.
