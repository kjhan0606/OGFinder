# ICL 200–300 kpc aperture covariance diagnostic (2026-11-03)

## Result

The offline diagnostic contracted the saved symmetric 9×9 coadd lag kernels
with each frame's actual valid-pixel weights in the 200–300 kpc ICL annulus.
It processed all 224 supported target/control pairs (71 `i`, 153 `i2`), or 448
frames. Reconstructed aperture coverage and mean surface brightness were
checked against the saved profiles before aggregation. The coadds were not
modified, no HSC query was made, and no covariance correction was applied.

The finite-lag kernel raises the estimated independent-frame pixel-noise
stack sigma from 0.06274 to 0.07869 nJy arcsec⁻² in `i` (factor 1.254), and
from 0.06006 to 0.07356 nJy arcsec⁻² in `i2` (factor 1.225). Across frames,
the target sigma-inflation median is 1.2588 in `i` (q16–q84 1.2328–1.2851)
and 1.2576 in `i2` (1.2244–1.2857). Matched controls have medians 1.2546
(`i`, 1.2179–1.2753) and 1.2441 (`i2`, 1.2104–1.2757).

These finite-lag pixel-noise terms are 7.56% (`i`) and 6.94% (`i2`) of the
half-widths of the existing spatial-block-bootstrap intervals. They provide
an additional approximate pixel-noise floor, not a replacement for those
intervals or a total uncertainty estimate. The similarity between target and
control inflation medians also gives no basis for a target-only scalar
correction.

## Limits and interpretation

The 9×9 kernels are medians over sampled patches and assume local stationarity
within each aperture. Lags beyond four pixels, upstream sky transfer, source
Poisson noise, extended PSF wings, mask uncertainty, and spatial covariance
between frames are omitted. Target and control pixel-noise contributions are
combined as independent. Larger-scale field covariance remains represented
only by the existing spatial block bootstrap. This is an exploratory
diagnostic; retain `science_measurement=false` and do not correct the measured
ICL profile from these factors.

## Reproducibility

Run `analyze_icl_annulus_aperture_covariance.py` against the saved annulus
profiles, masks, and independently audited coadd lag kernels. Output:

`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_aperture_covariance_200_300_2026-11-03_01/`

The summary SHA-256 is
`82b0f2ee2d364b204552933f60dfcc0853e5c9869d014957a07801dfe090fb03`.
The 448-frame table SHA-256 is
`d7632bf7b43a978d31045a1bb7919bddb7237580b7664e51243f598b7188605a`; the
224-pair table SHA-256 is
`b1b084de4c638da6d35b009834f48185913535c4a7ee616c91f6c3f97d9b47c6c7`.
The independent output-hash manifest agrees with all three values.
