# HSC warp-variance comparison for the supported 94-pair ICL sample

## Scope and acquisition

The user approved plan SHA-256 `bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a`, covering 94 matched ICL–satellite target/control pairs (26 `i`, 68 `i2`), or 188 fixed 2′×2′ HSC-I image/mask/variance warp cutouts. No catalog, shape, photo-z, other filter, or wider-area image request was included.

The downloader completed 184 archives and hit HTTP 404 twice for each of four request positions, leaving four center-level gaps. It downloaded 9,939,640,320 bytes (9.26 GiB), below the approved 12 GiB hard cap; the pilot-scaled estimate was 9,072,971,094 bytes (8.45 GiB). All 184 returned archives contain every metadata-expected visit. The endpoint returned 161 additional visit files within the same approved cutouts; these were retained and audited. No additional positions or filters were requested.

Four unavailable requests:

- `i2:target:3010`
- `i2:control:control_7281_43773774104897710`
- `i:control:control_418_36411590274189484`
- `i:control:control_6229_41567062857876122`

## Method

`analyze_icl_supported94_warp_noise.py` checks the approved plan and download summary hashes, each archive checksum, expected visit inventory, FITS plane alignment, mask definitions, and WCS. Warp cutouts have fractional-pixel offsets and small scale/rotation differences relative to coadds, so each visit WCS is mapped to its coadd with a least-squares affine fit on a 7×7 grid. Fits with maximum residual above 0.05 coadd pixel are rejected. Warp variance and masks are nearest-neighbor sampled onto coadd pixel centers. For each coadd pixel, the diagnostic predicts variance as `1 / sum(1 / visit_variance)` using finite positive visit variances and the no-data mask. It compares this to the delivered coadd VARIANCE.

The primary support requires at least three valid visit samples and a valid coadd variance pixel. The strict mask sensitivity additionally excludes mask planes shared by all returned warps and the coadd, and requires at least three clean visit samples. Radial bins are in proper kpc at the matched target redshift. Pair differences in the 200–300 kpc annulus use a 5,000-replicate bootstrap that resamples connected 5-degree sky blocks, not pixels.

## Results

- 184/188 centers passed archive, expected-visit, FITS, mask-bit, and affine-WCS checks; the other four are the HTTP 404 gaps above. This leaves 88 complete target/control pairs (24 `i`, 64 `i2`) for paired comparisons.
- Across the 184 centers, the median of per-center pixel-median coadd/predicted variance ratios is 1.00669 in 200–300 kpc. The central 68% range of those center medians is 1.00475–1.01109. These are descriptive center distributions, not confidence intervals.
- The strict common-mask median ratio in 200–300 kpc is 1.00658. Median strict coverage is 81.9% of the full geometric annulus (primary coverage is 100%). At 300–400 kpc, median strict coverage falls to 53.1%; this outer bin should not be treated as uniformly supported.
- In 200–300 kpc, median target-minus-control center ratio is −0.000905 in `i` (24 pairs; block-bootstrap 95% interval −0.00306 to +0.00115) and +0.000559 in `i2` (64 pairs; −0.000177 to +0.00175). Both intervals include zero.

Reproducible output:

- Summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported94_warp_noise_analysis_affine_03_2026-11-06/summary.json`
- Per-center radial table: same directory, `per_center_radial.tsv`
- Summary SHA-256: `991274bb24525442b44d457540eaa7bd0b0e479f540b43e66d58cbd16d3fd031`
- Approved download summary SHA-256: `0039321462f8158d1c006a1df64bb3e96ae988744ed780c845f74a60c79b95b4`
- Analysis runner SHA-256: `36560f729ace577fb1904ec7da52207d67596b91bee09e1d39ce82daeb722590`

## Interpretation and limits

The result supports consistency between the delivered coadd variance plane and a simple inverse-variance reconstruction from the returned visit warps, with a small roughly 0.7% ratio offset. Target and control variance ratios are statistically compatible under the descriptive sky-block resampling. This is a product/variance consistency diagnostic only. Nearest-neighbor sampling does not reproduce HSC's exact resampling kernel or covariance, and the comparison does not measure source-Poisson/effective gain, correlated resampling noise, residual-sky covariance, extended PSF wings, or the total ICL-estimator uncertainty. Pixel quantiles do not treat pixels as independent. This does not establish or refute an ICL–satellite association.

`science_measurement=false` remains. No GPU or Slurm job was used. No commit or push was made.

## Pipeline context

The HSC pipeline describes image coaddition as resampling onto the tract grid and combining exposures; the LSST `ScaleVarianceTask` documentation states that warping/convolution leaves covariance that a variance-plane rescaling does not track explicitly. These references support treating the present nearest-neighbor variance comparison as an intermediate diagnostic and retaining the kernel-aware covariance gate:

- [HSC SSP pipeline overview](https://hsc.mtk.nao.ac.jp/ssp/pipeline/)
- [LSST Science Pipelines: ScaleVarianceTask](https://pipelines.lsst.io/api/lsst.meas.algorithms.ScaleVarianceTask.html)
