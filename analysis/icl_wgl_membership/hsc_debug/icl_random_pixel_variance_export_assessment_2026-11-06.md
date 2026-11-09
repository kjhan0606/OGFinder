# HSC random-point variance export audit (2026-11-06)

## Approved scope and transfer

The exact approved plan SHA-256 was
`eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52`. It queried
the HSC PDR3 Wide `random` table at 800 target and 798 matched-control centers,
within the 200–300 proper-kpc annulus. The export contains random-point
coordinates/IDs, patch identity, `i_pix_variance`, patch sky statistics, input
count, local PSF moments, and selected pixel flags. It returned 40,048 rows
(estimated 46,424) in 16/16 completed batches, totaling 11,138,863 CSV bytes.
No image, shape, photo-z, object-catalog, or credential data were returned;
credentials were not saved.

## Coverage and local conditions

The independent post-download auditor verified every batch hash, row contract,
center identity, and 200–300 kpc angular boundary. It found 40,048 unique
random-point rows at 1,587/1,598 requested centers: 791/800 targets and
796/798 controls. Both endpoints have a valid random-point `i_pix_variance`
summary for 787/798 matched pairs. The point-count distribution over all
planned centers is q16/median/q84 = 11/22/39 points per center.

Across the center-level summaries, the median of the per-center median
`i_pix_variance` is 0.003649 (q16–q84: 0.002738–0.005267 in the catalog's
native units). Corresponding median per-center `i_sky_std` is 2.384, median
`i_inputcount_value` is 7, and median local PSF `shape11` is 0.06255. Null
counts are zero for `i_pix_variance`, `i_sky_std`, and `i_inputcount_value` in
returned rows.

These are sparse random-point local statistics, not target-minus-control ICL
fluxes. They do not provide a coadd weight image, pixel-to-frame/amplifier
mapping, exact resampling covariance, source-Poisson noise, or sky-transfer
covariance. No variance correction or science promotion follows; the result
remains `science_measurement=false`.

## Reproducibility

- Approved plan SHA-256: `eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52`.
- Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_pixel_variance_results_2026-11-06_01/`.
- Run checkpoint SHA-256: `c6d34ba4b3798ea55a26ef7736db32c2be70ae10e571454c909e55ada117fd12`.
- Upstream summary SHA-256: `0546a7a289bee5858324c01060f99c63b4e50e444014953606de137dd74daa08`.
- Independent audit: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_pixel_variance_audit_2026-11-06_01/summary.json`, SHA-256 `96245d565375635a994c93ec3294fc581d007161c4479232d5091bbd3cf52389`.
- Submitter SHA-256: `17f1de9d1f0ca7ce347e46df821dc558964ecdb1568c5dc4e50764131502e4ac`.
- Auditor SHA-256: `05352d72efbd95d97ce540e5008970b4ff75fcb7fe121efb467b59e913bb80cc`.
