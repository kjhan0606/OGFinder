# HSC sky-subtraction product provenance for the ICL stack — 2026-11-03

## Full-sample audit

Astrafek checked the existing target and matched-control download manifests and
download-status records for all 1,598 images. All 800 target images and all 798
matched controls are declared `sky=global-sky`, `product=coadd/bg`, with reviewed
product provenance. The 1,598 completed image paths exist and their on-disk byte
sizes match the recorded download sizes. This metadata pass did not reread pixels or
recompute the 90.61 GB of image-content hashes.

The target collection has 309 `i` and 491 `i2` coadds; controls have 312 `i` and
486 `i2`. This is consistent with the existing 783 same-filter target/control pairs
(302 `i`, 481 `i2`), leaving 15 filter-mismatched controls and two targets without a
control.

The reusable audit is
[`audit_icl_sky_product_provenance.py`](audit_icl_sky_product_provenance.py); three
focused tests cover accepted global-sky products, rejection of local-sky products,
and file-size mismatches. The result is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_sky_product_provenance_2026-11-03_01/summary.json`
(SHA-256 `235c4f5cedd9506eb874c9cfe0ac09507a71b1a5f729b1816779b63da4f2c839`). The
audit used no network and did not download data.

## Scientific interpretation

The HSC PDR3 documentation distinguishes the global-sky patch coadd under
`deepCoadd/` from the local-sky `calexp` under `deepCoadd-results/`; the latter is
used for source detection and measurements. HSC's cutout manual maps the saved
`coadd/bg` product to the global-sky coadd before aggressive background subtraction.
Thus the current ICL images use the intended extended-light product, rather than the
locally subtracted measurement coadd. See the primary HSC [PDR3 available-data
description](https://hsc-release.mtk.nao.ac.jp/doc/index.php/available-data__pdr3/),
[PDR3 processing description](https://hsc-release.mtk.nao.ac.jp/doc/index.php/processing__pdr3/),
and the archived product interpretation in the target manifests.

“Global sky” does not mean that the images are sky-free. PDR3 still subtracts
exposure-level background models and a residual local-sky pattern before coaddition.
The product audit therefore closes the product-selection question, but it does not
measure loss or residual transfer of diffuse cluster light in upstream processing.
The existing 600–900 kpc plane subtraction, PSF-wing, matched-control, and injection
tests remain separate diagnostics and do not recover that transfer function.

Keep `science_measurement=false`. The next sky gate is to quantify the spatially
correlated target-minus-control residual on the existing global-sky images and
propagate its uncertainty through the fixed ICL and satellite weights. Do not use
the local-sky `calexp` as the ICL science image or interpret the present provenance
audit as a sky-transfer calibration.
