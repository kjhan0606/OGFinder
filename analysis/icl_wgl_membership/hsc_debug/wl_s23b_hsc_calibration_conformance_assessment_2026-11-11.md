# HSC S19A estimator calibration-conformance audit — 2026-11-11

## Scope

This was a static review of Astrafek's S19A pair estimator and the source-selection correction, against the official HSC PDR3 shape-catalog documentation and the calibration utility linked from that documentation. No source CSVs were reread, no estimator was rerun, and no query or download was made.

## Field and equation mapping

The HSC documentation maps the exported fields to sky-coordinate ellipticity, per-component shape error and RMS, inverse-variance weight, per-object multiplicative and additive bias, resolution R2, 1-arcsec aperture magnitude, HSC Y3 tomographic bin, B-mode mask, and PSF second moments. Astrafek uses those same fields in s19a_shape_stack.py.

For each radial bin, the implementation uses the HSC responsivity

R = 1 - sum(w_i beta_i^2 rms_e_i^2) / sum(w_i beta_i^2),

where w_i is the catalog shape weight and beta_i is the pair-averaged inverse critical surface density. It applies per-object m, c1, and c2 corrections, the radial-bin edge-selection correction, and then stacks the tangential and cross components with calibrated inverse-variance weights. This matches the S19A estimator structure described by the HSC catalog documentation.

The edge-selection helper was compared term by term with get_sel_bias in the HSC-linked hsc-y3-calib utility:

- The magnitude-edge density uses mag_A10 >= 25.5 - 0.025 and width 0.025.
- The resolution-edge density uses R2 <= 0.3 + 0.01 and width 0.01.
- m_sel = -0.059 p_mag + 0.019 p_R2.
- a_sel = 0.0064 p_mag + 0.0063 p_R2.
- The additive selection term is formed from a_sel times the weighted mean PSF ellipticity derived from the PSF second moments.

Astrafek uses the same thresholds, coefficients, quoted coefficient uncertainties, PSF-moment conversion, and lensing weights for these corrections. No formula discrepancy was found.

The completed source-row QA also shows that the retained sample does not sit on the strict/inclusive boundaries: maximum ellipticity amplitude is below 2, maximum sigma_e is below 0.4, minimum R2 is above 0.3, and maximum mag_A10 is below 25.5. Thus the boundary convention does not change this export's retained rows.

## Remaining calibration and science gates

This checks implementation conformity only. It does not validate the source catalog, the survey footprint, or the physical lensing result. Keep science_measurement=false.

The HSC documentation warns that photo-z bias and source dilution need correction. Astrafek uses the pinned HSC Y3 tomographic n(z) for sample-level lensing efficiency, but has not propagated its posterior uncertainty or applied an independent photo-z bias calibration. The released edge-selection coefficient uncertainties are also not propagated. Reduced shear, magnification, member dilution/boost, target-control depth/mask/PSF/sky matching, and validated contiguous-footprint covariance remain open.

The next spatial-covariance step remains exact HSC tract-polygon assignment after the four-file geometry download receives separate approval.

## Provenance

- Local implementation: s19a_shape_stack.py.
- Source-row QA assessment: wl_s23b_wl_common_shape_row_qa_assessment_2026-11-11.md.
- No data, query, download, or new numeric result was produced by this review.

## Primary sources

- [Official HSC S19A Shape Catalog (PDR3) documentation](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/).
- [HSC-linked S19A/Y3 calibration utility, utils.py](https://github.com/mr-superonion/hsc-y3-calib/blob/main/utils.py).
