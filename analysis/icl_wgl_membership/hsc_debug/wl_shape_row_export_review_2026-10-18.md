# HSC shape-row export review (2026-10-18)

## Current evidence

The approved count-only preflight completed for all 800 cluster lenses. Its 32
CSV batches have verified SHA-256 hashes, 25 rows per batch, and an exact match
to the 800 unique lens IDs. The sum is 3,069,544 lens-cone source candidates;
168 lenses have zero candidates, and the per-lens median and maximum are 4,438
and 9,767. These are raw shape-catalog counts before redshift-bin, magnitude,
resolution, B-mode, and radial analysis cuts. A source in two cluster cones is
counted twice because the science data unit is a lens-source pair.

The result and audit are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_preflight_2026-10-18_01/cas_counts/`.
Both its checkpoint and summary say `source_rows_downloaded=false`,
`science_measurement=false`, and `credentials_saved=false`.

## Proposed row export

The prepared plan is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_rows_plan_2026-10-18_02/`.
It contains 32 count-matched SQL batches, each covering 25 lenses. The SQL
uses the same PDR3 joins and 3 Mpc cone geometry as the completed count-only
preflight, but returns only these 22 explicit columns:

`lens_id`, lens RA/Dec/redshift, `object_id`, source RA/Dec, e1/e2, shape
uncertainty/RMS/weight, per-object m/c1/c2, 1-arcsec aperture magnitude,
resolution, three PSF second moments, raw B-mode mask value, and HSC Y3 source
tomographic bin.

Each result must have exactly the count-only row total for every lens and no
duplicate `(lens_id, object_id)` pair. Repeated `object_id` values across
different lenses are retained, because each occurrence is a distinct
lens-source pair needed by the stack and its covariance analysis.

The expected total is exactly 3,069,544 pair rows. A rough uncompressed CSV
planning range is 256–480 bytes per row, or **0.73–1.37 GiB**. This is a
serialization estimate, not a measured transfer size; the exact bytes can be
reported only after an approved export. No images or full per-object photo-z
PDFs are in this query.

## Redshift calibration and scientific limits

The HSC S19A shape table provides `hsc_y3_zbin`; zero denotes sources excluded
from the four HSC Y3 tomographic bins (including sources with bimodal P(z)), and
1–4 identify the assigned bins. The shape-stack adapter now accepts those bins
with the published HSC Y3 calibrated source `n(z)` distributions, so the planned
row export does not need per-object full-PDF rows. The HSC Y3 redshift-calibration
analysis combines photometry with clustering-redshift information and treats
lensing weights consistently. The released SACC `wl_0`--`wl_3` tracers consumed
by the adapter provide fiducial `n(z)` curves; they are not posterior samples.
Posterior/systematic uncertainty must be propagated separately through the
`DeltaSigma` profile. HSC's documentation says shape rows
join to positions by `object_id`, and point estimates or full PDFs need
dilution/photo-z correction; this is why calibrated tomographic `n(z)` is the
preferred row-export path.

The full DNNz-Wide PDF tar has HTTP `Content-Length` 144,312,176,640 bytes
(134.4 GiB); this was measured with a HEAD request only. It is intentionally
outside this export plan. The separate `n(z)` calibration product has not been
downloaded or size-checked. The row export alone is therefore not a completed
weak-lensing dataset, and any resulting profile must keep
`science_measurement=false` until the calibrated `n(z)`, source selection,
photo-z/dilution response, random-lens null, and covariance gates are complete.

The official shape-catalog page identifies `b_mode_mask` as the mask for
removing the large-B-mode region, but does not state whether stored 0 or 1
means masked. The query therefore preserves the raw value. The adapter now has
`decode_hsc_b_mode_mask(raw, masked_value=...)`, which requires the stored-value
polarity to be verified explicitly and rejects non-binary values. Do not decode
or stack the exported rows until that polarity is documented or independently
verified from the release schema.

References:

- HSC S19A shape catalog and required columns/calibration notes: <https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/>
- HSC PDR3 photo-z products and per-tract full-PDF description: <https://hsc-release.mtk.nao.ac.jp/doc/index.php/photometric-redshifts__pdr3/>
- HSC Y3 weak-lensing redshift-distribution results: <https://hsc-release.mtk.nao.ac.jp/doc/index.php/wly3/>

## Execution status

`prepare_wl_shape_rows.py` builds the query plan locally and makes no network
request. `submit_wl_shape_rows.py --plan` validates the plan without contacting
CAS. Submission requires `--submit --approved-source-row-export` in an
interactive `lageunha` session. No source-row SQL has been submitted, and no
source row has been downloaded. This document records a request ready for
review; it is not a data-release approval.
