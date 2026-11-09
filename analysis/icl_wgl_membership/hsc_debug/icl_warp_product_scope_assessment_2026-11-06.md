# HSC warp-product scope assessment — 2026-11-06

## Finding

The approved CCD-frame metadata export is complete: 67,071 rows for 1,598
centers (800 targets, 798 controls), with no image, shape, object-catalog, or
photo-z rows. It supplies per-frame exposure, amplifier gains, sky scatter,
calibration, and visit/CCD identifiers. It does not contain the per-pixel
coadd weights needed to propagate a synthetic source's Poisson variance.

The authenticated PDR3 DAS Cutout manual cached at
`/gpfs/kjhan/Astrafek/data/logs/cutout_global_manual.html` documents a `warp`
file type. It returns a tar archive of the single-visit warped CCD products
contributing to a coadd. Each FITS cutout has fixed HDUs for image, mask, and
variance; the request can select which layers are returned. This is the most
direct next product for recovering local input variance/weight behavior, but
the variance planes and frame metadata still require a calibrated mapping to
the synthetic-source Poisson term. Their existence alone does not pass that
gate.

## Scope counted from the approved metadata

An offline deduplication of the 67,071 returned rows found:

- 1,662 unique `(tract, patch, filter)` tuples;
- 12,839 unique `(tract, patch, filter, visit)` tuples;
- 56,731 unique `(tract, patch, filter, visit, CCD)` tuples;
- 15,194 center-specific `(center, tract, patch, filter, visit)` tuples.

These are metadata-derived counts, not an HSC request count or a byte estimate.
The cutout endpoint returns all warp products intersecting a requested
rectangle, so a full-sample request could return overlapping products many
times. Do not launch a full-sample warp transfer until a deduplicated request
manifest and measured pilot sizes are available.

## Bounded pilot proposed for separate approval

One 2′ × 2′ `pdr3_wide`, `HSC-I`, `type=warp` cutout at a median-coverage
target center:

- center ID: `target:6966`
- RA, Dec: `346.7915030000°, +1.0399770000°`
- tract/patch: `9703/306`
- metadata footprint: 9 distinct visits in this patch/filter
- requested layers: image, mask, variance
- semi-width/semi-height: `1arcmin` each, giving a 2′ × 2′ tangent-plane size;
- expected uncompressed planes: 55,212,300 bytes (52.65 MiB) if all nine visits
  fill the square at 0.168 arcsec/pixel with 32-bit image, mask and variance
  planes; the downloader stops at a 70 MiB ceiling. Actual tar size depends on
  overlap, FITS headers and file representation.

The pilot would verify archive contents, visit correspondence, FITS headers,
variance/mask alignment, and actual byte size. It would not be used to estimate
an ICL signal. If successful, construct a hash-bound deduplicated manifest and
present the full transfer scope and storage estimate for separate approval.

The offline, immutable one-position request plan is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_warp_pilot_plan_2026-11-06_02/plan.json`
(SHA-256 `a436df245a33bda09cd5599e564b2b9302e5d0d4c006b59fcf6fab1c2c33a65e`).
The Astrafek planner and approval-gated downloader are
`hsc_debug/prepare_icl_warp_pilot_plan.py` and
`hsc_debug/download_icl_warp_pilot.py`. The latter is plan-only by default,
checks the exact plan hash and 70 MiB cap before submission, and reuses the
existing authenticated PDR3 cutout helper. This plan has not been submitted;
no warp pixels have been requested or downloaded.

- Planner SHA-256: `5151f003a0c01bce3128e84a3d9d266cff1670825418fc57732bdfc90673f04a`.
- Approval-gated downloader SHA-256: `766dbedb2b76951307b6c10a9d4be309f96ba1f509d6e09073c0716a6c71e291`.

The downloader's offline FITS audit also records finite image pixels, mask
occupancy, finite/positive variance fractions, per-HDU WCS alignment, whether
the requested sky center falls inside each image, and whether the WCS pixel
scale is consistent with 0.168 arcsec/pixel (5% tolerance). A response with
the right visit set but failed layout/WCS/variance checks is labeled
`complete_with_qa_warnings`; this gate does not promote `science_measurement`.

## Science status

Current coadd cutouts already have image, mask, and variance planes. The
amplitude-grid tests show a stable relative transfer response but only about
23% absolute diffuse-light response; they do not justify a flux correction.
Source Poisson noise, per-pixel coadd weights, upstream sky transfer, extended
PSF wings, and validated covariance remain unresolved. Keep
`science_measurement=false`.

## References

- [HSC PDR3 DAS Cutout manual](https://hsc-release.mtk.nao.ac.jp/das_cutout/pdr3/doc/)
- [HSC PDR3 available data](https://hsc-release.mtk.nao.ac.jp/doc/index.php/available-data__pdr3/)
- [HSC PDR3 processing](https://hsc-release.mtk.nao.ac.jp/doc/index.php/processing__pdr3/)
