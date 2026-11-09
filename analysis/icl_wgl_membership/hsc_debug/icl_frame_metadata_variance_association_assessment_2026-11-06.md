# Frame metadata and delivered coadd-variance balance (2026-11-06)

## Question and scope

Using only the already approved HSC frame-metadata export and existing
target/control coadd variance products, this offline analysis asks whether
patch-level frame metadata differences predict the delivered target/control
variance ratio in the 200–300 proper-kpc ICL annulus.

The join retains 774 supported same-filter pairs: 299 `i` and 475 `i2`. Each
pair has more than the existing 50-pixel variance floor at both endpoints.
The variance outcome is `log(target/control median delivered VARIANCE)`.
Five exploratory predictors are paired log differences in unique frame count,
unique visit count, median sky RMS converted by median FLUXMAG0, median
FLUXMAG0, and the geometric mean of the four median amplifier gains. Every
association adjusts for target redshift and log richness.

The analysis forms a shared 5-degree connected-block partition from both
endpoints of all 774 pairs and uses 5,000 common block-bootstrap resamples.
The joint family contains all ten filter-by-predictor associations. A centered
single-step max-T interval provides an approximate familywise sensitivity
check; it is not a validated HSC-field or finite-sample guarantee.

## Results

No metadata association excludes zero after the joint max-T adjustment.
The closest result is the `i`-band frame-count difference (`r=-0.147`,
simultaneous interval `[-0.312, 0.017]`, adjusted `p=0.088`). In `i2`, the
largest point estimate is the FLUXMAG0 difference (`r=0.182`, simultaneous
interval `[-0.247, 0.611]`, adjusted `p=0.750`). All remaining adjusted p-values
are at least 0.60. This provides no familywise evidence that the exported
coarse frame summaries explain the observed coadd-variance imbalance.

The connected partition has 38 blocks and its largest block contains 106 of
the 774 pairs. The broad simultaneous intervals and approximate spatial
partition limit the strength of a null conclusion. The earlier pointwise
`i` frame-count interval excludes zero, but its simultaneous interval includes
zero; it is not an independent detection.

## Limits and interpretation

The frame table identifies contributing frame/visit/CCD memberships at patch
level. It does not map each science-image pixel to its contributing frames or
amplifiers. Unique frame counts are therefore coarse patch summaries, not
effective exposure counts at the measured pixel. The sum of exposure times
was deliberately not treated as a pixel-level depth measure because rows can
repeat exposure metadata across CCDs. The normalized sky RMS is also a proxy,
not a coadd noise estimate. The delivered VARIANCE plane omits source-Poisson
and correlated-noise validation.

HSC PDR3 also provides a coadd `nImage` layer recording the number of visits
contributing to each pixel ([official product description](https://hsc-release.mtk.nao.ac.jp/doc/index.php/available-data__pdr3/)).
It could map spatial coverage in a small, separately approved pilot, but visit
counts are not inverse-variance weights and cannot replace the pending random
point variance audit or a source-Poisson model.

This result does not correct the ICL or satellite profiles, does not explain
the ICL–satellite association, and does not change `science_measurement=false`.
The next useful noise input remains the separately prepared random-point
pixel-variance export; it has not been submitted. Per-pixel warp weights or
equivalent coadd provenance would still be required for a source-Poisson model.

## Reproducibility

- Astrafek runner: [`analyze_icl_frame_metadata_variance_association.py`](analyze_icl_frame_metadata_variance_association.py).
- Result directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_variance_association_2026-11-06_02/`.
- Summary SHA-256: `bf387b4f3ae5b7af6ecf51bc2bcb47daac16623ec6bdd353bdc013c06dea6aee`.
- Joined-pair table SHA-256: `fccb97357222abcda48eb9de3dd024adf3c8bffd975dd49397e04bc3d2c60666`.
- Runner SHA-256: `6ce31dd4de10c2f4d0a70ed86851d19dce92f770d5e7fb96d4f1fb15f8ded5a1`.
- The summary records the approved metadata-audit, aggregate-table,
  variance-summary, and pair-comparison input hashes, plus the script hash.
- No HSC query or download was made; `science_measurement=false` remains.
