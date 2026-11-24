# S23B member-mask component pilot — 2026-11-17

A deterministic pilot reprocessed two i and two i2 target frames using the
previously approved images. It instrumented the original mask function during
the usual source-aware measurement and separated the quality/invalid base,
neighbor segmentation, neighbor Kron/fallback source shapes, and circles
around unmatched candidates. The base includes the selected HSC quality bits,
bright-object bit, nonfinite pixels and nonpositive variance. It does not
separate those causes individually.

For every stamp, the union of the four reconstructed masks was exactly equal
to the original mask. The 20 retained profile rows for the selected frames
reproduced the saved full-sample fluxes and coverage to an absolute tolerance
of 1e-9. The pilot recorded 297 target/control stamps. No network request or
data download occurred; `science_measurement=false`.

## Pixel-weighted area in the 50–75 kpc annulus

Percentages are pooled pixel counts divided by the sum of **in-image** annular
pixels across stamps, not an average over clusters. The same image contributes
multiple correlated target and rotated-control stamps. Components overlap, so
their footprints must not be added. Exclusive area means covered by that one
component and none of the other three; it is not a causal estimate of the
scientific benefit of removing that mask.

| Filter / role | Stamps | Unmasked | Base footprint | Neighbor-shape footprint | Base exclusive | Neighbor-shape exclusive |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| i target | 17 | 35.58% | 45.43% | 40.99% | 21.57% | 13.01% |
| i controls | 136 | 21.53% | 66.80% | 32.35% | 44.85% | 7.90% |
| i2 target | 16 | 43.30% | 29.35% | 34.06% | 21.92% | 20.82% |
| i2 controls | 128 | 25.70% | 56.19% | 29.57% | 43.49% | 12.27% |

The base and neighbor shapes both cover substantial annular area in this
pilot, and the base's exclusive footprint is larger in rotated controls than
in target stamps. This small, deterministic pilot does **not** estimate the
full 979-frame population, prove why the quality/invalid base is large, or
justify relaxing any mask. Segmentation and unmatched-circle counts, including
overlap and exclusive counts, are preserved in the output table. In-image
fractions do not include missing area beyond the image boundary.

## Reproduction

Runner: `analysis/icl_wgl_membership/hsc_debug/audit_s23b_member_mask_components_pilot.py`.
Output root:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_mask_components_pilot_2026-11-17_01/`.

- `summary.json` SHA-256: `f83d1919f3c7607296e62e31bd243f786afa9e0e53ef585f60c367c2295bdf7e`
- `stamp_component_counts.tsv` SHA-256: `c0bcd9081b1a88de46d2abd96ca6d52f4765d39e40c096398d0c1a6cb150dee3`
- The deterministic sample IDs and script/input hashes are in the summary.

Next, separate base mask bits from nonfinite/nonpositive-variance pixels,
quantify boundary loss, and examine a larger predeclared sample with cluster
as the uncertainty unit. Keep current masks and all PSF, sky, Poisson,
coadd-covariance and membership gates until validated.
