# S23B CAMIRA member-file raw-field schema consistency audit

## Published catalog context

The S23B CAMIRA paper reports 12,162 Wide clusters with the bright-star mask
and more than 600,000 member galaxies with membership weight
`w_mem >= 10^-3`. It defines richness through stellar-mass and spatial
filters plus a compensated background filter; richness is not expected to
equal an unqualified sum of exported member weights. The paper also says no
`w_mem` cutoff is applied when calculating richness. See [Oguri et al. 2026,
Sections 2.3 and 3.1](https://arxiv.org/html/2512.13954v1#S2.SS3) (especially
the richness/filter description and member-weight threshold). The [official
S23B catalog README](https://raw.githubusercontent.com/oguri/cluster_catalogs/main/hsc_s23b_camira/README.md)
lists the cluster catalog files but does not document the headerless member
file's column order.

## Offline file audit

`audit_s23b_camira_member_raw_field_schema.py` rejoined all rows of the
already approved member file to the pinned cluster catalog. No network request,
new download, or image-pixel read occurred.

- All **637,146 rows** map uniquely to the **12,162** clusters in the pinned
  S23B Wide star-mask catalog.
- Raw field 8 ranges from **0.001000069 to 1.0**. Every row is at or above
  `10^-3`, matching the published member-catalog threshold. Together with the
  published member count and field's range, this is strong evidence that
  field 8 is `w_mem`.
- Raw field 7 ranges from **9.984 to 12.546**, with median **10.701**. This
  resembles a log stellar-mass field, but this audit does not establish its
  meaning.
- Summing field 8 does not reproduce catalog richness within rounding for
  12,159 of 12,162 clusters. That is consistent with the published use of
  stellar-mass and spatial filters and compensated background subtraction;
  it is not evidence against a `w_mem` interpretation.

The strongest supported interpretation is **field 8 is probably CAMIRA
`w_mem`**. The file is headerless and the official public README does not state
its column order, so keep the source field as `field8_raw` in data products
until a file-specific schema or author source confirms the mapping. The
previous raw-field-8 profile weighting stays a sensitivity analysis and is not
upgraded to a science detection.

## Products

Runner:
`analysis/icl_wgl_membership/hsc_debug/audit_s23b_camira_member_raw_field_schema.py`.

Output directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_schema_audit_2026-11-07_04/`.

- `summary.json` SHA-256: `f57de9630c101d5eaf686b29d4a6a5c8a73de6b54028b2bc444ca4f202c06364`
- `cluster_field_sums.tsv` SHA-256: `a5dc8fa57047f00874dc5fee3cdd5fc760e7d66f4e73832cc61f2adc1e5c3609`

Keep `science_measurement=false`. No commit or push has been made.
