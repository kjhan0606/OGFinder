# ICL–satellite association: connected HSC-region influence — 2026-10-20

## Scope and validation

This leave-one-region check uses the previously validated fixed 94 cluster/
filter rows: 26 in `i` and 68 in `i2`. It assigns both the cluster target and
its matched-control coordinate to connected components of HSC patch-QA tract
centers, then omits a pair if either endpoint belongs to the removed region.
After each deletion, it recomputes the partial Pearson correlation controlling
for redshift and richness, plus the unadjusted Spearman correlation.

The 840-tract graph has four components, three occupied by the 800-cluster
manifest (499, 246, and 55 clusters). The target/control endpoint pair counts
and grouping are unchanged for tract-center adjacency thresholds from 1.6 to
3.0 degrees. The full-sample partial Pearson values exactly reproduce the
archived 5-degree baseline to `1e-12`: `r=0.29264` for `i2` and `r=0.07898`
for `i`. The fixed sample, selected radii, and mask factor are unchanged.

## Results

| Filter | Full pairs | Full partial Pearson `r` | Full Spearman `rho` | Omitted region | Pairs omitted | Partial `r` after deletion | Spearman `rho` after deletion |
|---|---:|---:|---:|---:|---:|---:|---:|
| `i2` | 68 | 0.293 | 0.111 | 0 | 34 | 0.504 | 0.361 |
| `i2` | 68 | 0.293 | 0.111 | 1 | 29 | 0.048 | −0.048 |
| `i2` | 68 | 0.293 | 0.111 | 2 | 5 | 0.286 | 0.082 |
| `i` | 26 | 0.079 | 0.275 | 0 | 8 | 0.178 | 0.397 |
| `i` | 26 | 0.079 | 0.275 | 1 | 18 | −0.378 | −0.214 |

Removing region 1 reduces the `i2` partial correlation from `0.293` to `0.048`
and changes Spearman `rho` from `+0.111` to `−0.048`. That deletion removes
29 of the 68 `i2` rows (43%). The relation is therefore materially sensitive
to which connected HSC region is included. The `i` deletion that changes both
statistics' signs leaves only eight pairs, so it is too small to interpret as
evidence for a negative population relation.

## Interpretation and limits

This check shows that the positive `i2` partial correlation is not robust to
the connected-region deletion, even though earlier 2.5–10-degree rectangular
cell deletions left its partial-Pearson sign positive. The coarser HSC tract
components expose dependence on a large geographic region, while the unadjusted
rank statistic also changes sign after that region is removed. The result
does not identify whether this dependence comes from astrophysical variation,
selection, sky/PSF transfer, or another regional systematic.

Three occupied regions are insufficient for a covariance or interval
estimate. These are influence diagnostics only, not a significance test or
population detection. The association remains conditional on the selected
94-row sample, limited ICL pixel support, and unresolved image-level
systematics. The weak-lensing random subtraction and paired field covariance
are also still missing. `science_measurement=false` remains in force.

## Reproduction

- Astrafek script: [`audit_satellite_icl_tract_region_influence.py`](audit_satellite_icl_tract_region_influence.py)
- Tests: [`test_audit_satellite_icl_tract_region_influence.py`](test_audit_satellite_icl_tract_region_influence.py), 3 passed
- Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_tract_region_influence_2026-10-20_02/summary.json`
- Output SHA-256: `7acc3063ddd9ebae30d28c1856cc3bc985d2a5fa4bc399e4a05fc537295586dc`
- Runner SHA-256: `ad1c1e755995d53ec3d4a2ab25f9d06521beaffc2ae3aba74ee4c9e45fd9d0e7`
- Provenance: exact association, ICL annulus, extraction, and baseline checksums are recorded in the output JSON.
