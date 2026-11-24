# S23B member outer-profile retention audit — 2026-11-17

## Scope

This local analysis uses the existing source-aware member stack's selected
frames, candidate diagnostics, retained radial bins, and cluster profiles.
It reads no image pixels and makes no network requests. The five upstream
file hashes are checked against the existing outer-relation runner's pinned
hashes. All 979 selected frames, 16,503 candidate positions and 13,269 retained
candidate/bin rows reconcile. Every retained bin has the original >=60%
target/control-mean coverage and >=3 rotated-control support. Cluster counts
also reproduce the saved cluster-profile table.

The denominator is the already selected frame/candidate population. Retention
is not a detection-completeness estimate for all satellites or all 1,043 targets.
CAMIRA positions remain photometric member candidates.

## Results

At 50–75 proper kpc around each candidate:

| Filter | Clusters retained / selected | Candidate positions retained / selected | No unique SEP match | Matched but no retained outer bin |
| --- | ---: | ---: | ---: | ---: |
| i | 62 / 226 (27.43%) | 160 / 3,504 (4.57%) | 1,281 | 2,063 |
| i2 | 365 / 753 (48.47%) | 1,197 / 12,999 (9.21%) | 4,269 | 7,533 |

The three candidate categories partition each filter's denominator. The saved
records cannot distinguish insufficient target area from rotated-control
failure within the final category, so neither cause is assigned individually.
At 5–10 kpc the cluster retention is 136/226 (60.18%) and 586/753 (77.82%);
profiles at different radii therefore describe different selected subsets.

Retained-minus-lost standardized mean differences use
`sqrt((sample_var_retained + sample_var_lost)/2)` as their denominator:

| Covariate | i | i2 |
| --- | ---: | ---: |
| Redshift | 0.052 | 0.053 |
| log10 CAMIRA richness | 0.336 | 0.072 |
| Relative ICL contrast | −0.161 | −0.030 |
| Number of candidate positions | 0.190 | 0.090 |

The i retained sample has higher mean log richness than the lost sample.
Its lowest and highest within-filter richness quartiles retain 9/57 (15.79%)
and 22/57 (38.60%) clusters. The middle quartiles retain 32.14% and 23.21%,
so this is not a strictly monotonic relation. For i2 the lowest/highest
quartiles retain 44.15%/50.26%. Quartile cut points are computed on each
filter's full selected population, with ties assigned to the upper bin;
quartile sizes can differ. These are descriptive comparisons without spatial
uncertainty estimates or detection significance.

## Decision and next analysis

The present outer-light association remains conditional on substantial
profile-area/control selection, particularly in i. The new audit does not
establish missingness at random or justify extrapolation to the full selected
sample. It does not support inverse-selection correction or imputation.
Keep `science_measurement=false` and retain the previous unresolved physical
association conclusion.

Before rerunning images, add per-bin target coverage, valid-control counts,
and explicit failure reasons to the member-profile measurement diagnostics.
Then test coverage/control-threshold sensitivity on common candidate and
cluster support, preserving separate i/i2 analyses. New PSF requests and
S16A random-center queries still require their separate exact-plan approvals.
Existing PSF/Poisson, sky, coadd covariance, centering and membership gates
remain open. The official [S23B author README](https://github.com/oguri/cluster_catalogs/tree/main/hsc_s23b_camira)
was checked again; it still supplies no file-specific member column dictionary,
so raw fields 7 and 8 retain their existing interpretation limits.

## Reproduce

From `/home/kjhan/BACKUP/Astrafek/OGFinder`:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/audit_s23b_member_profile_retention.py
```

The runner refuses to overwrite an existing output. Use `--output` with a new
path for another run. The default system Python lacks matplotlib; the existing
Miniconda environment was used without installing packages.

Outputs are outside Git at:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_profile_retention_2026-11-17_01/`.

- `cluster_radial_retention.tsv`: full frame/bin attrition ledger.
- `radial_retention.tsv`: filter/bin retention totals.
- `outer_retention_balance.tsv`: retained/lost covariate balance.
- `outer_retention_quartiles.tsv`: covariate-quartile conditional retention.
- `retention.png`: retention and descriptive balance plot.
- `summary.json`: input/output/script hashes, scope and limitations.

Summary SHA-256:
`1a12dd5bdde8adf4a67a9c72f0a8352d6ae074fbadc47d12be905b0cdbcc4c0d`.
Runner SHA-256:
`cecd3ddd61da47f0a06ea63b28a394321fea464b96272e4cab9e6c3dd0a00a76`.
All recorded output hashes and candidate attrition partitions were separately
checked after writing. The plot was visually inspected. This run adds no CAS
submission, data transfer or scientific detection claim.
