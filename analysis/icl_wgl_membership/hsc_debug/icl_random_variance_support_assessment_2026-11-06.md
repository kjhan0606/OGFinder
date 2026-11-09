# Existing HSC random-point support and next ICL noise-data plan (2026-11-06)

## Existing-data audit

`audit_icl_existing_random_noise_support.py` examined the 160,000 HSC PDR3
random-point rows already saved for the 800 target centers. The original query
selected 200 random points per center from a 1,800-arcsec cone. The HSC random
table provides `pix_variance`, `sky_std`, `inputcount_value`, and local PSF
moments; the official documentation describes `pix_variance` as pixel variance
at a random point and notes that an object's reported variance sums over its
footprint ([HSC PDR3 random points](https://hsc-release.mtk.nao.ac.jp/doc/index.php/random-points__pdr3/)).

We converted the saved points to proper cluster-centric radius with the sample
cosmology (`H0=70 km s^-1 Mpc^-1`, `Omega_m=0.3`) and each target's spectroscopic
redshift. The existing sparse sample has virtually no support in the required
ICL 200--300 kpc annulus:

| Proper annulus | Points across 800 targets | Centers with any point | Centers with >=10 points |
|---|---:|---:|---:|
| 50--75 kpc | 3 | 3 | 0 |
| 75--150 kpc | 23 | 21 | 0 |
| 150--200 kpc | 28 | 27 | 0 |
| 200--300 kpc | 109 | 98 | 0 |
| 300--500 kpc | 306 | 245 | 0 |

The 200--300 kpc sampled-point `i_pix_variance` median is 0.00355 in the
catalog's native units, but only 109 positions contribute and no target has ten
points. The random-point rows cover targets only, with no corresponding rows
for the 798 matched controls. Thus they cannot provide target/control-matched
variance or a robust per-center annular noise estimate.

## Prepared next query

Astrafek now prepares a narrowly scoped PDR3 `pdr3_wide.random` query for the
200--300 proper-kpc annulus around 800 targets and 798 matched controls. It
uses the target redshift to set the same physical annulus at both endpoints,
and `adjust_density < 0.2` to retain a nominal 20 random points per square
arcminute. The planned fields are random-point coordinates/IDs, patch identity,
`i_pix_variance`, patch sky statistics, input count, local PSF moments, and
selected pixel flags. Images, galaxy shape rows, photo-z values, and galaxy
object-catalog rows are excluded.

The prepared plan estimates 46,424 rows (roughly 10 MB) in 16 SQL batches. It
has not been submitted and still requires explicit approval of its exact plan
hash. Even a successful export would be a local noise/coverage diagnostic, not
a full coadd weight map or complete source-Poisson/correlated-noise correction.
It does not provide sufficient sampling for the much narrower 50--75 kpc
satellite annulus.

## Reproducibility and state

- Existing-data audit runner:
  `audit_icl_existing_random_noise_support.py`, SHA-256
  `0e2b87449a5557ae94d6cdc7b78c9776e56d6918977d91f3d6bcdc38a7be7b73`.
- Audit summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_existing_random_noise_support_2026-11-06_02/summary.json`,
  SHA-256 `42218211153128870a1e387447a9d7705e552a575bf955e69afc6c2eba0a7c73`.
- Per-center radial support:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_existing_random_noise_support_2026-11-06_02/random_point_radial_support_by_cluster.tsv`.
- SQL-plan builder:
  `prepare_icl_random_pixel_variance_plan.py`, SHA-256
- Approval-gated submitter:
  `submit_icl_random_pixel_variance.py`, SHA-256
  `17f1de9d1f0ca7ce347e46df821dc558964ecdb1568c5dc4e50764131502e4ac`.
- Prepared plan:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_pixel_variance_plan_2026-11-06_03/plan.json`,
  SHA-256 `eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52`.
- SQL aggregate SHA-256:
  `5d82df638435ee39029dc51acdcf580ab6e242d21562bc5523ef9c807368f585`.
- No HSC request or image pixel read was made by the offline audit or plan builder.
- `science_measurement=false` remains in force.

An initial audit output directory ending in `_01` was invalidated because its
angular-to-proper-radius conversion divided by 60 twice. The corrected runner
uses the `_02` output above; do not use `_01` results.
