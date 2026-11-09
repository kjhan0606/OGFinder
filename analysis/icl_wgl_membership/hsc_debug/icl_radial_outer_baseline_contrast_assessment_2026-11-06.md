# ICL radial contrast relative to the 300–500 kpc profile

## Result

The existing HSC-I matched-control stack measures six annuli from 50 to 500
proper kpc. The exploratory sample requires at least 40% common target/control
coverage in every annulus with a 2.5 Kron mask. The fixed complete-case support
contains 20 `i` pairs in 10 connected sky blocks and 48 `i2` pairs in 30 blocks.

The primary contrast is the target-minus-control surface brightness at 200–300
kpc minus the area-weighted target-minus-control surface brightness at 300–500
kpc. The outer-annulus weights are 0.4375 for 300–400 kpc and 0.5625 for
400–500 kpc. The weights follow the projected annular areas.

| Filter | Extraction | Pairs / blocks | ΔSB 200–300 | Outer weighted ΔSB | Contrast | 95% block-bootstrap interval |
|---|---:|---:|---:|---:|---:|---:|
| `i` | raw | 20 / 10 | 3.620 | 0.840 | 2.780 | [0.825, 5.115] |
| `i` | outer-plane corrected | 20 / 10 | 4.123 | 1.629 | 2.494 | [0.640, 4.774] |
| `i2` | raw | 48 / 30 | 3.644 | 1.284 | 2.359 | [−0.110, 4.219] |
| `i2` | outer-plane corrected | 48 / 30 | 3.151 | 0.992 | 2.159 | [−0.296, 4.047] |

Surface-brightness units are nJy arcsec⁻². Each interval comes from 20,000
resamples of connected paired target/control 5-degree sky blocks. The bootstrap
resamples blocks and weights each draw by its number of cluster pairs.

## Interpretation

The measured profile declines between the 200–300 and 300–500 proper-kpc
annuli in both filters. The raw `i` and outer-plane-corrected `i` intervals
exclude zero. Both `i2` intervals include zero. The `i` and `i2` samples are
small and strongly selected by full-aperture coverage. The four rows are not
independent tests because raw and corrected measurements share the same images,
masks, cluster pairs, and outer-region sky fits.

The 300–500 kpc annuli can contain genuine extended cluster light. The outer
profile is not blank sky, so the contrast measures radial concentration rather
than an absolute ICL component. The exploratory 40% coverage floor is not an
adopted science threshold. The profile does not separate BCG light from diffuse
intracluster light or establish an ICL–satellite causal relation.

The intervals describe block-sampling variability in the selected pairs.
Source completeness, source-Poisson variance, exact coadd covariance, upstream
sky transfer, and extended PSF wings remain unresolved. No measurement
correction was applied. `science_measurement=false` remains in force.

## Reproducibility

Runner: [`analyze_icl_radial_outer_baseline_contrast.py`](analyze_icl_radial_outer_baseline_contrast.py)
(SHA-256 `33fac24a5ca352e43214ed2e89f4cf937d0c1a99751a263cc0738f4a4e46a991`).

Output directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_radial_outer_baseline_contrast_2026-11-06_01/`.

- `summary.json` SHA-256 `96382058878820f1d9882745b7029aae36cb96fe6a6474d96373cc2b407f2eb2`
- `outer_baseline_contrast.tsv` SHA-256 `86fe3ae589a9bffc4f1a449c2d5a5b10efdb124b3d88e6973dde4dc9cd3a4f72`
- `block_bootstrap_outer_contrast_draws.npz` SHA-256 `b51cbd130f12a784d4a9ebef6c1564cf41cc8ddd01603ee99d6e2b997b7a1a5e`

The runner verified the archived target profiles, 798 control baselines, saved
sample indices, annular means, and input hashes. The runner read profile
archives only. No image pixels or HSC network data were accessed.
