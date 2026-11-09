# Variance-adjusted ICL–satellite association sensitivity (2026-11-06)

## Purpose and scope

This offline diagnostic asks whether the delivered HSC coadd VARIANCE-plane
balance changes the saved ICL–satellite association in the fixed supported
sample. It joins the 94 supported matched target/control pairs to the
full-sample radial variance table at 200–300 proper kpc, checking pair IDs,
filter, control ID, and agreement of the saved ICL, satellite, and redshift
values. The sample has 26 `i` and 68 `i2` pairs.

For each filter, the outcome association is the partial Pearson correlation
between target-minus-control ICL surface brightness at 200–300 kpc and the
satellite 50–75 kpc contrast. Three models control (1) redshift and richness,
(2) those covariates plus log target/control delivered variance ratio, and
(3) those covariates plus variance ratio and mean log target/control variance
as a coarse pair-depth proxy. The models use identical 5,000 resamples of
connected 5-degree target/control sky blocks, so the change from baseline is
paired. No image pixels were read and no HSC query was submitted.

## Results

| Filter | Model | Partial r | Paired spatial-block 95% interval | Change from baseline | Paired change 95% interval |
|---|---|---:|---:|---:|---:|
| `i` (N=26) | redshift + richness | 0.079 | −0.191–0.500 | — | — |
| `i` | + log variance ratio | 0.059 | −0.277–0.482 | −0.020 | −0.192–0.036 |
| `i` | + variance ratio + mean depth | 0.075 | −0.225–0.506 | −0.004 | −0.182–0.200 |
| `i2` (N=68) | redshift + richness | 0.293 | 0.043–0.498 | — | — |
| `i2` | + log variance ratio | 0.285 | 0.031–0.499 | −0.008 | −0.078–0.034 |
| `i2` | + variance ratio + mean depth | 0.260 | −0.018–0.493 | −0.032 | −0.148–0.027 |

For context, the median 200–300 kpc target/control variance ratio in this
supported subset is 1.024 in `i` and 0.982 in `i2`. It is closer to one than
the full `i2` sample median of 0.961. Adding variance balance shifts the
`i2` point estimate by only −0.008; adding the coarse mean-depth term shifts
it by −0.032. Both paired change intervals include zero. The wider
uncertainty in `i` also includes zero for all models. This sensitivity test
does not show that the full-sample `i2` delivered-variance difference
explains the supported-sample ICL–satellite association.

## Interpretation limits

These are selected-support, low-power exploratory correlations, not a
causal adjustment or a detection. The coadd VARIANCE plane is not a calibrated
total-noise estimate: the effective-gain/source-Poisson contribution and
per-pixel coadd weights are unavailable. Mean log variance is only a coarse
depth proxy. The analysis does not estimate the covariance of the masked ICL
estimator and does not address sky transfer, extended PSF wings, completeness,
or astrophysical covariance. The connected 5-degree components are
approximate blocks rather than validated HSC fields. Keep
`science_measurement=false`.

## Reproduction and integrity

Runner: [`analyze_icl_variance_adjusted_satellite_assoc.py`](analyze_icl_variance_adjusted_satellite_assoc.py),
SHA-256 `a51e03db90fcbf1eed7f37e3c44adf2e3713433cf5c2c3eee4a1b62b2e976f09`.

The runner pins the 94-pair input SHA-256
`79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da` and the
full-sample comparison SHA-256
`db77c9a32797aafa79efaac7f7e18d7e14a8f79c96b18e8fe296aceb6f35830d`.
Results are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_variance_adjusted_satellite_assoc_2026-11-06_01/`:
`summary.json` SHA-256
`8731105a3d5a9c0430ae2f688442beea640277d06af7d70ff03b1f48e9eb0644`, and
`joined_supported_pairs.tsv` SHA-256
`a2eb2d628912dbb029c1116d8eb0819d1ec0f0658b489be9690d53cfb9eb4b34`.
The 5,000 joint replicates were valid for every filter/model. No new HSC data
request or download occurred.
