# Independent-support ICL estimator sensitivity — 2026-10-19

## Input and method

This assessment compares the independent-support estimator—target and matched-control means on each frame's own valid pixels—with the strict common-support estimator on identical target/control pairs and identical 25 kpc annuli. The three source runs use mask factors K=2.5, 3.5, and 5.0.

Each input run passed `verify_icl_independent_support.py`. All outputs remain exploratory and keep `science_measurement=false`. The tabulated pair SEM is descriptive; it is not a spatial-block covariance interval.

## 50–200 kpc radial bins

The support-gain columns count pairs with a valid independent target-minus-control mean versus pairs that also have a valid common-support mean. Coverage columns are median geometric pixel fractions across available pairs.

| Mask | Filter | R (kpc) | N independent | N common | Added | Independent Δ ± pair SEM | Independent Δ on common pairs | Common Δ on same pairs | Shift (%) | Median target / control / common coverage |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| K=2.5 | i | 50–75 | 218 | 212 | 6 | 60.41 ± 4.946 | 60.45 | 56.69 | 6.64% | 0.153 / 0.598 / 0.0778 |
| K=2.5 | i | 75–100 | 243 | 234 | 9 | 31.16 ± 2.109 | 30.26 | 32.76 | -7.62% | 0.246 / 0.564 / 0.129 |
| K=2.5 | i | 100–125 | 260 | 248 | 12 | 20.12 ± 1.555 | 19.43 | 19.2 | 1.21% | 0.305 / 0.517 / 0.135 |
| K=2.5 | i | 125–150 | 271 | 260 | 11 | 12.95 ± 1.497 | 12.59 | 15.35 | -18% | 0.325 / 0.492 / 0.141 |
| K=2.5 | i | 150–175 | 278 | 264 | 14 | 9.074 ± 1.148 | 8.932 | 10.21 | -12.5% | 0.317 / 0.464 / 0.132 |
| K=2.5 | i | 175–200 | 284 | 270 | 14 | 7.012 ± 0.9104 | 6.952 | 6.403 | 8.57% | 0.308 / 0.434 / 0.111 |
| K=2.5 | i2 | 50–75 | 347 | 333 | 14 | 62.13 ± 3.147 | 63.25 | 66.67 | -5.13% | 0.189 / 0.63 / 0.102 |
| K=2.5 | i2 | 75–100 | 387 | 376 | 11 | 33 ± 1.8 | 33.92 | 34.38 | -1.34% | 0.323 / 0.589 / 0.155 |
| K=2.5 | i2 | 100–125 | 406 | 394 | 12 | 19.43 ± 1.019 | 19.67 | 17.6 | 11.7% | 0.348 / 0.56 / 0.162 |
| K=2.5 | i2 | 125–150 | 438 | 425 | 13 | 13.44 ± 0.9567 | 13.85 | 13.14 | 5.42% | 0.363 / 0.532 / 0.162 |
| K=2.5 | i2 | 150–175 | 447 | 436 | 11 | 9.462 ± 0.7964 | 9.538 | 8.86 | 7.64% | 0.373 / 0.51 / 0.161 |
| K=2.5 | i2 | 175–200 | 456 | 444 | 12 | 6.518 ± 0.8367 | 6.786 | 7.593 | -10.6% | 0.364 / 0.503 / 0.16 |
| K=3.5 | i | 50–75 | 190 | 172 | 18 | 54.53 ± 4.774 | 55.04 | 52.01 | 5.82% | 0.0332 / 0.495 / 0.00893 |
| K=3.5 | i | 75–100 | 232 | 215 | 17 | 25.51 ± 1.97 | 25.74 | 27.45 | -6.24% | 0.147 / 0.473 / 0.0511 |
| K=3.5 | i | 100–125 | 257 | 242 | 15 | 16.32 ± 1.674 | 15.99 | 15.49 | 3.21% | 0.213 / 0.439 / 0.0776 |
| K=3.5 | i | 125–150 | 269 | 257 | 12 | 11.32 ± 1.492 | 10.89 | 12.01 | -9.27% | 0.235 / 0.421 / 0.076 |
| K=3.5 | i | 150–175 | 277 | 262 | 15 | 7.916 ± 1.176 | 8.125 | 8.993 | -9.65% | 0.247 / 0.388 / 0.0853 |
| K=3.5 | i | 175–200 | 282 | 262 | 20 | 6.216 ± 0.8608 | 6.245 | 4.902 | 27.4% | 0.235 / 0.364 / 0.0695 |
| K=3.5 | i2 | 50–75 | 301 | 277 | 24 | 61.51 ± 4.055 | 63.53 | 64.76 | -1.91% | 0.0485 / 0.515 / 0.0128 |
| K=3.5 | i2 | 75–100 | 371 | 355 | 16 | 29 ± 1.995 | 30.1 | 30 | 0.324% | 0.172 / 0.489 / 0.062 |
| K=3.5 | i2 | 100–125 | 401 | 379 | 22 | 16.47 ± 1.101 | 16.9 | 16.34 | 3.42% | 0.231 / 0.468 / 0.0862 |
| K=3.5 | i2 | 125–150 | 436 | 416 | 20 | 11.63 ± 0.9749 | 12.17 | 10.19 | 19.5% | 0.259 / 0.449 / 0.0933 |
| K=3.5 | i2 | 150–175 | 446 | 430 | 16 | 8.366 ± 0.8411 | 8.549 | 8.713 | -1.88% | 0.276 / 0.429 / 0.0979 |
| K=3.5 | i2 | 175–200 | 456 | 439 | 17 | 5.816 ± 0.8377 | 6.203 | 6.892 | -9.99% | 0.285 / 0.423 / 0.107 |
| K=5 | i | 50–75 | 111 | 92 | 19 | 50.99 ± 8.6 | 55.25 | 52.61 | 5.03% | 0 / 0.33 / 0 |
| K=5 | i | 75–100 | 180 | 149 | 31 | 20.14 ± 2.363 | 19.67 | 18.38 | 7.03% | 0.0237 / 0.304 / 0 |
| K=5 | i | 100–125 | 229 | 186 | 43 | 13.17 ± 1.75 | 12.68 | 12.81 | -0.978% | 0.0718 / 0.293 / 0.0101 |
| K=5 | i | 125–150 | 254 | 215 | 39 | 6.789 ± 1.489 | 7.546 | 8.418 | -10.4% | 0.101 / 0.274 / 0.0173 |
| K=5 | i | 150–175 | 270 | 238 | 32 | 5.425 ± 1.137 | 5.29 | 4.773 | 10.8% | 0.123 / 0.253 / 0.0252 |
| K=5 | i | 175–200 | 280 | 249 | 31 | 5.624 ± 0.9295 | 5.24 | 4.041 | 29.7% | 0.133 / 0.253 / 0.026 |
| K=5 | i2 | 50–75 | 178 | 148 | 30 | 66.41 ± 6.084 | 65.97 | 63.68 | 3.6% | 0 / 0.347 / 0 |
| K=5 | i2 | 75–100 | 289 | 233 | 56 | 25.28 ± 2.267 | 24.5 | 28.59 | -14.3% | 0.0255 / 0.338 / 0 |
| K=5 | i2 | 100–125 | 358 | 312 | 46 | 15.05 ± 1.492 | 15.23 | 17.66 | -13.7% | 0.0747 / 0.328 / 0.0129 |
| K=5 | i2 | 125–150 | 409 | 360 | 49 | 9.217 ± 1.24 | 9.756 | 9.138 | 6.76% | 0.116 / 0.312 / 0.0253 |
| K=5 | i2 | 150–175 | 435 | 401 | 34 | 6.294 ± 0.9406 | 6.517 | 6.451 | 1.03% | 0.139 / 0.298 / 0.031 |
| K=5 | i2 | 175–200 | 451 | 421 | 30 | 4.442 ± 0.8797 | 4.677 | 6.971 | -32.9% | 0.162 / 0.288 / 0.0367 |

## Reading the comparison

- /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_2d_independent_support_k2p5_2026-10-19_02: 783 matched pairs (i=302, i2=481); all source-input paths exist and the saved annular aggregates reproduce from per-pair profiles.
- /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_2d_independent_support_k3p5_2026-10-19_02: 783 matched pairs (i=302, i2=481); all source-input paths exist and the saved annular aggregates reproduce from per-pair profiles.
- /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_2d_independent_support_k5p0_2026-10-19_02: 783 matched pairs (i=302, i2=481); all source-input paths exist and the saved annular aggregates reproduce from per-pair profiles.
- K=2.5, i: maximum absolute same-pair estimator shift over 50–200 kpc is 17.97% across the six annuli.
- K=2.5, i2: maximum absolute same-pair estimator shift over 50–200 kpc is 11.75% across the six annuli.
- K=3.5, i: maximum absolute same-pair estimator shift over 50–200 kpc is 27.39% across the six annuli.
- K=3.5, i2: maximum absolute same-pair estimator shift over 50–200 kpc is 19.46% across the six annuli.
- K=5, i: maximum absolute same-pair estimator shift over 50–200 kpc is 29.67% across the six annuli.
- K=5, i2: maximum absolute same-pair estimator shift over 50–200 kpc is 32.90% across the six annuli.

The independent estimator recovers more pairs or valid pixels only where the target and control have separate support. The shift can therefore reflect both extra area and differing azimuthal/source-mask environments. It is a sensitivity diagnostic, not a correction to the common-support stack.

These runs do not resolve residual sky transfer, position-dependent extended PSF wings, extinction gradients, coadd noise covariance, or satellite-selection completeness. They do not measure ICL–satellite association significance. Keep the common-support estimator as the reference and retain `science_measurement=false` pending those gates and covariance-aware analysis.

## Reproducibility

- Astrafek runner: [`stack_hsc_icl_independent_support.py`](stack_hsc_icl_independent_support.py)
- Independent verifier: [`verify_icl_independent_support.py`](verify_icl_independent_support.py)
- Assessment generator: [`summarize_icl_independent_support.py`](summarize_icl_independent_support.py)
- icl_2d_independent_support_k2p5_2026-10-19_02: `summary.json` SHA-256 `6411eb9851a393ce0426ecd058a13665454ef5f7645ce1b16466a7ea3c9f0370`; verifier SHA-256 `bf8971c9ed18e71b36b35ab4382399ddfe66058d9eecf703e755a6b433d4ae2f`.
- icl_2d_independent_support_k3p5_2026-10-19_02: `summary.json` SHA-256 `39fabfd55e5353c065bc57914470bfb1a2d2b4ee8ef6de8d9a68fc9454eb1ff6`; verifier SHA-256 `4e117897dd11657fc1d2396e62509c6ef203e349566d385364c45706a6d65582`.
- icl_2d_independent_support_k5p0_2026-10-19_02: `summary.json` SHA-256 `0be8094301c0335c7da33dac3e6b85e4535d60fdc055c419e10b5e127893747a`; verifier SHA-256 `c0b95d94754f573c0d243c34b84434eef965db4ed7fb3b2f530b24f595448818`.
