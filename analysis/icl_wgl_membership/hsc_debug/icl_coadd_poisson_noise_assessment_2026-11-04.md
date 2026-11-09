# ICL coadd Poisson-noise input audit (2026-11-04)

## Question and method

The existing ICL–satellite injection pipeline adds synthetic flux to an
in-memory science plane but does not change its variance plane. Before adding a
noise term, this audit checked whether the downloaded target and matched-control
coadds contain the gain/exposure metadata needed to translate synthetic flux
into an effective coadd shot-noise variance.

Astrafek tool:
[`audit_icl_poisson_noise_inputs.py`](audit_icl_poisson_noise_inputs.py).
It read the FITS headers for the 800 target images and 798 matched controls.
It did not read pixel arrays, make network requests, or modify any image. The
audit checked the saved byte size and the saved IMAGE/MASK/VARIANCE HDU mapping.
It searched every primary and image-plane header for `GAIN`, `GAIN1`, `GAIN2`,
`EFFGAIN`, `COADDGAIN`, `EXPTIME`, `EXPOSURE`, `EFFEXPT`, `SUMEXPT`, `NEXP`, and
`NINPUT`.

## Result

All 1,598 files exist, match their recorded byte sizes, and have IMAGE, MASK,
and VARIANCE HDUs at the indices recorded in the download reviews. This holds
for all four subsets: target `i` 309/309, target `i2` 491/491, control `i`
312/312, and control `i2` 486/486.

All 1,598 contain a variance plane. None has any of the checked explicit gain
or exposure header keys. They carry `FLUXMAG0`, which sets the photometric
conversion, but that alone does not supply the effective gain for propagating
the Poisson variance of a synthetic source through a multi-visit coadd.
Therefore the existing coadd bundle is insufficient to update the injection
variance with a physically grounded synthetic-source shot-noise term. This is
a limitation of the present inputs, not evidence that source Poisson noise is
negligible.

HSC PDR3 describes the coadd as a combination of overlapping visit images and
lists per-visit `warp` products and a per-pixel visit-count `nImage` product.
Those products suggest a route to investigate effective-noise propagation, but
the visit count alone does not determine the coadd gain or its pixel weights.
The PDR3 database documentation lists `mosaic`, `frame`, and `warped` as image
metadata tables; the Wide `warped` table contains metadata for about 1.58 million
warped images. A follow-up should first map each target/control coadd patch to
its contributing visits through metadata, then determine whether the available
per-visit variance and coadd weights suffice for shot-noise propagation. This
metadata query has not been submitted.
See the HSC [PDR3 processing description](https://hsc-release.mtk.nao.ac.jp/doc/index.php/processing__pdr3/)
and [available data products](https://hsc-release.mtk.nao.ac.jp/doc/index.php/available-data__pdr3/),
and [PDR3 database tables](https://hsc-release.mtk.nao.ac.jp/doc/index.php/database__pdr3/).

## Reproduction and status

Run on `lageunha` from the Astrafek repository:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/audit_icl_poisson_noise_inputs.py
```

The hash-bound output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_poisson_noise_input_audit_2026-11-04_01/summary.json`,
SHA-256 `f75e704de3ae8dbdfc6c2f52700bebb9407e13532d0f1dd10eb6b71990fd533c`.
The tool SHA-256 is
`bd42fc675aaa7c4fd2c33c22d9bbbcb35309e6518fd38351ac301e7770c04e2b`.
Download-status SHA-256 values are embedded in the summary. No new HSC data were
requested or received. The ICL injection model remains noise-incomplete and
`science_measurement=false`; no correction was applied.
