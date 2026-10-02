# Extended morphology (plugin `morph_ext`, step "Extended Morphology")

Adds what `morphometry/` (CAS, Gini, M20, one Petrosian radius) lacked, in the new module `morphometry/extended.py` (pure functions) and the driver `ds9/library/ds9_morph_ext.py`.
The existing step and its columns (`CONC ASYM GINI M20 R_PETRO`) are unchanged.

## Measurements (columns `MX_*`)
| column | definition |
|---|---|
| `MX_RP`, `MX_RP_LO`, `MX_RP_HI` | Petrosian radii (semi-major axis, px) for eta = 0.2 / 0.1 / 0.3 (parameters): eta(r) = I(r) / <I>(<r) with I(r) from the 0.8 r - 1.25 r annulus (Bershady et al. 2000), fine log-spaced growth curve in elliptical apertures, linear interpolation of the crossing |
| `MX_FLUX_P`, `MX_MAG_P` | Petrosian flux = flux inside 2 R_P(eta) and its magnitude |
| `MX_R20 MX_R50 MX_R80 MX_R90` | radii containing 20/50/80/90 % of the Petrosian flux (from the growth curve) |
| `MX_CONC` | 5 log10(R80 / R20) |
| `MX_KRON_R`, `MX_KRON_MAG` | Kron radius (first moment, iterated once) and the magnitude in `kron-scale` x r_k (>= `kron-min`) elliptical aperture |
| `MX_SMOOTH` | smoothness / clumpiness S (Conselice 2003): sum(abs(I - I_boxcar) - abs(B - B_boxcar)) / sum(I) inside 1.5 R_P, boxcar width 0.25 R_P, central box excluded, B = Gaussian noise patch of the measured rms (noise correction) |
| `MX_GINI_P`, `MX_M20_P` | Gini and M20 on the Lotz et al. (2004) segmentation: pixels inside 1.5 R_P brighter than the surface brightness at R_P |
| `MX_FLAG` | 1 object too close to the image edge / too small, 2 no Petrosian radius found, 4 2 R_P beyond the search radius or image, 8 neighbours masked, 1024 no position |

Apertures follow the catalog ellipse (A_IMAGE, B_IMAGE, THETA_IMAGE) unless *Circular apertures*. Other catalog objects (> 3 % of the target flux) are masked with ellipses of `neighbour-radius` x A_IMAGE;
masked pixels are replaced by the mean of the unmasked pixels at the same elliptical radius (sep's area-ratio correction is not used: it is wrong when the mask sits in the outskirts).
Growth curves (`morph_ext_growth.tsv`: NUMBER R FLUX ETA) and a growth-curve / eta plot of the brightest objects (`morph_ext_curves.png`) are written to `<work>/morph_ext/`
(menu: *Growth Curves Plot…*, *Growth Curve Table…*). Parameters live in the dialog of the "Morph+" chip (Measure tab); catalog keys `morphext,*`.

    python ds9/library/ds9_morph_ext.py IMAGE --catalog TSV [--mask FITS] [--work DIR] [--eta 0.2] [--kron-scale 2.5] [--mask-neighbours] [--circular] ...

## Validation (analytic Sersic references: exact Petrosian eta from the incomplete gamma function, noise-free renderings, 401² images, F_tot = 1e5; `plugins/morph_ext/tests`)
| case | R_P (eta 0.1/0.2/0.3) rel. error | R20/R50/R80/R90 max rel. error | C error | F_P/F_tot (true) |
|---|---|---|---|---|
| n=1, re 15, q 1 | +0.0001 / 0.0000 / -0.0004 | 0.0010 | 0.000 | 0.9931 (0.9933) |
| n=1, re 15, q 0.6, PA 40 | +0.0004 / +0.0002 / +0.0001 | 0.0029 | -0.005 | 0.9930 (0.9933) |
| n=4, re 12, q 1 | +0.0004 / +0.0017 / +0.0024 | 0.0072 | -0.013 | 0.8169 (0.8167) |
| n=2, re 12, q 0.7, PA 100 | +0.0005 / +0.0009 / +0.0026 | 0.0062 | -0.012 | 0.9427 (0.9431) |

* Kron: the flux in 2.5 r_k equals the analytic Sersic aperture fraction at the measured r_k to 0.0005 (n=1: 0.9214 vs 0.9216; n=4: 0.7429 vs 0.7430; q 0.6 and 0.7 cases the same); r_k = 14.99 px (n=1, re 15).
* Noise (12 realisations, n=1, re 15, q 0.7): R_P bias 0.000 / scatter 0.001 at rms 0.5, bias 0.000 / scatter 0.013 at rms 8; smoothness with the noise correction 0.022 (rms 0.5) and 0.005 ± 0.006 (rms 8), without it 0.045 and 0.376.
* Smoothness: smooth disk 0.007, the same disk with 10 Gaussian clumps (20 % of the flux) 0.208.
* Gini_P / M20_P: n=1 0.437 / -1.77, n=4 0.575 / -2.35, n=4 with noise 0.577 / -2.40 (peak S/N 989).
* Companion in the outskirts (30 % of the flux, 14 px mask): R_P error +77.5 % unmasked, +0.1 % masked; F_P/F_tot 1.299 → 0.994 (true 0.993).
* Driver on a synthetic 400² field (4 galaxies, noise 1, F = 4e4): R_P errors +0.5 / -1.4 / +0.2 / +0.1 %, F_P -0.1…-1.5 %, R50 ≤ 2.2 %, S 0.01-0.04, no flags.

## Limitations
Apertures use the catalog ellipse (no iterative centre / shape refinement); the asymmetry index of the existing step is not changed (its background correction and centre search are the old ones);
smoothness uses a simulated Gaussian noise patch (correlated noise, e.g. drizzled data, under-corrects: use a real empty-sky patch via the Python API `noise_patch=`); R_P needs the profile to
reach eta < level inside the cutout (flag 2/4 otherwise); the growth curve is not PSF-corrected.
