# Spectroscopy (plugin `spectra`, tab Measure, chip "Spectra")

Spectra linked to catalog objects: 1D spectra, 2D spectra and IFU cubes; Gaussian line fits and redshifts from emission lines; viewer windows. Engine `ogfkit/spectra.py`
(pure functions, JSON-serialisable results), CLI `plugins/spectra/spectra.py`, synthetic truth `ogfkit/spectrasynth.py`, demo data `plugins/spectra/make_demo.py`.

## Object -> spectrum links

*Link Spectra to Catalog* builds `spectra_links.tsv` (NUMBER FILE KIND X Y ROW OK) from either a file-name pattern (`spec-dir` + `spec-pattern`, `{NUMBER}` replaced) or a user
link table (`link-file`: NUMBER FILE [KIND X Y ROW]; relative FILE names are relative to the table or to `spec-dir`) and adds `SP_FILE SP_KIND SP_OK`. KIND is detected: FITS table or
1D image = `1d`, 2D image = `2d`, 3D = `cube`. 2D spectra: dispersion along axis 1 (WCS CRVAL1/CDELT1/CRPIX1, or CTYPE2 wavelength: transposed), object row from the link table or
the catalog column `row-column`. Cubes: wavelength on axis 3; the catalog `X_IMAGE/Y_IMAGE` are read as spaxel coordinates (`cube-xy`) or X, Y of the link table. Readers: FITS binary tables
(wave/wavelength/lambda/loglam, flux, err/ivar), 1D images with linear or log WCS, ASCII/CSV; wavelength units from CUNIT/TUNIT or `wave-unit` (micron-sized values without unit are read as micron);
inverse variance 0 -> no error.

## Fit Lines + Redshift

Per object: extraction (2D: boxcar with fractional edges or Horne-optimal with a Gaussian profile fitted to the collapsed spectrum, optional sky rows; cube: circular aperture with 5x5 sub-pixel
weights), running-median continuum with iterative line masking, noise from the MAD (or the file's errors), Gaussian matched-filter line detection (`snr-min`), blind redshift by matching
detections to a 22-line vacuum rest list (score = sum of line weight x S/N over matched lines, minus a small penalty for strong lines in range that are missing), Gaussian + local linear
continuum fits of the matched lines (weighted least squares, covariance errors scaled by sqrt(chi2_red) if > 1) and a weighted-mean redshift of the lines with S/N > 2.5.
Quality `SP_ZQ`: 0 = none or **ambiguous** (best score < 1.2 x the best alternative redshift, or a single weak line), 1 = single line, 2 = two lines, 3 = three or more, 4 = fixed at the
catalog redshift (`z-column`). Columns: `SP_Z SP_ZERR SP_ZQ SP_NLINES SP_SNR` (continuum S/N per pixel) and, for the chosen `fit-line`, `SP_LINE_FLUX SP_LINE_FLUXERR SP_LINE_FWHM_KMS`
(observed, includes the instrument) `SP_LINE_EW` (rest frame). `spectra_results.json` keeps all detections and fitted lines per object. Failures of a single file are reported on stderr and
leave empty cells.

Cube helpers in the engine (`extract_cube`, `white_light`, `line_moments`: flux / velocity / dispersion maps from continuum-subtracted moments) are used by the viewer and tests; the plugin steps
extract spectra, they do not write maps.

## Viewer

*Spectrum Viewer...* opens a window with a plot of the object (spectrum, continuum, error, fitted lines; for 2D spectra the 2D image with the extraction window, for cubes the white-light image
with the aperture), a NUMBER combobox, Prev / Next buttons; it starts at the selected catalog object. *Link / Result Table...* shows the link table. Catalog keys `spectra,links_file`, `spectra,results_file`.

## Validation (`plugins/spectra/tests/test_spectra.py`, 9 tests, ~45 s)

| Test | Result |
| --- | --- |
| line fit, 300 noise realisations (flux S/N ~ 16) | flux pull mean +0.00 / std 1.00, centroid pull mean -0.07 / std 0.99; continuum-subtracted noise estimate / truth 0.989 |
| continuum under strong lines | median abs deviation 0.009 (noise 0.1), max 0.090 |
| blind redshift vs line strength (40 spectra each, z 0.2-1.4, 3 pixel tolerance) | strongest line S/N 7: 34/40 correct, 0 wrong, 6 no answer; S/N 13: 40/40; S/N 21: 40/40 (scatter 3.4e-5); lowest strength: 2/40 correct, 0 wrong, 38 no answer |
| 2D extraction (Gaussian profile, sky residual) | flux/truth boxcar 0.984 (rms 0.017), optimal 0.997 (rms 0.015); pixel rms error 1.41 vs 1.13 (x1.24 better) |
| IFU cube 21x21x240, rotating disc (150 km/s) | aperture H-alpha flux 691.5 +- 2.1 vs 696.8 truth; velocity map: 429/441 spaxels, rms 25.7 km/s (bright spaxels: 5.2 km/s), median +0.7; dispersion 80 km/s (truth 80) |
| CLI: 12 x 1D + 1 x 2D + 1 x cube + 1 missing file | 14/14 redshifts, max abs dz 1.2e-4; missing file -> SP_OK 0; deterministic output; viewer PNGs 640x288 (1D), 640x416 (2D, cube) |
| 10 real SDSS DR17 spectra (plate 266; 7 galaxies, 1 QSO, 1 star, ...) | quality >= 2 for 6, all within 94 km/s of the SDSS pipeline redshift (rms 44 km/s); the 4 others (absorption-line galaxies, star, weak emission) get quality 0 |

## Kinematics: 2D rotation curves and IFU velocity / dispersion maps (`--task kin`, step "Fit Kinematics", viewer "Kinematics Viewer")

`ogfkit/kinematics.py` (engine, no GUI dependencies), driven by `spectra.py --task kin` for objects linked to 2D spectra (long slit) or cubes.

**2D (long slit).** The emission line (default Halpha, `--kin-line`; redshift from the catalog column `--z-column`, `--kin-z`, or the line-fit redshift) is fitted per spatial row with a Gaussian + local
continuum in a +-`--kin-window-kms` window; rows are binned along the slit until the line S/N reaches `--kin-bin-snr` (1D adaptive binning). Velocity v(y), error and intrinsic sigma
(instrument FWHM `--kin-inst-fwhm` A subtracted in quadrature) give an arctan rotation curve v(y) = vsys + v_obs (2/pi) arctan((y-y0)/r_t) by Levenberg-Marquardt with multi-start. With `--kin-inc` the
observed amplitude is deprojected, v_c = v_obs / (sin i cos psi) (`--kin-slit-psi` = angle between slit and major axis; **no correction is made for psi unless it is given**).

**Cubes.** Every spaxel (or accretion bin, `--kin-bin-snr`, otherwise spaxels with S/N > `--kin-snr-pix`) gets a Gaussian+continuum fit -> maps FLUX, VEL, VELERR, SIGMA, SIGERR, SNR, plus the best-fit
disc MODEL, RESID and the BIN id (`kin_<N>_maps.fits`; `kin_<N>_vel/sigma/flux.fits` are single-HDU copies the viewer opens in frames). The velocity map is then fitted with a thin-disc model
(vsys, v_c, r_t, PA, inclination, centre; inclination optionally fixed from `--kin-inc`), giving kinematic PA and inclination with errors from the Jacobian.
Conventions: **PA = direction of the receding side, degrees counter-clockwise from +x** (image axes), v_c is the *deprojected* asymptotic speed, `SP_KIN_VSINI` = v_c sin i is what is observed. The flux-weighted mean intrinsic
dispersion goes into `SP_KIN_SIGMA`. Columns: `SP_KIN_VSYS, SP_KIN_VSINI, SP_KIN_VC, SP_KIN_VC_ERR, SP_KIN_RT, SP_KIN_PA, SP_KIN_INC, SP_KIN_SIGMA, SP_KIN_CHI2R, SP_KIN_N`; details in `spectra_kin.json`.

**Validation (all synthetic, `plugins/spectra/validation/run_kin_validation.py`, 20 random discs per configuration, report `kin_validation_report.json`; unit/CLI tests `tests/test_kinematics.py`, 9 tests).** Discs: vc 120-300 km/s,
r_t 2-6 px, i 30-75 deg, random PA/centre, sigma 25-70 km/s, 41x41 spaxels, [NII] doublet included, instrument FWHM 2 A. Errors are fit - truth (median / std / max |.|):

| configuration | PA (deg) | inclination (deg) | v_c rel. | r_t rel. | sigma ratio | notes |
|---|---|---|---|---|---|---|
| high S/N, inclination free | 0.02 / 0.11 / 0.23 | -0.02 / 0.86 / 2.9 | 0.0001 / 0.016 / 0.058 | 0.001 / 0.009 / 0.020 | 1.006 | corr(v_c, i) median -0.74 |
| high S/N, inclination fixed | 0.01 / 0.12 / 0.35 | - | -0.0003 / 0.009 / 0.034 | -0.003 / 0.015 / 0.052 | 1.004 | |
| low S/N, per-spaxel | -0.29 / 21 / 86 | 0.9 / 25 / 69 | 0.063 / 4.8 / 20 | 0.005 / 0.47 / 1.0 | 1.13 | **heavy tails: unusable** |
| low S/N, binned to S/N 8 | -0.12 / 1.1 / 2.5 | -3.1 / 13.8 / 47 | -0.008 / 1.7 / 6.5 | -0.06 / 0.18 / 0.47 | 1.10 | PA good; inclination, v_c poorly constrained (corr -0.88) |
| seeing 2.5 px FWHM | -0.01 / 0.13 / 0.35 | -2.4 / 2.0 / 5.3 | +0.053 / 0.045 / 0.19 | +0.20 / 0.07 / 0.33 | 1.05 | beam smearing, see below |
| slit, psi = 0 | v_obs +0.0014 / 0.017 / 0.063 | - | - | 0.007 / 0.024 / 0.085 | 0.993 | centre 0.04 px |
| slit, psi = 30 deg, uncorrected | v_obs -0.159 / 0.11 / 0.38 | - | - | - | 0.997 | expected cos 30 = 0.87 loss |
| slit, noisy (psi 0) | v_obs -0.049 / 0.12 / 0.28 | - | - | -0.05 / 0.16 / 0.41 | 1.04 | |
| slit, seeing 2.5 px | +0.055 / 0.042 / 0.19 | - | - | +0.20 / 0.09 / 0.44 | 1.04 | |

What this says: with adequate S/N (noise-free-like) PA is recovered to 0.1 deg and v_c to ~1-2 %, the inclination to ~1 deg (it is degenerate with v_c; fixing it from imaging halves the v_c scatter).
At low S/N per-spaxel fitting produces catastrophic outliers; binning fixes the PA but the inclination/v_c pair stays poorly determined, so quote `SP_KIN_VSINI` (well-constrained) rather than `SP_KIN_VC` in that regime,
and look at `chi2r` and the errors. Seeing is **not corrected**: the velocity gradient inside the beam inflates sigma (by 5 % here, 4.6 km/s for a 40 km/s disc in the unit test), flattens the curve
(r_t overestimated by 20 %) and biases v_c by ~5 %. The high-S/N cube fits have chi2r median 1.2 (mean 1.5, max 3.5), i.e. slightly above 1; the cause (line-profile mismatch from [NII] blending / seeing / error underestimate) was not investigated. Keep the window below ~700 km/s for Halpha so [NII] stays outside it.

Limits specific to kinematics: validated **only on synthetic discs from a generator that shares the disc model with the fitter** (a fit of its own model family; no real IFU cube or long-slit data was fitted,
so line-profile asymmetries, non-circular motions, warps, bars and real PSFs are untested); thin disc only; one Gaussian per spaxel; no PSF deconvolution / beam-smearing correction; the 2D path assumes a straight trace and
a slit through the centre (psi must be supplied for the deprojection); instrument FWHM is a single number; air/vacuum wavelengths are not converted; adaptive/accretion binning is S/N-based, not Voronoi.

## Limits

* Redshifts come from **emission lines only**: absorption-line galaxies, stars and featureless spectra give no answer (quality 0); a single emission line is never accepted unless the redshift
  is given. Without the ambiguity cut the first version returned wrong z for 4 of the 10 SDSS spectra (confidence 1.04-1.1); the cut removes them but also rejects marginal true detections.
* The 22-line list is vacuum wavelengths; air-wavelength data are not converted. Line widths are not deconvolved from the instrument. Blended doublets are fitted line by line (the [OII] doublet as one line).
* 2D extraction assumes a straight trace along a row (no curvature / tilt), no cosmic-ray rejection; the redshift fit does not use cubes (kinematics below do).
* SDSS spectra are in 1e-17 erg/s/cm^2/A: fluxes are reported in the file's units.
