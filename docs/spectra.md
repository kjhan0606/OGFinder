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

## Limits

* Redshifts come from **emission lines only**: absorption-line galaxies, stars and featureless spectra give no answer (quality 0); a single emission line is never accepted unless the redshift
  is given. Without the ambiguity cut the first version returned wrong z for 4 of the 10 SDSS spectra (confidence 1.04-1.1); the cut removes them but also rejects marginal true detections.
* The 22-line list is vacuum wavelengths; air-wavelength data are not converted. Line widths are not deconvolved from the instrument. Blended doublets are fitted line by line (the [OII] doublet as one line).
* 2D extraction assumes a straight trace along a row (no curvature / tilt), no cosmic-ray rejection; cubes: no spatial smoothing or Voronoi binning; line maps are not a plugin output.
* SDSS spectra are in 1e-17 erg/s/cm^2/A: fluxes are reported in the file's units.
