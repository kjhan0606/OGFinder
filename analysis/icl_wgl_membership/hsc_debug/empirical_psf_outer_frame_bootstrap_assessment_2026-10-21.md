# Filter-separated outer-field PSF-wing support — 2026-10-21

## Scope and correction

The existing matched target/control PSF pilot used 18 Gaia stars on 18 frames above 800 kpc. A pixel-free support sweep with saved HSC-I images and Gaia/mask ledgers found 91 eligible star/control pairs above 600 kpc, on 78 independent target frames. They split into `i`: 29 stars/25 frames and `i2`: 62 stars/53 frames. Twelve frames contain multiple stars, at most three, so the stack averages stars within each coadd frame and resamples frames.

The first execution (`_01`) mistakenly pooled `i` and `i2` in the stack key. **That stack and its pooled filter counts are superseded and must not be used.** The corrected `_02` execution stores separate cells for each subfilter and its independent validator reproduces all 48 cells from the saved profiles. No earlier HSC data were changed or downloaded.

The 600 kpc outer-star cut keeps the star-centered 0–100 kpc satellite-scale profile outside the measured 200–300 kpc cluster annulus. It does not make the 200–300 kpc PSF wings safe: those outer stellar rings can approach the cluster ICL and retain cluster-correlated background.

## Filter-specific support

The strict per-profile rule requires at least 50% unmasked ring area and all eight azimuth sectors. Independent matched-frame counts are:

| Filter | 50–75 kpc | 75–100 kpc | 200–225 kpc | 225–250 kpc | 250–275 kpc | 275–300 kpc |
|---|---:|---:|---:|---:|---:|---:|
| `i` | 0 | 1 | 1 | 1 | 0 | 0 |
| `i2` | 6 (120–150 px) / 4 (160–190 px) | 6 / 4 | 0 | 0 | 0 | 0 |

The local-sky setting does not change the decisive lack of `i` support at satellite scales or the lack of either-filter support over the ICL annulus. Under a diagnostic relaxation to at least 20% ring area and six of eight sectors, `i` has only 1–2 frames per 200–300 kpc bin; `i2` has 3–6. Requiring all eight sectors leaves zero frames for `i2` throughout 200–300 kpc. Partial coverage is not a mask-flux correction and cannot be reported as a calibrated PSF.

Geometry further limits an outer background estimate: among all 91 pairs, 63 stars/57 frames have matched target/control radial coverage to 350 kpc, 37/37 to 400 kpc, only 4/4 to 500 kpc, and none to 600 kpc. The existing local-sky annuli (120–150 and 160–190 pixels) are much closer to each star than 300 kpc; their subtraction can remove extended stellar light. The current images therefore do not give a well-supported background annulus beyond the measured wing profile.

## Decision

The expanded, correctly filter-separated sample does not calibrate the PSF at either the satellite-outskirts scale in `i` or the 200–300 kpc ICL scale in either filter. The `i2` 50–75 kpc pilot has only 4–6 independent frames and remains a diagnostic, not a correction. `psf_approved=false` and `science_measurement=false` remain in force. A validated extended PSF/sky-transfer reference or substantially better mask and frame support is still required before interpreting the ICL–satellite association.

## Reproduction and provenance

- Pixel-free threshold sweep: [`audit_empirical_psf_outer_star_support.py`](audit_empirical_psf_outer_star_support.py); summary SHA-256 `0fb88bf6b7507807276e75b282f62a3803cff1b01495d3d87b821d0a2f9fe920`.
- Filter-separated runner: [`stack_empirical_psf_outer_field_wings_frame_bootstrap.py`](stack_empirical_psf_outer_field_wings_frame_bootstrap.py); run on lageunha CPU with 10,000 frame-bootstrap replicates.
- Independent validator: [`validate_empirical_psf_outer_field_wings_frame_bootstrap.py`](validate_empirical_psf_outer_field_wings_frame_bootstrap.py). It reproduced 4,368 profile rows, all 48 per-filter cells, and the 10,000-replicate frame bootstrap.
- Tests: eight focused support/runner/validator tests passed, including a regression test that keeps `i` and `i2` separate; `git diff --check` passed.
- Corrected output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_outer_field_wings_frame_bootstrap_2026-10-20_02/`.
- `summary.json` SHA-256: `556fe0e828a5f2b0ff1600c4e4f2449147e9b383a115d3757eb4b5ea18faa388`.
- `validation.json` SHA-256: `844337eb4fa89f72648ec1d660f57cc89208df2fa3ff9b44f15b406a3116da6f`.
- Per-filter stack SHA-256: `8882c9be1f825108b69e9b04db68e07d6c2e8aea45580c2252ddb6ded241406c`.
- Profile table SHA-256: `f256ef84eb610560045994bf0ada51a81061b5162b710a2fd4d5871449127785`.
