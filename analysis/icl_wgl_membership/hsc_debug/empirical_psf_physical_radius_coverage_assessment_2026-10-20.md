# Empirical-star PSF radial coverage versus ICL scales — 2026-10-20

## Method

The audit joined the 333 hash-recorded Gaia/HSC star stamps (227 cluster
frames) to the already downloaded 800-frame manifest. It used each frame's
reviewed HSC-I pixel scale, lens redshift, and the same flat cosmology as the
ICL profile pipeline (`H0=70 km s⁻¹ Mpc⁻¹`, `Omega_m=0.3`, `Tcmb0=0`). The
provenance check matched every stamp hash in the star ledger to both the
existing frame-weighted PSF stack and the star-stamp quality audit.

The existing empirical profile stack has radial edges only through 96 pixels;
the input stamp is 401×401 pixels, so its inscribed circular radius is about
200 pixels. The ICL target annulus is 200–300 kpc and the outer satellite
profile extends to 75 kpc. This calculation tests only physical-radius
coverage. It reads no image pixels and does not evaluate PSF signal or mask
correctness.

## Results

The 48–96 pixel empirical-wing interval maps to a proper radius of 38.3–50.9
kpc at its inner edge (10th–90th percentile) and 76.6–101.7 kpc at its outer
edge; the median interval is 44.1–88.2 kpc. The stack's 96-pixel endpoint
corresponds to 70.6–104.6 kpc across the 333 stamps, with a median of 88.4 kpc.
It reaches 75 kpc in 316/333 stamps, geometrically overlapping the outer
satellite-profile scale, although the separate mask and sky audit still
prevents interpreting that overlap as a calibrated wing.

No empirical radial profile reaches the ICL annulus inner edge of 200 kpc in
any of the 227 frames. The full stamp's inscribed-circle radius has median
184.1 kpc (10th–90th percentile 159.8–210.8 kpc). A complete circular stamp
footprint reaches 200 kpc in only 52/227 frames and reaches 300 kpc in 0/227.
Thus even the saved raw stamps do not span the full 200–300 kpc ICL annulus;
the measured profile stack stops much farther inward.

This establishes a direct coverage gap: the current empirical-star stack can
inform inner PSF behavior at scales overlapping satellite outskirts, but it
cannot empirically constrain the PSF contribution across the ICL measurement
annulus. Existing injection exercises do not close that gap because they test
the downstream analysis transfer and/or synthetic PSF inputs, not an observed
stellar profile through 200–300 kpc. The current ICL–satellite result therefore
still needs a validated wide-radius PSF/sky-transfer model before it can be
promoted to a science measurement.

## Reproduction and provenance

- Astrafek audit: [`audit_empirical_psf_physical_radius_coverage.py`](audit_empirical_psf_physical_radius_coverage.py)
- Tests: [`test_audit_empirical_psf_physical_radius_coverage.py`](test_audit_empirical_psf_physical_radius_coverage.py), 2 passed.
- Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_physical_radius_audit_2026-10-20_01/`
- `summary.json` SHA-256: `6d88022648930be6fd6f5f6785d66271cb399b1f6c4107dd7305453fead5d872`
- Per-stamp table SHA-256: `b1b4cfee369a20939dbf3b6b040ff4a08b7efb4239e155829f65fc1bd99a1f32`
- Audit script SHA-256: `6338c8df135f7f2204fd981bfd407e8c2d4e86c55ccda6305b0a60134c9f551d`

The per-stamp table is an audit of geometry only. `science_measurement=false`
and `psf_approved=false` remain in force.
