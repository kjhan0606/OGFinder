# Weak-lensing random-center candidate-count audit — 2026-10-20

## Scope and validation

The audit revalidated the existing approved count-only HSC CAS result for all
800 lens centers. It checked the count plan, all 32 returned count CSV hashes,
the aggregate 800-lens summary, and the prepared coordinate-only SQL plan. It
made no CAS request and did not read random coordinates, source rows, shapes,
photo-z values, or images.

Each count is the number of quality-filtered `pdr3_wide.random` points in the
annulus from twice that lens's angular proper-3-Mpc aperture to 1,800 arcsec.
The area normalization uses a spherical annulus. The query also applies
`adjust_density < 0.01`, primary-source selection, and the HSC-I edge and bright
object center masks.

## Results

The 800 approved count rows sum to 780,191 random-point candidates. Counts range
from 364 to 1,443 per lens, with median 1,000. Every lens exceeds the planned
200-candidate target before the broader all-cluster geometry screen. Thus the
count-only result does not indicate a shortage of raw random-table points.

The available annulus area has median 1,842 arcmin² (range 1,186–2,080 arcmin²).
The observed candidate surface density after the query's thinning and quality
cuts has median 0.566 arcmin⁻² (range 0.214–0.791 arcmin⁻²). Raw counts correlate
with lens redshift in this fixed sample (Spearman ρ=0.595), while dividing by
the per-lens available annulus area reduces that descriptive association to
ρ=0.080. This indicates that much of the raw count trend follows the varying
angular exclusion area; neither coefficient is an inferential significance
test.

| Lens redshift | Lenses | Median candidates | Median candidates / arcmin² |
|---|---:|---:|---:|
| 0.25–0.30 | 9 | 580 | 0.476 |
| 0.30–0.35 | 138 | 776.5 | 0.557 |
| 0.35–0.40 | 117 | 905 | 0.554 |
| 0.40–0.45 | 204 | 1,022.5 | 0.569 |
| 0.45–0.50 | 162 | 1,107 | 0.576 |
| 0.50–0.55 | 159 | 1,160 | 0.571 |
| 0.55–0.70 | 11 | 1,218 | 0.588 |

This is availability before checking overlap with any of the 800 CAMIRA
apertures, coordinate duplication, random-aperture overlap, image footprint,
mask/depth/PSF/sky matching, or source-shape quality. The first 50-lens geometry
pilot retained only 2,129 of 33,514 candidates (6.35%) and found at least one
geometry-clean point for 31/50 lenses. The full 800-lens geometry screen is
therefore the next limiting gate; the first-50 survival fraction must not be
extrapolated to the full sample.

## Reproduction and provenance

- Astrafek audit: [`audit_wl_random_coordinate_candidate_counts.py`](audit_wl_random_coordinate_candidate_counts.py)
- Tests: [`test_audit_wl_random_coordinate_candidate_counts.py`](test_audit_wl_random_coordinate_candidate_counts.py), 2 passed.
- Result: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_candidate_count_audit_2026-10-20_02/`
- `summary.json` SHA-256: `9be7c73d515a5ce42edd7dad4684469a15409969ac546c3a727407ffa3d32070`
- Per-lens table SHA-256: `0fc89751cfd78e0e7656c09ac47dfd2c1e714fdd04fc0e890cb91e63fcbd8267`
- Audit script SHA-256: `fb2f733e93f2fffb7c52f1d4e7c47e88e8345c316ab10db26cea566ee7bb2e09`
- Coordinate-plan SHA-256: `e8672433abdfe97c3f612c9327d9c449aa5c966b888702313bc15d4a95004827`
- Count-plan SHA-256: `70bac1c067a36fb13d06b39bfb9728fbc4bff7c61e45c907ac01fb74d5715a70`

`science_measurement=false`; no random positions have been selected or used in
a weak-lensing stack.
