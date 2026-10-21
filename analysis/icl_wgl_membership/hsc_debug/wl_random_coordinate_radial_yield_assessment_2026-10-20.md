# Radial support in the partial HSC random-coordinate export — 2026-10-20

## Question and inputs

The current random-coordinate plan queries from twice each lens's proper-3-Mpc
angular radius to 30 arcmin. This audit asks how much radial width that leaves
for the first two completed batches, and whether the geometry-clean candidates
are concentrated within that interval. It uses only the already approved
33,514 random coordinates and their existing CAMIRA-overlap/duplicate screen;
it makes no HSC request and reads no shapes, photo-z values, or images.

The first 50 plan-order lenses are explicitly a partial, nonrepresentative
pilot. The script verifies the partial geometry output and coordinate-plan
hashes, recomputes great-circle distance from each lens center, and bins points
by offset from the inner 2×3-Mpc exclusion boundary. It also computes the ideal
spherical annulus area; this is not corrected for the actual HSC footprint.

## Results

The first 50 lens exclusion boundaries leave only 7.14–8.76 arcmin of radial
width before the 30-arcmin outer limit (median 8.10 arcmin). Thus the absence of
points in offset bins beyond 10 arcmin is imposed by the current query geometry,
not evidence of a catalog or HSC-coverage hole.

| Offset from inner boundary | Raw points | Geometry-clean points | Clean fraction | Lenses with ≥1 clean point | Median clean points/lens |
|---|---:|---:|---:|---:|---:|
| 0–5 arcmin | 19,920 | 1,199 | 6.02% | 27/50 | 2 |
| 5–10 arcmin, truncated by each lens's 30′ outer limit | 13,594 | 930 | 6.84% | 30/50 | 2 |

In total, 2,129 of 33,514 points (6.35%) pass the existing CAMIRA-aperture
overlap and exact-coordinate-duplicate cuts. Counts are highly uneven: the
first shell has 23 lenses and the outer shell 20 lenses with no geometry-clean
point. For these 50 lenses, the ideal area of a separate 30–60 arcmin annulus is
6.02–7.15 times the current 2×3-Mpc–30′ annulus (median 6.42). This ratio is a
geometric comparison only; real HSC random-catalog coverage and depth may lower
the yield.

The evidence strengthens the case for measuring the 30–60 arcmin candidate
counts, but it does not establish that those positions are valid random lenses.
The prepared count-only plan remains unsubmitted pending its separate approval.
Even a successful count query would need full CAMIRA overlap screening, HSC
image/mask/depth/PSF matching, random-aperture overlap checks, and covariance
validation before shape-row export or lensing subtraction.

## Provenance and reproduction

Astrafek implementation: [`audit_wl_random_coordinate_radial_yield.py`](audit_wl_random_coordinate_radial_yield.py).
Two focused tests pass in [`test_audit_wl_random_coordinate_radial_yield.py`](test_audit_wl_random_coordinate_radial_yield.py).
The audit script SHA-256 is `65b6e761226dce1ce300681cb781476d68d85a75f99529c884248e813ec207df`.

Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_coordinate_radial_yield_audit_2026-10-20_03/`.
The summary SHA-256 is `cd5a30122a2beb2ad875dec3f59bbca35286ee0862cbc1473d021cbb240402ad`.
The radial table SHA-256 is `d44df93ec47f2f230f5d4aaa5c93fce9f69944ad703e2734fdcbce9c14378164`;
the per-lens table SHA-256 is `9e1be3126071926d946a8afd450f78b29c42a61ce4dce4c42845db266031177d`.
The audit summary records source-plan and partial-pilot hashes. Both
`science_ready` and `science_measurement` remain false.
