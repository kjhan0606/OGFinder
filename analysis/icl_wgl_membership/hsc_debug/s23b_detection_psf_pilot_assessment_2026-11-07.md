# S23B detection-PSF candidate pilot (2026-11-07)

## Purpose and scope

This CPU pilot checks whether the already received S23B target cutouts contain
enough isolated, compact sources to build a local PSF-core sample. It reads
only the approved member-position table and the 50 selected HSC cutouts. It
does not contact HSC, CAS, Gaia, or any other service, and it downloads no
data. The output is an engineering candidate inventory, not a PSF model or a
science measurement.

The 50-image subset contains 25 frames with primary coadd filter `i` and 25
with `i2`. Within each filter, five equal-count RA strata were crossed with
five within-stratum redshift strata; the middle row in each cell was selected.
The frames include 2,044 CAMIRA member positions used to avoid selecting
known member locations as field-star candidates. Input FITS sizes and SHA-256
values match the verified image receipt for all 50/50 frames.

The existing `prepare_detection_psf.py` runner used SEP detections at 1.5σ,
its pre-existing isolation/shape/SNR cuts, and at most 25 candidates per
frame. It saved raw 101×101-pixel stamps. Stamps were not recentered,
normalized, or background-subtracted. The primary measurement flag remains
`science_measurement=false`.

## Results

| Primary filter | Prepared / requested | Candidate count, median (range) | Frames with ≥4 candidates | Frames with ≥10 | Unique central-position detection |
|---|---:|---:|---:|---:|---:|
| `i` | 25 / 25 | 8 (0–24) | 22 / 25 (88.0%) | 12 / 25 (48.0%) | 8 / 25 |
| `i2` | 24 / 25 | 10.5 (2–24) | 22 / 24 (91.7%) | 13 / 24 (54.2%) | 6 / 24 |
| Combined | 49 / 50 | 10 (0–24) | 44 / 49 (89.8%) | 25 / 49 (51.0%) | 14 / 49 |

`S23B_10647` failed with `Too few valid pixels`; it is the one unprepared
frame, not a zero-candidate result. Across the 49 prepared frames, the median
detected-source count is 3,004. The per-frame FWHM-locus median is 5.125
pixels (range 4.125–7.375); this is a provisional moment-based selection
diagnostic, not a seeing calibration. Candidate-stamp valid-pixel fraction
has median 0.863 (10th–90th percentiles 0.687–0.957); 166/492 stamps have at
least 90% valid pixels.

The legacy runner labels a unique detection within 1 arcsec of the supplied
center as a BCG. Here the supplied centers are S23B CAMIRA centroids, explicitly
not BCG positions, so the table's central-position counts are not BCG
identifications. The red cross in the review images marks the CAMIRA centroid.
Visual spot checks of the 24-candidate frame `S23B_6162` and zero-candidate
frame `S23B_2643` show why candidate counts cannot be treated as a complete
star census: both fields contain compact-looking sources, while the bright
object masks and isolation cuts reject many. No Gaussian or S19A PSF fallback
was substituted. The local S19A Gaia result subset contains 944 sources; its
nearest source is 523 arcsec from any of the 492 S23B morphology candidates,
so it provides no S23B candidate identification. This zero-overlap screen is
not evidence that the candidates are non-stellar.

## Assessment and next gate

The pilot shows that many S23B images can support a useful *candidate* PSF
core sample, but it does not establish candidate purity, subpixel centering,
coadd-epoch astrometry, target/control PSF matching, or PSF wings. The raw
stamps extend only 8.4 arcsec from the center; they cannot calibrate the broad
wings relevant to diffuse-light measurements. Keep the member-light profiles
exploratory and keep `science_measurement=false`.

A Gaia DR3 positional cross-match plan is prepared for these 492 positions:
2-arcsec cones, two batches of at most 250, and 18 returned columns comprising
source ID, astrometry/epoch, quality flags, magnitude/variability/classifier
fields, and match separation. The original rough estimate was 62
candidate-source rows (about 12 kB), extrapolated from the previous S19A
candidate-match rate. After explicit approval, the query completed with 77
rows. It requested no images, shape rows, photo-z values, HSC data, or CAS
data. The follow-up audit found 48 conservative diagnostic candidates; the
wide-stamp profile assessment is in
`s23b_gaia_empirical_psf_wing_assessment_2026-11-07.md`.

## Reproducibility

- Runner: `prepare_detection_psf.py`, SHA-256
  `1dec4fa9e95b7977f7584003e991c8b8c0ab61773fd18863a0d17ab3e5c7a55f`.
- Pilot manifest:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_detection_psf_pilot_2026-11-07_01/manifest.json`
  (SHA-256 `be3ed0bea7706aef37cced52e56ef11e57a4d50aac914da4bc2b72c9d5ce0fa1`).
- Selection and input hashes:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_detection_psf_pilot_2026-11-07_01/plan.json`
  (SHA-256 `37b4b8558018ba5639a05a8bbe616df45eb29a13f1084fe7a081ac6849e330fc`).
- Preparation summary:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_detection_psf_pilot_2026-11-07_01/preparation/summary.json`
  (SHA-256 `5a8646689173bc084e44a5722d26887f0a4db30cbca408e8d41abdfc6f5745db`).
- Approved Gaia plan:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_gaia_plan_2026-11-07_01/plan.json`
  (SHA-256 `394083a8cdaf0770376959e71c8bd46bd004c762dcb3e77a6256bac656ee52f0`).
- Gaia request positions:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_gaia_plan_2026-11-07_01/candidate_positions.tsv`
  (492 rows; the JSON requests file is hash-pinned by the plan).

The Gaia response SHA-256 is
`1f1f8b947a3b23e1dad4675447247137d66d5f2586d6aea3a73989df3cb88e6b`.
The query and subsequent profile stacks remain diagnostic only;
`science_measurement=false` and `psf_approved=false`.

No commit or push has been made.
