# HSC coadd epoch metadata audit for empirical PSF stars — 2026-10-20

## Finding

The Gaia-star ledger contains 333 usable stamps across 227 distinct HSC target
coadd frames. An Astrafek header-only audit searched every HDU in each existing
`cluster_<id>.fits` file for standard observation date/MJD cards. None of the
227 coadds contains `DATE-OBS`, `DATE-AVG`, `DATE-BEG`, `DATE-END`, `MJD-OBS`,
`MJD-AVG`, `MJD-MID`, `MJD-BEG`, or `MJD-END`.

This does **not** prevent image-space recentering: [`diagnose_empirical_psf.py`](diagnose_empirical_psf.py)
sets each stamp center from the HSC detection's native `X_IMAGE`/`Y_IMAGE`, not
from Gaia coordinates. Missing epoch metadata only leaves the Gaia-to-HSC
cross-match identity unverified at the exact observation epoch. A separate
existing sensitivity test propagated all 480 Gaia candidates through ±10-year
scenarios while holding the HSC positions fixed; every candidate remained
within the 2-arcsec match radius. An independent recalculation from its saved
table confirms 480/480 at both ±10 years, with maximum separations of 1.536
arcsec at −10 years and 1.572 arcsec at +10 years. This supports match-radius
stability under those scenarios, but it does not prove stellar identity in
crowded fields or approve any PSF profile. Input and output hashes for that
sensitivity are recorded in [`extended_psf_README.md`](extended_psf_README.md).

## Method and reproducibility

[`audit_empirical_psf_hsc_coadd_epochs.py`](audit_empirical_psf_hsc_coadd_epochs.py)
reads the 333-row Gaia stamp ledger, opens each of the 227 referenced coadd
FITS files lazily, and hashes the serialized HDU headers. It does not access
pixel arrays or make HSC requests. The output table records each frame's file
size, HDU count, supported epoch cards, and header digest. Three focused tests
cover epoch cards outside the primary HDU, reject a non-epoch exposure-time
card, and validate the full header-only run output.

The clarified run is stored at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_coadd_epoch_header_audit_2026-10-20_02/`.
The initial `_01` audit is preserved; `_02` adds the stamp-centering and Gaia
identity scope explicitly.

| Artifact | SHA-256 |
|---|---|
| `summary.json` | `fb5ff0c7b87d1d17bdeb42ef824980b3a50eb7c7f0a487d7449b1fd9fe88d93c` |
| `coadd_epoch_headers.tsv` | `f22959a0a0299dcfb329bfa0eca58ddbfbaed2516c11eb16f7921b389d4215a5` |
| Input `star_ledger.tsv` | `7dbe3cd079b91c78648f8f90132a9de9d96911e6e389819659a3da5c4fc01952` |
| Audit script | `7b152e56ddae0e6344add6872120437b3f824f57c1ac7419cd38cd97ae6c4c09` |

The summary reports 333 stamps, 227 unique frames, zero standard epoch cards,
`stamp_center_source=HSC detection X_IMAGE/Y_IMAGE`,
`gaia_epoch_propagation_required_to_recenter_stamp=false`,
`image_pixels_read=false`, `hsc_requests_made=false`,
`science_measurement=false`, and `psf_approved=false`.

## Consequence for ICL analysis

Continue measuring empirical profiles around the observed HSC detection centers;
the missing Gaia epoch is not a centering gate. Treat Gaia matches as
astrometrically plausible under the tested ±10-year scenarios, while keeping
crowding, stellar classification, unknown variability, saturation, mask
support, local-background contamination, and extended-wing transfer as separate
validation questions. Do not apply the profiles as an ICL PSF correction until
those empirical-wing and background-transfer checks pass. Keep
`psf_approved=false` and `science_measurement=false`.
