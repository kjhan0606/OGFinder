# S23B member-profile injection: position-specific PSF request plan

## Purpose

The exploratory source-aware member-halo injection recovery used 427
S23B CAMIRA member positions, one in each target/filter frame (62 HSC `i`,
365 `i2`). The injected halo was not convolved with the local coadd PSF.
Before interpreting that recovery as a realistic measurement-transfer
calibration, prepare a position-specific PSF request for those same positions.

## Coordinate and input audit

The plan builder rechecked the selected-member table against the imported
CAMIRA member-position table, the saved injection-recovery receipts, the
target manifest, image-mask records, and all 427 target coadd FITS headers and
checksums. The selected coordinates are unique and remain inside the approved
image footprints. The headers split into 62 `i` and 365 `i2` frames; the PSF
Picker request uses the established `HSC-I` filter designation for this image
family.

Three previously received and locally audited PSF banks contain 4,291 files
in total (744, 1,542, and 2,005). None of the 427 requested member coordinates
is an exact match to those banks; none lies within 60 arcsec, and the nearest
existing PSF coordinate is 69.6869 arcsec away. Reusing a nearby model would
not provide a position-specific PSF for these members.

## Prepared request

- HSC PDR3 `pdr3_wide`, centered coadd model-PSF FITS cores, `HSC-I`.
- 427 exact CAMIRA member coordinates from the audited injection sample.
- Estimated transfer size: 8,608,320 bytes if every file matches the largest
  file (20,160 bytes) observed in the existing audited PSF banks. This is a
  review cap, not a guaranteed server-side size; stop if the returned data
  would exceed it.
- No coadd images, shape rows, photo-z data, or catalog exports are requested.

Plan directory (outside Git):
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_profile_psf_plan_2026-11-07_03/`.

- Plan SHA-256: `03f974ad24ead32eabd103ef7a6e0d6b52317aa2647023a5f7402003a2e3b25d`
- Coordinate-list SHA-256: `7327d411bacab15954538c64cb44deae9141955377325bc459a1e6a162db7b31`
- Position-mapping SHA-256: `3a25fd2e7180a8380bcfe902a0d3157715361256a331a5f65490b35269af8735`
- Planner SHA-256: `b35bad5464d7b25f423fd30aa864328dc0f6d1c924c61321f6b4724b5ebc209e`

The plan is `prepared_not_submitted`. No network request was made, and
`science_measurement=false` remains in force. This request needs its own
explicit approval.

## Next analysis after receipt

After an approved request is received, validate every FITS file, byte count,
coordinate match, filter, and checksum before using it. Convolve the injected
member-halo model with the PSF at its own position and rerun the saved
injection/recovery pipeline. Report PSF-convolved recovery separately from
the earlier unconvolved result. This step still does not add source-Poisson
noise or establish completeness; those limitations and the science flag must
remain explicit.

The scoped batch auditor accepts this plan schema and requires the exact
approved plan hash. The injection runner now has an explicit PSF-bank mode
that requires the same plan hash, uses each selected member's own PSF, keeps
the signed PSF pixels, and records the per-file and plan provenance. Neither
the data request nor the PSF-mode calculation has been submitted or run.
