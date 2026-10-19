# HSC PDR3 PSF Picker batch audit — 2026-10-19

## Approved request and provenance

The user approved the remaining 744 model-PSF requests for the existing
injection sample: 56 additional cluster-center positions and 688 satellite
injection positions. Combined with the already received center test, this
covers 57 cluster centers and 688 satellite positions (745 unique positions).
The batch used the official HSC PDR3 PSF Picker with `rerun=pdr3_wide`,
`type=coadd`, `filter=HSC-I`, and centered stamps. No science images, shape
catalog rows, or photo-z data were requested. Credentials were not written to
the repository or audit manifest.

The approved input plan is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_injection_plan_2026-10-19_02/`.
The coordinate-list SHA-256 is
`762f048188ef13717d4ae5fe25ccf6e2c14d88bcbe95df3c82b0c77c124194f0`; the
position-mapping SHA-256 is
`0d679aa6999fef99b28b380b423618b8b738b3256f47d0b53d8af625f0eecb5d`.
The plan is retained unchanged as the approved request record.

## Retrieval and audit

The 744 returned FITS files are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_injection_2026-10-19_01/`.
The audit manifest is `psf_picker_audit.json` in that same directory. Its
SHA-256 is
`6375bc58df6d5a15d837cf5b13d3db79ebdc6944dbf7510164f465b6cbe72f7e`; the
manifest records a SHA-256 for every FITS file.

All 744 names preserved a unique planned RA/Dec, tract, patch, and `HSC-I`
filter, with no missing or extra coordinate. All files contained one finite,
two-dimensional primary image and passed unit-flux normalization. Their total
size is 13,947,840 bytes. Stamp shapes were 43×43 (379), 43×41 (125), 41×41
(67), and 41×43 (173) pixels. Sums ranged from 0.9999999999999982 to
1.0000000000000018. The maximum positive-flux centroid offset from the stamp
center was 0.2146 pixel. Negative interpolation pixels were present (372,203
in total); preserve the signed PSF values in convolution rather than clipping
them. The audit output marks `science_measurement=false` and confirms that no
image, shape, or photo-z rows were returned.

The audit is reproducible with
`audit_hsc_psf_picker_batch.py --plan-dir <approved-plan> --download-dir <received-files> --output <new-audit-json>`.
It verifies the immutable coordinate/mapping checksums, one-to-one coordinate
coverage, FITS structure, finite pixels, normalization, and centroid offset.
It refuses to overwrite an existing audit output.

## Scientific limits and next gate

This completes retrieval of position-dependent *model PSF cores* for the
existing injection positions; it does not validate extended PSF wings, image
noise, or the injection-and-recovery result. The model stamps are only about
7 arcsec across, while ICL radii span much larger angular scales. Retain the
previously obtained extended-PSF reference for wing sensitivity and separately
test its applicability. Before using these PSFs in injection/recovery, verify
the exact request-to-file mapping on a small set, preserve negative pixels,
and compare recovery across PSF and noise treatments. No further HSC data
request is included in this authorization. No ICL or weak-lensing science
measurement was made (`science_measurement=false`).
