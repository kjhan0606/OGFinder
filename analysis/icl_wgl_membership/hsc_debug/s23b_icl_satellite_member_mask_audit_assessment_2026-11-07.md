# S23B CAMIRA member-position mask audit — 2026-11-07

## Scope and validation

The audit compared all 45,634 imported S23B CAMIRA member positions with the
already-created target mask products. It did not download data, contact HSC,
or read science-image pixel values. It read only the saved mask planes.

The member-position table and 1,039 successful mask records were pinned by
SHA-256. The audit verified every one of the 1,039 mask-file hashes and
matched its stored image SHA-256 identifier to the member import. Reprojecting
each member coordinate with the saved mask WCS reproduced the imported
position to a maximum difference of 1.08e-7 pixel.

Four target images have no successful mask product because the preparation
runner reported “Too few valid pixels.” Those targets contain 292 member
rows; 255 are geometrically inside their image pixel-center bounds but cannot
be sampled against a saved mask. Among the other 1,039 targets, 38,911 member
centers were sampled and 6,431 member positions fall outside the image bounds.
Across all targets, 6,468 of 45,634 positions are outside the footprint.

## Member-center mask overlap

A member position is sampled at its nearest pixel center. Fractions below are
among in-footprint member positions with a successful mask product.

| Filter | Positions | BASE_MASK | Kron 2.5 | Kron 3.5 | Kron 5.0 |
| --- | ---: | ---: | ---: | ---: | ---: |
| HSC-i | 8,679 | 48.9% | 98.5% | 98.8% | 99.1% |
| HSC-i2 | 30,232 | 41.4% | 98.4% | 98.7% | 99.1% |

The generic automatic source masks therefore exclude almost every CAMIRA
candidate center. This is expected when the mask includes detected galaxy
sources, but it means those combined masks cannot be applied unchanged to
member-centered satellite-light profiles.

The existing mask products contain PRIMARY, BASE_MASK, and the three combined
Kron masks; they do not retain SEP segmentation labels or per-detection mask
ownership. Consequently, this audit cannot distinguish a candidate's own
galaxy mask from a neighboring-source mask, nor can it safely restore the
candidate's own pixels while keeping neighboring objects excluded.

## Image-area context

The member-center BASE_MASK fractions are comparable to or below the median
whole-image BASE_MASK area fraction for the corresponding target sample. The
median per-image masked area and q10–q90 range are:

| Filter | BASE_MASK | Kron 2.5 | Kron 3.5 | Kron 5.0 |
| --- | --- | --- | --- | --- |
| HSC-i | 55.2% (35.4–72.0%) | 64.2% (47.2–78.3%) | 70.1% (55.3–81.7%) | 79.8% (68.5–87.7%) |
| HSC-i2 | 46.5% (32.9–64.3%) | 56.6% (44.7–71.9%) | 63.8% (52.4–76.3%) | 75.0% (65.2–83.9%) |

These whole-image fractions are a coarse context, not a spatially matched
random-position test: the masks are spatially structured and member positions
concentrate near cluster centers. The large masked fractions reinforce the
need to model and validate usable area in any future profile stack.

## Assessment and next gate

The mask-coordinate join is sound, but the saved combined source masks are not
fit for direct use in a satellite outer-light stack. The next preparation
product should retain object segmentation labels and a detection table, match
CAMIRA candidate coordinates to detections, and construct per-candidate masks
that preserve the candidate galaxy while excluding unrelated sources. It
should then measure masked-area support in member-centered annuli and matched
control positions. Injection/recovery and the sky, PSF-wing, and correlated
noise gates remain necessary afterward.

This audit establishes neither satellite detections nor a physical ICL
relation. CAMIRA candidates remain photometric; the mask system is automatic;
and `science_measurement=false` and `images_ready_for_science=false` remain.

## Reproducibility

Runner:
`analysis/icl_wgl_membership/hsc_debug/audit_s23b_camira_member_mask_overlap.py`
(SHA-256 `c494bda32a606d17bee12780c7271f916a201e0267e29f9353f4a24d18657750`).

Output directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_member_mask_audit_2026-11-07_02/`

- Summary SHA-256: `98f4d4e5abe043457a19cba229c2b7c0a7f31e3e07ceeb25d7b5e84d1d713710`
- Per-member table SHA-256: `b9a807a01d254954c7014e9d47f20507050049daadc39858eec052162a08419b`
- Per-target table SHA-256: `360378b88942bb2c0d76171d1fe09185d9e4078fddddbad159c4feb6075a6164`

All products are audits only and retain `science_measurement=false`.
