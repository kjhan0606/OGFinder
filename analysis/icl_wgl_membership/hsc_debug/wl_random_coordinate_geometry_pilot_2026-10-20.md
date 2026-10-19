# Partial HSC random-coordinate geometry pilot — 2026-10-20

## Scope

The approved 800-lens coordinate export is still running. This pilot consumes
only the two completed batches already present in the immutable export plan;
it does not resubmit a CAS query or wait for later batches. It verifies the
per-lens candidate counts and file hashes, then applies the same 3-proper-Mpc
aperture-overlap and exact-coordinate-duplicate rules as the full geometry
screen. The 50 lenses are the first two plan batches and are not a
representative subsample of all 800 lenses.

The Astrafek implementation is
[`screen_wl_random_completed_batches.py`](screen_wl_random_completed_batches.py),
run by [`run_wl_random_partial_geometry.sbatch`](run_wl_random_partial_geometry.sbatch).
Slurm job 428780 completed successfully in 3 seconds. Ten focused tests passed,
including the completed-batch hash/count checks and source-manifest cosmology
checksum. The pilot output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_geometry_partial_2026-10-20_04/`.

## Results

| Quantity | Pilot result |
|---|---:|
| Completed coordinate batches | 2 of 32 |
| Target lenses represented | 50 |
| Candidate coordinate rows | 33,514 |
| Geometry candidates after overlap and duplicate cuts | 2,129 (6.35%) |
| Rows overlapping at least one CAMIRA 3-Mpc aperture | 31,254 |
| Rows duplicating an exact coordinate elsewhere in the pilot | 2,307 |
| Rows meeting both rejection conditions | 2,176 |
| Lenses with at least one geometry candidate | 31 of 50 |
| Lenses with at least 8 geometry candidates | 24 of 50 |
| Per-lens candidate count, median [min, max] | 5.5 [0, 436] |

The overlap and duplicate counts intersect and must not be added. Most rows
were rejected by the CAMIRA-aperture overlap rule. Candidate availability is
highly uneven across the 50 represented lenses, and 19 have no geometry-clean
coordinate in this partial result. This is evidence that the local 30-arcmin
candidate annulus may not provide the same random-control support for every
target; the remaining batches are required to measure that across the full
sample.

## Limits and next gate

This is a pipeline pilot, not a final random-lens pool. The full coordinate
export has not reached a terminal state, and the first 50 plan lenses may not
represent the remaining 750. No random centers were selected. Full-aperture
image coverage, mask fraction, HSC-I depth, PSF, sky, extinction, overlapping
random apertures, and spatial-block covariance have not been checked. No
random-source shape or photo-z rows and no images were requested. Both
`science_ready` and `science_measurement` remain false.

The authenticated submitter remains alive. On 2026-10-20 04:45 KST, PID
1049265 had been running for about 17 h 37 m. The checkpoint still showed
2/32 batches complete and CAS job 1062312 submitted; its mtime remained at
2026-10-19 11:07:53 KST. The live process has the older submitter code loaded:
it does not persist each successful CAS status poll, so the current CAS status
is unknown. Keep the live process and do not resubmit the batch. After a
verified complete export, run the existing 780k-coordinate geometry screen
and assess per-lens support before selecting controls. If the complete sample
remains sparse, prepare a new exact coordinate plan with a wider search region
and obtain separate approval before submitting it.

## Provenance

- Coordinate-plan SHA-256: `e8672433abdfe97c3f612c9327d9c449aa5c966b888702313bc15d4a95004827`.
- Completed batch 000: 15,879 rows; SHA-256 `a14b74ea174e790cd57dcaa4516f8dd95fe446df227f0b1f5323f8839fffd8d9`.
- Completed batch 001: 17,635 rows; SHA-256 `3119eacb04de443ad3c8b99cd57f5064c545e99685b58142e916caf2ae59628e`.
- Summary SHA-256: `28ae29da34dc71a8f966561eff537664ea8ac53d63707fd1831b7e390e2a6075`.
- Candidate audit SHA-256: `ef0cafc844df5af6db2db870640cb43eed0d05837f483caab41dd6bc8c220e46`.
- Per-lens table SHA-256: `1e84407b73d8d6587e672108670fbbd14049d7867c12b6505fc14ff7d5cfbc51`.
- Pilot screen SHA-256: `c9c8ad953163af43d6b39ed825980452ce9d18f2f88a8932a5531027b34a6901`.
- Geometry function SHA-256: `ca67cc7227994499835da49dd432223260473da019a9937693bbd93d009e4198`.
