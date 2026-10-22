# Weak-lensing random-coordinate export recheck — 2026-10-22

## Current evidence

The local, hash-checked checkpoint audit at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_coordinate_checkpoint_audit_2026-10-22_01/`
verifies 2 of 32 approved coordinate batches: 33,514 of 780,191 planned
lens–candidate rows. Batches `000` and `001` match their recorded sizes,
per-lens counts, and SHA-256 values. Batch `002` remains recorded as submitted
under HSC CAS job `1062312`; batches `003`–`031` remain pending.

The last successful persisted CAS status poll for job `1062312` says `running`
at `2026-10-22T07:34:55.736381Z`. A one-time process check on `lageunha` found
that the former authenticated submitter PID `1049265` no longer exists. This
only establishes that the submitter is absent; it does **not** establish
whether HSC CAS job `1062312` is still running, complete, or failed.

## Decision

Do not resubmit batch `002` or start batches `003`–`031` until the existing CAS
job is checked through the authenticated HSC status endpoint. No HSC query or
data transfer was made during this recheck. The original user-approved scope
remains limited to the 800-lens random-coordinate columns
`lens_id, ra_deg, dec_deg`; no source rows, shapes, photo-z values, or images
are included. Random-lens subtraction and weak-lensing science measurement
remain incomplete.

## Reproduction

- Astrafek audit:
  [`audit_wl_random_coordinate_checkpoint.py`](audit_wl_random_coordinate_checkpoint.py)
- Output JSON SHA-256:
  `3950611b12e217aa6aedbef0c464970bb1f6c7b37eb45e7be9b39789931af9d1`
- Input checkpoint SHA-256:
  `d9bf46d1e016596950dcdd6afa1b6d8c07fe13eae59979a6989a3c036cb495f0`
- Approved coordinate plan SHA-256:
  `e8672433abdfe97c3f612c9327d9c449aa5c966b888702313bc15d4a95004827`
- Status: `verified_partial`; `science_measurement=false`.
