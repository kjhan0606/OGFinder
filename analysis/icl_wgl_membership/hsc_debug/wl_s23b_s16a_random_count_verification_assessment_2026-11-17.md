# S16A random-center count verification and approval boundary

The existing approved count-only run completed normally on 17 November 2026.
The scope was 570 random centers for 114 selected lenses, five per lens,
3 proper-Mpc cones and 23 sequential SQL batches. The exact approved plan
SHA-256 is `dfd2fc75650fc08ce4acfe1eef6d0ed476a7bdf8c84b5d083b607f4009678de0`.
The interrupted batch was resumed by its recorded CAS job ID. No source rows,
source coordinates, shapes, photo-z values, images or PSFs were returned by
this count-only run. `science_measurement=false`.

## Verification result

An independent standard-library audit checked the canonical approved plan,
all SQL hashes and sizes, the exact completed-job set, all 23 CSV hashes and
sizes, the two returned columns, nonnegative integer counts, exact batch IDs,
cross-batch uniqueness and complete coverage of all 570 planned centers. It
reconstructed the combined TSV byte for byte from the original CSVs and checked
the combined table, summary and checkpoint hashes and arithmetic. The audit
made no network requests and preserved its own code beside its JSON result.

| Returned quantity | Value |
|---|---:|
| Completed SQL batches | 23 / 23 |
| Unique returned random centers | 570 / 570 |
| Lenses, with exactly five planned centers each | 114 |
| Centers with positive raw counts | 525 |
| Centers with zero raw counts | 45 |
| Minimum / maximum raw cone count | 0 / 9,920 |
| Sum of raw cone counts | 2,442,267 |

The sum can count the same source in overlapping cones. It is **not** a count
of unique sources across the survey or a count surviving the scientific cuts.

| Positive-count random centers out of five | Number of lenses |
|---|---:|
| 5 | 96 |
| 4 | 8 |
| 3 | 1 |
| 2 | 5 |
| 1 | 0 |
| 0 | 4 |

Thus 96 lenses have nonzero raw counts at all five planned random centers.
The other 18 need a support diagnosis; four have zero raw counts at all five.
These are diagnostics, not a revised selected sample. No center was replaced,
excluded or promoted by this verification. Zero catalog counts do not by
themselves identify an imaging defect, and positive counts do not establish
aperture coverage, matched depth, retained shape support or calibration.
S16A measurements cannot be promoted on these counts or on the existence of
the separate fixed-52 analysis; release-specific cuts and calibration still
need to be established.

## Preserved evidence

- Original plan/results directory:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_expansion_random_count_plan_2026-11-11_01/`
- Independent audit directory, containing `audit.py` and `verification.json`:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_random_count_verification_2026-11-17_01/`
- Audit JSON SHA-256:
  `8394aeeda9b1d4f885b2a884bd7be14fde9278ed2dcf750aebe8d601d937578f`
- Audit code SHA-256:
  `d1c826be19d4dad619152fc3e9b4eed86718996cc45a28b41d1707356e5be90d`
- Count summary SHA-256:
  `c3a92b0420c5c8af8cb6191ec1264305e3fc38e503365b456f5dcad3f07332ef`
- Combined TSV SHA-256:
  `91f7521d49dc8c734fadf375f80fe58d610ff72b851569c2c38e3a5dc6eb2b77`

The submitter exited zero, the checkpoint status is `complete`, every job is
complete and `current_query` is null. The checkpoint retains the historical
`error_type=KeyboardInterrupt` marker from the deliberate earlier pause; that
marker is not a failure of the completed resumed run. No credentials were
saved.

The same final verification also matched six pinned existing PSF/injection/WL
summary hashes and all 438 listed output-file hashes, including 427 PSF cores.
This verifies preservation of those products, not additional physical
calibration. PSF wings, sky transfer, source Poisson/coadd weights and spatial
covariance remain unresolved. The separate web replay verification corrected
missing-required-output handling; its ten focused checks and real browser ZIP
to native-Python E2E passed. Details are in the respective project handovers.

## Next stage: await user approval

Proposed next work uses existing local data only: compare the raw target and
random counts, existing footprint/selection provenance and per-lens support;
diagnose the 18 incomplete five-random groups; and check the available S16A
calibration path before drafting any later source-row transfer plan. Preserve
the original 570-center ledger and report any proposed changes separately.
No new CAS query, image/PSF download or source-row export is included in this
proposed local stage. A later transfer needs its own exact scope and approval.

The user requested completion of the current verification followed by an
approval boundary. The proposed next work has **not** started. Await approval;
do not automatically launch a new analysis, selection or transfer.
