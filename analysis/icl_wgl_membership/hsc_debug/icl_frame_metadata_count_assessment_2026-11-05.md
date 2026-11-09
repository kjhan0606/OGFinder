# HSC PDR3 coadd frame-count audit (2026-11-05)

## Verified query and result

The user approved count-only plan v3, SHA-256
`3dac4db631d0dcbee5c5f4509e639f2bd90ca2af99fa7cf5f66c2a9f60bd9014`.
Astrafek submitted its 16 batches to HSC PDR3 and verified exactly one returned
row for each of the 800 target and 798 matched-control image centers. All 1,598
centers map to a containing HSC-I patch. The result returns patch identifiers
and frame/visit/CCD counts only; no image, shape, photo-z, gain, exposure, sky
value, or credential was returned or saved.

| Per-center count | Targets: median (p10–p90), mean | Controls: median (p10–p90), mean |
|---|---:|---:|
| Unique frames | 42 (30–54), 42.05 | 42 (31–53), 41.89 |
| Unique visits | 9 (7–12), 9.52 | 9 (7–12), 9.50 |
| Unique CCDs | 38 (28–45), 36.88 | 37 (28–44), 36.67 |

Among the 798 matched target-control pairs, the median target-minus-control
difference is zero for all three counts. Frame counts are equal in 36.8% of
pairs; targets have more in 33.8% and controls in 29.3%. Visit counts are
equal in 49.6% of pairs. This indicates no gross central coadd-coverage offset
in these integer counts. It does not establish matched pixel depth or Poisson
variance: per-frame exposure, gain, sky variance, coadd weights, masks, and
pixel-level overlap are absent from this query.

## Scientific use

The counts can guide the next variance audit and identify unusual low-coverage
centers. They cannot be converted directly to an ICL shot-noise map or used to
rescale the stack. Any row export of frame-level calibration/variance values
is a separate data request and needs its own exact plan and approval. The
count result remains exploratory with `science_measurement=false`.

Reproduce with [`audit_icl_frame_metadata_count_results.py`](audit_icl_frame_metadata_count_results.py).
Verified row output and summary are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_count_audit_2026-11-05_02/`.
The summary SHA-256 is
`63e7424cb6abe3f89b7138bb99ea1072fbb367ab1e18b014974981334522810b`.
