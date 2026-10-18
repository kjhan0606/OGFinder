# Paired-field ICL plane-estimator injection audit — 2026-10-19

## Scope

Astrafek tool `audit_icl_plane_injection.py` injects four elliptical exponential
diffuse profiles into the existing HSC coadd/bg science arrays for the fixed
74-pair ICL/satellite support sample. Each profile uses scale length 100 or
250 kpc and axis ratio 1.0 or 0.65, with 1 nJy/arcsec² at elliptical radius
100 kpc. The same K2.5 mask, 200–300 kpc measurement annulus, delivered target
and matched random-control fields, and connected target/control sky blocks are
used throughout.

For each original and injected image, the tool recomputes Astrafek's 600–900
kpc tile-plane estimate and measures both raw and plane-corrected annular
means. It reports target and control recovery fractions relative to the
injected light actually visible through the fixed annular mask, plus paired
target-minus-control block-bootstrap intervals. Source detections and masks
are held fixed; no source Poisson noise or PSF convolution is added.

This tests downstream local-plane transfer on the delivered global-sky
coadd/bg images. It does not test HSC's upstream sky subtraction before these
images were delivered, nor does it validate PSF wings, mask growth,
completeness, extinction, or correlated coadd noise. Results remain
exploratory and cannot be used as a science correction or final covariance.

## Execution record

The run uses only already downloaded target/control images and manifest-linked
masks. No HSC query or download is part of this job. The first Slurm attempt,
427839, measured all 74 pairs but failed in final aggregation because the
shape-group check compared model rows with pair rows. The failure checkpoint
is preserved at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_plane_injection_2026-10-19_01/`;
its progress file records `stopped_on_error`. Astrafek now tests shape grouping
against the fixed per-filter pair count. The corrected rerun uses output folder
`icl_plane_injection_2026-10-19_02`. Job 427845 completed successfully in
4m50s (peak RSS 571,292 KiB). It produced 296 paired rows, representing 74
target/control pairs times four injected shapes. All six recorded code hashes
match the current Astrafek tools; the output-row hash and small manifest/input
hashes were independently verified. The analysis passed its four targeted
tests before submission.

## Results

Raw annulus recovery is 1.0000 for all shapes, as expected from the direct
unmasked-annulus measurement. Recomputing the local plane reduces recovery in
both cluster and control fields. Values below are mean recovery fractions;
the final column is the 95% connected-block bootstrap interval for the paired
mean difference, target minus control.

| Filter | Profile | Pairs | Cluster plane recovery | Random-field plane recovery | Paired difference (95% interval) |
|---|---|---:|---:|---:|---:|
| i2 | 100 kpc, q=1.0 | 53 | 0.9911 | 0.9911 | 0.0001 (−0.0002, 0.0004) |
| i2 | 100 kpc, q=0.65 | 53 | 0.9954 | 0.9950 | 0.0004 (−0.0001, 0.0010) |
| i2 | 250 kpc, q=1.0 | 53 | 0.8598 | 0.8593 | 0.0005 (−0.0012, 0.0027) |
| i2 | 250 kpc, q=0.65 | 53 | 0.9077 | 0.9055 | 0.0022 (−0.0024, 0.0075) |
| i | 100 kpc, q=1.0 | 21 | 0.9911 | 0.9912 | −0.0001 (−0.0006, 0.0007) |
| i | 100 kpc, q=0.65 | 21 | 0.9951 | 0.9948 | 0.0004 (−0.0006, 0.0015) |
| i | 250 kpc, q=1.0 | 21 | 0.8596 | 0.8600 | −0.0004 (−0.0044, 0.0048) |
| i | 250 kpc, q=0.65 | 21 | 0.9058 | 0.9020 | 0.0038 (−0.0047, 0.0139) |

All paired target-minus-control intervals include zero. Under these four
analytic profiles, the downstream plane estimator removes about 0.5–0.9% of
the 100 kpc-scale signal and about 9–14% of the 250 kpc-scale signal. The
similar target and random-field transfer suggests little differential bias
for these matched fields and shapes; it does not justify applying a universal
correction to observed ICL or satellite associations.

The authoritative summary and row-level results are
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_plane_injection_2026-10-19_02/summary.json`
and `transfer_rows.tsv`. `science_measurement=false` remains in force. The next
gate is full-pipeline source re-detection/mask response, PSF-wing and source
noise injection, broader ICL morphologies, and complete spatial covariance.
