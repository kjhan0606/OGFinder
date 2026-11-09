# Diffuse ICL injection amplitude response (2026-11-06)

## Design

The same 137 position-dependent-PSF-convolved HSC frames (44 `i`, 93 `i2`)
were rerun at diffuse amplitudes 0.5, 1.0, and 2.0 nJy arcsec⁻² at 100 kpc.
Satellite injections stayed at 1.0 nJy arcsec⁻² per galaxy at 50 kpc. SEP
detection and source masks were regenerated for each run. The paired Astrafek
analysis verified identical frame IDs, input images, position-dependent PSFs,
baseline detections, and all non-amplitude injection settings.

## Results

| Filter | Input amplitude | ICL annulus response per unit input (median, q16–q84) | Paired change vs 1.0 (5° block 95% interval) |
|---|---:|---:|---:|
| `i` | 0.5 | 0.240 (0.189–0.275) | 0.019 [−0.010, 0.026] |
| `i` | 1.0 | 0.226 (0.178–0.249) | — |
| `i` | 2.0 | 0.224 (0.204–0.244) | 0.007 [−0.007, 0.014] |
| `i2` | 0.5 | 0.225 (0.178–0.278) | −0.002 [−0.013, 0.006] |
| `i2` | 1.0 | 0.226 (0.194–0.257) | — |
| `i2` | 2.0 | 0.226 (0.207–0.252) | −0.002 [−0.004, 0.006] |

The per-unit response is stable across the tested two-fold amplitude range:
all paired intervals include zero. The median absolute recovery remains about
0.23 of input. Regenerated source masks changed in most frames, so similar
response coefficients do not mean that the full measurement is mask-invariant.

## Interpretation and limits

This is a pipeline-transfer sensitivity on the existing coverage-selected
sample, not an astrophysical measurement or a completeness correction. The
experiment omits injected-source Poisson noise, per-pixel coadd weights,
upstream HSC sky processing and extended PSF wings. Its 5° connected cells are
an approximate partition, not validated independent HSC fields. Keep
`science_measurement=false` and do not divide the observed ICL stack by this
response as a calibrated correction.

## Reproduction and provenance

- Runner: [`analyze_icl_diffuse_injection_amplitude_grid.py`](analyze_icl_diffuse_injection_amplitude_grid.py), SHA-256 `8eef6823016e35538a059b0c125d12e88b98e9f8307c52cc21f88ceb04af7448`.
- Paired summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_diffuse_injection_amplitude_grid_supported137_2026-11-06_01/summary.json`, SHA-256 `538656e4b3c5141282ddd24c16ad8ce0d20454692228f38b437ad996b29fcd96`.
- Paired frame table SHA-256 `0920a85c932f0fc6addf746ee857bd44b45f5044b31935676169155d1aa0706d`.
- 0.5, 1.0, and 2.0 source run summary SHA-256 values, respectively: `ec5a48917a3397077717fd1bb6a4b2224b95e6a8295864ad1b33f6c095c58385`, `be3c2e8cb22185c06768d9022ae3e6192072f180fa54d55e21b85920e5572750`, and `e4ac20691975ed9e968417ae4628221333b600ff3cc4d08844c7333d1703483a`. Each run completed 137 frames with zero errors.

## HSC frame metadata export

The separately approved query plan SHA-256 `49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3` returned all 67,071 approved rows for 800 target and 798 matched-control centers. The independent offline audit found zero nulls across exposure, four gains, FLUXMAG0, sky RMS and MJD; matched-pair median differences were zero for frame, visit, CCD, summed exposure, sky RMS, FLUXMAG0 and all four gains. This is only a coarse patch-membership balance diagnostic: no per-pixel coadd weights or amplifier maps were returned, and it does not provide a source-Poisson correction.

- Submitter: [`submit_icl_frame_metadata_rows.py`](submit_icl_frame_metadata_rows.py), updated to accept the CAS original-frame filter codes `i` and `i2`.
- Independent offline auditor: [`audit_icl_frame_metadata_row_results.py`](audit_icl_frame_metadata_row_results.py), SHA-256 `d0340160ed0f7c2ed7c8e7570b2687ed1009917380d8c8dda1625607d5922336`.
- Audit summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_row_audit_2026-11-06_01/summary.json`, SHA-256 `0f361fcb496fbd7be7c1059f638414f6900535f73a9ac48696030480f400a237`.
- CAS summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_row_results_2026-11-06_01/summary.json`, SHA-256 `3925f74be34b19cf16ea0183055ba4d74d91568f66150c569e2d48b6bd2f7b53`.

The coadd and frame products differ in filter notation: the coadd selection is
`HSC-I`, while individual contributing rows report `i` or `i2`. The export did
not contain images, shape rows, photo-z, object rows or credentials. Keep
`science_measurement=false`.
