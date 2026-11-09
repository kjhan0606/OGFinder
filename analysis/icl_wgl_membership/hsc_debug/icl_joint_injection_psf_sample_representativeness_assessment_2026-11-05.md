# Representativeness of the PSF-supported joint-injection subset

## Scope

This offline Astrafek audit compares the 57 cluster/filter pairs with existing
position-dependent PSF-convolved joint-injection results against the 224-pair
parent with ICL 200–300 kpc coverage support. It reads saved profile and
summary tables only, makes no HSC request, reads no image pixels, and retains
`science_measurement=false`.

## Results

| Filter | PSF subset / parent | Δ mean redshift (SMD) | Δ mean richness (SMD) | Δ mean ICL coverage (SMD) | Satellite support 60%, subset / parent |
|---|---:|---:|---:|---:|---:|
| `i` | 13/71 | +0.0128 (+0.179) | −1.372 (−0.128) | +0.0812 (+1.061) | 12/13 vs 24/71 |
| `i2` | 44/153 | +0.0272 (+0.378) | −1.752 (−0.260) | +0.0602 (+0.736) | 42/44 vs 74/153 |
| All | 57/224 | +0.0244 (+0.340) | −1.879 (−0.229) | +0.0649 (+0.800) | 54/57 vs 98/224 |

The PSF-injection subset has substantially higher ICL-annulus coverage and
outer-satellite support than the parent. Its measured injection transfer
therefore cannot be generalized to all 224 pairs. The 57-frame amplitude
response check is a useful local transfer diagnostic; it is not a calibrated
selection or completeness correction.

## Reproduction and provenance

- Astrafek script: [`audit_joint_injection_psf_sample_representativeness.py`](audit_joint_injection_psf_sample_representativeness.py)
- Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_joint_injection_psf_sample_representativeness_2026-11-05_01/`
- Summary SHA-256: `e9bd84a6ccb309fc57b405edaca18c0c0c71f91702fda88364987a6e4fb6bbf3`
- Pair selection table SHA-256: `1d25da036cfdeeecb75337738dde2881aa515a06c56fcb8520a097141c59c20d`
- Analysis script SHA-256: `76486d5da014c13c48ceeb1a966aeb8db01a8bdea2fc08f4e600df4c7a17a298`

The run used no HSC network and did not read image pixels. It does not correct
selection and does not establish population representativeness.
