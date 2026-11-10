# S23B fixed-52 ICL–satellite–weak-lensing familywise audit

This analysis applies a restricted Webb six-point wild-cluster bootstrap with a shared 5-degree target/control block multiplier to the 116 partial-correlation tests in the saved fixed-52 relation. Four tests are direct ICL–satellite relations; 112 test the WL profile against ICL or satellite structure. The analysis reads only saved tables and makes no HSC/CAS request, image read, raw source-row CSV read, or download. `science_measurement=false`.

The fixed 52-lens sample forms 30 connected 5-degree blocks (largest block: 7 lenses). The full family uses 10,000 shared-multiplier replicates; 10,000 were valid. CR1 cluster-studentized OLS slopes are tested under the corresponding restricted zero-slope model. This is a complementary null-based familywise procedure; it does not adjust the percentile intervals in the relation report.

The earlier percentile-bootstrap report has 8 of 112 WL partial-correlation intervals that exclude zero pointwise. The wild-cluster max-T family has 10 unadjusted p-values below 0.05 and 0 max-T adjusted p-values below 0.05. Its smallest adjusted p-value is 0.5154 for i / support_ge_0p10 / wl_vs_satellite / tangential / radial bin 3 (partial r = -0.549).

Interpretation remains limited to this selected 52-lens sample. The 5-degree connected-cell partition is an approximate spatial dependence model, not a validated HSC survey covariance. The analysis does not resolve source photo-z, shear response, additive/multiplicative residuals, PSF, magnification, member dilution, masks, or source reuse. It does not establish a physical ICL–satellite–halo relation or a detection.

## Provenance

- Runner SHA-256: `52aa8140fdb433ce28b544efa918c78765ff28dabd9568025b9179dec0dda57f`
- Relation input runner SHA-256: `27839cf7bfa1683d94e6bf3dd30e2a7eac9a0f71d56c4bedeaa574beb3a41c8a`
- Fixed lens ID SHA-256: `18312f27139b3750abdbf71543f22d9b7b18b2d3b24cbd29ed686a377524bb7d`
- Family size: 116 (4 direct ICL–satellite, 112 WL association tests)
- Result directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_wl_fixed52_maxt_2026-11-11_01`
- Full per-test estimates and p-values: `max_t_metrics.tsv`
