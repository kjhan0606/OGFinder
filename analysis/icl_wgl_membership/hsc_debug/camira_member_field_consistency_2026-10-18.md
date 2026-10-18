# CAMIRA S20A member-field consistency audit

**Date:** 2026-10-18 (Asia/Seoul)<br>
**Status:** schema consistency evidence only; `science_measurement=false`<br>
**Conclusion:** column 7 is *likely* log stellar mass and column 8 behaves like a membership weight, but the exact column order for this headerless S20A file remains unverified.

## Why the interpretation is plausible, but not definitive

The local file `camira_s20a_wide_sm_member.dat.gz` has no header. The acquisition record matches its SHA-256 (`bf0c1fec5e11c6378ea69e3a65e49e9b32beecc83e5e68000bfa997e4ef2059e`) and labels its final columns `field7` and `weight`; that label is a parser choice, not a source schema.

CAMIRA literature describes member catalogs containing stellar mass `M*` and membership probability `w`; it explains that `w` combines the number parameter, stellar-mass filter, and spatial filter. The richness is a weighted member count above roughly `10^10.2 M_sun` within `1 h^-1` physical Mpc. The detailed description located is for HSC S21A, not this exact S20A file, so it cannot establish the column order here. S20A analyses likewise describe cluster stellar mass as the member stellar masses summed with membership-probability weights. See the [CAMIRA S21A catalog description](https://doi.org/10.3847/1538-4357/ad76ad) and [S20A BCG-mass study](https://arxiv.org/abs/2108.11288).

The selected-file values are consistent with that pair of quantities:

| Check | Result |
|---|---:|
| Selected clusters / member rows | 800 / 36,513 |
| Column 7 1st / 50th / 99th percentiles | 10.031 / 10.732 / 11.794 |
| Column 7 at or above 10.2 | 89.58% |
| Column 8 1st / 50th / 99th percentiles | 0.00122 / 0.4071 / 1.000 |
| Column 8 exactly 1 | 27.28% |
| Median `sum(column 8) / N_mem` across clusters | 0.923 |
| 10th–90th percentile of that ratio | 0.810–1.059 |
| Median absolute relative difference from `N_mem` | 8.97% |
| Clusters within 20% of `N_mem` | 91.0% |
| Median `sum(column 8 where column 7 >= 10.2) / N_mem` | 0.910 |
| Projected member radius median / 95th / 99th percentile / max | 0.536 / 1.313 / 1.413 / 1.446 Mpc |
| Members beyond `1 h^-1 Mpc` (using `h=0.7`) | 0.46% |

These checks are consistent with a log-mass field and a probability-like weight, and with members already selected near the CAMIRA aperture. The near-unity weighted-count comparison is not an identity test: masking corrections and catalog construction can change `N_mem`. Numeric range and behavioral agreement cannot rule out another field with similar values.

## Implication for the satellite–ICL analysis

Keep the existing labels `candidate_log10_mstar` and `field7` provisional. Do not present mass-bin trends as a physical stellar-mass result until the publisher's dictionary or companion source code maps column 7 for this exact S20A member file. The CAMIRA-weighted sensitivity remains useful as a robustness check because column 8 is explicitly named `weight` in the acquisition record and has probability-like values, but it is still not an independently calibrated satellite sample. Completeness, projection contamination, centering, masking, PSF/sky transfer, and profile covariance remain unresolved.

## Reproduction

Run from the OGFinder checkout on a host with the existing data mounted:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/audit_camira_member_field_consistency.py
```

The script verifies the catalog and acquisition-summary hashes, 800 cluster identities, each selected cluster's `N_mem` and `z_cl`, the weight range, and the expected 36,513 selected rows before printing the metrics. It reads existing files only; it does not contact HSC or export source rows.
