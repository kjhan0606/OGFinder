# Independent audit: hybrid-PSF satellite-shell exclusion — 2026-11-06

## Scope and reproducibility

Re-audited the existing synthetic-injection outputs for the fixed ICL–satellite
association sample. The offline summarizer verified the pinned 74-row
association table (SHA-256
`483f81540543fc46ec1cf4333bbd3e238799b493e48fe03235da3038f5714917`), 12/12
`i` frames, 36/36 `i2` frames, all per-frame JSON hashes, and 48 unique
target/filter keys. No HSC request, download, or observation-image read took
place during this audit.

The fixed synthetic model used a 2.4-arcsec splice radius, 0.6-arcsec blend,
and wing amplitude 0.65. The public HSC-I wing is not position-, epoch-, or
seeing-matched; its use for `i2` is only a proxy sensitivity. These are not
approved observational PSFs.

## Results

| Filter and injected component | Fixed-mask median response | Shell-excluded median | Paired q16–q84 change |
|---|---:|---:|---:|
| `i`, diffuse → ICL | 0.22538 | 0.22472 | −0.00167 to +0.00102 |
| `i`, satellite envelope → ICL | 0.46817 | 0.46687 | −0.01379 to +0.00274 |
| `i2`, diffuse → ICL | 0.22568 | 0.22602 | −0.00110 to +0.00132 |
| `i2`, satellite envelope → ICL | 0.68183 | 0.65932 | −0.06592 to +0.01161 |

The diffuse-response changes are small in the paired medians and their q16–q84
sample ranges include zero. Satellite-envelope transfer is highly dispersed
between frames; its paired ranges also include zero. No scalar transfer
correction is justified.

Excluding supported satellite-shell pixels removed 62,858 `i` and 257,890
`i2` pixels. Median 200–300 kpc mask coverage changed from 0.681 to 0.618 in
`i` and from 0.652 to 0.611 in `i2`; the support loss must accompany any
interpretation of the reduced model point-source wing fraction. The median
unit-point-source wing fraction changed from 0.000994 to 0.000902 (`i`) and
0.000970 to 0.000896 (`i2`).

## Interpretation

The audit verifies the saved model-sensitivity calculations and their sample
membership. It does not validate observed ICL fluxes. The data still lack
source-Poisson calibration, complete coadd and sky covariance, matched-control
uncertainty, a filter-specific validated `i2` extended PSF, and uncertainty
from reselecting masks. Keep `science_measurement=false`,
`psf_approved=false`, `observational_validation_passed=false`, and
`flux_covariance_available=false`.

## Artifacts

- Runner: [`summarize_icl_hybrid_shell_exclusion.py`](summarize_icl_hybrid_shell_exclusion.py)
- Machine summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_hybrid_shell_exclusion_audit_2026-11-06_01/summary.json`
- Summary SHA-256: `599006d4f9d2d0f114df413defe3f9040a0c8bfa96fadf08ccf72509534c9aa9`
- Runner SHA-256: `f94191ddee4b9ad163a0ba02587cb67d885e65d359a2f9cb024baad57a2865b8`
- Source batch summaries: `icl_hybrid_shell_exclusion_association_i12_2026-11-03_01/summary.json` and `icl_hybrid_shell_exclusion_association_i2_36_2026-11-03_01/summary.json` under the same data root.

Reproduce from the OGFinder project root:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/summarize_icl_hybrid_shell_exclusion.py \
  --association-table /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_outer_sky_association_reblocked_2026-11-03_01/association_rows.tsv \
  --i-summary /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_hybrid_shell_exclusion_association_i12_2026-11-03_01/summary.json \
  --i2-summary /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_hybrid_shell_exclusion_association_i2_36_2026-11-03_01/summary.json
```
