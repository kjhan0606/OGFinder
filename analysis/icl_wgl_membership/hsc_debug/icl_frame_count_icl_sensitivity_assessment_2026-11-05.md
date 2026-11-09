# HSC frame-count sensitivity of the ICL target–control contrast (2026-11-05)

## Question

The approved HSC PDR3 count-only query found similar median frame, visit, and
CCD counts at 800 target and 798 matched-control centers. This offline test
asks whether the *paired count offsets* nevertheless track the saved ICL
surface-brightness contrast in the 200–300 kpc aperture-supported sample.

## Method

The Astrafek runner joins the verified per-center frame-count table to the
existing 200–300 kpc target-minus-control ICL profile for 224 cluster/filter
pairs (71 `i`, 153 `i2`). It computes the partial Pearson correlation between
the ICL contrast and each target-minus-control count difference (unique frames,
visits, and CCDs), controlling for cluster redshift, log richness, and common
aperture coverage. Each filter is analyzed separately. Uncertainty intervals
come from 5,000 resamples of connected 5-degree target/control endpoint blocks.
The runner reads saved tables only; it made no HSC request and read no image
pixels.

## Results

All six point estimates are near zero and all block-bootstrap intervals include
zero:

| Filter | Count offset | Partial *r* | 95% connected-block interval | Pairs / blocks |
|---|---|---:|---:|---:|
| `i` | Unique frames | −0.085 | [−0.312, 0.233] | 71 / 18 |
| `i` | Unique visits | −0.061 | [−0.296, 0.148] | 71 / 18 |
| `i` | Unique CCDs | +0.078 | [−0.145, 0.311] | 71 / 18 |
| `i2` | Unique frames | +0.017 | [−0.122, 0.172] | 153 / 40 |
| `i2` | Unique visits | −0.041 | [−0.150, 0.083] | 153 / 40 |
| `i2` | Unique CCDs | −0.008 | [−0.126, 0.116] | 153 / 40 |

This provides no evidence in this supported subset that *integer frame-count
imbalance* drives the ICL contrast after the listed adjustments. The `i`
intervals remain broad. A null count correlation does not show that target and
control pixels have equal depth, sky variance, coadd weighting, or source-noise
properties.

## Limits

The input count is the number of CCD-frame memberships in the nearest HSC-I
mosaic patch. It is not a pixel-level exposure or inverse-variance map. The
profile sample is selected by ICL-aperture coverage and is not the full 800
cluster parent. The intervals are conditional on the saved connected 5-degree
endpoint partition. This diagnostic does not establish absence of imaging
systematics, supply a variance correction, or support a causal interpretation.
Keep `science_measurement=false`.

## Reproducibility

Run `analyze_icl_frame_count_icl_sensitivity.py` on the saved manifests,
profile table, and approved count-only audit. The joined input/output table and
machine summary are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_count_sensitivity_2026-11-05_01/`.

- Runner SHA-256: `de0ae9afbe5c6261b24009a9b9b6cdc0f9b24c6d4832e47f26d1eea9968b6bd4`
- Joined table SHA-256: `470b30cdd0d398bc7f80c4739a2e3c0f367e8e899056e02809a789b12e5a9dcb`
- Summary SHA-256: `dd00410649e230359df6a6ea6582be2df6881e247bf81d29347bdca33d569cca`
- Count-only parent plan SHA-256: `3dac4db631d0dcbee5c5f4509e639f2bd90ca2af99fa7cf5f66c2a9f60bd9014`
