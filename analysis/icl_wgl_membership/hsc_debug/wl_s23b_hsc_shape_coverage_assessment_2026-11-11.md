# S23B HSC weak-lensing shape coverage assessment (2026-11-11)

## Measured S19A coverage

The approved count-only audit evaluated 427 S23B CAMIRA centers that all occur in the 1,043-target ICL image manifest. Within a 3 proper-Mpc cone, 329 centers had zero HSC S19A shape-catalog rows and 98 had at least one. Only 82 reached 1,000 raw rows. Of those, 80 had five unique matched-control candidates under the current separation rule; two failed the five-control requirement.

| Selection stage | Centers | Raw lens-source rows |
|---|---:|---:|
| Count-only S23B candidate set | 427 | 356,191 |
| Positive S19A shape count | 98 | — |
| At least 1,000 raw shapes | 82 | — |
| Five-control candidate requirement passed | 80 | 340,256 |
| After shape and integrated support cuts | 74 | 239,565 retained lens pairs |
| Fixed support across seven annuli, ≥4/5 controls | 52 | 138,849 lens-source pairs after bin/calibration-input selection |

The approved 80-lens row export contains 340,256 of the 356,191 raw S19A lens-source count rows in the 427-center set (95.5%). The remaining candidate centers therefore do not provide a path to hundreds of S19A weak-lensing lenses: most have no S19A source rows, and the remaining positive-count subset is small. The fixed 52-lens profile is a useful exploratory stack, but its sample cannot support a hundreds-cluster S23B weak-lensing relation under this shape catalog.

## Data-release implications

HSC lists S16A and S19A as its released weak-lensing shape-catalog products. The dedicated S19A page describes the S19A/Y3 catalog and its use for weak-lensing analyses. These documents do not provide a newer general PDR3 shape catalog that fills the missing S19A lens positions. An S16A overlap would require a separate position and calibration audit before combining it with S19A; simply concatenating the catalog rows would not be a valid expansion.

For the stated ICL–satellite objective, the existing 1,043 HSC ICL targets and CAMIRA member-position catalog can continue as a separate analysis sample. For weak lensing, retain the 52-lens S19A stack as an exploratory pilot. A several-hundred-lens S23B weak-lensing sample would require either a validated additional shape catalog with compatible calibration or an image-based shape-measurement and simulation-calibration pipeline. Any new remote query or data transfer remains subject to the user's explicit approval.

## Reproducibility and scope

The coverage assessment was prepared by [prepare_s23b_wl_common_support_plan.py](prepare_s23b_wl_common_support_plan.py) from the approved S23B count-only audit and the saved ICL target manifest. It made no new CAS request or download. All shape-count, profile, and overlap products retain `science_measurement=false`.

- Count-only plan SHA-256: `b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984`
- Count-only results summary SHA-256: `748fea398f144165a25210071a5f6ab1866792df3be231ccfdb75282829935b8`
- Local S23B support assessment SHA-256: `15b30e2436e0abb4d7b23013f83313b1e434b5dab89a441dc1dd5dea231369c6`
- ICL target manifest SHA-256: `e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44`
- Approved shape-row plan SHA-256: `7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`

## HSC references

- [HSC weak-lensing results and data products](https://hsc-release.mtk.nao.ac.jp/doc/index.php/weak-lensing-results-and-data-products/)
- [HSC S19A Shape Catalog (PDR3)](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/)
- [HSC PDR3 available data](https://hsc-release.mtk.nao.ac.jp/doc/index.php/available-data__pdr3/)
