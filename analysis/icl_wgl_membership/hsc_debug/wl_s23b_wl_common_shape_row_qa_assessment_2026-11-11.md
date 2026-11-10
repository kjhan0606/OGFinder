# S23B weak-lensing source-row quality and common-support audit

Approved row-plan SHA-256: `7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`  
Input HSC Y3 n(z) SHA-256: `a28f9e2e088e92d96d2d87083a6eeeac4c0c84c6b48e033bd3d1dc1bf0d8958f`  
Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/quality_audit_2026-11-11_01`  
No network request, additional download, or lensing stack was run.

## Source cuts and semantics

The same source cuts were applied to target lenses and matched controls: `|e|<2`, `0<=sigma_e<=0.4`, `mag_a10<=25.5`, `resolution>=0.3`, HSC Y3 `source_zbin=1..4`, exclude the verified raw B-mode value 0, positive fiducial HSC Y3 `n(z)` lensing efficiency, and proper radius 0.1–3 Mpc. The official HSC S19A page defines `b_mode_mask` as a large-B-mode region cut but does not document the stored 0/1 direction. Prior coordinate validation of the same release identified raw 0 as masked; the new sample's coordinate comparison to the published GAMA09H rectangle found **0 mismatches** (exact_with_published_GAMA09H_rectangle).
The exact coordinate cross-check found 48,765 raw-0 pairs inside the rectangle, 0 raw-1 pairs inside, 0 raw-0 pairs outside, and 1,863,944 raw-1 pairs outside.

Official definitions: [HSC S19A Shape Catalog (PDR3)](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/). Prior independent polarity audit: `wl_bmode_coordinate_validation_2026-10-19.md`.

## Retained pair support

Lens apertures: 340,256 input pairs; 240,068 pass the HSC shape, magnitude, resolution, tomographic-bin, and B-mode cuts; 240,068 also have positive fiducial n(z) efficiency; 239,565 remain after 0.1–3 Mpc and calibration-input-validity checks.
Five-control apertures: 1,572,453 input pairs; 1,107,161 pass the same HSC cuts; 1,107,161 also have positive fiducial n(z) efficiency; 1,105,779 remain after the same radial and calibration-input-validity checks.
Retained HSC Y3 z-bin counts (1/2/3/4): lenses 60,561/76,131/65,905/36,968; controls 254,424/366,090/314,763/170,502.
At the paired-lens level, 2 lenses have zero retained lens rows; 9 have at least one zero-support control; 0 have all five controls at zero.
The q16/median/q84 ratios of mean-five-control to lens are 0.690/0.925/1.311 by count, 0.686/0.929/1.311 by HSC weight sum, and 0.690/0.956/1.387 by geometry-weight proxy sum.
The table below shows how many nonzero-lens sets remain under a minimum number of positive-support controls. Ratios use only positive-support controls for each set; the table is a support sensitivity, not a selected estimator.

| Minimum positive controls | Eligible lens sets | Lens sets retained | Control/lens count ratio q16 / median / q84 |
|---:|---:|---:|---:|
| 3/5 | 74 | 92.5% | 0.726/0.940/1.315 |
| 4/5 | 74 | 92.5% | 0.726/0.940/1.315 |
| 5/5 | 71 | 88.8% | 0.740/0.950/1.319 |

Center-level `N_eff` q16/median/q84 is 1854.246/2964.564/3769.294 (lens) and 1556.506/2840.735/3676.239 (control) for catalog weight; geometry-weight proxy `N_eff` is 1248.329/1910.568/2802.954 and 1080.821/1970.442/2766.651.

Retained shape and weight distributions (q16 / median / q84):

| Retained quantity | Lenses q16 / median / q84 | Controls q16 / median / q84 |
|---|---:|---:|
| |e| | 0.236/0.529/0.818 | 0.238/0.532/0.819 |
| sigma_e | 0.037/0.106/0.216 | 0.037/0.105/0.214 |
| rms_e | 0.375/0.402/0.422 | 0.375/0.402/0.422 |
| HSC weight | 4.860/5.379/6.667 | 4.866/5.384/6.655 |
| multiplicative m | -0.203/-0.102/-0.066 | -0.199/-0.101/-0.065 |
| mag_a10 | 23.048/24.207/24.947 | 23.097/24.235/24.957 |
| resolution | 0.429/0.622/0.795 | 0.428/0.620/0.794 |
| PSF moment trace | 0.084/0.113/0.151 | 0.083/0.112/0.156 |

Radial/tomographic support is pair-counted, with the number of apertures having at least one retained source in each cell:

| Role | R proper (Mpc) | Pair rows | z-bin pair rows (1/2/3/4) | Centers with support (1/2/3/4) |
|---|---:|---:|---:|---:|
| lens | 0.1–0.2 | 1,224 | 546/310/225/143 | 65/63/59/52 |
| lens | 0.2–0.4 | 4,281 | 1,603/1,205/944/529 | 70/70/70/68 |
| lens | 0.4–0.7 | 10,011 | 3,186/2,947/2,454/1,424 | 74/73/72/71 |
| lens | 0.7–1.0 | 14,652 | 4,036/4,454/4,002/2,160 | 77/77/76/76 |
| lens | 1.0–1.5 | 33,922 | 8,572/10,766/9,321/5,263 | 78/78/78/78 |
| lens | 1.5–2.0 | 46,088 | 11,153/14,869/12,721/7,345 | 78/78/78/78 |
| lens | 2.0–3.0 | 129,387 | 31,465/41,580/36,238/20,104 | 78/78/78/78 |
| control | 0.1–0.2 | 4,088 | 895/1,331/1,192/670 | 285/303/293/262 |
| control | 0.2–0.4 | 15,671 | 3,660/5,313/4,358/2,340 | 332/335/335/330 |
| control | 0.4–0.7 | 42,666 | 9,695/14,295/12,051/6,625 | 347/348/347/348 |
| control | 0.7–1.0 | 64,917 | 14,848/21,746/18,415/9,908 | 353/353/354/353 |
| control | 1.0–1.5 | 156,772 | 36,335/51,686/44,733/24,018 | 360/360/361/361 |
| control | 1.5–2.0 | 214,064 | 49,747/70,863/60,842/32,612 | 366/365/366/365 |
| control | 2.0–3.0 | 607,601 | 139,244/200,856/173,172/94,329 | 376/376/376/376 |

Per-aperture retained counts and effective-support measures are in `per_center.tsv`; paired five-control comparisons are in `per_lens_five_control_support.tsv`; 0.1–3 Mpc annular counts and weights are in `radial_support.tsv`. The 13,440-row `per_center_radial_tomographic_support.tsv` and its 56-row `radial_tomographic_support.tsv` aggregate the same support by radial annulus and HSC Y3 source bin.

## Interpretation limits

`N_eff_shape_weight` uses only the HSC catalog weight. The reported `shape_weight*(D_l beta)^2` is a geometry-weighted support proxy from fiducial HSC Y3 n(z); it omits per-object shear calibration, responsivity, HSC selection-bias correction, and n(z) uncertainty. These summaries do not constitute a shear profile or measurement. Repeated source IDs across apertures are retained pair rows and are not independent samples. PSF spatial matching, control contamination, source-reuse/field covariance, and injection/recovery remain unvalidated. `science_measurement=false`.
