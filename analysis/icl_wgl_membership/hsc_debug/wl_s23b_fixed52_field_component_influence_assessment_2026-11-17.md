# Fixed-52 WL field influence: lens and random components

Using only the saved per-lens calibrated profiles and existing official HSC
footprint polygons, the 52 lens centers were assigned to AUTUMN (21), HECTOMAP
(3), and SPRING (28). The pinned input and geometry hashes, seven bins per lens,
field counts, and complete-sample profile were checked. For each lens-center
field and radial bin, omitting the lens sets in that field changes the mean
random-subtracted profile by `(change in mean lens term) - (change in mean
equal-center random term)`. This is an arithmetic decomposition, not an
identification of causal systematics.

The result was independently compared to the saved exact field deletion, which
removes a pair if either endpoint is in the field and recalculates support.
The maximum difference over 42 field/bin/component rows is 0.0625 M_sun/Mpc²,
negligible relative to the 10^13 M_sun/Mpc² scale shown below. Their close
numerical equality is a check on these aggregate profiles; it does not prove
that individual support membership is identical. No image or source-row
reading, network request, or new transfer was used.

## Quantitative result

Units throughout the table are 10^13 M_sun/Mpc²; changes are relative to the
full 52-lens equal-center profile. `− random` is the negative change of the
mean matched-random term, so the two contributions sum to the total shift.

| Component and removed field | Lens contribution | − random contribution | Total shift |
|---|---:|---:|---:|
| Tangential, 1.5–2.0 Mpc, AUTUMN | +0.6603 | +0.1093 | +0.7696 |
| Tangential, 1.5–2.0 Mpc, HECTOMAP | −0.0617 | −0.0425 | −0.1042 |
| Tangential, 1.5–2.0 Mpc, SPRING | −0.7269 | −0.0544 | −0.7813 |
| Cross, 0.1–0.2 Mpc, AUTUMN | +1.6450 | +1.7742 | +3.4192 |
| Cross, 0.1–0.2 Mpc, HECTOMAP | +2.4203 | +0.6310 | +3.0513 |
| Cross, 0.1–0.2 Mpc, SPRING | −7.0662 | −3.5800 | −10.6462 |

In this selected sample, most of the two large-field shifts in the outer
bin comes from the lens term; the inner cross shift also has a material random
component. That does not prove that the former is astrophysical or that the
latter is caused by the control selection. The lens and random terms retain
shape noise, calibration, source reuse and observational systematics. With only
three occupied fields, their changes do not furnish independent uncertainties
or statistical significance; no mass or detection is inferred. Geometry labels
are footprint labels, not a solution to spatial covariance, photo-z/shape
calibration, PSF/mask matching or boost uncertainties. Science remains false.

- Runner: `analyze_s23b_wl_field_component_influence.py`
- Tests: `test_s23b_wl_field_component_influence.py`
- Products: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed52_field_component_influence_2026-11-17_01/`
- Summary SHA-256: `4fbb668d73ad6866d24f0e438a8d356f198f67cebf874239176b6179b444a026`
