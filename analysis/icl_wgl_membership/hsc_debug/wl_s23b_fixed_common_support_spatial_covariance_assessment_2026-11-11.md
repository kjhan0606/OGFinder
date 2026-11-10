# S23B fixed 52-lens matched-random spatial covariance (2026-11-11)

## Scope and estimator

This Astrafek CPU analysis uses only the already approved S23B HSC S19A shape-row
export, the pinned HSC Y3 tomographic `n(z)`, and existing Astrafek support and
profile products. It made no CAS query, network request, download, image read,
or halo-mass fit. `science_measurement=false` remains in force.

The fixed sample contains 52 matched lens sets. Each must have positive lens
support and at least four positive-support matched random centers in every one
of seven annuli from 0.1 to 3 proper Mpc, using source bins 2–4. For each lens
set the calibrated lens stack is computed separately, the five matched random
centers are pooled by source-pair lensing weights, and the random-subtracted
profile is then averaged equally over lens sets.

Each delete-one replicate removes a lens/random-center–source pair if either
endpoint falls in the deleted equal-area RA/sin(Dec) cell. Lens and random
profiles use the same deleted cell. The full S19A calibration is recomputed for
each affected set. A lens set contributes to a replicate only if both profiles
are present in all seven annuli and at least four random centers remain in each
annulus; the same set intersection is used for all profile components. The
grid widths are 2.5° and 5°.

The full-data point estimate reproduces the pinned fixed-profile result, with a
maximum absolute difference of 0.03125 Msun Mpc⁻² across the seven bins. Across
795,289 selected target/control pair rows there are 493,402 unique source IDs;
2,099 occur in both target and control rows. Repeated source IDs have no
coordinate or tomographic-bin mismatch.

## Profile and delete-one errors

Values and errors are in `10^13 Msun Mpc^-2`. Errors are the diagonal standard
deviations from the paired endpoint spatial covariance. The full point profile
is common to both block widths.

| R proper (Mpc) | ΔΣt | Error 2.5° | Error 5° | Cross ΔΣ | Cross error 2.5° | Cross error 5° |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1–0.2 | 23.5538 | 5.2856 | 4.4101 | 7.1589 | 6.5972 | 5.5325 |
| 0.2–0.4 | 6.0478 | 2.3783 | 2.0390 | −0.6018 | 2.9736 | 2.7464 |
| 0.4–0.7 | 4.0030 | 1.4803 | 1.3399 | −0.4770 | 1.5286 | 1.6330 |
| 0.7–1.0 | 3.5211 | 1.5141 | 1.2548 | −1.4481 | 1.2778 | 1.1719 |
| 1.0–1.5 | 2.5206 | 1.0271 | 0.9687 | 0.5832 | 0.8681 | 0.7651 |
| 1.5–2.0 | 0.6750 | 0.6897 | 0.5796 | 0.8598 | 0.7573 | 0.6633 |
| 2.0–3.0 | 1.2551 | 0.5563 | 0.5352 | 0.1079 | 0.5900 | 0.6013 |

| Cell width | Occupied cells | Covariance rank | Minimum common lens sets after deletion | PSD | Exploratory numerical gate |
|---:|---:|---:|---:|:---:|:---:|
| 2.5° | 64 | 14/14 | 49/52 | Yes | Pass |
| 5° | 28 | 14/14 | 47/52 | Yes | Pass |

The gate checks finite delete-one replicates, at least 15 occupied cells, a
full-rank positive-semidefinite 14×14 matrix, and at least 80% of the fixed
lens sets retained by every deletion. Passing is a numerical stability check
for this approximate grid, not validation of the HSC survey covariance. No
profile significance or mass is inferred. Errors change with block width,
particularly in the inner annuli.

## Remaining limitations

- Equal-area cells are not validated contiguous HSC regions and can combine
  disconnected survey footprints.
- Fiducial `n(z)` uncertainty and photo-z dilution calibration are not
  propagated; reduced shear and magnification are not modeled.
- Random apertures are not yet demonstrated to match target depth, mask, PSF,
  sky, boost, and cluster-member contamination.
- Selection-response uncertainty and injection/recovery remain to be assessed.

Consequently, this is an exploratory covariance diagnostic, not a publishable
weak-lensing measurement. `science_measurement=false` remains in force.

## Provenance

- Runner: [`analyze_s23b_wl_fixed_common_support_spatial_covariance.py`](analyze_s23b_wl_fixed_common_support_spatial_covariance.py)
- Runner SHA-256: `662b17470b342f0f2486a2fd23472b7105ba6a76e5354fa37b1aa679a0010256`
- Approved row-plan SHA-256: `7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`
- Fixed lens-set membership SHA-256: `18312f27139b3750abdbf71543f22d9b7b18b2d3b24cbd29ed686a377524bb7d`
- Radial-support summary / per-lens table SHA-256: `c51a4d9268935a0919d0abd3da2f10b6cd16ee7e35f1e6e23a8f3fe8c9e51234` / `6744ece545781faffebef65fc7a0143284da3a234425b6bed540bd1e9d19625c`
- Fixed-profile / QA summary SHA-256: `4d270b81737e5bcb6188941f9f26cd6255614fd7ae8d1e6f2977cf51a4686930` / `eb73aa0365813d880d3edf60b12e4039f2575f8a97eab8868fc918121a34d445`
- HSC Y3 `n(z)` SHA-256: `a28f9e2e088e92d96d2d87083a6eeeac4c0c84c6b48e033bd3d1dc1bf0d8958f`
- Final results (outside Git): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed_common_support_spatial_covariance_2026-11-11_02/`
- Final summary SHA-256: `acdd4da69e4d7fe3dbb52d8584b3996bd836a41196da59aecfe3fcdf3e1fc6e9`
- Assessment / deletion-support / replicate / covariance-matrix SHA-256:
  `d82ce35c95cb7ba1184f185ceb5b1a23aa751397cb94f3718ed4a1abc1184a05` /
  `810e80bdf9cd28bad42c456968e885ab61d417f68a5c823c5fad5ca98ba6b790` /
  `158bbaf5518b91589c6477673b259a79de9deeab997e8069abd0995bd46a8f27` /
  `bc718cb4c5734e7a981048b735ad6239486aeb99cbc6ee26155228906ceed9b4`
- A first `_01` output was superseded after its block-deletion summary retained
  only the 5° rows. `_02` preserves all 64 2.5° and 28 5° deletion rows; the
  covariance results are numerically unchanged.
