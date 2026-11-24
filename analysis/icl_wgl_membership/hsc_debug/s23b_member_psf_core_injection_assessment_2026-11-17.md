# Position-specific PSF-core member injection: exploratory assessment

The approved 427 HSC-I member-position model-PSF cores were applied to the
same seeded member halo in each of 62 i and 365 i2 frames, at three injected
amplitudes. The CPU run completed 427/427 frames, with zero failed frames and
1,281 member-amplitude rows. Its 1,219 supported PSF rows are **not** a
population completeness estimate. The comparison verifies exact matching of
member/amplitude keys and baseline replay against the pinned unconvolved run.
All listed output hashes were verified. `science_measurement=false`.

For the 1.0 nJy arcsec^-2 input amplitude at 50 kpc, 58/62 i and 351/365 i2
pairs passed the profile selection in both runs. On this conditional support,
the median of the paired PSF-core-minus-unconvolved **normalized response** was
-0.00146 (95% 5-degree-cell bootstrap interval [-0.00581, 0.00196], 20 cells)
for i and +0.00200 [0.00083, 0.00329] (61 cells) for i2. At the same amplitude,
the corresponding median paired **recovered surface-brightness increment**
was +0.00361 [0.00252, 0.00449] nJy arcsec^-2 for i and +0.00418
[0.00382, 0.00461] nJy arcsec^-2 for i2. The response ratios use different
nominal input denominators in the convolved and unconvolved cases, so the
physical increment is reported separately. These are exploratory marginal
intervals; they are not multiplicity-adjusted detections or population-wide
corrections.

At 1.0 amplitude, one i frame changed from supported to too few rotated
controls; in i2 one changed supported to below coverage and one changed in the
opposite direction. Selection is recomputed after injection. Results conditioned
on both being supported omit those changes and cannot by themselves establish
selection completeness. The circular exponential injected halo (15-kpc scale)
is one model, and 50–75 kpc is one chosen profile bin. The finite PSF cores do
not assess extended scattering wings. Source Poisson noise, effective coadd
weights, residual sky transfer, membership uncertainty and full covariance
remain unresolved; this is not a calibrated ICL–satellite science result.

- Input plan SHA-256: `03f974ad24ead32eabd103ef7a6e0d6b52317aa2647023a5f7402003a2e3b25d`
- Injection summary SHA-256: `7cdd34b10f3e976a1ecd94518c338e286aaa59fd4b42c310c3607da4c27431c4`
- Paired comparison summary SHA-256: `3ecb1ad6900ac6255b4143b461782a3978e020e3dbed9e22b9209e534bec19cf`
- Products: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_psf_injection_run_2026-11-17_01/`
- Comparison: `compare_s23b_member_psf_injection.py`; tests: `test_s23b_member_psf_comparison.py`.
- No additional data downloads or network queries were made in this analysis.
