# Full-sample local-core/public-wing overlap sensitivity — 2026-10-21

## Scope and decision

Astrafek compared the approved HSC PDR3 `HSC-I` symmetric extended-PSF reference
with the position-specific coadd model-PSF cores at all 800 cluster centers and
798 matched-control centers. It retained signed PSF pixels and fit the local
profile amplitude against the public profile in four radial windows between
1.5 and 3.3 arcsec. The two overlapping fit windows provide a radial-shape
sensitivity check; they are not independent measurements. The result is a
core/reference overlap diagnostic. It constructs and applies no hybrid kernel,
does not read coadd science pixels, and does not correct the ICL stack.

The target filter came from the validated `cluster_profiles.tsv`; each control
filter came from that control's actual `coadd_header_filter` in the approved
control manifest. The output contains all 1,598 center fits. Target/control
comparisons include only the 783 exact same-filter pairs: 302 `i` and 481
`i2`. Fifteen pairs have different actual filters and were excluded from paired
statistics; two targets have no matched control. The first exploratory `_01`
output inherited target filters for controls and is superseded; use only the
corrected `_02` output.

## Results

The scale is the weighted least-squares factor `local profile / public profile`
over the named overlap window. Across the full 1.5–3.3 arcsec window, median
scale was 0.713 for `i` targets and 0.740 for `i` controls; for `i2` it was
0.787 and 0.797. The corresponding position-to-position q16–q84 ranges were
0.535–0.919, 0.557–0.950, 0.633–0.938, and 0.622–0.935. These are
distributions across positions, not confidence intervals.

The scale changes with radius. Median scales in the inner 1.5–2.4 arcsec and
outer 2.4–3.3 arcsec windows were 0.731 and 0.520 for `i` targets, 0.754 and
0.560 for `i` controls, 0.803 and 0.555 for `i2` targets, and 0.817 and 0.572
for `i2` controls. Full-window weighted shape-residual medians were 12.4%,
11.7%, 12.8%, and 12.1%, respectively. The radial change and residuals show
that one amplitude does not make the survey-level PSF profile match the local
core over this interval.

Within exact same-filter pairs, the median target-minus-control full-window
scale difference was −0.0239 for `i` (302 pairs; q16–q84 −0.1294 to 0.0644)
and −0.0033 for `i2` (481 pairs; −0.1181 to 0.0988). Median absolute paired
scale differences were 0.0625 and 0.0764. These are descriptive; no hypothesis
test or uncertainty model was applied. Similar medians do not validate outer
wings or prove matched sky transfer.

## Interpretation and next gate

The public PSF is useful as an extended-wing sensitivity template, but the
position dependence, radial scale drift, and 11.7–12.8% overlap shape residuals
do not justify a single global splice amplitude. The local model stamps end
near 3.6 arcsec, far inside the physical 200–300 kpc ICL annulus for this
sample. Existing outer-star support is sparse and does not bridge that gap.

The next gate is injection/recovery in the existing target/control image masks
over a predeclared grid of splice radii and amplitudes, with point-source and
diffuse-light profiles evaluated separately, and with sky-transfer sensitivity
kept independent. A candidate hybrid remains a sensitivity scenario only until
those recovery checks and local outer-wing evidence pass. Current flags remain
`psf_approved=false` and `science_measurement=false`.

## Reproduction and provenance

Analysis:
[`analyze_public_psf_target_control_overlap.py`](analyze_public_psf_target_control_overlap.py)
and synthetic regression tests. Independent FITS-level recalculation:
[`validate_public_psf_target_control_overlap.py`](validate_public_psf_target_control_overlap.py).
Seven focused tests passed; both scripts compile, and `git diff --check` passed.

Authoritative outputs are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/extended_psf_target_control_overlap_2026-10-21_02/`.
The summary SHA-256 is
`8ba338cf139ab469d1c21a10f4401ea3d6a86e3dd402ea07bc49ca512be2e6d8`; the
1,598-row position table SHA-256 is
`ffe0a557d47694f73156b8f084d22f5977d658bb1b6932fe1cb5400dba22c62a`; the
783-row exact-filter pair table SHA-256 is
`b6276e52678abc91499ce56a65acbaec0ceaa055599aa56fbecb2a6be1686556`; and the
independent validation JSON SHA-256 is
`add0ce56323fe32564147f6ab9ec91d05d8ee3d2810f62c69494ac1b38b607ec`.

Analysis-script SHA-256:
`374eade29c92022d75cf5ee415cf9651cbdae18e4e848bcf2e0c66a30be249e2`;
validator SHA-256:
`0ac5bcf29b733cb157a7f7b27bebbbdc00286d9ea5ca75181c851720e9c38903`.
