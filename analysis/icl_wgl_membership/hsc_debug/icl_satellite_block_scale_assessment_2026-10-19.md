# ICL–satellite spatial-block scale sensitivity — 2026-10-19

## Scope and reproducibility

The exploratory ICL–satellite association used connected target/control sky
blocks formed from 5-degree RA/Dec cells. This audit reruns the same selected
cluster/filter pairs with 2.5-, 5-, and 10-degree cells to measure how the
bootstrap interval depends on that spatial-dependence assumption. It does not
select a block width from the observed signal, calculate a new significance
test, query HSC, or download data.

Astrafek's [`audit_satellite_icl_block_scale_sensitivity.py`](audit_satellite_icl_block_scale_sensitivity.py)
uses the checksum-verified 74-row table from
`satellite_icl_geometric_strata_2026-10-18_01`, with 5,000 connected-block
bootstrap draws per supported filter/stratum and a fixed seed of 20261020.
The 5-degree run was compared against the archived calculation across 60
numeric results and matches to 1e-12. The authoritative output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_block_scale_sensitivity_2026-10-19_04/`.
Its input table SHA-256 is
`a0c71945f664688ee1f77eaf07ec91ceeb44edc49f0cf5ecc2c3e258ca3af471`, result
TSV SHA-256 is
`49bc38313e16d4e2c3000fdbd963870439859940ed7c5696b348b22dce2ea2ce`, and
summary SHA-256 is
`5652fb9276de901a1c0e6f322bbc669927b67f632da2a155a49fbc824bf1c4f0`.

## i2 results

| Stratum | Cell width | Pairs | Connected blocks | Partial r | Spatial-block 95% interval |
|---|---:|---:|---:|---:|---:|
| All supported pairs | 2.5° | 53 | 36 | 0.234 | [−0.013, 0.391] |
| All supported pairs | 5° | 53 | 28 | 0.234 | [0.014, 0.418] |
| All supported pairs | 10° | 53 | 20 | 0.234 | [0.016, 0.396] |
| Geometrically non-overlapping | 2.5° | 33 | 25 | 0.043 | [−0.274, 0.328] |
| Geometrically non-overlapping | 5° | 33 | 22 | 0.043 | [−0.318, 0.345] |
| Geometrically non-overlapping | 10° | 33 | 17 | 0.043 | [−0.330, 0.409] |
| Overlap possible | 2.5° | 45 | 35 | 0.280 | [0.005, 0.482] |
| Overlap possible | 5° | 45 | 28 | 0.280 | [0.042, 0.538] |
| Overlap possible | 10° | 45 | 20 | 0.280 | [0.030, 0.481] |

The point estimates do not change because the selected observations and
regression are fixed. The full-sample i2 interval crosses zero at 2.5 degrees
but not at 5 or 10 degrees, so whether that interval excludes zero depends on
the assumed spatial correlation scale. The geometrically non-overlapping
stratum includes zero at all three widths. The paired overlap-minus-nonoverlap
correlation difference is `0.143`; its 95% interval includes zero at every
width (`[−0.450, 0.661]`, `[−0.419, 0.759]`, and `[−0.574, 0.754]`).

This is not a detection. The cell widths are post hoc sensitivity choices, the
association sample was selected exploratorily, and no multiplicity correction
is applied here. The separate familywise block-wild analysis found no adjusted
significant result. These intervals also omit pixel/coadd covariance, PSF
wings, sky transfer, extinction, and membership-completeness uncertainties.
`science_measurement=false` remains in the machine-readable output.

## Pipeline and run record

The five-degree result reproduces the original seeded bootstrap exactly, which
checks the generalized connected-block code against the prior implementation.
Five local tests pass, including cell assignment, deterministic bootstrap,
invalid-coordinate rejection, output schema, and 5-degree reproduction.

Three preliminary output directories `_01` through `_03` were preserved after
successive report-writing/provenance fixes. `_04` is the complete authoritative
run; all output paths refuse overwrite. No source data or prior analysis output
was modified.
