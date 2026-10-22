# HSC independent-footprint expansion gate — 2026-10-22

## Why the current sample must grow by usable sky, not catalog rows alone

The existing HSC PDR3 analysis starts from 800 CAMIRA clusters and has 33
independent i2 spatial blocks for the ICL–satellite endpoint. The current
sample-size sensitivity estimated 92 blocks for nominal 80% power at a
post-hoc effect of `r=0.293`, and 180 blocks under a 30-comparison Bonferroni
scenario. Those are planning values, not a confirmed effect or a final power
claim. The local S23B feasibility audit found up to 1,333 additional
catalog-level candidates after its stated cuts, but it did not verify their
PDR3 image, mask, depth, PSF, satellite, or shape-catalog support.

The public [HSC PDR3 CAMIRA documentation](https://hsc-release.mtk.nao.ac.jp/doc/index.php/camira-cluster-catalog__pdr3/)
identifies its cluster catalog as constructed from the S20A dataset and
provides Wide and Deep/UltraDeep samples. The [HSCLA2020 release page](https://hscla.mtk.nao.ac.jp/doc/)
reports substantially broader PI-program coverage, about 7,500 square degrees
in at least one filter, but the archive's [data-access documentation](https://hscla.mtk.nao.ac.jp/doc/data-access/)
warns that multi-filter area is sparse and uses a separate account from HSC-SSP.
The [HSCLA product inventory](https://hscla.mtk.nao.ac.jp/doc/available-data/)
also describes patchy field coverage. This broader archive may help construct
an ICL imaging sample, but it is not by itself evidence of uniform CAMIRA
membership, matched controls, or a compatible weak-lensing shape catalog.
HSC's [PDR3 release page](https://hsc-release.mtk.nao.ac.jp/doc/index.php/sample-page/pdr3/)
separates the Year-3 shape catalog from its imaging products; the official
[S19A shape-catalog page](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/)
must remain the authority for lens-source availability and calibration.

## Offline work before any new request

1. Preserve the current PDR3/S19A analysis as the reference sample and freeze
   one primary ICL–satellite endpoint before using another catalog. Treat
   `r=0.293` as a post-hoc design scenario only; also report smaller effects.
2. Build an Astrafek candidate ledger from already available CAMIRA tables,
   recording catalog release, selection, RA/Dec, redshift, richness, overlap
   with the 800-cluster set, and coarse sky block. Do not count a catalog row
   as usable data.
3. For every candidate, separately record whether there is a calibrated
   single-filter coadd suitable for ICL, a matched control field, required
   masks/variance/PSF support, and a validated weak-lensing source-shape
   catalog. Only candidates with the full support needed for a specific
   analysis contribute to that analysis's effective sample.
4. Compare the independent sky-block yield—not just cluster count—with the
   design target. If HSCLA contributes imaging without compatible shape rows,
   keep it as a separate ICL-only extension; do not pool it into weak lensing.
5. Prepare exact count-only availability plans first. Any new catalog or
   image/shape transfer remains outside this plan and requires the user's
   explicit approval for its precise scope.

## Current status and next gate

No HSCLA or additional HSC query, download, or account operation was performed
for this assessment. Existing PDR3 images, the 800-cluster parent, approved
shape rows, and the S23B catalog audit remain the only inputs. The next
Astrafek implementation should join candidate coordinates to locally
available image/shape footprint manifests and compute unique independent sky
blocks, without querying CAS. If the local products cannot establish coverage,
the next step is an exact count-only request plan for review before execution.

The current weak-lensing random-control coordinate export is a separate
pending prerequisite: CAS job `1062312` was last observed `running`; its
checkpoint contains two verified batches (33,514 of 780,191 planned
coordinate rows). Do not start a second query set while this job is unresolved.
All reported stacks remain exploratory; `science_measurement=false`.
