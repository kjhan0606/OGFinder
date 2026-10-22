# Alternative HSC-PDR3 cluster selection for ICL expansion — 2026-10-22

## Finding

There is a plausible complementary candidate list for an **ICL-only** expansion:
Vicentin et al. (2025), *Selecting Clusters and Protoclusters via Stellar
Mass Density II*, identify 16,007 candidate structures over approximately
850 square degrees of the HSC-SSP Wide PDR3 footprint, across
`0.1 < z_phot < 2`. Their catalog is selected through massive candidate
central galaxies and includes probabilistic membership and richness estimates.
The authors report a simulation-predicted structure purity of about 90%, but
only about 65% correct dominant-galaxy assignment; they also identify groups
near `10^13.5 Msun` as a main contaminant. They report up to 40% unique
candidates relative to other methods at `0.1 < z_phot < 0.45`.

This is not an independent imaging footprint: both this catalog and the HSC
PDR3 CAMIRA catalog are based on HSC Wide PDR3/S20A imaging. It may identify
additional centers or fill cells absent from the current CAMIRA target list,
but it cannot add sky outside HSC PDR3. Candidate count alone therefore does
not resolve the current shortage of independent spatial blocks. The catalog
also mixes virialized clusters and protoclusters up to `z=2`; high-redshift
objects cannot be silently pooled with the current low-redshift ICL sample.

## Fit to the two research branches

- **ICL and satellite association:** worth an offline overlap audit if the
  machine-readable table can be obtained. Keep the selection as its own
  stratum, report CAMIRA overlap and unique 5-degree cells, and do not combine
  richness or satellite counts across finders without calibration. Freeze a
  cluster-only redshift interval and a centering/membership rule before
  examining measured ICL. The paper's structure-purity estimate is not a
  substitute for image, mask, depth, PSF, sky, or matched-control support.
- **Background-galaxy weak lensing:** this catalog by itself adds no shape
  measurements. HSC's PDR3 release page says a PDR3 weak-lensing release was
  not yet available there; the existing S19A shape sample remains a separate
  footprint/data product. Do not count these candidates toward the WL lens
  sample until shape-catalog overlap and random-lens support are established.
- **Membership data:** the paper describes probabilistic member selection,
  but the currently inspected manuscript evidence does not establish that a
  member-by-member table is included with the published cluster table. If only
  cluster rows are available, a consistent satellite catalog requires a
  separate, approved HSC source-catalog workflow.

## Decision and next gate

Treat this as a promising **candidate-source audit**, not as an approved sample
or replacement for CAMIRA. First identify the publisher's exact machine-
readable table artifact, column definitions, license, byte size, and checksum
without retrieving it. Then prepare a small, reviewable download request for
that single cluster table only. After user approval, Astrafek should compare
its selected low-redshift candidates against S20A CAMIRA, the 800-cluster
parent, and the existing HSC image manifest; report unique candidates and
new spatial cells before any image or object-row request. Any latter request
requires separate approval.

No catalog, image, shape, photo-z, or member data was downloaded or queried
for this assessment.

## Sources

- Vicentin et al. (2025), author manuscript / primary paper:
  <https://arxiv.org/abs/2510.10736>
- Official HSC PDR3 CAMIRA documentation (catalog constructed from PDR3/S20A):
  <https://hsc-release.mtk.nao.ac.jp/doc/index.php/camira-cluster-catalog__pdr3/>
- Official HSC PDR3 release page (weak-lensing release status):
  <https://hsc-release.mtk.nao.ac.jp/doc/index.php/sample-page/pdr3/>
