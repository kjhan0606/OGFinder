# S23B CAMIRA catalog expansion gate — 2026-10-20

## Why this catalog is being considered

The current ICL–satellite study uses the HSC PDR3/S20A CAMIRA pool. The
Astrafek design assessment found 33 current connected sky blocks and an
exploratory target of 180 blocks; even the broadest candidate expansion from
the local S20A pool was bounded by 84 target cells. A separate catalog release
is therefore worth screening before more image requests are designed. This is
a sample-expansion feasibility step, not evidence that S23B adds usable images
or independent blocks.

Oguri et al. (2026) publish the HSC-SSP final-year S23B CAMIRA catalogs. The
paper reports 12,162 clusters in the Wide catalog with the bright-star mask,
over 1,220 deg², for `0.1 <= z_cl <= 1.38` and richness `N_mem >= 15`. The
catalogs are distributed in the public
[`oguri/cluster_catalogs` repository](https://github.com/oguri/cluster_catalogs).

## Exact proposed transfer

The repository API directory listing was read as metadata only. It identifies
the single star-masked Wide catalog as:

- File: `camira_s23b_wide_sm.tbl`
- Repository commit: `d80734b10bd823d13b9cd5292dadc77c9c05f683`
- Pinned URL: <https://raw.githubusercontent.com/oguri/cluster_catalogs/d80734b10bd823d13b9cd5292dadc77c9c05f683/hsc_s23b_camira/camira_s23b_wide_sm.tbl>
- Size: 827,339 bytes
- Git blob SHA-1 from repository metadata: `066c2c8beca8f3a65186db7511fddca52d298908`
- Repository README SHA-1: `2213d144376bd00ef91483ee4ea72e63dfe3206b`

The README's displayed table uses a different filename spelling, so the API
directory listing above supplied the filename and object identifier. The user
approved this one-file catalog-only transfer. The file was downloaded on
`lageunha` to
`/gpfs/kjhan/Astrafek/data/catalogue_s23b_camira_2026-10-20_01/` and verified
at 827,339 bytes with the pinned Git blob SHA-1 above. Its SHA-256 is
`5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71`.
The adjacent manifest records the pinned URL, commit, hashes, timestamp, and
scope. No other CAMIRA table, HSC image, shape row, random-source row, or
photo-z distribution was fetched.

## Data-release compatibility check

The catalog paper identifies its parent as the HSC-SSP final-year S23B data
and says S23B will correspond to Public Data Release 4 (PDR4); it reports
roughly 1,200 deg² for the parent sample. The current project images and
shape-source data are from PDR3. HSC's official PDR3 availability page reports
over 600 deg² at nominal depth. Therefore the S23B catalog cannot be treated
as an automatic expansion of the current image sample: the approved catalog-only
screen can identify overlap with the existing PDR3 image inventory, but no
candidate outside that verified support should enter a stack. Any later image
request requires confirming the relevant release, product, and full-aperture
coverage as a separate gate.

The official HSC-SSP release schedule currently says PDR4 is the final release
but lists its date as TBD and says the timing is not yet known. Thus the
S23B/PDR4 paper catalog is only a candidate-discovery input for this project
until compatible PDR4 image and shape products are actually available; this
plan does not imply or request PDR4 imaging. The exact schedule page was
checked on 2026-10-20.

Sources: [Oguri et al., S23B CAMIRA catalog paper](https://arxiv.org/abs/2512.13954),
[official HSC PDR3 availability page](https://hsc-release.mtk.nao.ac.jp/doc/index.php/available-data__pdr3/),
and [official HSC-SSP release schedule](https://hsc.mtk.nao.ac.jp/ssp/data-release/).

## Catalog-only audit result

The table has seven whitespace-separated columns, as stated in its header:
ID, name, RA, Dec, cluster photometric redshift, richness, and BCG
spectroscopic redshift (`-1` means unavailable). It contains 12,162 unique IDs.
The audit uses unit-sphere nearest-neighbor separations, rather than assuming
that IDs are stable across releases.

- 4,492 S23B rows have an S20A CAMIRA position within 1 arcsec; 5,192 are
  within 1 arcmin. The 1-arcsec threshold is the high-confidence duplicate
  count; the wider radius is shown as a positional-neighbor sensitivity.
- 607 S23B rows match one of the current 800 target positions within 1 arcsec
  (661 within 1 arcmin). The existing target-image download manifest is
  complete for 800/800 targets. Of the 607 close matches, 559 also satisfy
  `0.30 <= z_cl < 0.55`, positive BCG spectroscopy, and
  `|z_BCG-z_cl|/(1+z_cl) < 0.01`.
- After the same redshift/BCG cut and 1-arcsec S20A duplicate removal, 1,440
  catalog rows remain. A richness-ordered greedy selection with at least
  6 arcmin from the existing 800 targets and previously added rows retains
  1,333 candidates. These occupy 22 new 5-degree RA/Dec cells, for 100 such
  cells including the current sample.
- The 1,333-candidate pool is a catalog-level upper bound. Current validated
  image manifests cover only the existing 800 targets. A catalog position
  outside those targets is not evidence either for or against usable PDR3
  image, mask, PSF, depth, or full-aperture support. No new candidate images
  or matched controls were checked or downloaded. The connected target/control
  block count cannot be computed before supported controls exist; 100 cells
  are not 100 independent blocks.

The one-time transfer is recorded in
[`acquire_s23b_camira_catalog.py`](acquire_s23b_camira_catalog.py). Reproduce
the offline audit from the saved table with
[`audit_s23b_camira_catalog.py`](audit_s23b_camira_catalog.py). The audit
summary and greedy candidate list are preserved outside Git at
`/gpfs/kjhan/Astrafek/data/catalogue_s23b_camira_audit_2026-10-20_02/`.
The scripts verify the exact acquisition manifest and refuse to overwrite
existing outputs.

All catalog-pool results must remain `science_measurement=false`. A positive
candidate count is not an ICL stack, a weak-lensing stack, or a measurement of
the ICL–satellite relation. The paper reports that the final-year catalogs are
larger than earlier CAMIRA releases; that does not establish new usable
coverage for this particular HSC PDR3 image set.

## References

- Oguri et al. (2026), PASJ 78, 416: <https://ui.adsabs.harvard.edu/abs/2026PASJ...78..416O>
- Catalog repository and S23B directory: <https://github.com/oguri/cluster_catalogs/tree/main/hsc_s23b_camira>
