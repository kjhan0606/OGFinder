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
- URL: <https://raw.githubusercontent.com/oguri/cluster_catalogs/main/hsc_s23b_camira/camira_s23b_wide_sm.tbl>
- Size: 827,339 bytes
- Git blob SHA-1 from repository metadata: `066c2c8beca8f3a65186db7511fddca52d298908`
- Repository README SHA-1: `2213d144376bd00ef91483ee4ea72e63dfe3206b`

The README's displayed table uses a different filename spelling, so the API
directory listing above is the filename and object identifier to verify if
transfer is approved. No catalog bytes have been fetched or saved. User
approval for this one file is pending. No other CAMIRA table, HSC image, shape
row, random-source row, or photo-z distribution is included in this proposed
transfer.

## Audit after approval

1. Retrieve only the one file above; verify exact byte count and Git blob
   identity before parsing.
2. Inspect its actual delimiter, header, units, null conventions, coordinate,
   redshift, richness, cluster-ID, and BCG fields. Preserve the raw table and
   its checksum outside Git; keep the Astrafek importer and audit script in
   this directory.
3. Cross-match against the existing S20A catalog and the current 800-lens
   sample using a spherical separation calculation. Report duplicate and new
   cluster counts before any selection.
4. Count candidates under the current fixed redshift, richness, BCG
   spectroscopic, and image-quality cuts; separately show candidates before
   those cuts so the loss is explicit.
5. Recompute target and matched-control 5-degree cells with the existing
   connected-pair-block definition. Report the resulting block count and a
   cell-level upper bound; do not equate a cell count with independent
   statistical information.
6. Cross-check candidate positions against the already-held HSC PDR3 image,
   mask, PSF, and depth inventory. The S23B catalog footprint and the S20A
   image footprint are not interchangeable. Present any additional image
   request as a separate, exact approval scope.

All catalog-pool results must remain `science_measurement=false`. A positive
candidate count is not an ICL stack, a weak-lensing stack, or a measurement of
the ICL–satellite relation. The paper reports that the final-year catalogs are
larger than earlier CAMIRA releases; that does not establish new usable
coverage for this particular HSC PDR3 image set.

## References

- Oguri et al. (2026), PASJ 78, 416: <https://ui.adsabs.harvard.edu/abs/2026PASJ...78..416O>
- Catalog repository and S23B directory: <https://github.com/oguri/cluster_catalogs/tree/main/hsc_s23b_camira>
