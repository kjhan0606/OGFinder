# S23B CAMIRA catalog audit — 2026-10-20

## Scope

The user approved one S23B CAMIRA Wide star-masked table for duplicate and
coverage feasibility checks. Astrafek fetched only
`camira_s23b_wide_sm.tbl` from `oguri/cluster_catalogs` commit
`d80734b10bd823d13b9cd5292dadc77c9c05f683`. It is 827,339 bytes; the verified
Git blob SHA-1 is
`066c2c8beca8f3a65186db7511fddca52d298908` and SHA-256 is
`5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71`.
The raw table and acquisition manifest are in
`/gpfs/kjhan/Astrafek/data/catalogue_s23b_camira_2026-10-20_01/`, outside Git.
No image, shape, random-source, or photo-z asset was requested or received.

The seven whitespace-separated columns match the file's header: cluster ID,
name, RA, Dec, photometric cluster redshift, richness, and BCG spectroscopic
redshift (`-1.0` denotes unavailable). The table has 12,162 unique IDs.

## Results

Nearest-neighbor matching on the unit sphere avoids assuming S20A and S23B IDs
are interchangeable:

| Audit | Result |
|---|---:|
| S23B rows within 1 arcsec of any S20A CAMIRA position | 4,492 |
| S23B rows within 1 arcmin of any S20A CAMIRA position | 5,192 |
| S23B rows within 1 arcsec of an existing 800-target position | 607 |
| Same, within 1 arcmin | 661 |
| Those close matches satisfying the current redshift/BCG cut | 559 |
| Eligible after the cut and 1-arcsec S20A deduplication | 1,440 |
| Greedy candidates after 6-arcmin separation from current and new targets | 1,333 |
| New occupied 5-degree RA/Dec cells / total after adding candidates | 22 / 100 |

The selection cut was `0.30 <= z_cl < 0.55`, positive BCG spectroscopy, and
`|z_BCG-z_cl|/(1+z_cl) < 0.01`. The greedy priority is descending richness,
then catalog ID. It is a transparent pool estimate, not a globally optimal
selection or a science sample. The 1-arcsec S20A match is treated as a
high-confidence duplicate; the 1-arcmin count is a wider positional-neighbor
sensitivity, not a duplicate count.

The existing PDR3 target-image manifest is complete for 800/800 targets. The
607 close S23B matches identify catalog positions already represented among
those target images; 559 also pass the redshift/BCG cut. The images were not
reprocessed in this audit. The 1,333 newly proposed catalog positions lie
outside the existing target list, so this audit does not establish whether
PDR3 coadds, masks, PSFs, depths, full apertures, or matched controls exist at
those locations. The S23B paper catalog is associated with final-year data / a
future PDR4 release, while the working images are PDR3. Image support and any
image transfer require a separate exact-scope review.

The 100 5-degree coordinate cells are only an upper bound on possible spatial
coverage. The current connected-block definition joins target and matched
control positions. Since new candidates have no verified controls, their
connected target/control block count cannot be computed from the catalog. The
new cells therefore do not imply 22 independent blocks.

## Reproduction and artifacts

The one-time acquisition was performed by the pinned script. To reproduce the
offline audit from that saved file on `lageunha`, use a fresh output path:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 analysis/icl_wgl_membership/hsc_debug/audit_s23b_camira_catalog.py \
  --s23b /gpfs/kjhan/Astrafek/data/catalogue_s23b_camira_2026-10-20_01/camira_s23b_wide_sm.tbl \
  --output /gpfs/kjhan/Astrafek/data/catalogue_s23b_camira_audit_NEW
```

The acquisition script is pinned to the exact authorized URL and refuses a
different filename, response size, hash, or existing output. The audit checks
the acquisition manifest, S20A parent table, 800-target manifest and complete
image-download manifest before writing any result. Its files are preserved in
`/gpfs/kjhan/Astrafek/data/catalogue_s23b_camira_audit_2026-10-20_02/`:

- `summary.json` (2,790 bytes), SHA-256
  `f98e9bdb0b7f5cfaf217132cafa6c8a0cf9187a0ec448cc9ee1c50826c800e20`
- `greedy_candidate_additions.tsv` (1,333 rows, 94,368 bytes), SHA-256
  `1b59ca9941849fa7970c43c6e0374db6587716fc4ae5d71b11f84a245ebe6b2c`

No HSC query was submitted. `science_measurement=false` remains in force.
