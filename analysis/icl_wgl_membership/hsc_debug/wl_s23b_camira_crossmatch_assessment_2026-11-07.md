# Existing weak-lensing sample vs. S23B CAMIRA member-profile sample

## Question

Can the existing HSC S19A weak-lensing stack be joined by cluster identity to
the new S23B CAMIRA member-centered outer-profile sample, so that one can test
a same-cluster weak-lensing–ICL–satellite relation?

## Method and provenance

`audit_wl_s23b_camira_crossmatch.py` compares the 800 saved WL lens centers
against all 1,043 S23B CAMIRA target centers using great-circle angular
separation. It also records each WL lens's nearest S23B target and separately
flags numeric catalog-ID collisions. The comparison covers the 620 lenses
with selected WL random centers and the 435 lenses in the saved common
cluster-shape/random-shape support.

Before matching, the runner reconstructs those 435 common-support lens IDs
from the saved count-only results and verifies their sorted-ID SHA-256 against
the existing weak-lensing stack. It confirms the 427 clusters represented in
the saved member outer-profile relation (62 HSC-i and 365 HSC-i2 rows). Inputs
are checksum-pinned in the output summary.

No network request, CAS submission, source-row transfer, science-image pixel
read, or new download occurred. The approved S23B member file is included in
the provenance hashes but was not reprocessed by this positional audit.

## Results

- The 800-lens WL sample, 620 random-center subset, and 435 common-support
  subset have **zero S23B CAMIRA centers within 60 arcsec**. The result is
  also zero when restricted to the 427 S23B clusters with saved member outer
  profiles.
- Across all 800 WL centers, the nearest S23B target is at least
  **364.846 arcsec** away (6.08 arcmin); the median nearest separation is
  1,934.995 arcsec.
- Numeric equality occurs for 66 WL and S23B catalog IDs, including 29 IDs
  among the member-profile targets. None of these same-number pairs is within
  60 arcsec. These are catalog-namespace collisions, not valid cluster
  matches; the IDs must not be used to join the samples.
- The reconstructed 435-lens support hash matches the saved WL stack exactly:
  `82c0e737ebdbab8027f94267c3f7c5d77e6ed9f3535c6b20b8e8ddf31394eed1`.

The current WL and S23B member-profile products therefore do not form a
same-cluster sample. Their three-way weak-lensing–ICL–satellite association
cannot be measured by joining these products.

## Output and status

Machine-readable summary and per-lens nearest-target table:

`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_camira_overlap_audit_2026-11-07_01/`

- `summary.json` SHA-256: `7483c538a470fbbccd04be081f060356f058a25c70529ea2c9bd3960d3d3d47e`
- `wl_lenses_nearest_s23b_target.tsv` SHA-256: `29c7fc1453d810a53a1eedcac5c4c21d627d4acdf0dad8956494d8e5739f2b50`

Keep `science_measurement=false`. The existing WL stack remains exploratory;
this audit resolves sample identity only and does not change its measurement
status.

## Next scientifically valid route

A local count-only support plan is now prepared for the same 427 S23B
member-profile cluster centers. It contains 18 HSC PDR3 CAS batches, each
counting distinct S19A shape-catalog object IDs inside a 3 proper-Mpc cone.
The result scope is 427 `(S23B lens ID, integer count)` rows. The queries do
not request source positions, shapes, photo-z, images, PSF, or variance data.
The total SQL payload is 122,198 bytes. The plan has not been submitted and no
CAS result has been retrieved.

Plan directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_shape_count_plan_2026-11-07_01/`.
Query-plan SHA-256:
`b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984`.

The plan is wired to
`submit_s23b_wl_shape_counts.py`. Its `--plan` mode verifies provenance and
all 18 SQL hashes without network access. The submit mode requires both an
interactive `lageunha` terminal and the exact explicitly approved plan hash;
it checkpoints each job and validates the returned IDs and count-only schema.
The plan mode has been run successfully. No submit-mode invocation has been
made.

After the count result, separately validate shape support and select matched
random centers before designing any source-row export. Each HSC CAS submission
or source-row transfer is a separate data operation and requires explicit
approval.

The alternative is to establish the authoritative CAMIRA member positions
and catalog version for the existing 800-lens sample, then use the already
available WL and ICL products on that common sample. The embedded member
coordinates in the old sample manifest must have their catalog provenance,
selection, and completeness verified before they are treated as scientific
membership data.
