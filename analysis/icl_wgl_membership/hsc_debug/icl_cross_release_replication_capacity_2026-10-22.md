# S20A–S23B candidate pools: independent-footprint ceiling — 2026-10-22

## Result

The current exploratory i2 ICL–satellite association occupies 33 target 5-degree cells. The 1,333 selection-compatible S23B additions occupy 100 cells, of which 67 are new against those 33 discovery cells. The already-approved but not submitted S20A broad-redshift coverage plan has 2,596 cluster-only candidates in 84 cells, 51 of them new against discovery. Its broader `0.20 ≤ z_cl < 0.70` selection does not match the current BCG spectroscopic selection.

Taking the union of those two pools gives 102 candidate cells total, only **69 new target cells** beyond discovery. The S20A pool contributes just two new cells beyond S23B. Under the single frozen-endpoint Fisher-z approximation, r=0.293 has 80% power at 92 blocks; the 69-cell upper bound has power 0.675 and remains 23 short. At r=0.20, the 80% target is 196 blocks and power at 69 cells is 0.368. For the exploratory family of 30 comparisons, r=0.293 requires 180 blocks; even adding the 33 discovery cells to the 69 new-cell ceiling leaves an optimistic combined count of 102, 78 short.

The 69 cells are not usable or independent connected blocks. Candidate image, mask, depth, PSF, satellite completeness and matched-control gates can reduce the number, and controls may merge sky components. The catalog union also changes the redshift/spectroscopic selection, so it cannot be treated as one homogeneous sample without a new selection and calibration design. The previously approved S20A metadata query can improve coverage knowledge but will not materially solve the independent-footprint shortage.

## Reproducibility

- Astrafek tool: [`assess_icl_cross_release_replication_capacity.py`](assess_icl_cross_release_replication_capacity.py)
- Unit tests: [`test_assess_icl_cross_release_replication_capacity.py`](test_assess_icl_cross_release_replication_capacity.py), 2 passed
- Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_cross_release_replication_capacity_2026-10-22_01/`
- Summary SHA-256: `ef410a4c62e6e3b6bd8131c85faf5edab531c1a14cf867f4a7d70e45ee7a98f2`
- Assessment SHA-256: `e7b5c91f8c6b224644426052bdd8ed150b757ed048b3ec2ff1c546e4d5aa654c`

The calculation reads only hash-pinned local candidate coordinates, parent manifest and association rows. It submits no HSC query, retrieves no files, reads no image pixels and remains `science_measurement=false`.
