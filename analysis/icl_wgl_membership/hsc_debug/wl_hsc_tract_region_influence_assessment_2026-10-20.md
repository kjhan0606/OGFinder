# HSC tract-region influence check — 2026-10-20

## Scope

This is a field-scale influence diagnostic for the existing 800-lens,
cluster-only HSC Y3 weak-lensing stack. It uses the approved shape rows,
CAMIRA lens manifest, and the HSC patch-QA tract centers already present in
the project data. It makes no HSC request and adds no random-lens, shape, or
photo-z data. `science_measurement=false` remains in force.

The Astrafek runner groups the 840 unique HSC tracts by connected components
of their spherical-mean patch centers, with a 2.0-degree adjacency threshold.
Each lens and source is assigned to its nearest tract center. Recomputing the
stack after deleting each occupied component removes every row whose lens or
source belongs to that component. Threshold checks from 1.6 to 3.0 degrees
give the same four tract components and the same three components occupied by
the lens sample, with 499, 246, and 55 lenses. The grouping is an approximate
footprint graph; it is not a polygon intersection or an official HSC field
label.

## Results

The full-stack counts and all seven tangential/cross profile values match the
existing 52-cell spatial-jackknife run exactly. As before, 3,069,544 input
rows yield 2,153,004 selected pairs, 609 lenses, and 2,051,991 unique sources.
The field-scale leave-one-region results are:

| Deleted region | Pair-table lenses | Deleted pair rows | Used pairs left | Tangential `DeltaSigma` change by radial bin |
|---|---:|---:|---:|---|
| 0 (499 lenses) | 388 | 1,907,737 | 844,767 | −15.8%, −7.0%, −16.4%, −3.7%, −2.4%, −6.9%, −8.0% |
| 1 (246 lenses) | 190 | 918,016 | 1,480,526 | +5.5%, +5.7%, +4.2%, +0.8%, +0.1%, +0.9%, +2.6% |
| 2 (55 lenses) | 54 | 243,791 | 1,980,715 | +2.3%, −1.4%, +3.8%, +0.9%, +0.9%, +2.2%, +1.4% |

Radial bins are 0.1–0.2, 0.2–0.4, 0.4–0.7, 0.7–1.0, 1.0–1.5, 1.5–2.0,
and 2.0–3.0 Mpc. Deleting the largest connected component moves the inner
profile by up to 16.4%. Its innermost cross-component value shifts upward by
`14.3 × 10^12 M_sun Mpc^-2`, from a full-stack value of about
`2.95 × 10^12 M_sun Mpc^-2`. These are material footprint-influence signals
for this sample.

## Interpretation and limits

The 52-cell jackknife remains useful for within-footprint spatial variation,
but it does not represent independent survey regions. The field-scale check
shows the profile depends measurably on the largest tract-connected region,
especially inside 0.7 Mpc. Three occupied regions are too few for a stable
full covariance, so these leave-one-out shifts are not uncertainties,
confidence intervals, or a detection test. Nearest-center source assignment
and approximate tract adjacency also leave boundary uncertainty unmodeled.

The next lensing gate remains the pending approved random-coordinate CAS
export. After that export completes, run the full geometry and footprint
screen, select matched random centers, and obtain the separately approved
random-source shape rows before forming a random-subtracted stack with paired
spatial covariance. The present field-influence diagnostic does not replace
those steps.

## Reproduction

- Runner: [`run_wl_hsc_tract_region_influence.py`](run_wl_hsc_tract_region_influence.py)
- Slurm wrapper: [`run_wl_hsc_tract_region_influence.sbatch`](run_wl_hsc_tract_region_influence.sbatch)
- Focused tests: [`test_run_wl_hsc_tract_region_influence.py`](test_run_wl_hsc_tract_region_influence.py), 4 passed
- Slurm job: 428935, completed, exit 0, 1 min 28 sec, peak RSS 1.82 GiB
- Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_hsc_tract_region_influence_2026-10-20_01/summary.json`
- Output SHA-256: `71d2ad8df5df023770105bdbedcccde4f04da2e6ea87637e3ecff2367bb312d4`
- Runner SHA-256: `9d204f7d62e032aa59eb484397f002026ef0d44015a4068120ca3f30d5c1ddc9`
- Inputs: 32 shape-row batch hashes, plan SHA-256 `3fc2f11105f956c5f4c940225c7b268bb70b22bcf17aaf0c62f3306a6a527557`, SACC SHA-256 `a28f9e2e088e92d96d2d87083a6eeeac4c0c84c6b48e033bd3d1dc1bf0d8958f`, patch-QA SHA-256 `5c982c954b1a4026e5fe4ab75e978d4734b67d41e04aeddf8725eb5351ea3c64`.
