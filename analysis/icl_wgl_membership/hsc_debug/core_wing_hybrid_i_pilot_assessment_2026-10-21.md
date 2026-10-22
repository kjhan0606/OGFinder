# Local-core/public-wing hybrid PSF i-band pilot — 2026-10-21

## Run and validation

The predeclared 12-model grid completed on the real target/control image and
mask for HSC `i` cluster `1056`: splice radii 1.8/2.4/3.0 arcsec, a fixed
0.6-arcsec blend, and wing amplitudes 0.50/0.65/0.80/1.00. The frame has 10
matched satellite injection centers, close to the 13-frame sample median.
Diffuse ICL used the cluster-center local core joined to the public extended
HSC-I reference; each satellite envelope used its own position-specific local
core joined to the same public wing.

All 12 configurations completed without errors. The grid driver independently
validated each frame record, all local-core file hashes and hybrid metadata,
and recomputed transfer-response arithmetic with maximum error 0. Slurm job
432792 completed on `syntax`/`syn04` in 13m58s using one CPU, 32 GB requested
memory, and no GPU. Reported maximum RSS was 2,866,960 KiB.

## Sensitivity results

The entries below are min–max across the 12 PSF models for this one cluster;
they are not confidence intervals or population ranges. Responses are
measured surface-brightness changes per unit injected surface brightness.

| Injection → estimator | Redetection/mask rerun | Fixed baseline masks | Redetection − fixed |
|---|---:|---:|---:|
| Diffuse → ICL 200–300 kpc | 0.22111–0.22128 | 0.23060–0.23077 | −0.00949 (constant to shown precision) |
| Satellite envelopes → ICL 200–300 kpc | 0.36896–0.39892 | 0.50564–0.50727 | −0.13798 to −0.10740 |
| Diffuse → satellite 50–75 kpc | −0.06952 to −0.06948 | −0.01251 to −0.01247 | −0.05701 (constant to shown precision) |
| Satellite envelopes → satellite 50–75 kpc | 0.09297–0.16836 | 0.31063–0.32470 | −0.21851 to −0.15634 |

Across this frame and grid, diffuse-to-ICL recovery is nearly insensitive to
the chosen splice radius and wing amplitude. Satellite-envelope recovery in
the satellite outer-light estimator is much lower when source detection and
masks are rerun than when baseline masks are fixed. This identifies
reselection/mask response as a larger pilot effect than the tested PSF-model
grid for that estimator. It does not show whether that effect is shared by
the population: the pilot is one cluster and only tests synthetic additions.

The run does not include isolated point-source injection, source photon noise,
coadd pixel covariance, upstream sky-transfer recovery, or an i2-specific PSF.
No observed image was changed. This does not calibrate the PSF, establish
completeness, or measure an ICL–satellite association. All scientific and PSF
approval flags remain false.

## Reproduction

The frozen parameter grid and execution are in
[`run_core_wing_hybrid_grid.py`](run_core_wing_hybrid_grid.py),
[`run_core_wing_hybrid_i_pilot.sbatch`](run_core_wing_hybrid_i_pilot.sbatch),
and [`core_wing_hybrid_i_pilot_plan_2026-10-21.md`](core_wing_hybrid_i_pilot_plan_2026-10-21.md).
The independent grid assembler is
[`summarize_core_wing_hybrid_grid.py`](summarize_core_wing_hybrid_grid.py).

Output root: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_core_wing_hybrid_i_pilot_cluster1056_2026-10-21_01/`.
The final audit JSON is `grid_audit.json`, SHA-256
`18606dcd49e8b918186eda6158e2209f2755a76b767806756e3571d7b4990ee1`; its
checkpoint SHA-256 is
`ff6ed388dd36b727b74fc03fcffb3db89bd00bc8da2f2134dfaac9607d3450a2`.
The audit script SHA-256 is
`9794723b2b56c9656630111cbce021bb3d362f6c6ddaf1552654b7a89425bb6f`.
