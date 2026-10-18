# HSC-I BCG+diffuse kpc pilot

Approved sequence: `implementation_plan.md`, 2026-10-12 section. Data and outputs
stay under `/gpfs/kjhan/Astrafek/data`; no catalog, FITS, or cluster log is committed.
Existing 01–08, `history.json`, and `run_analysis.py` are preserved.

## Definition and scope

Keep the BCG and mask all other detections, defects, and bright-star areas.
This is BCG+diffuse, not a pure ICL fraction. SEP background subtraction is for
detection and Kron radii only; profiles use the original calibrated plane.
Use explicitly calibrated, Galactic-extinction-corrected nJy/arcsec² in the
observed HSC-I band. There is no default zero point or dust correction.
Negative linear flux is retained. Native-pixel annular profiles are averaged
with equal cluster weights. Coverage includes cropped area in the full annulus.

Area-overlap 2D images supply a separate radial diagnostic. They integrate
input pixel overlaps; their axes follow the input image and do not measure
BCG-aligned morphology. Square orthogonal pixels and <0.1% checked WCS scale
variation are required. Larger distortions require a full reprojection.
Defaults: 0–500 kpc, 25 kpc annuli, 5 kpc image spacing, coverage threshold 0.5.
These are engineering choices, not validated completeness criteria.
Kron factors 2.5, 3.5, and 5 are separate variants. Missing BCG detection stops
the frame. Unmatched member positions get max(2 arcsec, 1.5 pixels) circles.
Blending and PSF wings still need inspection and realistic injection/recovery.

## Prepare and authenticate

From OGFinder, with a new output directory:

```bash
python3 analysis/icl_wgl_membership/hsc_pilot.py prepare \
  --work /gpfs/kjhan/Astrafek/data/pilot_2026-10-12
python3 analysis/icl_wgl_membership/hsc_debug/probe_login.py
```

The prepared manifest samples twelve equally spaced richness ranks from the
preserved selected-cluster table and joins members. It is an engineering
pilot, not a representative sample. Seeing/depth/area still require review.
Image/calibration/quality fields deliberately remain empty. Requested field
half-width corresponds to 1000 kpc, not a zero-ICL background boundary.

The probe prompts for account and hidden password. Enter them after the
prompts, not in the command line. Slash and shell metacharacters are read
literally. No credentials are written to arguments/files/logs. The default
probe requests one `type=coadd`; that does not verify global-sky provenance.
The old IMAGE-only header also omits inherited primary calibration. To retain
the authenticated official manual and all headers of one `coadd/bg` probe:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/probe_login.py --global-sky
```

This requires that the official manual documents `coadd/bg`, then requests
one 20 arcsec half-width image with mask/variance. It stores new, exclusive
files `logs/cutout_global_manual.html`, `cutout_global_headers.txt`, and
`cutout_global_status.json`. It does not retain the image or bulk-download.
The authenticated PDR3 manual maps `coadd/bg` to
`{rerun}/deepCoadd/{filter}/{tract}/{patch}.fits` and describes it as the
intermediate coadd before aggressive background subtraction. The probe's
primary header has `FILTER=i`, `FLUXMAG0=63095734448.01944` (AB zero point
27.0); HDUs 1/2/3 are IMAGE/MASK/VARIANCE, each 239×239 and inheriting the
primary header. The image is calibrated instrumental flux per pixel despite
its absent `BUNIT`. Keep the product string `coadd/bg`; the reviewed source
record is `/gpfs/kjhan/Astrafek/data/logs/cutout_global_review.json`.
PDR3 mask cards give bit positions (`MP_*`); convert a position `n` to its
bit value `2**n` when filling `exclude_bits`. Exclude reviewed bad-data planes
and `BRIGHT_OBJECT`; leave `DETECTED` to the source-mask stage so the BCG core
can remain. Inspect real mask pixels and mask growth before fixing the policy.

This probe saves headers only; it does not retain the 239×239 pixel planes.
The twelve pilot frames, random fields, Galactic-extinction values and
sky/PSF/injection checks remain pending. Product mapping is verified for this
one request, not its residual-background performance. Preserve the product
name `coadd/bg`.
For a repeat probe, supply a new basename, for example
`--global-sky --record-name cutout_global_02`; existing records are preserved.

## Download the pilot pixel files

The actual 12-frame download is separate from the header probe. It uses the
network and CPU only; it does not use a GPU. The prepared transfer is about
0.62 GiB, so start it in the foreground from the syntax SSH shell and use this
direct command only while it remains a short transfer (a few minutes):

```bash
python3 /home/kjhan/BACKUP/Astrafek/OGFinder/analysis/icl_wgl_membership/hsc_debug/download_pilot.py \
  --user kjhan0606 \
  --manifest /gpfs/kjhan/Astrafek/data/pilot_2026-10-12/pilot_manifest.json \
  --work /gpfs/kjhan/Astrafek/data/pilot_2026-10-12/download_01
```

The login-node rule is to move work that lasts more than a few minutes into
Slurm. The currently listed Slurm partitions are GPU-backed, and the site rule
requires a GPU request for Slurm jobs even when a task leaves it idle; use that
route only if the transfer is too long for the login node.

The command uses the account name already confirmed in the syntax session and
prompts once for its hidden password. The password is read literally, so `/`
needs no escaping; it is not put in the command line, environment, or files.
It downloads IMAGE, MASK, and VARIANCE FITS planes for all twelve clusters,
validates their FITS structure, WCS, calibration header, and mask-plane
definitions, then records each SHA-256. PDR3 coadds combine the similar but
not identical `i` and `i2` filter variants, so the validator accepts either
header value and records the actual `FILTER` for each frame
(<https://hsc-release.mtk.nao.ac.jp/doc/index.php/processing__pdr3/>).
Success prints `Downloaded 12/12 verified FITS files` and leaves the pixel
files in `download_01/images/`, with `download_manifest.json` and
`download_status.json` beside that directory. The manifest deliberately keeps
`science_measurement=false`; extinction, matched random fields, and the sky,
PSF, mask, and injection checks still need review. A stopped/interrupted run
preserves its status and files; rerun against the same directory with
`--resume` to verify completed responses and download only missing frames.

## Required frame records

Paths are relative to the manifest unless absolute. Each frame needs:

* `image`, explicit `science_hdu` (index or name), catalog coordinates/redshift.
* `sky=global`/`global-sky`, `product=deepCoadd`/`coadd/bg`,
  `product_reviewed=true`, and `product_evidence`.
* `calibration`: `input_unit=nJy/pixel` or `instrumental_flux/pixel`, `evidence`,
  numeric `a_i_mag`, and `extinction_evidence`. Instrumental flux also needs
  `ab_zeropoint` matched to the stored values. Review HSC inherited primary
  calibration. Exposure normalization is not guessed or applied twice.
* `bad_mask`, `bright_star_mask`: `file`, `hdu`, `evidence`,
  `same_grid_reviewed=true`, `kind=boolean` (0/1) or `kind=bits`.
  Bit masks also need `exclude_bits`, powers-of-two **values**, identified from
  that product's mask definitions. Do not blindly mask DETECTED or use DS9
  mask-manager definitions as HSC definitions. A reviewed empty star raster
  is acceptable only for a clean field.
* `variance`: `file`, `hdu`, `evidence`, `same_grid_reviewed=true`.
  Nonfinite/nonpositive variance pixels are excluded. It is not used to
  manufacture independent-pixel errors after resampling.

`randoms` holds matched random frame records with these same inputs and the
associated lens redshift. Centres are the random sky positions. Supply
`random_matching_evidence` on each cluster and spatial `region` labels.
Match footprint, depth, seeing, extinction, masks, area and selection; arbitrary
offsets are not a scientific random catalog. Regions and random positions
are not generated by this tool.

## Inspect, then run on Slurm

```bash
python3 analysis/icl_wgl_membership/hsc_pilot.py inspect \
  --manifest /gpfs/kjhan/Astrafek/data/pilot_2026-10-12/pilot_manifest.json \
  --work /gpfs/kjhan/Astrafek/data/pilot_2026-10-12/inspection_reviewed
sbatch analysis/icl_wgl_membership/pilot_job.sh \
  /gpfs/kjhan/Astrafek/data/pilot_2026-10-12/pilot_manifest.json \
  /gpfs/kjhan/Astrafek/data/pilot_2026-10-12/measurement_01
```

All output directories must be new. Inspection reads headers and reports
missing inputs. Measurement refuses to run outside Slurm and validates all
headers first. The job requests 16 GB, two CPUs, and one GPU per cluster
policy (the calculations use CPUs); modules are reset inside the job.

Per Kron factor: `kpc_stack.fits`, `cluster_count.fits`, `frames.json`,
`profile.tsv`, and `covariance.json`. Image axes are marked in kpc. Images
are raw diagnostics; background correction belongs to the profile only.
The TSV has raw/paired-raw/random/corrected SB, raw/paired cluster counts, jackknife sigma,
image SB and image-minus-native SB. Missing values are NA/JSON null.
`CORRECTED_SB = PAIRED_RAW_SB - RANDOM_SB`; `RAW_SB` may include additional
clusters lacking random coverage and is used for the raw image comparison.
Frame records retain native profiles, coverage, usable area and source masks.

Random profiles are averaged per lens, then subtracted with equal cluster
weights. A bin requires every random field in that set. Regional delete-one
jackknife excludes cluster and random terms together. Without three regions
there is no covariance; only bins supported by every replicate receive it.
This is a pilot estimate, not a final survey covariance or systematic budget.
Every result has `science_measurement=false`. Sky transfer, extended PSF,
mask completeness, observed profile agreement and realistic injection/recovery
remain pending. There is no mass or light fraction.

## Regression evidence

```bash
sbatch analysis/icl_wgl_membership/pilot_job.sh --test
```

Tests cover checkerboard aliasing, integrated impulse, masked bright flux,
cropped coverage, calibration, non-member masks/BCG retention and joint
regional noise cancellation. These are small algorithm checks, not survey
injection/recovery. Job 420654: **35 passed in 8.94 s**, including the existing
scaled/menu regressions, inherited primary calibration, an offline literal-slash
authentication probe, and complete synthetic CLI measurement inside Slurm.
Log is under `data/logs/hsc-icl-pilot-420654.out`.

Same-lens covariance is now `covariance_recorded=true` when supplied and
`covariance_propagated=false`. The legacy survey checklist checks non-empty
metadata, not scientific acceptance. Mass-error inputs are not yet forwarded
by the scaled menu; selection/Eddington and joint lensing inference are later
stages.
