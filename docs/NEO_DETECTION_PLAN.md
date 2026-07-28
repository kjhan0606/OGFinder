# NEO Detection Pipeline Plan

## Objective

Add a dedicated NEO detection workflow to the pull-down menus in the
OGFinder/DS9 right-hand catalog panel.

The production baseline is:

```text
FITS sequence
  -> astrometric and timing validation
  -> ZOGY difference imaging
  -> point-source and trail detection
  -> same-night tracklet linking
  -> known-object matching
  -> MPC digest2 scoring
  -> preliminary orbit and ADES export
```

A single image can identify only streak-like moving-object candidates. A
defensible NEO detection requires at least three, and preferably four or more,
exposures of the same field with accurate observation times.

## Existing Integration Points

- `ds9/library/ds9.tcl`
  - Creates the top-level horizontal paned window.
  - Creates the right-hand catalog panel.
- `ds9/library/layout.tcl`
  - `CreateCatalogPanel` defines the current pull-down menus.
  - `CatalogPanelGetScript` locates Python analysis commands.
  - `CatalogPanelGetFITS` obtains the current FITS filename.
  - `ds9(frames)` provides the list of loaded DS9 frames.
  - Existing TSV loading, row selection, marker, status, and export mechanisms
    should be reused.
- Existing architecture:
  - Tcl/Tk owns menu state, dialogs, status, frame navigation, and markers.
  - A standalone Python CLI owns scientific processing.

Implementation must preserve this separation.

## Proposed Right-Panel Menu

Add a top-level `NEO` pull-down menu:

```text
NEO
  1. Build Image Sequence...
  2. Validate Astrometry & Timing
  3. Build Reference / Difference Images
  4. Detect Moving Sources
      Point Sources
      Streaks
      Synthetic Tracking
  5. Build Tracklets
  6. Match Known Objects (MPC)
  7. NEO Score (digest2)
  8. Preliminary Orbit...
  Review Candidates...
  Export ADES...
  Settings...
```

Place `NEO` between `LSBG` and `Analysis`. If the additional menu overflows the
right panel, change the menu container to a responsive two-row layout instead
of shortening scientifically meaningful labels.

## Data Requirements

Each input exposure must provide or obtain:

- FITS image and valid celestial WCS.
- Mid-exposure time derived from `DATE-OBS`/`MJD-OBS` and `EXPTIME`.
- Exposure-time convention: start, midpoint, or end must be explicit.
- Gain, read noise, saturation mask, bad-pixel mask, and variance where
  available.
- Filter and photometric zeropoint where available.
- Ground observatory MPC code and location, or a space-observer ICRF state
  vector.

The sequence builder will:

1. Collect filenames from `ds9(frames)` or a user-selected directory/list.
2. Remove FITS HDU suffixes only for filesystem access while retaining the
   selected HDU in the manifest.
3. Sort exposures by mid-exposure UTC.
4. Show a validation table containing filename, time, exposure, filter, WCS,
   pixel scale, seeing, and quality flags.
5. Reject duplicate times, missing WCS, inconsistent fields, or insufficient
   temporal baseline before detection.

## Scientific Pipeline

## Implemented DS9--CODES bridge

The first integration path is now present.

- `ds9/library/ds9_neo.py` calls the CODES Python package.
- `layout.tcl` adds the `NEO` pull-down menu.
- `NEO -> Fit CODES Orbit from Astrometry CSV...` accepts explicit optical
  astrometry and loads the resulting FK5 regions into the active DS9 frame.
- CODES writes `preliminary_orbit.json` and a refined state table using the
  normal DE442 or DE441 force model.
- The result is labelled preliminary and is not treated as an MPC orbit.

The current solver is an initial-state fit. A short optical arc can have a
distance and radial-velocity degeneracy, and a least-squares fit can select a
formally hyperbolic branch even for a known asteroid. The production orbit
path therefore requires statistical ranging or OpenOrb before orbit elements,
close-approach distances, or impact probabilities are shown as authoritative.

### 1. Calibration and registration

- Verify or solve the WCS with Astrometry.net.
- Refine relative astrometry from stationary stars.
- Resample images and masks onto a common tangent plane.
- Estimate a spatially varying background, variance map, and PSF for each
  exposure.
- Record astrometric residuals rather than assuming FITS WCS accuracy.

Reference:

- Astrometry.net: https://astrometry.net/doc/readme.html

### 2. Reference and proper difference imaging

Preferred reference order:

1. User-supplied deep reference from the same instrument/filter.
2. Archived external reference with compatible bandpass and adequate depth.
3. Leave-one-out robust reference constructed from the input sequence.

Use ZOGY proper subtraction to produce:

- Proper difference image.
- Difference PSF.
- Corrected significance image.
- Bad-pixel/artifact mask.

Use a maintained implementation such as `properimage` behind a local adapter.
Do not spread third-party API calls through the rest of the pipeline.

References:

- Zackay, Ofek, and Gal-Yam:
  https://arxiv.org/abs/1601.02655
- Properimage:
  https://properimage.readthedocs.io/

### 3. Point-source and trail detection

Run two deterministic detectors on each significance image:

- PSF matched filter for sources whose motion during one exposure is less than
  approximately half a PSF FWHM.
- Linear trail matched-filter bank over configurable position angle and trail
  length for faster objects.

Fit every detection for:

- Pixel and ICRF centroid.
- Mid-exposure UTC.
- Flux, magnitude, S/N, and uncertainties.
- Trail length and position angle.
- PSF/trail goodness of fit.
- Mask, edge, saturation, dipole, and subtraction-residual flags.

### 4. Stationary-source and artifact rejection

- Reject detections that recur at a fixed ICRF position.
- Treat variables separately from movers.
- Reject cosmic rays, diffraction spikes, saturated bleed trails, image
  boundaries, bad columns, and subtraction dipoles using explicit flags.
- Preserve rejected detections for completeness and false-positive analysis.

Machine learning must not be the primary detector. A DeepStreaks-like model may
later be added as an artifact-veto score after deterministic completeness is
measured.

Reference:

- DeepStreaks: https://arxiv.org/abs/1904.05920

### 5. Same-night tracklet linking

Use a Pan-STARRS/Rubin MOPS-style same-night linker:

1. Build possible detection pairs in angular-rate bounds using a KD-tree.
2. Propagate each pair to other exposure times.
3. Search for compatible detections in the astrometric uncertainty ellipse.
4. Fit constant angular velocity.
5. Optionally fit angular acceleration for close, fast NEOs.
6. Require at least three detections by default; four is the accepted-cadence
   preset.
7. Resolve shared-detection conflicts by maximizing a tracklet score based on
   detection S/N, fit residual, number of epochs, and artifact penalties.

Store the full detection-to-tracklet association. Do not retain only a fitted
line.

Reference:

- Pan-STARRS MOPS: https://arxiv.org/abs/1302.7281

### 6. Synthetic tracking

Provide an optional deep mode for many short exposures:

- Search a configurable two-dimensional angular-velocity grid.
- Shift each image or score map to a common epoch and coadd.
- Use a coarse-to-fine velocity search.
- Correct the detection significance for the number of velocity trials.
- Refine position, velocity, flux, and astrometric uncertainty after detection.

Use Numba for the required CPU backend. Add CuPy only as an optional GPU
backend; the scientific result must not depend on GPU availability.

Reference:

- NEO observations using synthetic tracking:
  https://arxiv.org/abs/2401.03255

### 7. Multi-night linking

Same-night tracklets are the first production target. For sparse, irregular
multi-night data, integrate THOR as an optional backend instead of developing a
new cadence-independent linker.

Reference:

- THOR: https://b612.ai/opensource/thor/

### 8. Known-object identification

- Query the MPC MPChecker service with tracklet epoch, position, field radius,
  magnitude, and angular motion.
- Cache all queries and responses with timestamps.
- Provide an offline mode using a cached MPC orbit catalog and OpenOrb
  propagation.
- Display match separation and velocity residual, not only the object name.
- Never classify an unmatched object as new solely because a network query
  failed.

Reference:

- MPC services: https://docs.minorplanetcenter.net/services/

### 9. NEO prioritization

Run the official MPC `digest2` package on accepted short-arc tracklets:

- Store all class scores and the NEO score.
- Treat `D2 >= 65` as a follow-up priority threshold.
- State clearly that digest2 is a short-arc pseudo-probability and is not an
  orbit-based NEO confirmation.

An object is confirmed as an NEO only after a sufficiently constrained orbit
shows perihelion distance `q < 1.3 au`.

Reference:

- MPC software and digest2:
  https://docs.minorplanetcenter.net/software/

### 10. Preliminary orbit

- Use OpenOrb statistical ranging for short arcs and non-Gaussian uncertainty.
- Show orbit samples and ephemeris uncertainty, not only a nominal six-element
  solution.
- Permit CODES as an optional refinement/visualization backend after the
  short-arc orbit interface is validated.
- Do not report impact probability from a nominal tracklet orbit.

The current CODES bridge implements the state hand-off and full-force-model
refinement. It does not yet claim to solve the short-arc ranging problem.
This distinction is enforced in the JSON summary and in the DS9 status text.

### 11. Public-image validation

`/home/kjhan/BACKUP/CODES/validate_dad_tracklet.py` uses the public NOIRLab
DECam Asteroid Database DR2 and its SIA calibrated-image service. The default
case downloads five CCD cutouts for MPC object 14941 and reads the measured
RA, Dec, and exposure time from the DAD tables. The runner compares the
astrometric fit with NASA/JPL Horizons elements. DAD also provides two
additional five-exposure tracklets for the same object on adjacent nights.
The single-night result validates image access, timing, WCS-compatible
astrometry ingestion, tracklet handling, and CODES propagation. A multi-night
statistical-ranging result is required before the element comparison can be
called an independent orbit validation.

Reference:

- OpenOrb:
  https://onlinelibrary.wiley.com/doi/10.1111/j.1945-5100.2009.tb01994.x

### 11. ADES export

- Export accepted astrometry in MPC ADES PSV.
- Include per-observation timing, RA/Dec uncertainties, magnitude uncertainty,
  band, station code, and observer state for space-based observations.
- Validate the generated file against the ADES schema.
- Require explicit human approval before submission.

Reference:

- MPC ADES documentation:
  https://docs.minorplanetcenter.net/software/

## Candidate Table and Review UI

The catalog panel should switch to a tracklet-level table with columns such as:

```text
TRACKLET_ID
NDET
RA
DEC
MJD_MID
RATE_RA
RATE_DEC
RATE_TOTAL
PA
MAG
SNR
TRAIL_LENGTH
ASTROM_RMS
KNOWN_ID
MATCH_SEP
D2_NEO
ORBIT_Q
STATUS
```

Selecting a tracklet must:

- Place predicted-position markers in every contributing DS9 frame.
- Pan to the candidate.
- Blink the sequence.
- Show fixed-size cutouts, difference cutouts, and residuals.
- Plot measured positions against the fitted motion model.
- Show detection flags, known-object match, digest2 score, and orbit
  uncertainty.

Suggested marker colors:

- Cyan: known minor planet.
- Yellow: unmatched moving-object candidate.
- Orange: high-priority digest2 candidate.
- Red: accepted urgent follow-up candidate.
- Gray: rejected/artifact candidate.

## Proposed Python Layout

```text
neo/
  __init__.py
  config.py
  sequence.py
  astrometry.py
  reference.py
  subtraction.py
  detection.py
  trails.py
  linking.py
  synthetic_tracking.py
  known_objects.py
  digest.py
  orbit.py
  ades.py
  injection.py
  products.py

ds9/library/
  ds9_neo.py

tests/neo/
  test_sequence.py
  test_astrometry.py
  test_subtraction.py
  test_point_detection.py
  test_trail_detection.py
  test_linking.py
  test_synthetic_tracking.py
  test_known_objects.py
  test_digest.py
  test_orbit.py
  test_ades.py
  test_end_to_end.py
```

`ds9_neo.py` should expose subcommands rather than separate scripts:

```text
ds9_neo.py sequence
ds9_neo.py calibrate
ds9_neo.py subtract
ds9_neo.py detect
ds9_neo.py link
ds9_neo.py synthetic
ds9_neo.py identify
ds9_neo.py score
ds9_neo.py orbit
ds9_neo.py export-ades
ds9_neo.py run
```

Each stage must write machine-readable products so it can be resumed without
rerunning earlier stages.

## Output Layout

Store run products outside the source tree:

```text
~/.ds9/neo/<sequence-id>/
  manifest.ecsv
  config.yaml
  calibrated/
  references/
  differences/
  score_maps/
  detections.ecsv
  tracklets.ecsv
  associations.ecsv
  known_matches.ecsv
  orbit_samples.ecsv
  candidates/
  accepted.ades.psv
  run.json
  log.txt
```

The run manifest must record software version, configuration, input checksums,
reference source, catalog versions, and network-query timestamps.

## Dependencies

Already available in the current environment:

- NumPy
- SciPy
- Astropy
- SEP
- scikit-learn
- Numba
- Astroquery

Additional required or optional dependencies:

- `properimage`: required for production ZOGY adapter.
- `reproject`: required for robust WCS image/mask reprojection.
- `digest2`: required for MPC NEO scoring.
- OpenOrb/PyOrb: optional until preliminary-orbit stage.
- CuPy: optional synthetic-tracking acceleration.

Pin NEO-specific dependencies in a separate requirements or environment file.
Do not make optional orbit/GPU packages mandatory for opening DS9.

## Implementation Phases

### Phase 1: Sequence and deterministic tracklet baseline

- Add the NEO menu and settings.
- Build/validate an image sequence from DS9 frames.
- Extract per-frame detections.
- Remove stationary sources.
- Link same-night tracklets.
- Display tracklets and blink candidates.

This establishes coordinate, time, UI, and product contracts before image
subtraction is introduced.

### Phase 2: Production difference-image detector

- Add astrometric refinement and reprojection.
- Add reference construction.
- Integrate ZOGY.
- Add point-source and trail matched filters.
- Propagate variance and quality flags.

### Phase 3: External identification and reporting

- Add MPC known-object matching and cache.
- Add digest2 scoring.
- Add OpenOrb statistical ranging.
- Add candidate acceptance workflow and ADES export.

### Phase 4: Deep and sparse-cadence modes

- Add CPU synthetic tracking.
- Add optional GPU backend.
- Add THOR multi-night backend.
- Add optional ML artifact-veto model only after deterministic validation.

## Validation and Acceptance

### Synthetic injection

Inject synthetic moving sources into real images before subtraction over:

- Magnitude and S/N.
- Angular rate and direction.
- Exposure-time trail length.
- Seeing and background.
- Distance to stars, galaxies, bad pixels, edges, and saturation artifacts.
- Linear and accelerated apparent motion.

Measure:

- Detection completeness.
- Tracklet completeness.
- Candidate purity.
- Astrometric and photometric bias.
- Velocity bias.
- False candidates per image and per square degree.
- Processing time and memory.

### Real-data validation

- Recover known NEOs and main-belt asteroids from public multi-epoch sequences.
- Compare identification to MPC results.
- Compare fitted positions and rates to official ephemerides.
- Verify digest2 scores against official examples.
- Validate preliminary orbit distributions against longer-arc official orbits.

### Required behavior

- No silent fallback when timing, WCS, or observer state is missing.
- No automatic claim that an unmatched tracklet is a new NEO.
- No automatic MPC submission.
- No nominal-orbit impact probability.
- Every rejection must retain a reason flag.
- Every result must be reproducible from the run manifest.

## Branching and Worktree Note

The repository currently contains unrelated modified and untracked files.
Implementation must begin on a dedicated branch or worktree. Do not reset,
clean, or overwrite the existing working tree when creating the NEO branch.
