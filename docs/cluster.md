# Cluster and lensing tools (plugin `cluster`, tab Measure, chip "Cluster")

Red sequence from the colour-magnitude diagram, cluster membership, strong-lensing arc candidates and overdensity maps. Engine `ogfkit/cluster.py` (pure numpy/scipy
functions, JSON-serialisable results, usable from a service), CLI `plugins/cluster/cluster.py`, synthetic truth generators `ogfkit/clustersynth.py`.

## Steps (all follow the add_columns contract and replay from the exported session script)

1. **Red Sequence + Members** (`--task members`; columns `CL_RS_RES CL_RS_PULL CL_RS CL_R CL_PMEM CL_MEMBER`). Colour = `blue-column` - `red-column`, magnitude `mag-column`,
   errors from the optional error columns plus a floor. The centre is the peak of the galaxy density map (or `center-x/-y`). Objects within `fit-radius` of the centre and in the
   magnitude window enter the fit: a slope grid gives the starting line (a tilted sequence is narrowest at the true slope), then a maximum-likelihood fit of
   `f N(colour - line; 0, scatter^2 + err_i^2) + (1-f)/(2W)` gives colour at `m0`, slope, intrinsic scatter and their errors (Hessian). `CL_RS` = within `rs-nsig` of the line.
   `CL_PMEM`: the radial distribution of the sequence galaxies is fitted by `S_bg + S_0 (1+(r/a)^2)^-beta` (unbinned Poisson likelihood over the area really inside the image);
   `p = 1 - S_bg/S(r)`. `CL_MEMBER` = sequence & p >= `pmem-min` & optional redshift window (`z-column`, `z-center` or the median of the sequence, `z-width` x (1+z); objects
   without redshift are not rejected) & optional `max-radius`.
2. **Arc Candidates** (`--task arcs`; needs image; columns `ARC_LEN ARC_WID ARC_LW ARC_RCURV ARC_SPAN ARC_ALIGN ARC_CC ARC_SNR ARC_SCORE ARC_FLAG`). Catalog objects with A/B >=
   `arc-elong-min` are cut out (half size 4A+10 px), background-subtracted, smoothed and the connected region above `arc-nsig` x sigma (measured on the smoothed image, so
   correlated noise is included) nearest to the catalog position is taken. Measured: principal axes, a weighted circle fit (curved if the radial rms is < 0.85 x the minor-axis rms and
   the span < 300 deg: closed blobs are never "curved"), length `L = sqrt(12) x rms of the along-ridge coordinate`, width = FWHM across, `LW = L/W`, radius, span, S/N, the angle
   between the major axis and the tangential direction about the centre (`ARC_ALIGN`, 0 = tangential) and `ARC_CC` = distance of the circle centre from the cluster centre / radius.
   `ARC_FLAG` = LW >= `arc-lw-min`, S/N >= `arc-snr-min`, curved (unless `arc-allow-straight`), alignment <= `arc-align-max`, `ARC_CC` <= `arc-cc-max`. `ARC_SCORE` (0-1) is a ranking
   heuristic, not a probability. The centre comes from the manual setting or from an earlier cluster step; without a centre the alignment cuts are skipped.
3. **Overdensity Map** (`--task density`; column `CL_SIGMA`; files `cluster_density.fits`, `cluster_peaks.tsv`). Counts per cell, Gaussian smoothing (`dens-sigma`), background rate from
   the clipped mean outside excesses, variance `lambda sum(k^2)` for the discrete kernel with the image footprint as validity mask; significance = (smoothed - lambda)/sd. Samples: all /
   red-sequence / members.

Windows: Colour-Magnitude Plot, Cluster Summary, Show Overdensity Map (new frame). Catalog keys `cluster,*`; catalog meta key `cluster` (centre, sequence, counts).

## Validation (`plugins/cluster/tests/test_cluster.py`, 12 tests, ~22 s; synthetic truth + real HUDF)

| Test | Result |
| --- | --- |
| red sequence, 40 clusters (150 members Plummer a=150 px + 3000 field galaxies with overlapping colours, errors 0.01-0.07 mag) | colour at m0 1.2013 (truth 1.2000), slope -0.0430 (-0.0400), scatter 0.0480 (0.0500); rms errors 0.0175 / 0.0097 / 0.0125; pull std 1.39 / 1.23 (quoted errors ~25 % optimistic) |
| membership, seed 3 (65 true members in 19.5-24 mag) | colour + radial: 70 selected, purity 0.83, completeness 0.89; + redshift: 55 selected, purity 0.96, completeness 0.82; sum of p_mem over the sequence 70.5 vs 66 true members; 12 clusters: sum p / true 1.06 |
| auto centre | 37 px from truth (smoothing 120 px) |
| density significance, 12 uniform Poisson fields | map mean 0.018, std 0.969 (expect 0, 1), max 4.19 |
| density peak, 250 members + 2500 field, 10 % of the field masked | 9.2 px from the centre, 34.8 sigma; background per cell 0.00977 vs 0.01000 |
| arcs, 6 fields, 59 arcs (R 60-220 px, span 35-110 deg, FWHM 3-6 px, amp 3-8 sigma) + radial/straight/blob decoys (252) | 57/59 flagged, 0/252 false; length error median -0.4 % (rms 3.4 %), radius -0.0 % (rms 2.4 %), span +0.5 deg, alignment <= 0.7 deg |
| same, faint arcs (amp 1.2-2.5 sigma) | 19/59 flagged, 0/252 false; length low by 47 % (fragmented), radius +1 % (rms 9.7 %) |
| decoys | radial lines align 88.7-90 deg (0 flagged); tangential straight lines 0/6 flagged, 6/6 with `arc-allow-straight` |
| real HUDF crop 1800x1800 (F105W, F160W; 1827 objects S/N>8) | density max 6.1 sigma (galaxy clustering; no cluster); CMD fit finds a "sequence" colour 0.257 scatter 0.226 but the radial model gives 16.6 expected members of 400 and 8 selected (0.4 % of the objects) |
| HUDF F160W + 10 injected arcs (8 sigma, empty sky), real noise/neighbours | 6/10 flagged (the others merged with bright neighbours in the detection catalog), length error median -0.3 %, 0 other flagged among 1202 detections |

## Limits

* The sequence is a straight line; no red-sequence evolution, no blue population, no BCG treatment. The intrinsic scatter is poorly constrained when measurement errors dominate (it
  can collapse to the 0.005 floor with < 60 members). Fit radius and magnitude window matter: a field-dominated fit region gives a meaningless sequence (the HUDF test shows it).
* `CL_PMEM` assumes a roughly circular profile and a uniform field; clustering of the field (HUDF) and foreground structure bias it.
* Arc measurement works on the catalog's segmentation by position: arcs merged with neighbours by the detection (4 of 10 injected on real HUDF) or fragmented (faint arcs) are not recovered;
  no lens modelling, no photometric-redshift check of the arc; spiral arms and ring galaxies can mimic curved features (not tested beyond Gaussian blobs).
* Density maps use the catalog extent when no image is available; the validity mask comes from nonzero finite pixels of the current image.
