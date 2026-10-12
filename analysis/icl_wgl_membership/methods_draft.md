# Methods draft: scaled cluster light and a weak-lensing mass in one desktop program

This is a methods draft for Astrafex Desktop. It describes software that puts a scaled surface-brightness stack and a weak-lensing mass on the same coordinates, the same flat \(\Lambda\)CDM settings, and the same membership table. It is not a measurement on a named survey. The synthetic run in this directory is a closed self-consistency check. Section 5 of `implementation_plan.md`, the gate that asks for a named data set, is not part of this draft.

## Light

`icl.wgl.stack_clusters` warps each cluster image to \(u=\theta/\theta_\Delta\). With \(\Delta=200\), \(\theta_\Delta\) follows from \(M_{200c}\) and the angular-diameter distance. The product key is `bcg+diffuse`, because the mean still contains the BCG. A pure diffuse fraction names the centre cut, the surface-brightness limit, or the model decomposition that defined it. The same frames are averaged on one proper-kpc grid and returned beside the \(u\) stack. The comparison coordinate is that kpc stack. The warp is counts per square arcsecond. That unit does not put frames with different zero points on one scale. The default mean divides by the number of covered pixels. Pixel weights in the warp are 0 or 1. An inverse-variance mean is a separate array, written only when the caller passed a variance map.

A colour can be formed from surface-brightness profiles in more than one band. The desktop colour profile writes those surface-brightness columns and their differences. A metallicity, an [Fe/H], or a \(Z/H\) is not a column of that menu and not a product of the scaled stack.

\(M_{200c}\) is the mass inside the radius where the mean density is \(200\rho_{\rm crit}\). \(M_{200m}\) is the mass inside the radius where the mean density is 200 times the mean matter density. The stack and the lens fit use \(M_{200c}\). A mass quoted as \(M_{200m}\) or \(M_{500c}\) needs its conversion written beside the radius. With redshift and cosmology fixed, a mass error moves the alignment radius by one third of the logarithmic mass change. The forward model averages the caller's mass draws. The \(u\) stack itself uses the point mass.

The desktop one-frame profile (`ds9/library/ds9_icl.py --mode profile`) is a separate path. The scaled stack is the menu item Stack BCG+diffuse. The program is `plugins/icl/icl_stack.py`. It calls `stack_clusters` and writes `icl_stack_u.fits`, `icl_stack_kpc.fits`, `icl_stack_profile.tsv` and `icl_stack.json`. A frame marked local-sky, or a named product other than `deepCoadd`, is refused. `science_measurement` stays false while the survey checklist is empty.

## Membership

`icl.members.classify_members` applies the caller's spectroscopic window. The software's authority is that window. A fixed \(dz\) has a rough velocity width \(\Delta v \approx c\, dz/(1+z_{\rm cl})\). Whether the window is a good physical cut is a separate question. A finite redshift with \(|z-z_{\rm cl}|<dz\) is a satellite. The boundary \(|z-z_{\rm cl}|=dz\) is outside the window. Every other finite redshift, including a foreground galaxy, receives the label `background`. A missing redshift stays `unknown` until photometry is supplied.

A point photo-\(z\) or a \(P(z)\) fills only rows that still have no finite spectroscopic redshift. A row in one of the caller's windows takes that label. A row in both, or in neither, stays `unknown`. The red sequence is a straight line fitted to spectroscopic satellites only. It does not recover a blue satellite. A galaxy off the line stays `unknown`, and missing the line does not create a background label. An `unknown` row on the line, brighter than the caller's magnitude limit, becomes a satellite. The surface-brightness stack masks satellites, background, stars and `unknown`, and leaves the BCG. A light fraction whose denominator is BCG plus satellites therefore drops an `unknown` blue satellite, and that fraction can rise. A projected red non-member promoted to a satellite enters the denominator, and the fraction can fall. A photometric label does not replace the redshift used for \(\Sigma_{\rm crit}\). A chosen BCG index and a star flag stay ahead of the redshift and colour labels.

Weak lensing keeps a row only when the label is `background`, the redshift is finite, and the source sits behind the lens by the redshift gap.

## Mass

The mass is a fit to the tangential shear of each galaxy,

\[
g_t = -g_1\cos 2\phi - g_2\sin 2\phi,
\]

with \(\phi\) measured from \(+x\) toward \(+y\). A positive \(g_t\) is a tangential alignment. The model is a spherical NFW halo in \(M_{200c}\). When the reduced-shear switch is on, each source is predicted as \(\gamma/(1-\kappa)\).

A stack of cluster images, with shapes measured on that coadd, is not this mass. The intracluster-light image stack is surface brightness. `stack_shear_profiles` averages shear catalogues on proper Mpc or on \(R/R_{200c}\), then fits that average. An image stack is rejected. The plugin command is `--task stack`. The synthetic test runs, from the OGFinder tree:

```bash
python plugins/wlcluster/wlcluster.py --task stack --clusters clusters.json \
  --work wlstack --radius mpc --rmin 0.4 --rmax 2.5 --nbins 4 \
  --concentration fixed --c-fixed 4 --no-two-halo --power linear \
  --no-lss-cov --sigma-e 0.001 --samples 0 --min-per-bin 8
```

The catalogue and work paths stand in for the temporary files in `plugins/wlcluster/tests/test_wlcluster.py`. The flags are the flags of `test_cli_stack_writes_the_average_profile`. The noiseless one-cluster test in the same file runs:

```bash
python plugins/wlcluster/wlcluster.py --task mass --catalog shapes.tsv \
  --work wlrun --coord xy --center-x 1000 --center-y 1000 --pixel-scale 1 \
  --z-lens 0.3 --z-source 1 --z-column '' \
  --concentration fixed --c-fixed 4 --no-two-halo \
  --power linear --no-lss-cov \
  --sigma-e 0.001 --samples 20 --seed 9 --min-per-bin 8
```

Menu defaults for the one-cluster mass task stay Duffy concentration, Takahashi et al. 2012 Halofit with the two-halo term, 200 samples, and a shape-noise plus Gaussian large-scale-structure covariance. The listings above are the closed tests, and they turn those pieces off or hold them fixed.

`DELTASIGMA` is the mean tangential shear divided by the mean of \(\Sigma_{\rm crit}^{-1}\), in solar masses per square parsec. The mass model uses the reduced shear of each source. Where reduced shear matters, the column is not an exact \(\Delta\Sigma\). \(\langle\Sigma_{\rm crit}^{-1}\rangle\) and \(\Sigma_{\rm crit}^{-1}(\langle z\rangle)\) are different quantities. A photometric \(n(z)\) is returned separately and is not copied into the spectroscopic redshift column.

`science_measurement` stays false until the caller has supplied a non-empty spectroscopic intersection of the shape catalogue with a trusted redshift catalogue, a constrained fit, a sensitivity table, notes on non-sphericity and projection, and a BCG–mass offset. A filled gate is not a named-data mass and not a survey calibration. `measure_stamps` uses one round Gaussian PSF and a per-object diagonal response. That is not an S19A re-Gaussianization, a KiDS lensfit calibration, or a DES Y3 metacalibration.

`stack_tangential_stamps` is a postage-stamp diagnostic, not a mass. Each stamp is rotated by minus its position angle and averaged in annuli. The unrotated average is the null image when the tangential directions differ. The tangential ellipticity on the aligned stamp is \(-e_1\), the same sign as \(g_t\). There is no menu task for this function. A cube of cluster images is rejected.

Mean stellar light and mean mass describe one population when the light weights and the lensing weights are the same weights. If they differ, the stack result says the two averages do not represent the same population. The saved run in this directory still fits one lens. It does not set `science_measurement`, and it does not contain the postage-stamp stack.

The tests named in the user-guide sections 15.5 and 15.6 are:

```bash
python -m pytest plugins/wlcluster/tests/test_wlcluster.py icl/tests/test_wgl.py icl/tests/test_members.py -q
```

## What the synthetic run records

| Standing | Statement |
|---|---|
| Recorded | Coordinate, unit and mask self-consistency of the scaled stack. |
| Recorded | The membership rules above, executed on the rows the caller passed. |
| Recorded | The fit path that is given an input NFW shear. |
| Recorded | Inputs and outputs that can be run again from `run_analysis.py`. |
| Not from this run | A surface-brightness accuracy on the sky, or a detection limit. |
| Not from this run | Membership completeness or purity. |
| Not from this run | A shear or mass accuracy on the sky. |
| Not from this run | A light fraction, self-similarity, or evolution measured on the sky. |

Individual NFW fits belong in the same methods paper as the light stack because both paths can use one coordinate system, one cosmology, and one membership table. That is a software example. A sentence that relates mean stellar light to mean mass waits until both averages are made on one named sample, with the weight difference recorded when the weights differ.

## Data availability

The synthetic run is in the Astrafex Desktop repository, <https://github.com/kjhan0606/OGFinder>, under `analysis/icl_wgl_membership/`. Cite the commit that contains this directory, `e441eb7f5faac9d2c278dfe65c8dff1d2cfd5e32`. `history.json` records `parent_commit` `9a4adcb8f421838115acbff9c8612671d4fff723`, the commit that was HEAD when that file was written. The two hashes are different commits. The citation hash is not written into `history.json`.
