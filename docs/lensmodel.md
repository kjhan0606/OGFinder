# Strong-lens modelling (`plugins/lensmodel`, `ogfkit/lensmodel.py`)

Parametric lens model for multiply imaged sources (round 2, item 7).  numpy/scipy/astropy only (matplotlib's `contourpy` is used for the critical-curve contours);
**lenstronomy is optional** and used only as an independent cross-check (`--check-lenstronomy`, `lenstronomy_alpha`).

## What it does

| piece | implementation |
|---|---|
| mass model | **SIE** (θ_E, axis ratio q, PA φ, centre) + **external shear** (γ, φ_γ); optional fixed **NFW** halo (κ_s, θ_s) and **PIEMD** cluster-member perturbers (b0 ∝ L^0.5, cut radius ∝ L^0.5, from a catalog flag column + magnitude) |
| fit | **image-plane χ²**: for each trial model the lens equation is solved *locally* (Newton from each observed image) for a fitted source position (β_x, β_y); residual = observed − predicted position / σ_pos.  Initial convergence uses the linearised (A⁻¹Δβ) residual; 8 multi-starts.  Optional flux-ratio term ln(μ_i/μ_0) vs ln(F_i/F_0).  Free set by number of images (`auto`: 2 images → θ_E; 3 → θ_E,q,φ; ≥4 → + γ,φ_γ) or by name; the lens centre is held at the lens galaxy |
| Einstein radius, mass | θ_E (parameter) and θ_E,eff = √(area/π) of the tangential **critical curve** (differs from θ_E once shear/perturbers are present); M(<θ_E,eff) = π (D_l θ)² Σ_cr, Σ_cr from astropy flat ΛCDM (H0, Ωm settable); SIS-equivalent σ_v |
| critical curves / caustics | contours of det A = 0 (801² grid, vertices slid onto the curve by Newton), caustic = β(curve) |
| magnification map | signed μ = 1/det A on the pixel grid of the image, FITS with the image WCS (`lens_magnification.fits`, clipped at ±`max-mu`) |
| image prediction | full lens-equation solution (grid triangle mapping + Newton) at the fitted source; **counter-images** = predicted images that match none of the given ones (listed with position, μ, parity) |
| time delays | Fermat potential, analytic ψ for SIE and shear, tabulated for NFW/PIEMD; days, relative to the first-arriving image, for the given z_l, z_s |
| source-plane reconstruction | every pixel of a cutout is ray-traced to the source plane (surface brightness conserved) and averaged into source pixels (`lens_source.fits`); the reconstructed source is lensed back (`lens_model_image.fits`) and subtracted (`lens_residual.fits`) |

Conventions: the model lives in the **pixel frame** (x right, y up, arcsec = pixel × scale; scale from the WCS or `pixscale`).  φ is the position angle of the **major axis of the mass**, degrees counter-clockwise from +x; the sky PA (east of north) is added when the header has a WCS.  θ_E is the area-equivalent Einstein radius of the SIE (κ = θ_E / (2 √(q x'² + y'²/q))).  Source/predicted positions are reported in arcsec and in 1-based pixels.

## GUI

Chip **Lens** (Measure tab): ▶ = *Fit Lens Model*.  Gear = parameter dialog (groups Images, Redshifts, Model, Perturbers, Outputs).  Menu entries: *Use Selected Rows as Images*, *Use Selected Row as Lens* (fill the parameters from the catalog selection), *Fit Lens Model* (adds columns `LENS_ROLE LENS_MU LENS_RES LENS_DT LENS_PARITY` to the image rows and draws the overlay: critical curve red, caustic cyan dashed, predicted images green (matched) / yellow (counter-image), source magenta; a previous overlay is replaced, other regions are untouched), *Critical Curves and Caustics*, *Predict Images / Delays*, *Magnification Map* and *Source-plane Reconstruction* (open in new frames), *Lens Model Results...* (notebook: summary with errors/physical quantities/warnings, table of predicted images with μ, parity, delay, residual, table of critical curves).  Every step is a CLI template, so the session recorder logs it (`analysis.lens_*`, class AUTO) and the exported script replays it; `scripts/verify_lensmodel.sh` checks that replay reproduces the catalog after each of the 7 steps.

## CLI

```
python3 plugins/lensmodel/lensmodel.py --task simulate --work DIR --seed 3          # synthetic quad: sim_lens.fits, sim_catalog.tsv, sim_truth.json
python3 plugins/lensmodel/lensmodel.py --task fit --work DIR --image sim_lens.fits --catalog sim_catalog.tsv --image-numbers 2,3,4,5 --lens-number 1 \
        --sigma-pos 0.005 --z-lens 0.5 --z-source 2.0 [--flux-column FLUX] [--member-column CL_MEMBER ...] [--nfw-kappa-s 0.2]
python3 plugins/lensmodel/lensmodel.py --task curves|magmap|source|predict --work DIR --image sim_lens.fits [...]
```

Positions can be given instead of catalog numbers (`--positions "x,y;x,y;…" --lens-x --lens-y`, 1-based pixels).  `ogfkit.lensmodel` can be used directly (`LensModel`, `fit_sie_shear`, `predict`, `time_delays`, `ray_trace_image`).

## Validation against truth (`plugins/lensmodel/validation/validation_report.json`, `run_validation.py`)

All numbers are for SIE + shear lenses with θ_E ∈ [0.9, 2.2]″, q ∈ [0.55, 0.95], γ ∈ [0, 0.10], random PAs, a point source inside the caustic (quads) and the lens centre fixed at the truth; images carry Gaussian position noise.  A quad gives 8 constraints for 7 free parameters (θ_E, q, φ, γ, φ_γ, β_x, β_y), i.e. **dof = 1**.

| test | result |
|---|---|
| noise-free recovery, 20 quads | θ_E, q, γ, φ, source recovered to < 4e-11 (relative), 3e-10, 8e-11, 1e-8 deg, 1e-10″ - the optimiser reaches the exact solution |
| noise 3 mas, 32 quads (+ 8 doubles) | Δθ_E/θ_E: median +6e-5, σ 9.3e-4, max 2.6e-3; Δq: σ 0.011 (max 0.033); ΔPA (q<0.9): σ 0.95° (max 3.0°); Δγ: σ 0.0036 (max 0.0115); source position error median 1.4 mas (max 4.8 mas); position rms median 0.55 mas (σ_pos 3 mas; expectation for 8 data, 7 parameters ≈ σ·√(1/8) = 1.1 mas rms per coordinate-pair mean); χ² mean 0.83 for dof 1 |
| noise 10 mas, 25 quads | Δθ_E/θ_E σ 4.8e-3 (max 1.4e-2); Δq σ 0.041 (max 0.116); ΔPA σ 3.6° (max 7.4°); Δγ σ 0.016; source 6 mas; rms median 3.1 mas; χ² mean 0.89 |
| doubles (θ_E fitted, q = 1, γ = 0 fixed) | Δθ_E/θ_E median +1.1 %, σ 5-6 % (max 10 %): the unmodelled ellipticity/shear is absorbed - doubles constrain θ_E to a few per cent only |
| predicted images (quads, fitted model) | 4 images in 32 / 32 |
| time delays (fitted vs true model, quads) | max |Δt| / (delay span): median 1.0 %, 84th percentile 2.7 %, **worst 38 %** (near-degenerate configurations with 3 mas noise) |
| counter-image prediction | unit test: dropping the 4th image of a quad predicts it back to < 1e-6″ |
| unmodelled perturber (PIEMD b0 = 0.12-0.15″ within 1″ of an image; 14 lenses) | model **with** the perturber: rms 1.0 mas, χ² median 0.43, Δθ_E/θ_E 5e-4;  **without**: rms 14.6 mas, χ² median 97, Δθ_E/θ_E median 4.0 % (max 19.5 %), Δq max 0.41 - χ² flags the missing mass, and ignoring it biases the Einstein radius |
| source-plane reconstruction (8 random lenses, Sérsic source r_eff 0.1″, noise 0.01) | true model: correlation with the analytic source 0.996, flux ratio 1.006, centroid error 16 mas, forward-model residual 7.5 × noise (binning/interpolation floor); θ_E wrong by 3 %: correlation 0.86, residual 37 × noise; wrong by 10 %: correlation -0.13, residual 55 × noise - the reconstruction is a sensitive test of the model |
| independent truth: **lenstronomy 1.14.2** | deflection of the SIE (50 random lenses × 30 points): max difference **9.8e-10″**; NFW: 4e-14″ (κ_s = 0.4, θ_s = 8″); the four images of 15 lenses found by lenstronomy's solver and by ours differ by < 9e-10″ (same count 4/4).  Images generated by lenstronomy, noise 3 mas, fitted by ogfkit (30 quads, `validation/crosscheck_lenstronomy.txt`): Δθ_E/θ_E σ 2.3e-3 (max 8e-3), Δq σ 0.016, ΔPA σ 0.64°, Δγ σ 0.0034, χ² median 0.44 |
| analytic checks (tests) | SIE κ vs analytic (1e-6), SIS images at β ± θ_E and μ = 1/(1 ∓ θ_E/r), SIS delay Δφ = 2θ_E β, critical-curve area of a lone SIE = πθ_E² (2e-3), ∇ψ = α for every component, NFW/PIEMD κ from deflection (2e-3) |

## Limitations and what is not done

* **No real lens was modelled.**  The fits above use synthetic lenses (the generator shares the engine with the fitter, which is why the lenstronomy-generated test and the deflection cross-check were added); real systems have lens-light residuals, position errors larger than the formal ones, substructure, line-of-sight mass and mass-sheet/ellipticity degeneracies that these tests do not contain.
* Only SIE + shear (centre optional) are **fitted**.  NFW and PIEMD members are fixed inputs (their strengths come from the scaling relation, not from the fit); PIEMD members are circular; no core radius in the SIE; no multipole/satellite fitting; no multiple sources at different redshifts (cluster-scale models with several image systems are out of scope).
* A quad has dof 1, a double has none unless parameters are fixed: χ² = 0 is not evidence for the model (a warning is printed when dof ≤ 0).  The errors are from the Jacobian (rescaled when χ²/dof > 1), not MCMC; time-delay uncertainties are not propagated (worst-case 38 % above).
* The model is in the pixel frame; WCS rotation/parity enters only through the sky-PA conversion.
* The image finder grid (401²) can miss a very close or very demagnified image (e.g. a central image of a cored/NFW lens); raise `grid`.  Extended arcs are not fitted in the image plane - only image positions (and optionally flux ratios); the source-plane reconstruction is a binning (no regularised inversion, no PSF), meant for model diagnostics.
* Distances use flat ΛCDM from astropy (no H0 tension handling); the mass is the mass inside the critical curve of the model, not a total cluster mass.
* GUI verification: `scripts/verify_lensmodel.sh` (37 checks on a synthetic lens in real ds9 + session replay), run on Xvfb only.
