# Multi-band, image-table link and shared mask (OGFinder)

New Tcl modules (all in `ds9/library/`, embedded into `bin/ds9` by `make`):

| file | purpose |
|---|---|
| `ogf_util.tcl`  | `OGFPython` (python interpreter, override with `OGFINDER_PYTHON`), small modal forms |
| `ogf_bands.tcl` | band registry, forced photometry driver, colour columns, tile view |
| `ogf_link.tcl`  | image <-> table selection, hover readout, multi-select, key stepping |
| `ogf_mask.tcl`  | shared mask manager + ICL/LSBG integration |
| `ds9_mask.py`   | mask back end (numpy/scipy/astropy/sep) |

The catalog-panel menubar now has ten menus; `CatMenu.TMenubutton` padding was
reduced from 10 to 2 px so they fit in the panel (menubar needs 565 px; panel is 561-585 px wide).

## Python environment

`ds9_mask.py`, `ds9_icl.py`, `ds9_lsbg.py` need astropy and sep. The GUI now runs
`$OGFINDER_PYTHON` if set, otherwise `python3` (exactly as before). Example:

    OGFINDER_PYTHON=/workspace/ogf_venv/bin/python3 bin/ds9 image.fits

## Bands menu

| entry | action |
|---|---|
| Register Current Frame as Band... | name (default = FILTER keyword), AB zeropoint (default from PHOTFLAM/PHOTPLAM via `ds9_sextract --info`), optional PSF FWHM; pixel scale read from the WCS |
| Load Band... | open a FITS file in a new frame and register it |
| Detection Band > | choose the band in which sources are detected |
| Remove Band > | unregister |
| List Bands... | table of bands, ZP, pixel scale, grid (same / WCS) |
| Detect in Detection Band | goes to that band's frame, sets the ZP, runs the normal Extract |
| Measure in All Bands... | forced photometry in every band (asks the S/N below which MAG=99) |
| Tile Bands / Single Frame View | tile display with `panzoom(lock)=wcs`, `crosshair(lock)=wcs`, zscale per frame |
| Copy Mask to Other Bands | reproject the detection-band mask to every other band (`ds9_mask.py --mode reproject`) |

Catalog columns added (merged by NUMBER with `CatalogPanelAddColumnsFromTSV`, so
sort/filter/save work): `MAG_<band>`, `MAGERR_<band>` for every band (including the
detection band, measured through the same forced path) and adjacent-band colours
`<bandA>-<bandB>` in pivot-wavelength order. Non-detections are `99` (MAG and MAGERR);
a colour is `99.000` if either magnitude is 99. Registered band frames share one
catalog (switching frames does not swap/blank it).

### ds9_sextract command line (C, `ds9/library/ds9_sextract.c`)

    ds9_sextract --info FILE
        key=value: NAXIS1/2 FILTER PUPIL INSTRUME BUNIT ZP_AB PIVOT EXPTIME PIXSCALE CRPIX* CRVAL* CD* WCSVALID
    ds9_sextract DET.fits --forced-catalog CAT.tsv --measure-image B.fits --band NAME
                 --mag-zeropoint ZP [--snr-min 1.0] [--back-size N] [--phot-aperture D]
        -> NUMBER X_<b> Y_<b> SCALE_<b> FLUX_AUTO_<b> FLUXERR_AUTO_<b> MAG_AUTO_<b> MAGERR_AUTO_<b>
           FLUX_APER_<b> FLUXERR_APER_<b> MAG_APER_<b> MAGERR_APER_<b> FLAGS_<b>

Per-band sep background/RMS and no-data mask (0/NaN). Catalog positions and Kron ellipses
(2.5 x KRON_RADIUS x A, B, THETA) are taken from the detection catalog; if the two WCSs
are identical the mapping is the identity, otherwise pixel -> RA/Dec -> pixel of the band
with the local Jacobian applied to the ellipse. FLAGS 0x1000 = outside image, 0x100 = no shape.

## Image <-> table link

* Click a source in the image -> row selected, scrolled into view, summary block filled,
  cyan cross + green ellipse drawn (in every registered band frame). Works **without** Mark All
  (hit test against a cached parse of the catalog; smallest ellipse wins).
* Click empty sky (no drag) or press **Esc** -> selection, markers and summary cleared.
* Hover over a source -> `Hover: #ID  MAG_AUTO = x` line in the summary area (throttled to
  ~15 Hz, ~0.08 ms per evaluation on 207 rows).
* **Shift+click** on a source toggles it in a multi-selection (extra members: cyan crosses,
  blue rows in the table, summary "N sources selected").
* **Up/Down/n/p** in the image (only while a source is selected, otherwise Up/Down keep
  warping the cursor) and **n/p** in the table step through the rows and centre the image.
* Ctrl+click merge selection, PSF-star removal, AI-merge key bindings are unchanged.

## Mask menu

One bit-flag mask per image: `~/.ds9/mask_<base>.fits` (uint8; 1 source, 2 bright star,
4 manual add, 8 manual erase/protect (wins), 16 imported). `mask_<base>_bool.fits` is the
effective 0/1 mask that ICL/LSBG read through `--mask`. Undo/redo stacks (10 deep) live in
`mask_<base>.undo/` and `.redo/`. The interpolated image `mask_<base>_masked.fits` is only
written on demand (Show Masked Image, ICL Background, ...).

| entry | action |
|---|---|
| Auto Mask... | shared parameters (detect threshold, min area, expansion, max dilate radius, bright-star mag limit & radius scale, mag threshold, LSB protect); buttons "ICL defaults"/"LSBG defaults" load the per-method values. Manual add/erase/import bits are preserved on re-run |
| Show Mask Overlay | DS9 native mask layer (`mask color/transparency/mark nonzero/blend screen`), no new frame |
| Overlay Colour/Transparency... | |
| Add / Erase Regions to/from Mask | rasterises the DS9 circle/ellipse/box/polygon/annulus regions of the current frame (image coordinates, `-` excluded shapes honoured) |
| Grow... / Shrink... | euclidean distance transform by N px |
| Invert / Clear Mask | |
| Undo / Redo | |
| Save As... / Import... | flag FITS or plain 0/1 FITS; import of a foreign mask goes through `reproject_mask` (WCS) |
| Show Masked Image | interpolated image in a new frame |
| Statistics | masked pixels, fraction (of image and of valid pixels), undo/redo depth in the status bar |

ICL "1. Source Masking (shared mask)", "Show Mask Overlay", "Save Mask As...", "Import Mask..."
and LSBG "1. Mask Bright Sources (shared mask)", "Show Masked Image", "Save/Import" now call
the same code. Existing shared masks are reused (manual edits kept). ICL Background feeds the
interpolated shared-mask image and `--mask mask_<base>_bool.fits`; with iterative refinement the
refined mask goes to `mask_<base>_iclrefined.fits` (the shared mask is not overwritten).
LSBG Clean/Run Full Pipeline obtain the shared mask (`OGFMaskEnsure`); `ds9_lsbg.py --mode run`
has a new `--mask-input FILE` option to reuse it instead of regenerating.

### ds9_mask.py

    ds9_mask.py IMAGE --mode auto    --mask M [--detect-thresh 5 --minarea 5 --expand-factor 1.5
                   --max-dilate-radius 20 --bright-star-mag-limit 18 --bright-star-radius-scale 10
                   --mag-threshold 99 --lsb-protect --mag-zeropoint 25 --catalog cat.tsv --fresh]
    ds9_mask.py IMAGE --mode add|erase --mask M --regions file.reg
    ds9_mask.py IMAGE --mode grow|shrink --mask M --pixels N
    ds9_mask.py IMAGE --mode invert|clear|undo|redo|stats --mask M
    ds9_mask.py IMAGE --mode import|export --mask M --file F [--boolean]
    ds9_mask.py IMAGE --mode reproject --mask M --target-image BAND.fits --output BANDMASK.fits
    ds9_mask.py IMAGE --mode masked --mask M --masked-output OUT.fits [--interp-method linear]

Every mode prints one `#MASK_STATS N_MASKED=.. FRACTION=.. N_SRC=.. N_STAR=.. N_ADD=.. N_ERASE=.. UNDO=.. REDO=..` line.

Mask and band steps are also recorded by the session recorder; see `docs/session_python_script.md` (hand-drawn add/erase regions are *manual* steps and are skipped in pipeline mode).
