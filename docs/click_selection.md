# Click selection with overlapping markers

File: `ds9/library/ogf_pick.tcl` (hooks in `frame.tcl` `Button1Frame` and `ogf_tile.tcl` `OGFTileClick`).

A left click on the image (pointer mode `none`) first calls `OGFPickClick`:

| hits within the radius | behaviour |
|---|---|
| 0 | returns 0; the old behaviour (marker click / link press / ds9 default) runs unchanged |
| 1 | the nearest row is selected directly (no popup) |
| >= 2 | a small chooser `.ogfpick` lists the candidates (columns **ID**, **kind**, **mag**, **d (px)**), nearest first, at most 12 lines (then "... n more"); a click or Enter on a line selects that table row |

Dismiss: **Esc**, or a click anywhere outside the chooser (no grab: that click keeps its normal meaning, only the
chooser closes).  Dismissing never changes the selection.

Candidates: moving, transient and detection rows (also `_pos` track points) are matched by the sky distance from the click
within `ogfpick(radius)` = 8 canvas pixels; galaxies use the old "inside the ellipse" rule (`OGFHitAll`, semi-axes floored at 2 px).
Only rows currently present in the table are offered.  Works for every kind in the single view and in tile mode (the clicked
tile's frame is made current first, then the row is selected via `CatalogPanelLinkSelect`).

## Test

    OGF_CHK_OUT=/tmp/chk.txt DISPLAY=:77 bin/ds9 /workspace/fits/m51.fits -geometry 1300x950 -source scripts/verify_click_chooser.tcl

Synthetic overlapping markers (3 moving rows + 1 transient within 3 px) on a real m51 extraction (181 galaxies).  Last run:
30 checks PASS, 0 failures.  Measured: candidates sorted by distance 0.00 / 1.30 / 2.24 / 2.92 px; chooser has 4 lines; no
selection before picking; picking line 1 selects M1 (table row 184); Esc and click-away dismiss without changing the selection;
a single hit selects directly with the nearest row; two overlapping galaxies -> chooser, pick selects galaxy #5; tile mode (2 frames):
both tiles find 4 candidates and `OGFTileClick` on the non-current tile opens the chooser and selects M0.  Geometry before, with
the chooser open, in tile view and after: table y=181, h=769, info area 154, main window 1300x950.

## Limitations

* Only mouse button 1 in `none` mode is handled.
* The test calls the click procs with synthetic canvas coordinates; real X mouse events were **not** generated.
* Galaxies in the "all" view are offered only if they are on the detection grid.
* The hit radius for time-domain rows is fixed (8 px), it does not scale with the zoom.
