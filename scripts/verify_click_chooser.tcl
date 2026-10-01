# Measured test of the overlapping-marker click chooser (ds9/library/ogf_pick.tcl).
#   rm -rf ~/ds9.auto ~/ds9.auto.dir; DISPLAY=:77 OGF_CHK_OUT=/tmp/chk.txt bin/ds9 /workspace/fits/m51.fits -geometry 1300x950 -source scripts/verify_click_chooser.tcl
# Synthetic moving / transient rows are placed on top of each other around image points (no pipeline needed); the galaxy catalog
# comes from a real extraction of the image.  Clicks go through the proc Button1Frame calls (OGFPickClick / OGFTileClick).
global current catpanel ds9 ogftd ogfpick
set ::out $::env(OGF_CHK_OUT)
set ::fh [open $::out w]
set ::nfail 0
set ::T0 [clock milliseconds]
rename tk_messageBox ::orig_mb
proc tk_messageBox {args} {puts $::fh "MESSAGEBOX $args"; flush $::fh; return ok}
proc bgerror {m} {puts $::fh "BGERROR $m $::errorInfo"; flush $::fh}
proc R {tag ok {detail {}}} {
    puts $::fh "[expr {$ok ? {PASS} : {FAIL}}] $tag $detail  @[expr {[clock milliseconds]-$::T0}]ms"; flush $::fh
    if {!$ok} {incr ::nfail}
}
proc geom {tag} {
    update idletasks; update
    set tf [winfo parent $::catpanel(tbl)]
    set tl [winfo toplevel $tf]
    return "tbl_y=[winfo rooty $tf] tbl_h=[winfo height $tf] info_h=[winfo height $::catpanel(infoarea)] main=[winfo width .]x[winfo height .]"
}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}

# canvas <-> image calibration of a frame (linear at the current zoom/pan).  In tile view every frame occupies a rectangle of the
# one ds9 canvas: its origin is $ds9(canvas) coords $fr (canvas coordinates of the click are canvas-wide).
proc origin {fr} {
    global ds9
    if {[catch {$ds9(canvas) coords $fr} o]} {return {0 0}}
    return [list [expr {int(ceil([lindex $o 0]))}] [expr {int(ceil([lindex $o 1]))}]]
}
proc calib {fr} {
    global ds9
    lassign [origin $fr] ox oy
    # inside the image: the centre of the frame's rectangle (canvas-wide coordinates)
    set c [$ds9(canvas) bbox $fr]
    set x0 [expr {([lindex $c 0]+[lindex $c 2])/2 - 30}]; set y0 [expr {([lindex $c 1]+[lindex $c 3])/2 - 30}]
    lassign [$fr get coordinates $x0 $y0 image] ax ay
    lassign [$fr get coordinates [expr {$x0+60}] $y0 image] bx by
    lassign [$fr get coordinates $x0 [expr {$y0+60}] image] cx cy
    set sx [expr {($bx-$ax)/60.0}]; set sy [expr {($cy-$ay)/60.0}]
    return [list $ax $ay $sx $sy $x0 $y0]
}
proc img2canvas {fr ix iy} {
    lassign [calib $fr] ax ay sx sy x0 y0
    return [list [expr {$x0 + ($ix-$ax)/$sx}] [expr {$y0 + ($iy-$ay)/$sy}]]
}
proc img2sky {fr ix iy} {
    lassign [img2canvas $fr $ix $iy] cx cy
    return [$fr get coordinates $cx $cy wcs fk5 degrees]
}
proc selnum {} {global catpanel; return [lindex $catpanel(sel,nums) end]}

proc run {} {
    global current catpanel ds9 ogftd ogfpick
    set fr $current(frame)
    update; wait_idle 500
    R geom_start 1 [geom start]
    set g0 [geom start]
    # real galaxy catalog
    set catpanel(param,detect-thresh) 2.0
    CatalogPanelExtract
    set t0 [clock milliseconds]
    while {$catpanel(alldata) eq {} && [clock milliseconds]-$t0 < 120000} {update; after 100}
    set ngal [expr {[llength [split $catpanel(alldata) \n]]-2}]
    R extract [expr {$ngal > 50}] "galaxies=$ngal"
    ZoomToFit; wait_idle 300
    # isolated test point: pick 2 image points far from galaxies; use fixed positions near the image centre
    lassign [$fr get coordinates [expr {[winfo width $ds9(canvas)]/2}] [expr {[winfo height $ds9(canvas)]/2}] image] icx icy
    set zoom [lindex [$fr get zoom] 0]
    R zoom [expr {$zoom > 0}] "zoom=$zoom"
    # ---- synthetic rows: A (3 overlapping moving + 1 transient within 3 px), B (1 isolated moving, 40 px away)
    set ax [expr {$icx + 150}]; set ay [expr {$icy + 120}]
    set bx [expr {$icx - 200}]; set by [expr {$icy - 160}]
    set mrows {}
    set k 0
    foreach {dx dy mag} {0.0 0.0 22.1  1.2 0.5 23.4  -2.0 1.0 24.0} {
	lassign [img2sky $fr [expr {$ax+$dx}] [expr {$ay+$dy}]] ra de
	lappend mrows [list $k $ra $de $mag 4 5.0 90 0.1 3.0 candidate]
	incr k
    }
    lassign [img2sky $fr $bx $by] ra de
    lappend mrows [list $k $ra $de 21.7 4 5.0 90 0.1 3.0 candidate]
    ::ogf::td::set_rows moving {id ra dec mag n rate_arcsec_h pa_deg rms_arcsec score status} $mrows
    lassign [img2sky $fr [expr {$ax+2.5}] [expr {$ay-1.5}]] ra de
    ::ogf::td::set_rows transient {id ra dec mag n_det files snr_max host_id host_sep host_z offset_re heuristic} [list [list 7 $ra $de 24.5 2 {x} 9.0 {} {} {} {} {}]]
    ::ogf::td::show all
    wait_idle 800
    R rows_all [expr {[llength [::ogf::td::visible_kinds]] >= 2}] "kinds=[::ogf::td::visible_kinds]"
    set nrows [expr {[$catpanel(tbl) cget -rows]-1}]
    R table_rows [expr {$nrows >= $ngal + 5}] "table rows=$nrows"

    # ---- 1. click on A: several hits -> chooser with 4 time-domain rows (+ any galaxy under it)
    lassign [img2canvas $fr $ax $ay] cx cy
    set cands [OGFPickCandidates $fr $cx $cy]
    set kinds [lmap c $cands {dict get $c kind}]
    set keys [lmap c $cands {dict get $c key}]
    R cand_count_A [expr {[llength [lsearch -all -glob $keys M*]] == 3 && [lsearch $keys T7] >= 0}] "cands=[llength $cands] keys=$keys kinds=$kinds"
    set d [lmap c $cands {format %.2f [dict get $c dist]}]
    R sorted_by_distance [expr {$d eq [lsort -real $d]}] "dist(canvas px)=$d"
    CatalogPanelClearSelection
    set ok [OGFPickClick $fr $cx $cy 400 300]
    wait_idle 300
    R chooser_open [expr {$ok == 1 && [winfo exists .ogfpick]}] "ret=$ok"
    if {[winfo exists .ogfpick]} {
	set lb .ogfpick.l
	R chooser_lines [expr {[$lb size] == [llength $cands]}] "lines=[$lb size] header='[.ogfpick.h cget -text]'"
	R chooser_text1 1 "first line: [$lb get 0]"
	R chooser_no_selection_yet [expr {[llength $catpanel(sel,nums)] == 0}]
	R geom_chooser_open [expr {[geom x] eq $g0}] [geom x]
	# pick the 2nd entry via the same proc a click/Enter on the line calls
	set want [dict get [lindex $cands 1] key]
	OGFPickChoose 1
	wait_idle 300
	R chooser_closed_after_pick [expr {![winfo exists .ogfpick]}]
	R picked_row_selected [expr {[selnum] eq $want}] "want=$want got=[selnum] sel_row=[OGFRowOfNumber [selnum]]"
    }
    # ---- 2. Esc dismisses without changing the selection
    CatalogPanelClearSelection
    OGFPickClick $fr $cx $cy 400 300; wait_idle 200
    R esc_open [winfo exists .ogfpick]
    event generate .ogfpick.l <Escape>; wait_idle 200
    R esc_dismisses [expr {![winfo exists .ogfpick] && [llength $catpanel(sel,nums)] == 0}] "sel=$catpanel(sel,nums)"
    # ---- 3. click-away dismisses
    OGFPickClick $fr $cx $cy 400 300; wait_idle 200
    set rx [expr {[winfo rootx .ogfpick]+[winfo width .ogfpick]+60}]; set ry [expr {[winfo rooty .ogfpick]+5}]
    OGFPickOutside . $rx $ry; wait_idle 300
    R clickaway_dismisses [expr {![winfo exists .ogfpick] && [llength $catpanel(sel,nums)] == 0}]
    # ---- 4. single hit selects the nearest row directly (isolated B, 0 px offset and 3 px offset)
    lassign [img2canvas $fr [expr {$bx+1.0}] [expr {$by+1.0}]] sx sy
    CatalogPanelClearSelection
    set c [OGFPickCandidates $fr $sx $sy]
    set ok [OGFPickClick $fr $sx $sy 400 300]
    wait_idle 300
    R single_hit_direct [expr {$ok == 1 && ![winfo exists .ogfpick] && [selnum] eq {M3}}] "cands=[llength $c] sel=[selnum]"
    # ---- 5. nearest: click 0.5 px from the closest of A's rows: chooser lists it first
    lassign [img2canvas $fr [expr {$ax+1.2}] [expr {$ay+0.5}]] nx ny
    set c [OGFPickCandidates $fr $nx $ny]
    R nearest_first [expr {[dict get [lindex $c 0] key] eq {M1}}] "order=[lmap x $c {dict get $x key}] d0=[format %.2f [dict get [lindex $c 0] dist]]"
    # ---- 6. galaxy rows: Galaxies view lists overlapping ellipses with kind galaxy; find a click with >= 2 galaxies if any
    ::ogf::td::show galaxies; wait_idle 500
    set best {}; set bn 0
    OGFCacheBuild
    set ix 0
    foreach n $::catcache(num) x $::catcache(x) y $::catcache(y) {
	lassign [img2canvas $fr $x $y] gx gy
	set cc [OGFPickCandidates $fr $gx $gy]
	if {[llength $cc] > $bn} {set bn [llength $cc]; set best [list $gx $gy $cc]}
	if {[incr ix] > 400} break
    }
    R galaxy_candidates [expr {$bn >= 1}] "max galaxies under a click (first 400 rows)=$bn kinds=[lmap c [lindex $best 2] {dict get $c kind}]"
    if {$bn >= 1} {
	lassign $best gx gy cc
	CatalogPanelClearSelection
	OGFPickClick $fr $gx $gy 400 300; wait_idle 300
	if {$bn == 1} {R galaxy_single_direct [expr {![winfo exists .ogfpick] && [selnum] ne {}}] "sel=[selnum]"} \
	else {R galaxy_chooser [winfo exists .ogfpick] "n=$bn"; OGFPickChoose 0; wait_idle 200
	      R galaxy_pick_selects [expr {[selnum] eq [dict get [lindex $cc 0] key]}] "sel=[selnum]"}
    }
    OGFPickDismiss
    # ---- 7. tile mode: two frames, markers drawn in both; click in the NON-current tile
    ::ogf::td::show all; wait_idle 500
    CreateFrame; LoadFitsFile /workspace/fits/m51.fits {} {}; wait_idle 800
    OGFUIDisplay tile; wait_idle 1000
    R tile_on [OGFTileIsOn] "frames=$ds9(frames) current=$current(frame)"
    set g1 [geom x]
    R geom_tile [expr {[geom x] eq $g1}] $g1
    set cur $current(frame)
    set other [lindex [lsearch -all -inline -not $ds9(frames) $cur] 0]
    # on each tile, find the canvas point of A's centre
    foreach f $ds9(frames) {
	if {[catch {img2canvas $f $ax $ay} p]} {R tile_calib_$f 0 $p; continue}
	lassign $p tx ty
	set c [OGFPickCandidates $f $tx $ty]
	R tile_cands_$f [expr {[llength [lsearch -all -glob [lmap x $c {dict get $x key}] M*]] >= 1}] "frame=$f canvas=([format %.0f $tx],[format %.0f $ty]) n=[llength $c] keys=[lmap x $c {dict get $x key}]"
    }
    # tile click through OGFTileClick (what Button1Frame calls for a non-current tile)
    GotoFrame $cur; wait_idle 300
    if {$other ne {}} {
	lassign [img2canvas $other $ax $ay] tx ty
	CatalogPanelClearSelection
	set r [OGFTileClick $other $tx $ty]
	wait_idle 300
	set open [winfo exists .ogfpick]
	R tile_click_chooser [expr {$r == 1 && $open}] "ret=$r chooser=$open current=$current(frame)"
	if {$open} {OGFPickChoose 0; wait_idle 300; R tile_pick_selects [expr {[selnum] ne {}}] "sel=[selnum]"}
	lassign [img2canvas $other $bx $by] tx ty
	CatalogPanelClearSelection
	set r [OGFTileClick $other $tx $ty]; wait_idle 300
	R tile_single_direct [expr {$r == 1 && ![winfo exists .ogfpick] && [selnum] eq {M3}}] "sel=[selnum]"
    }
    OGFPickDismiss
    OGFUIDisplay single; wait_idle 800
    R geom_end [expr {[geom x] eq $g0}] "start: $g0 | end: [geom x]"
    puts $::fh "SUMMARY failures=$::nfail"
    close $::fh
    exit
}
proc geom {x} {
    update idletasks; update
    set tf [winfo parent $::catpanel(tbl)]
    return "tbl_y=[winfo rooty $tf] tbl_h=[winfo height $tf] info_h=[winfo height $::catpanel(infoarea)] main=[winfo width .]x[winfo height .]"
}
after 4000 {if {[catch run err]} {puts $::fh "ERROR $err\n$::errorInfo"; close $::fh; exit 1}}
