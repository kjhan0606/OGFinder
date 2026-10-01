# REAL X-event test of the click chooser / tile click: the mouse is moved and clicked with xdotool (XTEST), so the events
# travel X server -> Tk -> canvas binding -> Button1Frame -> OGFPickClick / OGFTileClick.  Needs xdotool and a DISPLAY.
#   rm -rf ~/ds9.auto ~/ds9.auto.dir; DISPLAY=:77 OGF_CHK_OUT=/tmp/chkx.txt bin/ds9 /workspace/fits/m51.fits -geometry 1300x950 -source scripts/verify_click_xevent.tcl
# (scripts/verify_click_chooser.tcl calls the procs directly; this one does not.)
# Original header:
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


# canvas point -> root screen point -> real click
proc xclick {cx cy {btn 1}} {
    global ds9
    set c $ds9(canvas)
    set rx [expr {int([winfo rootx $c] + $cx - [$c canvasx 0])}]
    set ry [expr {int([winfo rooty $c] + $cy - [$c canvasy 0])}]
    exec xdotool mousemove --sync $rx $ry
    wait_idle 150
    exec xdotool click $btn
    wait_idle 400
    return [list $rx $ry]
}
proc xkey {key} { exec xdotool key $key; wait_idle 300 }
proc run {} {
    global current catpanel ds9 ogftd ogfpick
    set fr $current(frame)
    update; wait_idle 500
    if {[catch {exec xdotool version} v]} {R xdotool 0 $v; puts $::fh "SUMMARY failures=$::nfail"; close $::fh; exit 1}
    R xdotool 1 $v
    set g0 [geom start]
    set catpanel(param,detect-thresh) 2.0
    CatalogPanelExtract
    set t0 [clock milliseconds]
    while {$catpanel(alldata) eq {} && [clock milliseconds]-$t0 < 120000} {update; after 100}
    ZoomToFit; wait_idle 300
    lassign [$fr get coordinates [expr {[winfo width $ds9(canvas)]/2}] [expr {[winfo height $ds9(canvas)]/2}] image] icx icy
    set ax [expr {$icx + 150}]; set ay [expr {$icy + 120}]
    set bx [expr {$icx - 200}]; set by [expr {$icy - 160}]
    set mrows {}; set k 0
    foreach {dx dy mag} {0.0 0.0 22.1  1.2 0.5 23.4  -2.0 1.0 24.0} {
	lassign [img2sky $fr [expr {$ax+$dx}] [expr {$ay+$dy}]] ra de
	lappend mrows [list $k $ra $de $mag 4 5.0 90 0.1 3.0 candidate]; incr k
    }
    lassign [img2sky $fr $bx $by] ra de
    lappend mrows [list $k $ra $de 21.7 4 5.0 90 0.1 3.0 candidate]
    ::ogf::td::set_rows moving {id ra dec mag n rate_arcsec_h pa_deg rms_arcsec score status} $mrows
    lassign [img2sky $fr [expr {$ax+2.5}] [expr {$ay-1.5}]] ra de
    ::ogf::td::set_rows transient {id ra dec mag n_det files snr_max host_id host_sep host_z offset_re heuristic} [list [list 7 $ra $de 24.5 2 {x} 9.0 {} {} {} {} {}]]
    ::ogf::td::show all; wait_idle 800
    R pointer_mode [expr {[info exists ::current(mode)] && $::current(mode) eq {none}}] "mode=[expr {[info exists ::current(mode)]?$::current(mode):{?}}]"
    # ---- 1. real click on the overlapping cluster A
    lassign [img2canvas $fr $ax $ay] cx cy
    CatalogPanelClearSelection
    set rp [xclick $cx $cy]
    R x_chooser_open [winfo exists .ogfpick] "root=$rp last=[llength $::ogfpick(last)] keys=[lmap c $::ogfpick(last) {dict get $c key}] cands_at_target=[llength [OGFPickCandidates $fr $cx $cy]]"
    if {[winfo exists .ogfpick]} {
	R x_chooser_lines [expr {[.ogfpick.l size] == 4}] "lines=[.ogfpick.l size]"
	R x_no_selection_yet [expr {[llength $catpanel(sel,nums)] == 0}]
	R x_geom_chooser [expr {[geom x] eq $g0}] [geom x]
	# real click on the 2nd line of the listbox (index 1)
	set bb [.ogfpick.l bbox 1]
	set lx [expr {[winfo rootx .ogfpick.l] + [lindex $bb 0] + 20}]
	set ly [expr {[winfo rooty .ogfpick.l] + [lindex $bb 1] + [lindex $bb 3]/2}]
	set want [lindex [.ogfpick.l get 1] 0]
	exec xdotool mousemove --sync $lx $ly; wait_idle 150
	exec xdotool click 1; wait_idle 500
	R x_line_click_selects [expr {![winfo exists .ogfpick] && [selnum] ne {}}] "line='$want' sel=[selnum]"
    }
    # ---- 2. Esc key (real key event) dismisses
    CatalogPanelClearSelection
    xclick $cx $cy
    R x_esc_open [winfo exists .ogfpick]
    if {[winfo exists .ogfpick]} {
	focus -force .ogfpick.l; wait_idle 200
	xkey Escape
	R x_esc_dismisses [expr {![winfo exists .ogfpick] && [llength $catpanel(sel,nums)] == 0}] "sel=$catpanel(sel,nums)"
    }
    # ---- 3. click away (real click on the empty part of the canvas) dismisses
    xclick $cx $cy
    set o [winfo exists .ogfpick]
    lassign [img2canvas $fr [expr {$icx-60}] [expr {$icy+300}]] ox oy
    CatalogPanelClearSelection
    xclick $ox $oy
    R x_clickaway [expr {$o && ![winfo exists .ogfpick]}] "was_open=$o now=[winfo exists .ogfpick]"
    # ---- 4. isolated marker: real click -> direct selection
    lassign [img2canvas $fr [expr {$bx+1.0}] [expr {$by+1.0}]] sx sy
    CatalogPanelClearSelection
    xclick $sx $sy
    R x_single_direct [expr {![winfo exists .ogfpick] && [selnum] eq {M3}}] "sel=[selnum]"
    # ---- 5. tile mode, click in the non-current tile
    CatalogPanelClearSelection
    ::ogf::td::show all; wait_idle 500
    CreateFrame; LoadFitsFile /workspace/fits/m51.fits {} {}; wait_idle 800
    OGFUIDisplay tile; wait_idle 1000
    R x_tile_on [OGFTileIsOn] "frames=$ds9(frames) current=$current(frame)"
    set cur $current(frame)
    set other [lindex [lsearch -all -inline -not $ds9(frames) $cur] 0]
    GotoFrame $cur; wait_idle 300
    if {$other ne {}} {
	lassign [img2canvas $other $ax $ay] tx ty
	CatalogPanelClearSelection
	xclick $tx $ty
	set open [winfo exists .ogfpick]
	R x_tile_click_chooser $open "chooser=$open current=$current(frame)"
	if {$open} {
	    set bb [.ogfpick.l bbox 0]
	    exec xdotool mousemove --sync [expr {[winfo rootx .ogfpick.l]+[lindex $bb 0]+20}] [expr {[winfo rooty .ogfpick.l]+[lindex $bb 1]+[lindex $bb 3]/2}]
	    wait_idle 150; exec xdotool click 1; wait_idle 500
	    R x_tile_pick_selects [expr {[selnum] ne {}}] "sel=[selnum]"
	}
	lassign [img2canvas $other $bx $by] tx ty
	CatalogPanelClearSelection
	xclick $tx $ty
	R x_tile_single_direct [expr {![winfo exists .ogfpick] && [selnum] eq {M3}}] "sel=[selnum]"
    }
    OGFPickDismiss
    OGFUIDisplay single; wait_idle 800
    R x_geom_end [expr {[geom x] eq $g0}] "start: $g0 | end: [geom x]"
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
