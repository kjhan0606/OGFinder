#  OGFinder: click selection with overlapping markers (galaxy / moving / transient / detection).
#
#  A click on the image collects EVERY row of the table whose marker lies within the hit radius of the pointer
#  (OGFPickCandidates):
#    * time-domain rows (moving / transient / detection, ogf_td.tcl): sky distance from the click to the row position (or to
#      the nearest of its track points) <= ogfpick(radius) canvas pixels;
#    * galaxies: every catalog ellipse (floor of 2 px semi-axes, as the old smallest-ellipse test) that contains the click.
#  Candidates are sorted by distance (image pixels of the clicked frame).  One candidate -> it is selected directly (this is the
#  NEAREST row); several -> a small chooser lists  ID  kind  mag  d(px)  and selecting an entry selects that table row.
#  Esc, a click outside the chooser, or the window losing focus dismisses it without changing the selection.
#  Works in single and tile view (the clicked tile becomes the current frame, as OGFTileClick did).  No widget of the main
#  window is touched, so the layout invariants are unaffected.  See docs/click_selection.md.

package provide DS9 1.0

proc OGFPickInit {} {
    global ogfpick
    if {[info exists ogfpick(init)]} return
    array set ogfpick {init 1 enable 1 radius 8 maxlist 12 win .ogfpick last {} }
}
OGFPickInit

# image coordinates of canvas point (x,y) in frame `which`, zoom and the sky scale (arcsec per canvas pixel)
proc OGFPickFrameInfo {which x y} {
    if {[catch {$which get coordinates $x $y image} ic]} {return {}}
    set z 1.0
    catch {set z [lindex [$which get zoom] 0]}
    if {![string is double -strict $z] || $z <= 0} {set z 1.0}
    set sky {}; set scale 0.0
    if {![catch {$which get coordinates $x $y wcs fk5 degrees} c0] && [llength $c0] == 2 &&
	![catch {$which get coordinates [expr {$x+20}] $y wcs fk5 degrees} c1] && [llength $c1] == 2} {
	set sky $c0
	set s [OGFPickSep [lindex $c0 0] [lindex $c0 1] [lindex $c1 0] [lindex $c1 1]]
	set scale [expr {$s / 20.0}]
    }
    return [list [lindex $ic 0] [lindex $ic 1] $z $sky $scale]
}

# angular separation in arcsec (small-angle with cos(dec))
proc OGFPickSep {ra1 dec1 ra2 dec2} {
    set d [expr {($ra2-$ra1)}]
    if {$d > 180} {set d [expr {$d-360}]} elseif {$d < -180} {set d [expr {$d+360}]}
    set dx [expr {$d*cos(($dec1+$dec2)*0.5*0.017453292519943295)*3600.0}]
    set dy [expr {($dec2-$dec1)*3600.0}]
    return [expr {sqrt($dx*$dx+$dy*$dy)}]
}

# candidates under canvas point (x,y) of frame `which`: list of dicts {key kind id mag dist}, nearest first
proc OGFPickCandidatesRaw {which x y} {
    global ogfpick ogftd
    OGFPickInit
    set fi [OGFPickFrameInfo $which $x $y]
    if {$fi eq {}} {return {}}
    lassign $fi ix iy zoom sky scale
    set R $ogfpick(radius)                      ;# canvas pixels
    set cands {}
    # ---- time-domain rows
    if {$sky ne {} && $scale > 0 && [info exists ogftd(kind)] && $ogftd(kind) ne "galaxies"} {
	lassign $sky cra cde
	set rmax [expr {$R*$scale}]                 ;# arcsec
	set dcut [expr {$rmax/3600.0*1.5}]
	foreach k [::ogf::td::visible_kinds] {
	    set def $ogftd(def,$k)
	    set pre [dict get $def prefix]
	    lassign $ogftd(raw,$k) cols rows
	    set ci [lsearch -exact $cols id]; set cr [lsearch -exact $cols ra]; set cd [lsearch -exact $cols dec]
	    set cp [lsearch -exact $cols _pos]; set cm [lsearch -exact $cols mag]
	    foreach r $rows {
		set pos [list [lindex $r $cr] [lindex $r $cd]]
		if {$cp >= 0 && [lindex $r $cp] ne {}} {set pos [lindex $r $cp]}
		set best 1e30
		foreach {ra de} $pos {
		    if {![string is double -strict $ra] || ![string is double -strict $de]} continue
		    if {abs($de-$cde) > $dcut} continue
		    set s [OGFPickSep $cra $cde $ra $de]
		    if {$s < $best} {set best $s}
		}
		if {$best > $rmax} continue
		set mag {}
		if {$cm >= 0} {set mag [lindex $r $cm]}
		lappend cands [dict create key $pre[lindex $r $ci] kind $k id [lindex $r $ci] mag $mag \
		    dist [expr {$best/$scale/$zoom}]]
	    }
	}
    }
    # ---- galaxies (kind galaxies or all): all ellipses containing the click on the detection grid
    if {(![info exists ogftd(kind)] || $ogftd(kind) in {galaxies all}) && [OGFFrameIsDetGrid $which]} {
	foreach g [OGFHitAll $ix $iy] {
	    lassign $g n d
	    lappend cands [dict create key $n kind galaxy id $n mag [OGFMagOfNumber $n] dist $d]
	}
    }
    # only rows that are in the table can be selected
    set out {}
    foreach c $cands {
	if {[OGFRowOfNumber [dict get $c key]] >= 0} {lappend out $c}
    }
    return $out
}

# sort candidate dicts by distance
proc OGFPickSorted {cands} {
    set l {}
    foreach c $cands {lappend l [list [dict get $c dist] $c]}
    return [lmap e [lsort -real -index 0 $l] {lindex $e 1}]
}

# every catalog ellipse of the galaxy table containing image point (ix,iy): list of {number dist_px}
proc OGFHitAll {ix iy} {
    global catcache ogftd
    OGFCacheBuild
    if {$catcache(n) == 0} {return {}}
    set out {}
    foreach n $catcache(num) x $catcache(x) y $catcache(y) sa $catcache(sa) sb $catcache(sb) th $catcache(th) {
	if {[info exists ogftd(kind)] && [::ogf::td::key_kind $n] ne "galaxies"} continue
	if {$sa < 2.0} {set sa 2.0}
	set dx [expr {$ix - $x}]
	if {$dx > $sa || $dx < -$sa} continue
	set dy [expr {$iy - $y}]
	if {$dy > $sa || $dy < -$sa} continue
	if {$sb < 2.0} {set sb 2.0}
	set c [expr {cos($th)}]; set s [expr {sin($th)}]
	set rx [expr {$c*$dx + $s*$dy}]; set ry [expr {-$s*$dx + $c*$dy}]
	if {($rx*$rx)/($sa*$sa) + ($ry*$ry)/($sb*$sb) <= 1.0} {
	    lappend out [list $n [expr {sqrt($dx*$dx+$dy*$dy)}]]
	}
    }
    return $out
}

proc OGFPickCandidates {which x y} {
    return [OGFPickSorted [OGFPickCandidatesRaw $which $x $y]]
}

# ---- entry point called from Button1Frame / OGFTileClick.  Returns 1 when the click was handled.
proc OGFPickClick {which x y {rootx {}} {rooty {}}} {
    global ogfpick current
    OGFPickInit
    if {!$ogfpick(enable) || ![::ogf::cat::exists tbl]} {return 0}
    OGFPickDismiss
    set c [OGFPickCandidates $which $x $y]
    set ogfpick(last) $c
    if {[llength $c] == 0} {return 0}
    if {[llength $c] == 1} {
	OGFPickSelect $which [dict get [lindex $c 0] key]
	return 1
    }
    if {$rootx eq {}} {
	set rootx [winfo pointerx .]; set rooty [winfo pointery .]
    }
    OGFPickChooser $which $c $rootx $rooty
    return 1
}

proc OGFPickSelect {which key} {
    global current
    if {$which ne $current(frame)} {catch {GotoFrame $which}}
    CatalogPanelLinkSelect $key replace 1
}

proc OGFPickDismiss {} {
    global ogfpick
    if {![info exists ogfpick(win)]} return
    set w $ogfpick(win)
    # restore the application-wide bindings installed while the chooser is up
    foreach {ev key} {<ButtonPress> old_bp <Escape> old_esc} {
	if {[info exists ogfpick($key)]} {
	    catch {bind all $ev $ogfpick($key)}
	    unset ogfpick($key)
	}
    }
    if {[winfo exists $w]} {destroy $w}
}

proc OGFPickFormat {c} {
    set m [dict get $c mag]
    if {[string is double -strict $m]} {set m [format %.2f $m]} else {set m -}
    return [format "%-8s %-10s %7s %8.2f" [dict get $c id] [dict get $c kind] $m [dict get $c dist]]
}

# the chooser: header + one line per candidate, nearest first.  Enter / double click / click on a line = select it;
# Esc, a click anywhere outside, or focus loss = dismiss.
proc OGFPickChooser {which cands rootx rooty} {
    global ogfpick
    set w $ogfpick(win)
    OGFPickDismiss
    toplevel $w -class OGFPick -relief solid -bd 1
    wm withdraw $w
    wm overrideredirect $w 1
    wm transient $w .
    label $w.h -text [format "%-8s %-10s %7s %8s" ID kind mag "d (px)"] -font TkFixedFont -anchor w -bg #e8e8e8 -padx 4
    set n [llength $cands]
    set shown [expr {min($n, $ogfpick(maxlist))}]
    listbox $w.l -font TkFixedFont -height $shown -width 36 -activestyle none -exportselection 0 \
	-selectmode browse -bd 0 -highlightthickness 0 -selectbackground #3c78d8 -selectforeground white
    foreach c [lrange $cands 0 [expr {$shown-1}]] {$w.l insert end [OGFPickFormat $c]}
    pack $w.h -fill x
    pack $w.l -fill both -expand 1
    if {$n > $shown} {
	label $w.m -text "... [expr {$n-$shown}] more (closer ones listed first)" -font TkSmallCaptionFont -anchor w
	pack $w.m -fill x
    }
    set ogfpick(cands) [lrange $cands 0 [expr {$shown-1}]]
    set ogfpick(which) $which
    $w.l selection set 0
    $w.l activate 0
    # keep the chooser on the screen
    update idletasks
    set sw [winfo screenwidth .]; set sh [winfo screenheight .]
    set ww [winfo reqwidth $w]; set wh [winfo reqheight $w]
    set px [expr {min($rootx+6, $sw-$ww)}]; set py [expr {min($rooty+6, $sh-$wh)}]
    wm geometry $w +[expr {max($px,0)}]+[expr {max($py,0)}]
    wm deiconify $w
    raise $w
    bind $w.l <Return> {OGFPickChoose [%W index active]}
    bind $w.l <KP_Enter> {OGFPickChoose [%W index active]}
    bind $w.l <ButtonRelease-1> {if {%x >= 0 && %y >= 0 && %x < [winfo width %W] && %y < [winfo height %W]} {OGFPickChoose [%W nearest %y]}; break}
    bind $w.l <Motion> {%W selection clear 0 end; %W selection set [%W nearest %y]; %W activate [%W nearest %y]}
    foreach b [list $w $w.l $w.h] {bind $b <Escape> {after idle OGFPickDismiss}}
    # Esc anywhere and a click anywhere outside the chooser dismiss it (no grab: the original press/release of the click
    # must reach ds9, and a click-away keeps its normal meaning).  The old "all" bindings are restored by OGFPickDismiss.
    # The bindings are armed only after the current event has been fully processed: the ButtonPress that opened the chooser
    # is still being dispatched (widget -> class -> toplevel -> all) and would otherwise reach the new "all" binding as a
    # click-away and close the chooser at once (found with real X events, scripts/verify_click_xevent.tcl).
    after idle [list OGFPickArm $w]
    catch {focus -force $w.l}
}

proc OGFPickArm {w} {
    global ogfpick
    if {![winfo exists $w] || [info exists ogfpick(old_bp)]} return
    set ogfpick(old_bp) [bind all <ButtonPress>]
    set ogfpick(old_esc) [bind all <Escape>]
    bind all <ButtonPress> "$ogfpick(old_bp)\nOGFPickOutside %W %X %Y"
    bind all <Escape> "$ogfpick(old_esc)\nafter idle OGFPickDismiss"
}

# a ButtonPress anywhere in the application: outside the chooser rectangle = click-away
proc OGFPickOutside {w X Y} {
    global ogfpick
    set win $ogfpick(win)
    if {![winfo exists $win]} {OGFPickDismiss; return}
    set x0 [winfo rootx $win]; set y0 [winfo rooty $win]
    if {$X < $x0 || $Y < $y0 || $X >= $x0 + [winfo width $win] || $Y >= $y0 + [winfo height $win]} {
	after idle OGFPickDismiss
    }
}

proc OGFPickChoose {idx} {
    global ogfpick
    if {![string is integer -strict $idx] || $idx < 0 || $idx >= [llength $ogfpick(cands)]} return
    set key [dict get [lindex $ogfpick(cands) $idx] key]
    set which $ogfpick(which)
    OGFPickDismiss
    OGFPickSelect $which $key
}
