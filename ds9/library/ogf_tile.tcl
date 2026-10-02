#  OGFinder: image tab strip + one-canvas tile (mosaic) mode.
#
#  The tab strip under the image canvas has one tab per frame (select = GotoFrame) and, on the right,
#  [Single] [Tile] [Layout v].  "Tile" is ds9's own one-canvas tile display (grid / column / row) with the
#  options "WCS lock" (pan/zoom/crosshair locked through the sky coordinates) and "Share scale"; "Single"
#  returns to one frame at a time.  The right-click "Tile all" item of the strip and Tools > Tile all frames
#  call the same procs (OGFUIDisplay).  Tabs never replace the mosaic; both views show the same frames.
#
#  Catalog graphics in tile mode:
#    * markers (sextract_all) and the selection (sextract_sel / msel) are drawn in EVERY tiled frame that shares
#      the catalog = the registered bands (Bands > Register...); positions are mapped through the WCS when the
#      pixel grids differ.  Frames that are not registered as bands keep their own per-frame catalog (ds9
#      swaps it in when the frame becomes current), so they are not drawn into.
#    * a click on a marker in any tile selects the same table row (OGFTileClick, called from Button1Frame).

package provide DS9 1.0

proc OGFTileInit {} {
    global ogftile
    array set ogftile {strip {} cur {} pending {} tabs {} width 0}
}

# ---------------------------------------------------------------- frames sharing the catalog
proc OGFTileIsOn {} {
    global ds9
    return [expr {[info exists ds9(display)] && $ds9(display) eq "tile"}]
}

# draw `sextract_all_reg` (built for the current frame) into the other tiles, mapping positions when the
# band's pixel grid differs from the detection grid
proc OGFMarkSend {frame} {
    global sextract_all_reg
    catch {$frame marker catalog command ds9 var sextract_all_reg}
    foreach fr [OGFSelFrames] {
	if {$fr eq $frame} continue
	if {[catch {$fr has fits} h] || !$h} continue
	set reg [OGFMapRegions $fr $sextract_all_reg]
	global ogf_mark_tmp
	set ogf_mark_tmp $reg
	catch {$fr marker catalog command ds9 var ogf_mark_tmp}
    }
}

proc OGFMarkDelete {frame} {
    catch {$frame marker catalog sextract_all delete}
    foreach fr [OGFSelFrames] {
	if {$fr ne $frame} {catch {$fr marker catalog sextract_all delete}}
    }
}

# pixel-scale ratio detection->band (axis lengths scale with it)
proc OGFBandAxisRatio {frame} {
    global ogfband
    set b [OGFBandOfFrame $frame]
    if {$b eq {} || $ogfband(detect) eq {} || $b eq $ogfband(detect) || $ogfband($b,same)} {return 1.0}
    set di $ogfband($ogfband(detect),info); set bi $ogfband($b,info)
    if {[catch {
	set dd [expr {abs([dict get $di CD1_1]*[dict get $di CD2_2] - [dict get $di CD1_2]*[dict get $di CD2_1])}]
	set db [expr {abs([dict get $bi CD1_1]*[dict get $bi CD2_2] - [dict get $bi CD1_2]*[dict get $bi CD2_1])}]
	set r [expr {sqrt($dd/$db)}]
    }]} {return 1.0}
    return $r
}

proc OGFMapRegions {frame reg} {
    if {[OGFBandOfFrame $frame] eq {}} {return $reg}
    lassign [OGFBandPos $frame {} 1 1] _x _y same
    if {$same} {return $reg}
    set ratio [OGFBandAxisRatio $frame]
    set out {}
    foreach line [split $reg \n] {
	if {[regexp {^ellipse\(([^ ]+) ([^ ]+) ([^ ]+)i ([^ ]+)i ([^)]+)\)(.*)$} $line -> x y a b th rest]} {
	    set num {}
	    regexp {tag=\{sextract_src\.([^\}]+)\}} $rest -> num
	    lassign [OGFBandPos $frame $num $x $y] bx by sm
	    append out "ellipse($bx $by [expr {$a*$ratio}]i [expr {$b*$ratio}]i $th)$rest\n"
	} elseif {$line ne {}} {
	    append out $line \n
	}
    }
    return $out
}

# ---------------------------------------------------------------- click on a marker in any tile
# Called from Button1Frame (mode none) for a frame that is not the current one.  Returns 1 when the click
# hit a catalog marker and selected its table row.
proc OGFTileClick {which x y} {
    global current
    if {![OGFTileIsOn] || ![::ogf::cat::exists tbl]} {return 0}
    # all rows under the click in this tile (chooser when several, the nearest row when one), see ogf_pick.tcl
    if {![catch {OGFPickClick $which $x $y} _pk] && $_pk} {return 1}
    if {[catch {$which get marker catalog id $x $y} id] || $id == 0} {return 0}
    set tags [$which get marker catalog $id tag]
    set num [OGFTDKeyFromTags $tags]
    if {$num eq {}} {
	if {$which ni [OGFSelFrames]} {return 0}
	foreach tag $tags {
	    if {[string match "sextract_src.*" $tag]} {set num [string range $tag 13 end]; break}
	}
    }
    if {$num eq {}} {return 0}
    # the catalog is shared by the tiles: make this tile current (no state swap), then select the row
    GotoFrame $which
    CatalogPanelLinkSelect $num replace 1
    return 1
}

# ---------------------------------------------------------------- display modes
proc OGFUIDisplay {mode} {
    global ogfui current ds9 tile panzoom crosshair scale
    if {$mode eq "single"} {
	set ogfui(tile) 0
	set current(display) single
	DisplayMode
	OGFTileAfter
	::ogf::cat::set status "Single frame view"
	return
    }
    if {[llength $ds9(frames)] < 2} {
	set ogfui(tile) 0
	::ogf::cat::set status "Tile: load at least two frames"
	OGFTileAfter
	return
    }
    set ogfui(tile) 1
    set current(display) tile
    set tile(mode) $ogfui(tilemode)
    DisplayMode
    if {$ogfui(lockwcs)} {
	set panzoom(lock) wcs
	LockFrameCurrent
	set crosshair(lock) wcs
	LockCrosshairCurrent
    } else {
	set panzoom(lock) none
	set crosshair(lock) none
    }
    if {$ogfui(sharescale)} {
	set scale(lock) 1
	LockScaleCurrent
    } else {
	set scale(lock) 0
    }
    ZoomToFit
    OGFTileAfter
    ::ogf::cat::set status "Tile: [llength $ds9(active)] frames, layout $ogfui(tilemode)[expr {$ogfui(lockwcs) ? {, WCS lock} : {}}][expr {$ogfui(sharescale) ? {, shared scale} : {}}]"
}

proc OGFUITileToggle {} {
    global ogfui
    if {$ogfui(tile)} {OGFUIDisplay tile} else {OGFUIDisplay single}
}

proc OGFUITileOptions {} {
    global ogfui
    if {$ogfui(tile)} {OGFUIDisplay tile}
}

# after every Single/Tile change: redraw catalog graphics in all frames of the view, refresh the strip
proc OGFTileAfter {} {
    global ogfui
    set ogfui(tile) [OGFTileIsOn]
    catch {CatalogPanelSyncInfoHeight}
    if {[::ogf::cat::exists markall,on] && [::ogf::cat::get markall,on]} {catch {CatalogPanelCreateAllMarkers}}
    if {[::ogf::cat::exists sel,nums] && [llength [::ogf::cat::get sel,nums]]} {catch {OGFApplySelection 0}}
    OGFTileRefresh
}

# ---------------------------------------------------------------- the strip
proc OGFTileBuildStrip {parent} {
    global ogftile ogfui
    if {![info exists ogftile]} {OGFTileInit}
    if {![info exists ogfui]} {OGFUIInit}
    set s [ttk::frame $parent.tabs -height 24]
    ttk::style configure OGFTab.Toolbutton -padding {6 1}
    ttk::frame $s.t
    ttk::frame $s.r
    ttk::radiobutton $s.r.single -text Single -style OGFTab.Toolbutton -variable ogfui(tile) -value 0 -command {OGFUIDisplay single}
    ttk::radiobutton $s.r.tile -text Tile -style OGFTab.Toolbutton -variable ogfui(tile) -value 1 -command {OGFUIDisplay tile}
    ttk::menubutton $s.r.lay -text "Layout \u25be" -style OGFTab.Toolbutton -menu $s.r.lay.m
    menu $s.r.lay.m -tearoff 0
    OGFTileFillOptions $s.r.lay.m
    ttk::menubutton $s.r.all -text "\u25be" -style OGFTab.Toolbutton -menu $s.r.all.m
    menu $s.r.all.m -tearoff 0 -postcommand [list OGFTileFrameMenu $s.r.all.m]
    pack $s.r.all $s.r.lay $s.r.tile $s.r.single -side right -padx 1
    pack $s.r -side right -fill y
    pack $s.t -side left -fill both -expand 1
    pack propagate $s 0
    set ogftile(strip) $s
    # right-click anywhere on the strip: same functions
    menu $s.pop -tearoff 0
    $s.pop add command -label "Tile all" -command {OGFUIDisplay tile}
    $s.pop add command -label "Single" -command {OGFUIDisplay single}
    $s.pop add separator
    OGFTileFillOptions $s.pop
    foreach w [list $s $s.t $s.r] {bind $w <Button-3> [list tk_popup $s.pop %X %Y]}
    grid $s -row 2 -column 0 -sticky ew
    foreach v {::ds9(frames) ::current(frame) ::ds9(display)} {
	trace add variable $v write OGFTileSoon
    }
    return $s
}

proc OGFTileFillOptions {m} {
    foreach {v l} {grid Grid column Column row Row} {
	$m add radiobutton -label "Layout: $l" -variable ogfui(tilemode) -value $v -command OGFUITileOptions
    }
    $m add separator
    $m add checkbutton -label "Lock pan/zoom (WCS)" -variable ogfui(lockwcs) -command OGFUITileOptions
    $m add checkbutton -label "Share scale limits" -variable ogfui(sharescale) -command OGFUITileOptions
}

proc OGFTileFrameMenu {m} {
    global ds9 current
    $m delete 0 end
    foreach fr $ds9(frames) {
	$m add radiobutton -label [OGFTileFrameLabel $fr] -variable current(frame) -value $fr -command [list GotoFrame $fr]
    }
}

proc OGFTileFrameLabel {fr} {
    set b {}
    catch {set b [OGFBandOfFrame $fr]}
    if {$b ne {}} {return $b}
    set fn {}
    catch {set fn [file tail [string trim [$fr get fits file name full] "{}"]]}
    regsub {\[.*\]$} $fn {} fn
    if {$fn eq {}} {set fn $fr}
    return $fn
}

proc OGFTileSoon {args} {
    global ogftile
    if {![info exists ogftile(pending)] || $ogftile(pending) ne {}} return
    set ogftile(pending) [after idle OGFTileRefresh]
}

proc OGFTileRefresh {} {
    global ogftile ds9 current ogfui
    set ogftile(pending) {}
    if {![info exists ogftile(strip)] || ![winfo exists $ogftile(strip)]} return
    set t $ogftile(strip).t
    set want [lmap fr $ds9(frames) {list $fr [OGFTileFrameLabel $fr]}]
    if {$want ne $ogftile(tabs)} {
	foreach w [winfo children $t] {destroy $w}
	set i 0
	foreach item $want {
	    lassign $item fr lab
	    ttk::radiobutton $t.f$i -text $lab -style OGFTab.Toolbutton -variable current(frame) -value $fr \
		-command [list GotoFrame $fr]
	    pack $t.f$i -side left -padx 1
	    bind $t.f$i <Button-3> [list tk_popup $ogftile(strip).pop %X %Y]
	    incr i
	}
	set ogftile(tabs) $want
    }
    set ogfui(tile) [OGFTileIsOn]
}
