#  OGFinder: image <-> catalog-table link (click, hover, multi-select, keys).
#  Part of the OGFinder multi-band / mask / link extension.
#  See docs/multiband_mask_link.md

package provide DS9 1.0

proc OGFLinkInit {} {
    global catpanel catcache
    ::ogf::cat::set sel,nums {}
    ::ogf::cat::set sel,base {No source selected}
    ::ogf::cat::set hover,text {}
    ::ogf::cat::set hover,last 0
    ::ogf::cat::set hover,pending {}
    ::ogf::cat::set hover,after {}
    ::ogf::cat::set hover,num {}
    ::ogf::cat::set cache,dirty 1
    set catcache(n) 0
    ::ogf::cat::trace add alldata OGFCacheDirty
    [::ogf::cat::get tbl] tag configure msel -bg #9cc7f5 -fg black
    bind [::ogf::cat::get tbl] <Key-n> {CatalogPanelStepKey 1; break}
    bind [::ogf::cat::get tbl] <Key-p> {CatalogPanelStepKey -1; break}
}

proc OGFCacheDirty {args} {
    ::ogf::cat::set cache,dirty 1
}

# ---- parsed catalog cache (used for hit tests; rebuilt lazily) ----
proc OGFCacheBuild {} {
    global catcache
    if {![::ogf::cat::get cache,dirty]} return
    ::ogf::cat::set cache,dirty 0
    foreach k {num x y sa sb th mag} {set catcache($k) {}}
    set catcache(n) 0
    set _tsv [OGFViewTSV]
    if {$_tsv eq {}} return
    set lines [split $_tsv \n]
    set headers [split [lindex $lines 0] "\t"]
    set cx -1; set cy -1; set ca -1; set cb -1; set ct -1; set ci -1
    set cn -1; set cm -1
    set c 0
    foreach h $headers {
	switch -- [string trim $h] {
	    NUMBER {set cn $c} X_IMAGE {set cx $c} Y_IMAGE {set cy $c}
	    A_IMAGE {set ca $c} B_IMAGE {set cb $c} THETA_IMAGE {set ct $c}
	    ISO_RADIUS {set ci $c} MAG_AUTO {set cm $c}
	}
	incr c
    }
    if {$cx < 0 || $cy < 0} return
    set nums {}; set xs {}; set ys {}; set sas {}; set sbs {}; set ths {}; set mags {}
    set i 0
    foreach line [lrange $lines 1 end] {
	incr i
	if {[string trim $line] eq {}} continue
	set f [split $line "\t"]
	set x [string trim [lindex $f $cx]]; set y [string trim [lindex $f $cy]]
	if {![string is double -strict $x] || ![string is double -strict $y]} continue
	set n $i
	if {$cn >= 0} {
	    set nv [string trim [lindex $f $cn]]
	    if {$nv ne {}} {set n $nv}
	}
	set iso 5.0; set a 0; set b 0; set th 0
	if {$ci >= 0} {set v [string trim [lindex $f $ci]]; if {[string is double -strict $v] && $v > 0} {set iso $v}}
	if {$ca >= 0} {set v [string trim [lindex $f $ca]]; if {[string is double -strict $v] && $v > 0} {set a $v}}
	if {$cb >= 0} {set v [string trim [lindex $f $cb]]; if {[string is double -strict $v] && $v > 0} {set b $v}}
	if {$ct >= 0} {set v [string trim [lindex $f $ct]]; if {[string is double -strict $v]} {set th $v}}
	set sb $iso
	if {$a > 0 && $b > 0} {set sb [expr {$iso * $b / $a}]}
	set m {}
	if {$cm >= 0} {set m [string trim [lindex $f $cm]]}
	lappend nums $n; lappend xs $x; lappend ys $y
	lappend sas $iso; lappend sbs $sb; lappend ths [expr {$th * 0.017453292519943295}]
	lappend mags $m
    }
    set catcache(num) $nums; set catcache(x) $xs; set catcache(y) $ys
    set catcache(sa) $sas; set catcache(sb) $sbs; set catcache(th) $ths
    set catcache(mag) $mags
    set catcache(n) [llength $nums]
}

# smallest ellipse (floor of 2 px semi-axes) containing image point (ix,iy)
proc OGFHitTest {ix iy} {
    global catcache
    OGFCacheBuild
    if {$catcache(n) == 0} {return {}}
    set best {}; set ba 1e30
    set k -1
    foreach n $catcache(num) x $catcache(x) y $catcache(y) sa $catcache(sa) \
	sb $catcache(sb) th $catcache(th) {
	incr k
	if {$sa < 2.0} {set sa 2.0}
	set dx [expr {$ix - $x}]
	if {$dx > $sa || $dx < -$sa} continue
	set dy [expr {$iy - $y}]
	if {$dy > $sa || $dy < -$sa} continue
	if {$sb < 2.0} {set sb 2.0}
	set c [expr {cos($th)}]; set s [expr {sin($th)}]
	set rx [expr {$c*$dx + $s*$dy}]; set ry [expr {-$s*$dx + $c*$dy}]
	if {($rx*$rx)/($sa*$sa) + ($ry*$ry)/($sb*$sb) <= 1.0} {
	    set area [expr {$sa*$sb}]
	    if {$area < $ba} {set ba $area; set best $n}
	}
    }
    return $best
}

proc OGFMagOfNumber {num} {
    global catcache
    OGFCacheBuild
    set k [lsearch -exact $catcache(num) $num]
    if {$k < 0} {return {}}
    return [lindex $catcache(mag) $k]
}

# Is `which` a frame on the detection pixel grid (so X_IMAGE/Y_IMAGE apply)?
proc OGFFrameIsDetGrid {which} {
    if {[info commands OGFBandPos] eq {}} {return 1}
    lassign [OGFBandPos $which {} 0 0] x y scale
    return [expr {$scale == 1.0}]
}

# ---- table helpers ----
proc OGFTableCol {name} {return [::ogf::cat::table_col $name]}
proc OGFRowOfNumber {num} {return [::ogf::cat::row_of $num]}
proc OGFNumberOfRow {row} {return [::ogf::cat::number_of $row]}

# ---- summary text (selected-source block + optional hover line) ----
proc OGFRefreshSelText {} {
    set t [::ogf::cat::get sel,base]
    if {[::ogf::cat::get hover,text] ne {}} {append t "\n" [::ogf::cat::get hover,text]}
    ::ogf::cat::set sel,text $t
}

proc OGFSetSelBase {text} {
    ::ogf::cat::set sel,base $text
    OGFRefreshSelText
}

proc OGFMultiSummary {} {
    global catcache
    OGFCacheBuild
    set nums [::ogf::cat::get sel,nums]
    set n [llength $nums]
    set mags {}
    foreach u $nums {
	set m [OGFMagOfNumber $u]
	if {[string is double -strict $m] && $m < 90} {lappend mags $m}
    }
    set ids [join [lrange $nums 0 7] ", #"]
    if {$n > 8} {append ids ", ..."}
    set txt [format "%d sources selected  (last #%s)\n#%s" $n [lindex $nums end] $ids]
    if {[llength $mags]} {
	set mags [lsort -real $mags]
	append txt [format "\nMAG_AUTO %.2f .. %.2f" [lindex $mags 0] [lindex $mags end]]
    }
    return $txt
}

# ---- selection core ----
# mode: replace | add | toggle.  pan: recentre the current frame on the source.
proc CatalogPanelLinkSelect {num {mode replace} {pan 1}} {
    set row [OGFRowOfNumber $num]
    if {$row < 0} {
	::ogf::cat::set status "Source #$num is not in the table (filtered out?)"
	return
    }
    set nums [::ogf::cat::get sel,nums]
    switch -- $mode {
	replace {set nums [list $num]}
	add {
	    set k [lsearch -exact $nums $num]
	    if {$k >= 0} {set nums [lreplace $nums $k $k]}
	    lappend nums $num
	}
	toggle {
	    set k [lsearch -exact $nums $num]
	    if {$k >= 0} {set nums [lreplace $nums $k $k]} else {lappend nums $num}
	}
    }
    ::ogf::cat::set sel,nums $nums
    if {[llength $nums] == 0} {
	CatalogPanelClearSelection
	return
    }
    OGFApplySelection $pan
}

proc OGFApplySelection {pan} {
    global current
    set tbl [::ogf::cat::get tbl]
    set nums [::ogf::cat::get sel,nums]
    set prim [lindex $nums end]
    set prow [OGFRowOfNumber $prim]
    if {$prow < 0} return
    $tbl selection clear all
    $tbl selection set $prow,1
    $tbl see $prow,1
    # extra (non-primary) rows get the msel tag
    $tbl tag delete msel
    $tbl tag configure msel -bg #9cc7f5 -fg black
    foreach u [lrange $nums 0 end-1] {
	set r [OGFRowOfNumber $u]
	if {$r >= 0} {$tbl tag row msel $r}
    }
    if {[llength $nums] == 1} {
	catch {CatalogPanelUpdateSelInfo $prow}
	OGFSetSelBase [::ogf::cat::get sel,text]
    } else {
	OGFSetSelBase [OGFMultiSummary]
    }
    after cancel CatalogPanelGotoSource
    CatalogPanelGotoSource $prow $pan
    OGFDrawExtraMarkers
}

proc OGFSelFrames {} {
    global current
    set fr {}
    if {[info commands OGFBandFrames] ne {}} {set fr [OGFBandFrames 1]}
    if {$current(frame) ne {} && $current(frame) ni $fr} {lappend fr $current(frame)}
    return $fr
}

# small cyan crosses for the non-primary members of a multi-selection
proc OGFDrawExtraMarkers {} {
    global catcache
    OGFCacheBuild
    foreach fr [OGFSelFrames] {
	catch {$fr marker catalog sextract_msel delete}
    }
    set extra [lrange [::ogf::cat::get sel,nums] 0 end-1]
    if {[llength $extra] == 0} return
    global ogf_msel_reg
    foreach fr [OGFSelFrames] {
	if {![$fr has fits]} continue
	set reg "image\n"
	foreach u $extra {
	    set k [lsearch -exact $catcache(num) $u]
	    if {$k < 0} continue
	    set x [lindex $catcache(x) $k]; set y [lindex $catcache(y) $k]
	    if {[info commands OGFBandPos] ne {}} {
		lassign [OGFBandPos $fr $u $x $y] x y
	    }
	    append reg "cross point($x $y) # color=cyan width=2 point=cross 12 tag={sextract_msel} select=0 edit=0 move=0 rotate=0 delete=1\n"
	}
	set ogf_msel_reg $reg
	catch {$fr marker catalog command ds9 var ogf_msel_reg}
    }
}

proc CatalogPanelClearSelection {} {
    ::ogf::cat::set sel,nums {}
    catch {[::ogf::cat::get tbl] selection clear all}
    catch {[::ogf::cat::get tbl] tag delete msel}
    catch {[::ogf::cat::get tbl] tag configure msel -bg #9cc7f5 -fg black}
    foreach fr [OGFSelFrames] {
	catch {$fr marker catalog sextract_sel delete}
	catch {$fr marker catalog sextract_msel delete}
    }
    catch {OGFTDClearSel}
    OGFSetSelBase {No source selected}
}

# ---- mouse hooks (called from frame.tcl) ----
proc OGFHitAt {which x y} {
    if {[catch {$which get coordinates $x $y image} c]} {return {}}
    if {![OGFFrameIsDetGrid $which]} {return {}}
    return [OGFHitTest [lindex $c 0] [lindex $c 1]]
}

# Button-1 press on a source without DS9 catalog markers -> select it.
proc CatalogPanelLinkPress {which x y} {
    if {![::ogf::cat::exists tbl]} {return 0}
    set n [OGFHitAt $which $x $y]
    if {$n eq {}} {return 0}
    CatalogPanelLinkSelect $n replace 1
    return 1
}

proc CatalogPanelLinkShiftClick {which x y} {
    if {![::ogf::cat::exists tbl]} return
    set n [OGFHitAt $which $x $y]
    if {$n eq {}} return
    CatalogPanelLinkSelect $n toggle 0
}

proc CatalogPanelLinkRelease {which x y} {
    global ds9
    if {![info exists ds9(none_press_x)]} return
    if {abs($x - $ds9(none_press_x)) + abs($y - $ds9(none_press_y)) > 3} return
    if {[llength [::ogf::cat::get sel,nums]] == 0} return
    if {[::ogf::cat::exists add_objects_mode] && [::ogf::cat::get add_objects_mode]} return
    CatalogPanelClearSelection
    ::ogf::cat::set status "Selection cleared"
}

# ---- hover readout (throttled) ----
proc CatalogPanelHover {which x y} {
    global catcache
    if {![::ogf::cat::exists hover,last]} return
    if {![::ogf::cat::get cache,dirty] && $catcache(n) == 0} return
    set now [clock milliseconds]
    if {$now - [::ogf::cat::get hover,last] < 60} {
	::ogf::cat::set hover,pending [list $which $x $y]
	if {[::ogf::cat::get hover,after] eq {}} {
	    ::ogf::cat::set hover,after [after 70 OGFHoverFlush]
	}
	return
    }
    OGFHoverDo $which $x $y
}

proc OGFHoverFlush {} {
    ::ogf::cat::set hover,after {}
    if {[::ogf::cat::get hover,pending] ne {}} {
	lassign [::ogf::cat::get hover,pending] w x y
	::ogf::cat::set hover,pending {}
	OGFHoverDo $w $x $y
    }
}

proc OGFHoverDo {which x y} {
    ::ogf::cat::set hover,last [clock milliseconds]
    set n [OGFHitAt $which $x $y]
    if {$n eq [::ogf::cat::get hover,num]} return
    ::ogf::cat::set hover,num $n
    if {$n eq {}} {
	::ogf::cat::set hover,text {}
    } else {
	set m [OGFMagOfNumber $n]
	if {[string is double -strict $m]} {set m [format %.2f $m]} else {set m -}
	::ogf::cat::set hover,text "Hover: #$n   MAG_AUTO = $m"
    }
    OGFRefreshSelText
}

# ---- keyboard stepping ----
# dir = +1 / -1.  require_sel=1: only act when something is selected
# (used by Up/Down in the image so the cursor-warp keys keep working).
proc CatalogPanelStepKey {dir {require_sel 0}} {
    if {![::ogf::cat::exists tbl]} {return 0}
    if {[::ogf::cat::exists ai,active] && [::ogf::cat::get ai,active]} {return 0}
    if {[OGFViewTSV] eq {}} {return 0}
    if {$require_sel && [llength [::ogf::cat::get sel,nums]] == 0} {return 0}
    set nr [expr {[[::ogf::cat::get tbl] cget -rows] - 1}]
    if {$nr < 1} {return 0}
    set cur -1
    if {[llength [::ogf::cat::get sel,nums]]} {
	set cur [OGFRowOfNumber [lindex [::ogf::cat::get sel,nums] end]]
    }
    if {$cur < 1} {
	set new [expr {$dir > 0 ? 1 : $nr}]
    } else {
	set new [expr {$cur + $dir}]
    }
    if {$new < 1} {set new 1}
    if {$new > $nr} {set new $nr}
    set num [OGFNumberOfRow $new]
    if {$num eq {}} {return 0}
    CatalogPanelLinkSelect $num replace 1
    return 1
}
