# catalog/catalog_table.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/catalog/plugin.json.

proc CatalogPanelFilter {} {
    global catpanel

    # time-domain views filter their own rows
    if {[catch {::ogf::td::active} _tda]} {set _tda 0}
    if {$_tda} {::ogf::td::fill; return}
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set pattern $catpanel(search_var)

    global $catpanel(tbldb)

    # Unbind table while modifying
    $catpanel(tbl) configure -variable {}
    unset -nocomplain $catpanel(tbldb)

    set data $catpanel(alldata)
    set lines [split $data \n]

    # Header
    set headers [split [lindex $lines 0] "\t"]
    set ncols [llength $headers]

    for {set c 0} {$c < $ncols} {incr c} {
	set ${catpanel(tbldb)}(0,[expr {$c+1}]) \
	    [string trim [lindex $headers $c]]
    }

    # Filter data rows
    set row 1
    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	if {$pattern ne {} && ![string match -nocase "*${pattern}*" $line]} continue
	set fields [split $line "\t"]
	for {set c 0} {$c < $ncols} {incr c} {
	    set ${catpanel(tbldb)}($row,[expr {$c+1}]) \
		[string trim [lindex $fields $c]]
	}
	incr row
    }

    # Rebind table to trigger full refresh
    $catpanel(tbl) configure -variable $catpanel(tbldb) \
	-cols $ncols -rows $row

    catch {OGFTDAppendKindColumn $ncols $row}
    set ndata [expr {$row - 1}]
    if {$pattern eq {}} {
	set catpanel(status) "Showing all $ndata sources"
    } else {
	set catpanel(status) "Filtered: $ndata sources matching '$pattern'"
    }
}

# Show key columns of the selected catalog row in the info area
proc CatalogPanelUpdateSelInfo {row} {
    global catpanel
    if {![catch {::ogf::td::selinfo $row} _tds] && $_tds} return
    global $catpanel(tbldb)

    set want {NUMBER X_IMAGE Y_IMAGE ALPHA_J2000 DELTA_J2000 MAG_AUTO
	FWHM_IMAGE ELLIPTICITY CLASS_STAR}
    set idx {}
    set ncols [$catpanel(tbl) cget -cols]
    for {set c 1} {$c <= $ncols} {incr c} {
	if {[info exists ${catpanel(tbldb)}(0,$c)]} {
	    dict set idx [set ${catpanel(tbldb)}(0,$c)] $c
	}
    }
    set v {}
    foreach k $want {
	set val {-}
	if {[dict exists $idx $k]} {
	    set c [dict get $idx $k]
	    if {[info exists ${catpanel(tbldb)}($row,$c)]} {
		set val [set ${catpanel(tbldb)}($row,$c)]
	    }
	}
	if {[string is double -strict $val] && $k ne {NUMBER}} {
	    set val [format %.4g $val]
	}
	dict set v $k $val
    }
    set catpanel(sel,text) [format \
	"Source #%s   x,y = %s, %s\nRA,Dec = %s, %s\nMAG_AUTO = %s   FWHM = %s   e = %s   Star = %s" \
	[dict get $v NUMBER] [dict get $v X_IMAGE] [dict get $v Y_IMAGE] \
	[dict get $v ALPHA_J2000] [dict get $v DELTA_J2000] \
	[dict get $v MAG_AUTO] [dict get $v FWHM_IMAGE] \
	[dict get $v ELLIPTICITY] [dict get $v CLASS_STAR]]
}

proc CatalogPanelSelectCmd {prev cur} {
    global catpanel

    # cur is "row,col" of current selection
    set row [lindex [split $cur ,] 0]
    if {![string is integer -strict $row] || $row <= 0} return

    catch {CatalogPanelUpdateSelInfo $row}
    # keep the link state in step with clicks in the table itself
    set num [OGFNumberOfRow $row]
    if {$num ne {}} {
	set catpanel(sel,nums) [list $num]
	catch {$catpanel(tbl) tag delete msel}
	catch {$catpanel(tbl) tag configure msel -bg #9cc7f5 -fg black}
	OGFSetSelBase $catpanel(sel,text)
    }
    after cancel CatalogPanelGotoSource
    after 100 [list CatalogPanelGotoSource $row]
}

proc CatalogPanelGotoSource {row {pan 1}} {
    global catpanel
    global current
    global ds9

    if {![catch {::ogf::td::goto $row $pan} _tdg] && $_tdg} return

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return

    global $catpanel(tbldb)

    # Find column indices from header row
    set ncols [$catpanel(tbl) cget -cols]
    set col_x -1
    set col_y -1
    set col_a -1
    set col_b -1
    set col_theta -1
    set col_ir -1
    set col_reff -1
    for {set c 1} {$c <= $ncols} {incr c} {
	if {[info exists ${catpanel(tbldb)}(0,$c)]} {
	    set hdr [set ${catpanel(tbldb)}(0,$c)]
	    switch -- $hdr {
		X_IMAGE     { set col_x $c }
		Y_IMAGE     { set col_y $c }
		A_IMAGE     { set col_a $c }
		B_IMAGE     { set col_b $c }
		THETA_IMAGE { set col_theta $c }
		ISO_RADIUS  { set col_ir $c }
		R_EFF_PIX   { set col_reff $c }
	    }
	}
    }

    if {$col_x < 0 || $col_y < 0} return

    # Get coordinates from selected row
    if {![info exists ${catpanel(tbldb)}($row,$col_x)]} return
    set x [set ${catpanel(tbldb)}($row,$col_x)]
    set y [set ${catpanel(tbldb)}($row,$col_y)]

    if {![string is double -strict $x] || ![string is double -strict $y]} return

    # Get ellipse parameters (with NaN/Inf safety via catch)
    set iso_radius 10.0
    set a_image 0
    set b_image 0
    set theta 0

    if {$col_ir >= 0 && [info exists ${catpanel(tbldb)}($row,$col_ir)]} {
	set val [set ${catpanel(tbldb)}($row,$col_ir)]
	if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} {
	    set iso_radius $v
	}
    }
    if {$col_a >= 0 && [info exists ${catpanel(tbldb)}($row,$col_a)]} {
	set val [set ${catpanel(tbldb)}($row,$col_a)]
	if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set a_image $v }
    }
    if {$col_b >= 0 && [info exists ${catpanel(tbldb)}($row,$col_b)]} {
	set val [set ${catpanel(tbldb)}($row,$col_b)]
	if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set b_image $v }
    }
    if {$col_theta >= 0 && [info exists ${catpanel(tbldb)}($row,$col_theta)]} {
	set val [set ${catpanel(tbldb)}($row,$col_theta)]
	if {[catch {set v [expr {$val + 0.0}]}] == 0} { set theta $v }
    }

    # Get R_EFF_PIX if available (for LSBG circle display)
    set r_eff_pix 0
    if {$col_reff >= 0 && [info exists ${catpanel(tbldb)}($row,$col_reff)]} {
	set val [set ${catpanel(tbldb)}($row,$col_reff)]
	if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} {
	    set r_eff_pix $v
	}
    }

    # Compute ellipse: ISO_RADIUS as semi-major, scaled by B/A for semi-minor
    set semi_a $iso_radius
    set semi_b $iso_radius
    if {$a_image > 0 && $b_image > 0} {
	set semi_b [expr {$iso_radius * $b_image / $a_image}]
    }

    # Delete previous selection markers (in every registered band frame)
    set frame $current(frame)
    set frames [OGFSelFrames]
    foreach fr $frames {
	catch {$fr marker catalog sextract_sel delete}
    }

    # Rebuild sextract_all markers from alldata to keep image in sync
    if {[info exists catpanel(markall,on)] && $catpanel(markall,on)} {
	CatalogPanelCreateAllMarkers
    }

    # Use global variable for marker creation (var form requires global access)
    global sextract_sel_reg

    # Marker in each band frame (positions mapped through the WCS when the
    # band has a different pixel grid; identity on a shared grid)
    set num [OGFNumberOfRow $row]
    foreach fr $frames {
	if {![$fr has fits]} continue
	lassign [OGFBandPos $fr $num $x $y] bx by same
	set isdet [expr {$same}]
	set sextract_sel_reg "image\ncross point($bx $by) # color=cyan width=2 point=cross 15 tag={sextract_sel} select=0 edit=0 move=0 rotate=0 delete=1\n"
	catch {$fr marker catalog command ds9 var sextract_sel_reg}
	if {$fr ne $frame && !$isdet} continue
	if {$r_eff_pix > 0} {
	    set sextract_sel_reg "image\ncircle($bx $by ${r_eff_pix}i) # color=green width=2 dash=1 tag={sextract_sel} select=0 edit=0 move=0 rotate=0 delete=1\n"
	    catch {$fr marker catalog command ds9 var sextract_sel_reg}
	}
	set sextract_sel_reg "image\nellipse($bx $by ${semi_a}i ${semi_b}i $theta) # color=green width=2 dash=1 tag={sextract_sel} select=0 edit=0 move=0 rotate=0 delete=1\n"
	catch {$fr marker catalog command ds9 var sextract_sel_reg}
    }

    # Pan to the object (locked band frames follow through the WCS lock)
    if {$pan} {
	PanToFrame $current(frame) $x $y image {}
    }

    set catpanel(status) "Source at image ($x, $y)"
}

proc CatalogPanelTableClick {x y} {
    global catpanel

    set tbl $catpanel(tbl)
    set idx [$tbl index @$x,$y]
    set row [lindex [split $idx ,] 0]

    # Only handle header row clicks
    if {$row != 0} return

    set col [lindex [split $idx ,] 1]

    global $catpanel(tbldb)
    if {![info exists ${catpanel(tbldb)}(0,$col)]} return
    set colname [set ${catpanel(tbldb)}(0,$col)]

    # Toggle direction if same column clicked again
    if {$catpanel(sort,col) eq $colname} {
	if {$catpanel(sort,dir) eq "ascending"} {
	    set catpanel(sort,dir) descending
	} else {
	    set catpanel(sort,dir) ascending
	}
    } else {
	set catpanel(sort,col) $colname
	set catpanel(sort,dir) ascending
    }

    CatalogPanelSort $colname $catpanel(sort,dir)
}

proc CatalogPanelSort {colname direction} {
    global catpanel

    if {[catch {::ogf::td::active} _tda]} {set _tda 0}
    if {$_tda} {::ogf::td::sort $colname $direction; return}
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set lines [split $catpanel(alldata) \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]

    # Find column index
    set colidx -1
    for {set i 0} {$i < [llength $headers]} {incr i} {
	if {[string trim [lindex $headers $i]] eq $colname} {
	    set colidx $i
	    break
	}
    }
    if {$colidx < 0} return

    OGFSessLog catalog.sort auto {} -tool native -title "Sort catalog by $colname $direction" \
	-payload [dict create col $colname dir $direction] -requires catalog
    # Collect data rows (skip header and empty lines)
    set datarows {}
    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	lappend datarows $line
    }

    # Determine sort type: check first non-empty value
    set isnumeric 1
    foreach drow $datarows {
	set val [string trim [lindex [split $drow "\t"] $colidx]]
	if {$val ne {}} {
	    if {![string is double -strict $val]} {
		set isnumeric 0
	    }
	    break
	}
    }

    # Sort
    if {$isnumeric} {
	set cmd [list CatalogPanelSortCmpNum $colidx]
    } else {
	set cmd [list CatalogPanelSortCmpStr $colidx]
    }
    if {$direction eq "descending"} {
	set sortedrows [lsort -decreasing -command $cmd $datarows]
    } else {
	set sortedrows [lsort -command $cmd $datarows]
    }

    # Rebuild alldata with sorted rows
    set newdata $header
    foreach drow $sortedrows {
	append newdata "\n$drow"
    }
    set catpanel(alldata) $newdata

    # Reload table
    CatalogPanelLoadTSV $catpanel(alldata) "sorted"

    set catpanel(status) "Sorted by $colname $direction"
}

proc CatalogPanelSortCmpNum {colidx a b} {
    set va [string trim [lindex [split $a "\t"] $colidx]]
    set vb [string trim [lindex [split $b "\t"] $colidx]]
    if {![string is double -strict $va]} { set va 0 }
    if {![string is double -strict $vb]} { set vb 0 }
    if {$va < $vb} { return -1 }
    if {$va > $vb} { return 1 }
    return 0
}

proc CatalogPanelSortCmpStr {colidx a b} {
    set va [string trim [lindex [split $a "\t"] $colidx]]
    set vb [string trim [lindex [split $b "\t"] $colidx]]
    return [string compare $va $vb]
}

