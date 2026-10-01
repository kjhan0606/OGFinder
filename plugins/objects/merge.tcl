# objects/merge.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/objects/plugin.json.

proc CatalogPanelMergeSources {} {
    global catpanel
    global current

    if {!$catpanel(merge,active)} return
    if {[llength $catpanel(merge,list)] < 2} {
	set catpanel(status) "Need at least 2 sources to merge"
	return
    }

    # Parse alldata
    set lines [split $catpanel(alldata) \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]
    set ncols [llength $headers]

    # Find column indices
    set idx_num -1
    set idx_x -1
    set idx_y -1
    set idx_a -1
    set idx_b -1
    set idx_theta -1
    set idx_ir -1
    set idx_flux -1
    set idx_mag -1
    set idx_npix -1
    for {set i 0} {$i < $ncols} {incr i} {
	set h [string trim [lindex $headers $i]]
	switch -- $h {
	    NUMBER      { set idx_num $i }
	    X_IMAGE     { set idx_x $i }
	    Y_IMAGE     { set idx_y $i }
	    A_IMAGE     { set idx_a $i }
	    B_IMAGE     { set idx_b $i }
	    THETA_IMAGE { set idx_theta $i }
	    ISO_RADIUS  { set idx_ir $i }
	    FLUX_AUTO   { set idx_flux $i }
	    MAG_AUTO    { set idx_mag $i }
	    NPIX_ISO    { set idx_npix $i }
	}
    }

    # Collect data for merge sources and find brightest
    set merge_rows {}
    set other_rows {}
    set max_number 0
    set brightest_idx -1
    set brightest_flux -1e30

    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	set fields [split $line "\t"]
	set num_val [string trim [lindex $fields $idx_num]]

	# Track max NUMBER
	if {[string is integer -strict $num_val] && $num_val > $max_number} {
	    set max_number $num_val
	}

	if {[lsearch -exact $catpanel(merge,list) $num_val] >= 0} {
	    lappend merge_rows $fields
	    if {$idx_flux >= 0} {
		set fval [string trim [lindex $fields $idx_flux]]
		if {[string is double -strict $fval] && $fval > $brightest_flux} {
		    set brightest_flux $fval
		    set brightest_idx [expr {[llength $merge_rows] - 1}]
		}
	    }
	} else {
	    lappend other_rows $line
	}
    }

    if {[llength $merge_rows] < 2} {
	set catpanel(status) "Merge error: sources not found in catalog"
	return
    }

    if {$brightest_idx < 0} { set brightest_idx 0 }
    set new_num [expr {$max_number + 1}]

    # Compute merged values
    # Pass 1: flux-weighted centroid for X,Y + total flux/npix
    set total_flux 0.0
    set wx 0.0
    set wy 0.0
    set total_npix 0

    foreach row $merge_rows {
	set flux 1.0
	if {$idx_flux >= 0} {
	    set fv [string trim [lindex $row $idx_flux]]
	    if {[string is double -strict $fv] && $fv > 0} { set flux $fv }
	}
	set x [string trim [lindex $row $idx_x]]
	set y [string trim [lindex $row $idx_y]]
	if {![string is double -strict $x]} { set x 0 }
	if {![string is double -strict $y]} { set y 0 }

	set total_flux [expr {$total_flux + $flux}]
	set wx [expr {$wx + $x * $flux}]
	set wy [expr {$wy + $y * $flux}]

	if {$idx_npix >= 0} {
	    set nv [string trim [lindex $row $idx_npix]]
	    if {[string is integer -strict $nv]} {
		set total_npix [expr {$total_npix + $nv}]
	    }
	}
    }

    if {$total_flux <= 0} { set total_flux 1.0 }

    set new_x [expr {$wx / $total_flux}]
    set new_y [expr {$wy / $total_flux}]

    # Pass 2: second-moment tensor for merged ellipse (A, B, THETA)
    set Ixx 0.0
    set Iyy 0.0
    set Ixy 0.0

    foreach row $merge_rows {
	set flux 1.0
	if {$idx_flux >= 0} {
	    set fv [string trim [lindex $row $idx_flux]]
	    if {[string is double -strict $fv] && $fv > 0} { set flux $fv }
	}
	set x [string trim [lindex $row $idx_x]]
	set y [string trim [lindex $row $idx_y]]
	if {![string is double -strict $x]} { set x 0 }
	if {![string is double -strict $y]} { set y 0 }

	set ak 0.0
	set bk 0.0
	set thetak 0.0
	if {$idx_a >= 0} {
	    set av [string trim [lindex $row $idx_a]]
	    if {[string is double -strict $av]} { set ak $av }
	}
	if {$idx_b >= 0} {
	    set bv [string trim [lindex $row $idx_b]]
	    if {[string is double -strict $bv]} { set bk $bv }
	}
	if {$idx_theta >= 0} {
	    set tv [string trim [lindex $row $idx_theta]]
	    if {[string is double -strict $tv]} { set thetak $tv }
	}

	# Intrinsic second moments of this source's ellipse
	set rad [expr {$thetak * 3.14159265358979 / 180.0}]
	set cosT [expr {cos($rad)}]
	set sinT [expr {sin($rad)}]
	set a2 [expr {$ak * $ak}]
	set b2 [expr {$bk * $bk}]
	set ixx_k [expr {$a2 * $cosT * $cosT + $b2 * $sinT * $sinT}]
	set iyy_k [expr {$a2 * $sinT * $sinT + $b2 * $cosT * $cosT}]
	set ixy_k [expr {($a2 - $b2) * $sinT * $cosT}]

	# Parallel axis theorem: add offset from merged centroid
	set dx [expr {$x - $new_x}]
	set dy [expr {$y - $new_y}]
	set Ixx [expr {$Ixx + $flux * ($ixx_k + $dx * $dx)}]
	set Iyy [expr {$Iyy + $flux * ($iyy_k + $dy * $dy)}]
	set Ixy [expr {$Ixy + $flux * ($ixy_k + $dx * $dy)}]
    }

    # Normalize by total flux
    set Ixx [expr {$Ixx / $total_flux}]
    set Iyy [expr {$Iyy / $total_flux}]
    set Ixy [expr {$Ixy / $total_flux}]

    # Eigenvalue decomposition → A, B, THETA
    set trace [expr {$Ixx + $Iyy}]
    set det [expr {$Ixx * $Iyy - $Ixy * $Ixy}]
    set disc [expr {sqrt(abs(($Ixx - $Iyy) * ($Ixx - $Iyy) + 4.0 * $Ixy * $Ixy))}]
    set lam1 [expr {($trace + $disc) / 2.0}]
    set lam2 [expr {($trace - $disc) / 2.0}]
    if {$lam1 < 0} { set lam1 0.0 }
    if {$lam2 < 0} { set lam2 0.0 }
    set new_a [expr {sqrt($lam1)}]
    set new_b [expr {sqrt($lam2)}]
    set new_theta [expr {0.5 * atan2(2.0 * $Ixy, $Ixx - $Iyy) * 180.0 / 3.14159265358979}]

    # MAG_AUTO from total flux
    set mag_zp 25.0
    if {[info exists catpanel(param,mag-zeropoint)]} {
	set mag_zp $catpanel(param,mag-zeropoint)
    }
    set new_mag [expr {-2.5 * log10($total_flux) + $mag_zp}]

    # ISO_RADIUS from total NPIX
    set new_ir 5.0
    if {$total_npix > 0 && $new_a > 0 && $new_b > 0} {
	set ratio [expr {$new_b / $new_a}]
	if {$ratio <= 0} { set ratio 1.0 }
	set new_ir [expr {sqrt($total_npix / (3.14159265 * $ratio))}]
    }

    # Build merged row: copy from brightest, override computed fields
    set base_row [lindex $merge_rows $brightest_idx]
    set new_fields {}
    for {set c 0} {$c < $ncols} {incr c} {
	set val [string trim [lindex $base_row $c]]
	if {$c == $idx_num} { set val $new_num }
	if {$c == $idx_x} { set val [format "%.4f" $new_x] }
	if {$c == $idx_y} { set val [format "%.4f" $new_y] }
	if {$c == $idx_flux && $idx_flux >= 0} { set val [format "%.6g" $total_flux] }
	if {$c == $idx_mag && $idx_mag >= 0} { set val [format "%.4f" $new_mag] }
	if {$c == $idx_npix && $idx_npix >= 0} { set val $total_npix }
	if {$c == $idx_ir && $idx_ir >= 0} { set val [format "%.4f" $new_ir] }
	if {$c == $idx_a && $idx_a >= 0} { set val [format "%.4f" $new_a] }
	if {$c == $idx_b && $idx_b >= 0} { set val [format "%.4f" $new_b] }
	if {$c == $idx_theta && $idx_theta >= 0} { set val [format "%.4f" $new_theta] }
	lappend new_fields $val
    }
    set new_line [join $new_fields "\t"]

    OGFSessLog catalog.merge manual {} -tool native -requires catalog \
	-title "Merge sources [join $catpanel(merge,list) ,]" \
	-payload [dict create nums_list $catpanel(merge,list) mag_zp $mag_zp]
    # Rebuild alldata: header + other rows + merged row
    set newdata $header
    foreach row $other_rows {
	append newdata "\n$row"
    }
    append newdata "\n$new_line"
    set catpanel(alldata) $newdata

    # Clear merge state
    set nmerged [llength $catpanel(merge,list)]
    set catpanel(merge,list) {}
    set catpanel(merge,active) 0

    # Delete merge markers
    if {$current(frame) != {}} {
	catch {$current(frame) marker catalog sextract_merge delete}
    }

    # Reload table and markers
    CatalogPanelLoadTSV $catpanel(alldata) "merged"
    CatalogPanelCreateAllMarkers

    # Find merged source row and auto-select/navigate
    global $catpanel(tbldb)
    set ncols [$catpanel(tbl) cget -cols]
    set nrows [$catpanel(tbl) cget -rows]
    set col_num -1
    for {set c 1} {$c <= $ncols} {incr c} {
	if {[info exists ${catpanel(tbldb)}(0,$c)]} {
	    set hdr [set ${catpanel(tbldb)}(0,$c)]
	    if {$hdr eq "NUMBER"} {
		set col_num $c
		break
	    }
	}
    }
    set merged_row -1
    if {$col_num >= 0} {
	for {set r 1} {$r < $nrows} {incr r} {
	    if {[info exists ${catpanel(tbldb)}($r,$col_num)]} {
		set val [set ${catpanel(tbldb)}($r,$col_num)]
		if {$val eq $new_num} {
		    set merged_row $r
		    break
		}
	    }
	}
    }
    if {$merged_row >= 0} {
	$catpanel(tbl) selection set $merged_row,1
	$catpanel(tbl) see $merged_row,1
	CatalogPanelGotoSource $merged_row
    }

    set catpanel(status) "Merged $nmerged sources into #$new_num: pos=([format %.2f $new_x],[format %.2f $new_y]) mag=[format %.3f $new_mag]"
}

proc CatalogPanelMergeCancel {} {
    global catpanel
    global current

    # Delete all merge markers
    if {$current(frame) != {}} {
	catch {$current(frame) marker catalog sextract_merge delete}
    }

    set catpanel(merge,list) {}
    set catpanel(merge,active) 0
    set catpanel(status) "Merge cancelled"
}

proc CatalogPanelEscapeKey {} {
    global catpanel

    if {$catpanel(ai,active)} {
	CatalogPanelAIDone
	return
    }
    if {$catpanel(merge,active)} {
	CatalogPanelMergeCancel
    }
    if {[llength $catpanel(sel,nums)] > 0} {
	CatalogPanelClearSelection
	set catpanel(status) "Selection cleared"
    }
}

proc CatalogPanelSetLogScale {} {
    global current
    global scale

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return

    # Set scale to log
    set scale(type) log
    $current(frame) colorscale log $scale(log)
    $current(frame) colorscale log

    # Get min/max pixel values
    $current(frame) clip mode minmax
    set limits [$current(frame) get clip]
    set pmin [lindex $limits 0]
    set pmax [lindex $limits 1]

    # Guard against NaN/Inf/non-numeric
    if {[catch {expr {$pmin + 0.0}}]} { set pmin 0.001 }
    if {[catch {expr {$pmax + 0.0}}]} { set pmax 1000.0 }

    # Ensure positive values for log
    if {$pmin <= 0} { set pmin 0.001 }
    if {$pmax <= $pmin} { set pmax [expr {$pmin * 1000}] }

    # Compute display max: 0.8*(log(max)-log(min)) + log(min)
    set log_min [expr {log10($pmin)}]
    set log_max [expr {log10($pmax)}]
    set log_disp [expr {0.8 * ($log_max - $log_min) + $log_min}]
    set disp_max [expr {pow(10.0, $log_disp)}]

    # Apply user-defined limits
    set scale(min) $pmin
    set scale(max) $disp_max
    set scale(mode) user
    $current(frame) clip user $pmin $disp_max
    $current(frame) clip mode user
    UpdateScale
}

