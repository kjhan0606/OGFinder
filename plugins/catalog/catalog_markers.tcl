# catalog/catalog_markers.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/catalog/plugin.json.

# Build region string and create sextract_all markers from catpanel(alldata).
# This is the single source of truth for marker creation.
# Called by: CatalogPanelMarkAll, CatalogPanelMergeSources, AI merge, GotoSource.
proc CatalogPanelCreateAllMarkers {} {
    global current

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return
    if {![::ogf::cat::has]} return

    set frame $current(frame)

    # Delete previous sextract_all markers
    OGFMarkDelete $frame

    # Parse directly from alldata (authoritative data source)
    set lines [split [::ogf::cat::tsv] \n]
    if {[llength $lines] < 2} return

    set headers [split [lindex $lines 0] "\t"]
    set ncols [llength $headers]

    # Find column indices (0-based in tab-split fields)
    set col_x -1
    set col_y -1
    set col_a -1
    set col_b -1
    set col_theta -1
    set col_ir -1
    set col_num -1
    for {set c 0} {$c < $ncols} {incr c} {
	set hdr [string trim [lindex $headers $c]]
	switch -- $hdr {
	    NUMBER      { set col_num $c }
	    X_IMAGE     { set col_x $c }
	    Y_IMAGE     { set col_y $c }
	    A_IMAGE     { set col_a $c }
	    B_IMAGE     { set col_b $c }
	    THETA_IMAGE { set col_theta $c }
	    ISO_RADIUS  { set col_ir $c }
	}
    }
    if {$col_x < 0 || $col_y < 0} return

    # Build region strings in batches to avoid DS9 marker command size limits
    set batch_size 500
    set reg "image\n"
    set count 0
    set batch_count 0
    global sextract_all_reg

    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	set fields [split $line "\t"]

	set x [string trim [lindex $fields $col_x]]
	set y [string trim [lindex $fields $col_y]]
	if {![string is double -strict $x] || ![string is double -strict $y]} continue

	# Get source NUMBER for individual tag
	set src_num [expr {$i}]
	if {$col_num >= 0} {
	    set nv [string trim [lindex $fields $col_num]]
	    if {$nv ne {}} { set src_num $nv }
	}

	# Get ellipse parameters (NaN/Inf safe via catch)
	set iso_radius 5.0
	set a_image 0
	set b_image 0
	set theta 0

	if {$col_ir >= 0} {
	    set val [string trim [lindex $fields $col_ir]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} {
		set iso_radius $v
	    }
	}
	if {$col_a >= 0} {
	    set val [string trim [lindex $fields $col_a]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set a_image $v }
	}
	if {$col_b >= 0} {
	    set val [string trim [lindex $fields $col_b]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set b_image $v }
	}
	if {$col_theta >= 0} {
	    set val [string trim [lindex $fields $col_theta]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0} { set theta $v }
	}

	# ISO_RADIUS as semi-major, scaled by B/A for semi-minor
	set semi_a $iso_radius
	set semi_b $iso_radius
	if {$a_image > 0 && $b_image > 0} {
	    set semi_b [expr {$iso_radius * $b_image / $a_image}]
	}

	append reg "ellipse($x $y ${semi_a}i ${semi_b}i $theta) # color=yellow width=1 tag={sextract_all} tag={sextract_src.$src_num} select=0 edit=0 move=0 rotate=0 delete=1 highlite=1 callback=highlite CatalogPanelMarkerCB {$src_num} callback=unhighlite CatalogPanelMarkerUnCB {$src_num}\n"
	incr count
	incr batch_count

	# Flush batch when limit reached
	if {$batch_count >= $batch_size} {
	    set sextract_all_reg $reg
	    OGFMarkSend $frame
	    set reg "image\n"
	    set batch_count 0
	}
    }

    # Flush remaining markers
    if {$batch_count > 0} {
	set sextract_all_reg $reg
	OGFMarkSend $frame
    }

    if {$count == 0} return

    ::ogf::cat::set markall,on 1
    ::ogf::cat::set status "Marked $count sources (yellow ellipses)"
}

proc CatalogPanelMarkAll {} {
    global current

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return
    if {![::ogf::cat::has]} return

    set frame $current(frame)

    # If already marked, clear first then re-mark
    OGFMarkDelete $frame
    CatalogPanelCreateAllMarkers
}

proc CatalogPanelClearMarkers {} {
    global current

    if {$current(frame) == {}} return

    set frame $current(frame)
    OGFMarkDelete $frame
    ::ogf::cat::set markall,on 0
    ::ogf::cat::set status "Markers cleared"
}

proc CatalogPanelMarkerCB {num_str id} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return

    global [::ogf::cat::get tbldb]

    # Find NUMBER column index
    set ncols [[::ogf::cat::get tbl] cget -cols]
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

    # Find table row matching this source NUMBER
    set nrows [[::ogf::cat::get tbl] cget -rows]
    set target_row -1

    if {$col_num >= 0} {
	for {set r 1} {$r < $nrows} {incr r} {
	    if {[info exists ${catpanel(tbldb)}($r,$col_num)]} {
		set val [set ${catpanel(tbldb)}($r,$col_num)]
		if {$val eq $num_str} {
		    set target_row $r
		    break
		}
	    }
	}
    }

    if {$target_row < 0} return

    # Select, scroll, summary, markers (shared link path)
    CatalogPanelLinkSelect $num_str replace 1
}

proc CatalogPanelMarkerUnCB {num_str id} {
    # no-op
}

# Click handler called from Button1Frame in none mode
# Find the smallest ellipse source at canvas coordinate (cx, cy).
# Converts canvas→image coords via the frame, then tests all ellipses.
# When ellipses overlap, returns the source NUMBER with the smallest area.
# Returns "" if no ellipse contains the point.
proc CatalogPanelSmallestEllipseAt {frame cx cy} {

    if {![::ogf::cat::has]} { return {} }

    # Convert canvas coords to 1-indexed image coords
    set imgcoord [$frame get coordinates $cx $cy image]
    set imgx [lindex $imgcoord 0]
    set imgy [lindex $imgcoord 1]

    set lines [split [::ogf::cat::tsv] \n]
    if {[llength $lines] < 2} { return {} }

    set headers [split [lindex $lines 0] "\t"]
    set col_x -1; set col_y -1; set col_a -1; set col_b -1
    set col_theta -1; set col_ir -1; set col_num -1
    for {set c 0} {$c < [llength $headers]} {incr c} {
	switch -- [string trim [lindex $headers $c]] {
	    NUMBER      { set col_num $c }
	    X_IMAGE     { set col_x $c }
	    Y_IMAGE     { set col_y $c }
	    A_IMAGE     { set col_a $c }
	    B_IMAGE     { set col_b $c }
	    THETA_IMAGE { set col_theta $c }
	    ISO_RADIUS  { set col_ir $c }
	}
    }
    if {$col_x < 0 || $col_y < 0} { return {} }

    set best_num {}
    set best_area 1e30
    set pi 3.141592653589793

    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	set fields [split $line "\t"]

	# X_IMAGE, Y_IMAGE are 1-indexed image coords
	set sx [string trim [lindex $fields $col_x]]
	set sy [string trim [lindex $fields $col_y]]
	if {![string is double -strict $sx] || ![string is double -strict $sy]} continue

	set src_num $i
	if {$col_num >= 0} {
	    set nv [string trim [lindex $fields $col_num]]
	    if {$nv ne {}} { set src_num $nv }
	}

	# Reconstruct marker ellipse (same logic as CreateAllMarkers)
	set iso_radius 5.0
	set a_image 0; set b_image 0; set theta_deg 0
	if {$col_ir >= 0} {
	    set val [string trim [lindex $fields $col_ir]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set iso_radius $v }
	}
	if {$col_a >= 0} {
	    set val [string trim [lindex $fields $col_a]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set a_image $v }
	}
	if {$col_b >= 0} {
	    set val [string trim [lindex $fields $col_b]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set b_image $v }
	}
	if {$col_theta >= 0} {
	    set val [string trim [lindex $fields $col_theta]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0} { set theta_deg $v }
	}

	set semi_a $iso_radius
	set semi_b $iso_radius
	if {$a_image > 0 && $b_image > 0} {
	    set semi_b [expr {$iso_radius * $b_image / $a_image}]
	}

	# Point-in-ellipse test: rotate (dx,dy) into ellipse frame
	set dx [expr {$imgx - $sx}]
	set dy [expr {$imgy - $sy}]
	set theta_rad [expr {$theta_deg * $pi / 180.0}]
	set cosT [expr {cos($theta_rad)}]
	set sinT [expr {sin($theta_rad)}]
	set rx [expr { $cosT * $dx + $sinT * $dy}]
	set ry [expr {-$sinT * $dx + $cosT * $dy}]

	if {$semi_a <= 0 || $semi_b <= 0} continue
	set t [expr {($rx * $rx) / ($semi_a * $semi_a) + ($ry * $ry) / ($semi_b * $semi_b)}]

	if {$t <= 1.0} {
	    set area [expr {$semi_a * $semi_b}]
	    if {$area < $best_area} {
		set best_area $area
		set best_num $src_num
	    }
	}
    }

    return $best_num
}

proc CatalogPanelMarkerClick {which x y} {

    if {![::ogf::cat::exists tbl]} return
    # time-domain markers (moving / transient / detection): select the table row of that object
    if {![catch {$which get marker catalog id $x $y} _mid] && $_mid != 0} {
	set _tdk [OGFTDKeyFromTags [$which get marker catalog $_mid tag]]
	if {$_tdk ne {}} {CatalogPanelLinkSelect $_tdk replace 1; return}
    }
    if {![::ogf::cat::has]} return
    if {![$which has fits]} return

    # Quick test: is there any marker at this canvas position?
    set id [$which get marker catalog id $x $y]
    if {$id == 0} return

    # Among all overlapping ellipses, pick the smallest one
    set src_num [CatalogPanelSmallestEllipseAt $which $x $y]
    if {$src_num eq {}} {
	# Fallback: use the marker DS9 picked (original behaviour)
	set tags [$which get marker catalog $id tag]
	foreach tag $tags {
	    if {[string match "sextract_src.*" $tag]} {
		set src_num [string range $tag 13 end]
		break
	    }
	}
	if {$src_num eq {}} return
    }

    CatalogPanelMarkerCB $src_num $id
}

# Ctrl+Click handler called from ControlButton1Frame in none mode
proc CatalogPanelMarkerCtrlClick {which x y} {

    if {![$which has fits]} return

    # Quick test: is there any marker at this canvas position?
    set id [$which get marker catalog id $x $y]
    if {$id == 0} return

    # Check if this is a psf_star marker — if so, remove it
    set tags [$which get marker catalog $id tag]
    foreach tag $tags {
	if {[string match "psf_star.*" $tag]} {
	    set star_num [string range $tag 9 end]
	    # Delete the marker
	    catch {$which marker catalog tag $tag delete}
	    # Remove from star_indices list
	    if {[::ogf::cat::exists psf,star_indices]} {
		set idx [lsearch -exact [::ogf::cat::get psf,star_indices] $star_num]
		if {$idx >= 0} {
		    ::ogf::cat::set psf,star_indices [lreplace [::ogf::cat::get psf,star_indices] $idx $idx]
		}
		::ogf::cat::set status "Removed star $star_num ([llength [::ogf::cat::get psf,star_indices]] stars remaining)"
	    }
	    return
	}
    }

    # Not a star marker — proceed with source merge selection
    if {![::ogf::cat::exists tbl]} return
    if {![::ogf::cat::has]} return

    # Among all overlapping ellipses, pick the smallest one
    set src_num [CatalogPanelSmallestEllipseAt $which $x $y]
    if {$src_num eq {}} {
	# Fallback: use the marker DS9 picked (original behaviour)
	foreach tag $tags {
	    if {[string match "sextract_src.*" $tag]} {
		set src_num [string range $tag 13 end]
		break
	    }
	}
	if {$src_num eq {}} return
    }

    CatalogPanelCtrlSelect $src_num
}

proc CatalogPanelShowVisible {} {
    global catpanel
    global current
    global ds9

    if {![::ogf::cat::has]} return
    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return

    # checkbutton already toggled catpanel(visible_mode) before calling us
    if {![::ogf::cat::get visible_mode]} {
	CatalogPanelLoadTSV [::ogf::cat::tsv] "all"
	::ogf::cat::set status "Showing all sources"
	return
    }

    set frame $current(frame)

    # Get viewport center in image coordinates
    set cursor [$frame get cursor image]
    set cx [lindex $cursor 0]
    set cy [lindex $cursor 1]

    # Get zoom level
    set zoom [$frame get zoom]
    set zx [lindex $zoom 0]
    set zy [lindex $zoom 1]

    # Get canvas size
    set cw [winfo width $ds9(canvas)]
    set ch [winfo height $ds9(canvas)]

    # Compute viewport bounds in image coordinates
    set x_min [expr {$cx - $cw / 2.0 / $zx}]
    set x_max [expr {$cx + $cw / 2.0 / $zx}]
    set y_min [expr {$cy - $ch / 2.0 / $zy}]
    set y_max [expr {$cy + $ch / 2.0 / $zy}]

    # Parse alldata, find X_IMAGE/Y_IMAGE columns
    set lines [split [::ogf::cat::tsv] \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]
    set ncols [llength $headers]

    set idx_x -1
    set idx_y -1
    for {set i 0} {$i < $ncols} {incr i} {
	set h [string trim [lindex $headers $i]]
	if {$h eq "X_IMAGE"} { set idx_x $i }
	if {$h eq "Y_IMAGE"} { set idx_y $i }
    }
    if {$idx_x < 0 || $idx_y < 0} return

    # Filter rows within viewport
    set filtered $header
    set count 0
    set total 0
    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	incr total
	set fields [split $line "\t"]
	set x [string trim [lindex $fields $idx_x]]
	set y [string trim [lindex $fields $idx_y]]
	if {![string is double -strict $x] || ![string is double -strict $y]} continue
	if {$x >= $x_min && $x <= $x_max && $y >= $y_min && $y <= $y_max} {
	    append filtered "\n$line"
	    incr count
	}
    }

    CatalogPanelLoadTSV $filtered "visible"
    ::ogf::cat::set status "Visible: $count of $total sources in current view"
}

proc CatalogPanelCtrlSelect {src_num} {
    global current

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return

    set frame $current(frame)

    # Toggle: if already in list, remove; otherwise add
    set idx [lsearch -exact [::ogf::cat::get merge,list] $src_num]
    if {$idx >= 0} {
	# Remove from merge list
	::ogf::cat::set merge,list [lreplace [::ogf::cat::get merge,list] $idx $idx]
	# Delete this source's merge marker
	catch {$frame marker catalog sextract_merge.$src_num delete}
    } else {
	# Add to merge list
	::ogf::cat::lappend merge,list $src_num

	# Find source position from alldata
	set lines [split [::ogf::cat::tsv] \n]
	set header [lindex $lines 0]
	set headers [split $header "\t"]
	set ncols [llength $headers]

	set idx_num -1
	set idx_x -1
	set idx_y -1
	set idx_a -1
	set idx_b -1
	set idx_theta -1
	set idx_ir -1
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
	    }
	}

	# Find the matching line
	for {set i 1} {$i < [llength $lines]} {incr i} {
	    set line [lindex $lines $i]
	    if {[string trim $line] eq {}} continue
	    set fields [split $line "\t"]
	    set num_val [string trim [lindex $fields $idx_num]]
	    if {$num_val eq $src_num} {
		set x [string trim [lindex $fields $idx_x]]
		set y [string trim [lindex $fields $idx_y]]

		# Get ellipse params
		set iso_radius 5.0
		set a_image 0
		set b_image 0
		set theta 0
		if {$idx_ir >= 0} {
		    set val [string trim [lindex $fields $idx_ir]]
		    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} {
			set iso_radius $v
		    }
		}
		if {$idx_a >= 0} {
		    set val [string trim [lindex $fields $idx_a]]
		    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set a_image $v }
		}
		if {$idx_b >= 0} {
		    set val [string trim [lindex $fields $idx_b]]
		    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set b_image $v }
		}
		if {$idx_theta >= 0} {
		    set val [string trim [lindex $fields $idx_theta]]
		    if {[catch {set v [expr {$val + 0.0}]}] == 0} { set theta $v }
		}

		set semi_a $iso_radius
		set semi_b $iso_radius
		if {$a_image > 0 && $b_image > 0} {
		    set semi_b [expr {$iso_radius * $b_image / $a_image}]
		}

		# Create red thick merge marker
		global sextract_merge_reg
		set sextract_merge_reg "image\nellipse($x $y ${semi_a}i ${semi_b}i $theta) # color=red width=3 tag={sextract_merge} tag={sextract_merge.$src_num} select=0 edit=0 move=0 rotate=0 delete=1\n"
		catch {$frame marker catalog command ds9 var sextract_merge_reg}
		break
	    }
	}
    }

    ::ogf::cat::set merge,active 1
    set n [llength [::ogf::cat::get merge,list]]
    if {$n == 0} {
	::ogf::cat::set merge,active 0
	::ogf::cat::set status "Merge selection cleared"
    } else {
	::ogf::cat::set status "Merge: $n sources selected (Ctrl+M to merge, Esc to cancel)"
    }
}

