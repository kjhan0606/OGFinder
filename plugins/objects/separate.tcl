# objects/separate.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/objects/plugin.json.

proc CatalogPanelSeparateGetScript {} {
    set bindir [file dirname [info nameofexecutable]]
    set script [file join $bindir ds9_separate.py]
    if {![file exists $script]} {
	set libdir [file join [file dirname $bindir] ds9 library]
	set script [file join $libdir ds9_separate.py]
    }
    return $script
}

# Get the currently selected source NUMBER from the table
proc CatalogPanelGetSelectedSource {} {
    global catpanel

    if {![::ogf::cat::exists tbl]} { return {} }
    if {![::ogf::cat::has]} { return {} }

    # Get selected row
    set sel [[::ogf::cat::get tbl] curselection]
    if {$sel eq {}} { return {} }
    # sel is "row,col" — extract row
    set row [lindex [split [lindex $sel 0] ","] 0]
    if {$row <= 0} { return {} }

    global [::ogf::cat::get tbldb]
    set ncols [[::ogf::cat::get tbl] cget -cols]

    # Find column indices
    set col_num -1; set col_x -1; set col_y -1
    set col_a -1; set col_b -1; set col_theta -1; set col_ir -1
    for {set c 1} {$c <= $ncols} {incr c} {
	if {[info exists ${catpanel(tbldb)}(0,$c)]} {
	    switch -- [set ${catpanel(tbldb)}(0,$c)] {
		NUMBER      { set col_num $c }
		X_IMAGE     { set col_x $c }
		Y_IMAGE     { set col_y $c }
		A_IMAGE     { set col_a $c }
		B_IMAGE     { set col_b $c }
		THETA_IMAGE { set col_theta $c }
		ISO_RADIUS  { set col_ir $c }
	    }
	}
    }
    if {$col_x < 0 || $col_y < 0} { return {} }

    set result [dict create]
    dict set result row $row

    if {$col_num >= 0 && [info exists ${catpanel(tbldb)}($row,$col_num)]} {
	dict set result number [set ${catpanel(tbldb)}($row,$col_num)]
    } else {
	dict set result number $row
    }
    dict set result x [set ${catpanel(tbldb)}($row,$col_x)]
    dict set result y [set ${catpanel(tbldb)}($row,$col_y)]

    set a 10.0; set b 10.0; set theta 0.0; set ir 10.0
    if {$col_a >= 0 && [info exists ${catpanel(tbldb)}($row,$col_a)]} {
	set val [set ${catpanel(tbldb)}($row,$col_a)]
	if {[string is double -strict $val] && $val > 0} { set a $val }
    }
    if {$col_b >= 0 && [info exists ${catpanel(tbldb)}($row,$col_b)]} {
	set val [set ${catpanel(tbldb)}($row,$col_b)]
	if {[string is double -strict $val] && $val > 0} { set b $val }
    }
    if {$col_theta >= 0 && [info exists ${catpanel(tbldb)}($row,$col_theta)]} {
	set val [set ${catpanel(tbldb)}($row,$col_theta)]
	if {[string is double -strict $val]} { set theta $val }
    }
    if {$col_ir >= 0 && [info exists ${catpanel(tbldb)}($row,$col_ir)]} {
	set val [set ${catpanel(tbldb)}($row,$col_ir)]
	if {[string is double -strict $val] && $val > 0} { set ir $val }
    }
    dict set result a $a
    dict set result b $b
    dict set result theta $theta
    dict set result iso_radius $ir

    return $result
}

proc CatalogPanelDeleteSelected {} {

    if {![::ogf::cat::has]} {
	::ogf::cat::set status "No catalog data"
	return
    }

    # Get selected source info
    set src [CatalogPanelGetSelectedSource]
    if {$src eq {}} {
	::ogf::cat::set status "No source selected — click a source first"
	return
    }

    set src_num [dict get $src number]

    # Remove the row from alldata
    set lines [split [::ogf::cat::tsv] \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]

    # Find NUMBER column index
    set num_idx -1
    for {set i 0} {$i < [llength $headers]} {incr i} {
	if {[string trim [lindex $headers $i]] eq "NUMBER"} {
	    set num_idx $i
	    break
	}
    }

    set new_lines [list $header]
    set deleted 0
    foreach line [lrange $lines 1 end] {
	if {$line eq {}} continue
	set fields [split $line "\t"]

	set this_num ""
	if {$num_idx >= 0 && [llength $fields] > $num_idx} {
	    set this_num [string trim [lindex $fields $num_idx]]
	}

	if {$this_num eq [string trim $src_num]} {
	    set deleted 1
	} else {
	    lappend new_lines $line
	}
    }

    if {!$deleted} {
	::ogf::cat::set status "Source $src_num not found in catalog"
	return
    }

    OGFSessLog catalog.delete manual {} -tool native -requires catalog \
	-title "Delete source $src_num" -payload [dict create number $src_num]
    ::ogf::cat::set alldata [join $new_lines \n]

    # Delete the source marker
    global current
    if {$current(frame) ne {}} {
	catch {$current(frame) marker catalog sextract_src.$src_num delete}
	catch {$current(frame) marker catalog sextract_sel delete}
    }

    # Reload table and refresh markers
    CatalogPanelLoadTSV [::ogf::cat::tsv] "deleted"
    if {[::ogf::cat::exists markall,on] && [::ogf::cat::get markall,on]} {
	CatalogPanelCreateAllMarkers
    }

    ::ogf::cat::set status "Deleted source $src_num"
}

proc CatalogPanelSeparateSelected {} {
    global current

    if {$current(frame) eq {}} return

    set src [CatalogPanelGetSelectedSource]
    if {$src eq {}} {
	::ogf::cat::set status "Select a source first (click on an ellipse)"
	return
    }

    # Get FITS filename
    set fn {}
    if {[$current(frame) has fits]} {
	set fn [$current(frame) get fits file name full]
    }
    if {$fn eq {} || ![file exists $fn]} {
	::ogf::cat::set status "No FITS image loaded"
	return
    }

    set script [CatalogPanelSeparateGetScript]
    if {![file exists $script]} {
	::ogf::cat::set status "ERROR: ds9_separate.py not found"
	return
    }

    set src_num [dict get $src number]
    set src_x [dict get $src x]
    set src_y [dict get $src y]
    set src_a [dict get $src a]
    set src_b [dict get $src b]
    set src_theta [dict get $src theta]
    set src_ir [dict get $src iso_radius]

    ::ogf::cat::set status "Separating source $src_num ..."
    update idletasks

    set paramargs {}
    lappend paramargs "--x" $src_x
    lappend paramargs "--y" $src_y
    lappend paramargs "--a" $src_a
    lappend paramargs "--b" $src_b
    lappend paramargs "--theta" $src_theta
    lappend paramargs "--iso-radius" $src_ir
    lappend paramargs "--parent-number" $src_num
    lappend paramargs "--deblend-nthresh" [::ogf::cat::get param,sep-deblend-nthresh]
    lappend paramargs "--deblend-mincont" [::ogf::cat::get param,sep-deblend-mincont]
    lappend paramargs "--detect-thresh" [::ogf::cat::get param,sep-detect-thresh]
    lappend paramargs "--detect-minarea" [::ogf::cat::get param,sep-detect-minarea]
    lappend paramargs "--radius-factor" [::ogf::cat::get param,sep-radius-factor]
    lappend paramargs "--back-size" [::ogf::cat::get param,sep-back-size]

    OGFSessLog catalog.separate manual {} -tool native -title "Separate source $src_num (manual deblend)" \
	-note "hand-picked source; result is merged into the catalog by GUI code that is not replayable from the shell"
    set errfile [file join [file normalize ~] .ds9 separate_stderr.txt]
    if {[catch {set data [exec [OGFPython] $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	if {$stderr_msg ne ""} {
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    ::ogf::cat::set status "Separate error: $last_err"
	} else {
	    ::ogf::cat::set status "Separate error: $err"
	}
	return
    }
    catch {file delete $errfile}

    # Parse output
    set lines [split $data \n]
    set n_sub 0
    set sub_lines {}

    foreach line $lines {
	if {[string match "#SEPARATE*" $line]} {
	    foreach field [split $line "\t"] {
		if {[string match "N_SUB=*" $field]} {
		    set n_sub [string range $field 6 end]
		}
	    }
	    continue
	}
	if {[string match "NUMBER*" $line]} continue
	if {$line eq {}} continue
	lappend sub_lines $line
    }

    if {$n_sub == 0 || [llength $sub_lines] < 2} {
	::ogf::cat::set status "No sub-components found for source $src_num"
	return
    }

    # Replace the parent source in alldata with sub-sources
    CatalogPanelSeparateReplace $src_num $sub_lines

    ::ogf::cat::set status "Source $src_num separated into $n_sub sub-components"
}

proc CatalogPanelSeparateReplace {parent_num sub_lines} {

    if {![::ogf::cat::has]} return

    set all_lines [split [::ogf::cat::tsv] \n]
    set header [lindex $all_lines 0]
    set headers [split $header "\t"]

    # Find NUMBER column index
    set num_idx -1
    for {set i 0} {$i < [llength $headers]} {incr i} {
	if {[string trim [lindex $headers $i]] eq "NUMBER"} {
	    set num_idx $i
	    break
	}
    }

    # Sub-source output has fixed columns; map to parent catalog columns
    set sub_cols {NUMBER X_IMAGE Y_IMAGE A_IMAGE B_IMAGE THETA_IMAGE ISO_RADIUS FLUX_AUTO MAG_AUTO FLAGS}

    # Build column index map: sub_col_name → parent_col_index
    set col_map {}
    for {set si 0} {$si < [llength $sub_cols]} {incr si} {
	set scol [lindex $sub_cols $si]
	for {set pi 0} {$pi < [llength $headers]} {incr pi} {
	    if {[string trim [lindex $headers $pi]] eq $scol} {
		lappend col_map [list $si $pi]
		break
	    }
	}
    }

    # Build new alldata lines
    set new_lines [list $header]
    set ncols [llength $headers]

    foreach line [lrange $all_lines 1 end] {
	if {$line eq {}} continue
	set fields [split $line "\t"]

	set this_num ""
	if {$num_idx >= 0 && [llength $fields] > $num_idx} {
	    set this_num [string trim [lindex $fields $num_idx]]
	}

	if {$this_num eq [string trim $parent_num]} {
	    # Replace parent with sub-sources
	    foreach sub_line $sub_lines {
		set sub_fields [split $sub_line "\t"]
		# Start with empty row matching parent column count
		set new_row {}
		for {set c 0} {$c < $ncols} {incr c} {
		    lappend new_row ""
		}
		# Fill in mapped columns from sub-source
		foreach mapping $col_map {
		    set si [lindex $mapping 0]
		    set pi [lindex $mapping 1]
		    if {$si < [llength $sub_fields]} {
			lset new_row $pi [lindex $sub_fields $si]
		    }
		}
		lappend new_lines [join $new_row "\t"]
	    }
	} else {
	    lappend new_lines $line
	}
    }

    ::ogf::cat::set alldata [join $new_lines \n]

    # Reload table and markers
    CatalogPanelLoadTSV [::ogf::cat::tsv] "separated"
    if {[::ogf::cat::exists markall,on] && [::ogf::cat::get markall,on]} {
	CatalogPanelCreateAllMarkers
    }
}

proc CatalogPanelSeparateSave {} {

    if {![::ogf::cat::has]} {
	::ogf::cat::set status "No catalog to save"
	return
    }

    set types {
	{{TSV Files} {.tsv .txt}}
	{{All Files} *}
    }
    set outfile [tk_getSaveFile -filetypes $types \
		     -title "Save Catalog" \
		     -initialfile "catalog_separated.tsv"]
    if {$outfile eq {}} return

    if {[catch {
	set fd [open $outfile w]
	puts -nonewline $fd [::ogf::cat::tsv]
	close $fd
    } err]} {
	::ogf::cat::set status "Save error: $err"
	return
    }
    ::ogf::cat::set status "Catalog saved to $outfile"
}

proc CatalogPanelSeparateLoad {} {

    set types {
	{{TSV Files} {.tsv .txt}}
	{{All Files} *}
    }
    set infile [tk_getOpenFile -filetypes $types \
		    -title "Load Catalog"]
    if {$infile eq {}} return

    if {[catch {
	set fd [open $infile r]
	set data [read $fd]
	close $fd
    } err]} {
	::ogf::cat::set status "Load error: $err"
	return
    }

    ::ogf::cat::set alldata [string trim $data]
    CatalogPanelLoadTSV [::ogf::cat::tsv] "loaded"
    if {[::ogf::cat::exists markall,on] && [::ogf::cat::get markall,on]} {
	CatalogPanelCreateAllMarkers
    }
    ::ogf::cat::set status "Catalog loaded from $infile"
}

proc CatalogPanelAddObjectGetScript {} {
    set bindir [file dirname [info nameofexecutable]]
    set script [file join $bindir ds9_add_source.py]
    if {![file exists $script]} {
	set libdir [file join [file dirname $bindir] ds9 library]
	set script [file join $libdir ds9_add_source.py]
    }
    return $script
}

proc CatalogPanelAddObjectAtPosition {which imgx imgy} {
    global current

    # Check if add objects mode is enabled
    if {![::ogf::cat::exists add_objects_mode] || ![::ogf::cat::get add_objects_mode]} {
	::ogf::cat::set status "Enable Add Objects in Display menu first"
	return
    }

    if {$which ne $current(frame) || $which eq {}} return
    if {![$which has fits]} return
    if {![::ogf::cat::has]} {
	::ogf::cat::set status "Extract sources first, then enable Add Objects"
	return
    }

    # Get FITS filename
    set fn [$which get fits file name full]
    if {$fn eq {} || ![file exists $fn]} {
	::ogf::cat::set status "No FITS image loaded"
	return
    }

    set script [CatalogPanelAddObjectGetScript]
    if {![file exists $script]} {
	::ogf::cat::set status "ERROR: ds9_add_source.py not found"
	return
    }

    # Find max NUMBER in current catalog
    set max_num 0
    set lines [split [::ogf::cat::tsv] \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]
    set num_idx -1
    for {set i 0} {$i < [llength $headers]} {incr i} {
	if {[string trim [lindex $headers $i]] eq "NUMBER"} {
	    set num_idx $i
	    break
	}
    }
    if {$num_idx >= 0} {
	foreach line [lrange $lines 1 end] {
	    if {$line eq {}} continue
	    set fields [split $line "\t"]
	    set nv [string trim [lindex $fields $num_idx]]
	    # Handle fractional numbers like "42.1" from separate
	    set int_part [lindex [split $nv "."] 0]
	    if {[string is integer -strict $int_part] && $int_part > $max_num} {
		set max_num $int_part
	    }
	}
    }
    set new_num [expr {$max_num + 1}]

    ::ogf::cat::set status "Detecting source at ($imgx, $imgy) ..."
    update idletasks

    set paramargs {}
    lappend paramargs "--x" $imgx
    lappend paramargs "--y" $imgy
    lappend paramargs "--number" $new_num
    lappend paramargs "--detect-thresh" [::ogf::cat::get param,detect-thresh]
    lappend paramargs "--detect-minarea" [::ogf::cat::get param,detect-minarea]
    lappend paramargs "--deblend-nthresh" [::ogf::cat::get param,deblend-nthresh]
    lappend paramargs "--deblend-mincont" [::ogf::cat::get param,deblend-mincont]
    lappend paramargs "--back-size" [::ogf::cat::get param,back-size]
    lappend paramargs "--mag-zeropoint" [::ogf::cat::get param,mag-zeropoint]
    lappend paramargs "--phot-aperture" [::ogf::cat::get param,phot-aperture]

    OGFSessLog catalog.add_object manual {} -tool native -title "Add object at ($imgx, $imgy) by hand" \
	-note "hand-picked position; not replayable from the shell"
    set errfile [file join [file normalize ~] .ds9 add_source_stderr.txt]
    if {[catch {set data [exec [OGFPython] $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	if {$stderr_msg ne ""} {
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    ::ogf::cat::set status "Add source error: $last_err"
	} else {
	    ::ogf::cat::set status "Add source error: $err"
	}
	return
    }
    catch {file delete $errfile}

    # Parse output
    set result_lines [split $data \n]
    set found 0
    set source_line {}

    foreach line $result_lines {
	if {[string match "#ADD_SOURCE*" $line]} {
	    foreach field [split $line "\t"] {
		if {[string match "FOUND=*" $field]} {
		    set found [string range $field 6 end]
		}
	    }
	    continue
	}
	# Skip column header
	if {[string match "NUMBER*" $line]} continue
	if {$line eq {}} continue
	set source_line $line
    }

    if {$found == 0 || $source_line eq {}} {
	::ogf::cat::set status "No source detected at ($imgx, $imgy)"
	return
    }

    # Add the new source to alldata
    # Map output columns to existing catalog columns
    set out_cols {NUMBER X_IMAGE Y_IMAGE A_IMAGE B_IMAGE THETA_IMAGE ELLIPTICITY KRON_RADIUS ISO_RADIUS FLUX_AUTO FLUX_APER MAG_AUTO MAG_APER PEAK CLASS_STAR FLAGS FWHM_IMAGE}
    set src_fields [split $source_line "\t"]

    set ncols [llength $headers]
    set new_row {}
    for {set c 0} {$c < $ncols} {incr c} {
	lappend new_row ""
    }

    # Map each output column to the catalog column
    for {set si 0} {$si < [llength $out_cols]} {incr si} {
	set scol [lindex $out_cols $si]
	for {set pi 0} {$pi < $ncols} {incr pi} {
	    if {[string trim [lindex $headers $pi]] eq $scol} {
		if {$si < [llength $src_fields]} {
		    lset new_row $pi [lindex $src_fields $si]
		}
		break
	    }
	}
    }

    # Append to alldata
    ::ogf::cat::append alldata "\n" [join $new_row "\t"]

    # Reload table and markers
    CatalogPanelLoadTSV [::ogf::cat::tsv] "added"
    if {[::ogf::cat::exists markall,on] && [::ogf::cat::get markall,on]} {
	CatalogPanelCreateAllMarkers
    }

    ::ogf::cat::set status "Added source #$new_num at ($imgx, $imgy)"
}

