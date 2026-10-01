# galaxy_model/galaxy_model.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/galaxy_model/plugin.json.

proc CatalogPanelGalaxyFit {model} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "No sources — run Extract first"
	return
    }

    set catpanel(status) "Galaxy $model fitting — not yet implemented"
}

proc CatalogPanelGalaxyParams {} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "No sources — run Extract first"
	return
    }

    set catpanel(status) "Galaxy parameter extraction — not yet implemented"
}

proc CatalogPanelGalaxyMorphology {} {
    global catpanel
    global current
    global ds9
    if {[info commands OGFAIBackendHook] ne {} && [OGFAIBackendHook morphology]} return  ;# ogf_ai.tcl: backend local|external

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "No sources — run Extract first"
	return
    }

    # Get current FITS filename
    set fn {}
    if {$current(frame) != {}} {
	catch {set fn [$current(frame) get fits file name full]}
    }
    set fn [string trim $fn "{}"]
    regsub {\[.*\]$} $fn {} fn
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    # Find ds9_galaxy_morph.py
    set bindir [file dirname [info nameofexecutable]]
    set script [file join $bindir ds9_galaxy_morph.py]
    if {![file exists $script]} {
	set libdir [file join [file dirname $bindir] ds9 library]
	set script [file join $libdir ds9_galaxy_morph.py]
    }
    if {![file exists $script]} {
	set catpanel(status) "ERROR: ds9_galaxy_morph.py not found"
	return
    }

    # Save catalog to temp TSV
    set catfile [file join [file normalize ~] .ds9 morph_catalog.tsv]
    catch {file mkdir [file dirname $catfile]}
    if {[catch {
	set fd [open $catfile w]
	puts $fd $catpanel(alldata)
	close $fd
    } err]} {
	set catpanel(status) "Morphology error: cannot write catalog: $err"
	return
    }

    # Build arguments
    set paramargs {}
    lappend paramargs "--catalog" $catfile

    # Find checkpoint
    set ckpt [file join [file dirname $bindir] galaxy_morph data checkpoints cnn_morph_best.pt]
    if {[file exists $ckpt]} {
	lappend paramargs "--checkpoint" $ckpt
    }

    set catpanel(status) "Morphology: classifying sources on [file tail $fn] ..."
    update idletasks

    # Run classification
    set errfile [file join [file normalize ~] .ds9 morph_stderr.txt]
    OGFSessLog galaxy.morphology auto [list [OGFPython] $script $fn {*}$paramargs] -title {Galaxy morphology (CNN)} -post [dict create kind morph]
    if {[catch {set data [exec [OGFPython] $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	catch {file delete $catfile}
	if {$stderr_msg ne ""} {
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    set catpanel(status) "Morphology error: $last_err"
	    puts "Morphology full stderr:\n$stderr_msg"
	} else {
	    set catpanel(status) "Morphology error: $err"
	}
	return
    }
    catch {file delete $errfile}
    catch {file delete $catfile}

    # Parse results
    CatalogPanelMorphParseResults $data
}

proc CatalogPanelMorphParseResults {data} {
    global catpanel
    global current

    set lines [split $data \n]
    set n_classified 0

    # Parse header: #GALAXY_MORPH	N_CLASSIFIED=245	N_SOURCES=300
    foreach line $lines {
	if {[string match "#GALAXY_MORPH*" $line]} {
	    foreach field [split $line "\t"] {
		if {[string match "N_CLASSIFIED=*" $field]} {
		    set n_classified [string range $field 13 end]
		}
	    }
	    continue
	}
    }

    # Build morph_map: NUMBER -> {morph_type morph_conf color}
    array unset catpanel morph,*
    set catpanel(morph,map) {}

    foreach line $lines {
	if {[string match "#*" $line]} continue
	if {[string match "NUMBER*" $line]} continue
	if {[string trim $line] eq {}} continue

	# NUMBER MORPH_TYPE MORPH_DESC MORPH_CONF TOP1_CLASS TOP1_PROB ... MORPH_COLOR
	set fields [split $line "\t"]
	if {[llength $fields] < 11} continue

	set src_num [lindex $fields 0]
	set morph_type [lindex $fields 1]
	set morph_desc [lindex $fields 2]
	set morph_conf [lindex $fields 3]
	set color [lindex $fields 10]

	set catpanel(morph,$src_num) [list $morph_type $morph_desc $morph_conf $color]
	lappend catpanel(morph,map) $src_num
    }

    if {[llength $catpanel(morph,map)] == 0} {
	set catpanel(status) "Morphology: no galaxies classified"
	return
    }

    # Add MORPH_TYPE and MORPH_CONF columns to alldata
    CatalogPanelMorphAddColumns

    # Recolor markers by morphology
    CatalogPanelMorphColorMarkers

    set catpanel(status) "Morphology: $n_classified galaxies classified"
}

proc CatalogPanelMorphAddColumns {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set lines [split $catpanel(alldata) \n]
    if {[llength $lines] < 2} return

    # Parse header
    set headers [split [lindex $lines 0] "\t"]
    set ncols [llength $headers]

    # Find NUMBER column
    set col_num -1
    for {set c 0} {$c < $ncols} {incr c} {
	if {[string trim [lindex $headers $c]] eq "NUMBER"} {
	    set col_num $c
	    break
	}
    }
    if {$col_num < 0} return

    # Check if columns already exist (avoid duplicates)
    set has_morph 0
    for {set c 0} {$c < $ncols} {incr c} {
	if {[string trim [lindex $headers $c]] eq "MORPH_TYPE"} {
	    set has_morph 1
	    break
	}
    }

    # Build new data with added columns
    set newdata {}

    if {$has_morph} {
	# Find existing MORPH_TYPE, MORPH_DESC, MORPH_CONF column indices
	set col_mt -1
	set col_md -1
	set col_mc -1
	for {set c 0} {$c < $ncols} {incr c} {
	    set h [string trim [lindex $headers $c]]
	    if {$h eq "MORPH_TYPE"} { set col_mt $c }
	    if {$h eq "MORPH_DESC"} { set col_md $c }
	    if {$h eq "MORPH_CONF"} { set col_mc $c }
	}

	# Update existing columns
	append newdata [lindex $lines 0]
	for {set i 1} {$i < [llength $lines]} {incr i} {
	    set line [lindex $lines $i]
	    if {[string trim $line] eq {}} continue
	    set fields [split $line "\t"]
	    set src_num [string trim [lindex $fields $col_num]]

	    if {[info exists catpanel(morph,$src_num)]} {
		set info $catpanel(morph,$src_num)
		if {$col_mt >= 0} {
		    lset fields $col_mt [lindex $info 0]
		}
		if {$col_md >= 0} {
		    lset fields $col_md [lindex $info 1]
		}
		if {$col_mc >= 0} {
		    lset fields $col_mc [lindex $info 2]
		}
	    }
	    append newdata "\n" [join $fields "\t"]
	}
    } else {
	# Add new columns
	append newdata [lindex $lines 0] "\tMORPH_TYPE\tMORPH_DESC\tMORPH_CONF"
	for {set i 1} {$i < [llength $lines]} {incr i} {
	    set line [lindex $lines $i]
	    if {[string trim $line] eq {}} continue
	    set fields [split $line "\t"]
	    set src_num [string trim [lindex $fields $col_num]]

	    set mt ""
	    set md ""
	    set mc ""
	    if {[info exists catpanel(morph,$src_num)]} {
		set info $catpanel(morph,$src_num)
		set mt [lindex $info 0]
		set md [lindex $info 1]
		set mc [lindex $info 2]
	    }
	    append newdata "\n" $line "\t" $mt "\t" $md "\t" $mc
	}
    }

    set catpanel(alldata) $newdata
    CatalogPanelLoadTSV $catpanel(alldata) "morphology"
}

proc CatalogPanelMorphColorMarkers {} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    set frame $current(frame)

    # Delete previous sextract_all markers and recreate with morph colors
    OGFMarkDelete $frame

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set lines [split $catpanel(alldata) \n]
    if {[llength $lines] < 2} return

    set headers [split [lindex $lines 0] "\t"]
    set ncols [llength $headers]

    # Find needed columns
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

    # Build region strings with morph-specific colors
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

	set src_num [expr {$i}]
	if {$col_num >= 0} {
	    set nv [string trim [lindex $fields $col_num]]
	    if {$nv ne {}} { set src_num $nv }
	}

	# Ellipse parameters
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

	set semi_a $iso_radius
	set semi_b $iso_radius
	if {$a_image > 0 && $b_image > 0} {
	    set semi_b [expr {$iso_radius * $b_image / $a_image}]
	}

	# Color: use morph color if classified, else yellow
	set color yellow
	if {[info exists catpanel(morph,$src_num)]} {
	    set color [lindex $catpanel(morph,$src_num) 3]
	}

	append reg "ellipse($x $y ${semi_a}i ${semi_b}i $theta) # color=$color width=1 tag={sextract_all} tag={sextract_src.$src_num} select=0 edit=0 move=0 rotate=0 delete=1 highlite=1 callback=highlite CatalogPanelMarkerCB {$src_num} callback=unhighlite CatalogPanelMarkerUnCB {$src_num}\n"
	incr count
	incr batch_count

	if {$batch_count >= $batch_size} {
	    set sextract_all_reg $reg
	    OGFMarkSend $frame
	    set reg "image\n"
	    set batch_count 0
	}
    }

    if {$batch_count > 0} {
	set sextract_all_reg $reg
	OGFMarkSend $frame
    }

    set catpanel(markall,on) 1
}

