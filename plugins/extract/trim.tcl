# extract/trim.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/extract/plugin.json.

proc CatalogPanelTrimDialog {} {
    global ed

    if {![::ogf::cat::has]} {
	::ogf::cat::set status "No catalog data to trim"
	return
    }

    set w {.sextracttrim}

    set ed(ok) 0

    # Get column names from alldata header
    set lines [split [::ogf::cat::tsv] \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]
    set ncols [llength $headers]

    set ed(trim,cols) {}
    for {set i 0} {$i < $ncols} {incr i} {
	set colname [string trim [lindex $headers $i]]
	lappend ed(trim,cols) $colname
	# Initialize from existing trim values or empty
	if {[::ogf::cat::exists trim,$colname,min]} {
	    set ed(trim,$colname,min) [::ogf::cat::get trim,$colname,min]
	} else {
	    set ed(trim,$colname,min) {}
	}
	if {[::ogf::cat::exists trim,$colname,max]} {
	    set ed(trim,$colname,max) [::ogf::cat::get trim,$colname,max]
	} else {
	    set ed(trim,$colname,max) {}
	}
    }

    DialogCreate $w {Trim - Column Filter} ed(ok)

    # Scrollable frame for columns
    set sf [ttk::frame $w.param]
    set canvas_w [canvas $sf.c -width 400 -height 300 \
		      -yscrollcommand [list $sf.vs set]]
    ttk::scrollbar $sf.vs -orient vertical -command [list $canvas_w yview]
    set inner [ttk::frame $canvas_w.inner]
    $canvas_w create window 0 0 -anchor nw -window $inner

    # Header labels
    ttk::label $inner.hcol -text "Column" -font {Helvetica 10 bold} -width 15
    ttk::label $inner.hmin -text "Min" -font {Helvetica 10 bold} -width 12
    ttk::label $inner.htilde -text "" -width 2
    ttk::label $inner.hmax -text "Max" -font {Helvetica 10 bold} -width 12
    grid $inner.hcol $inner.hmin $inner.htilde $inner.hmax \
	-padx 2 -pady 2 -sticky w

    set row 1
    foreach colname $ed(trim,cols) {
	ttk::label $inner.l$row -text "$colname:" -anchor w -width 15
	ttk::entry $inner.emin$row -textvariable ed(trim,$colname,min) -width 12
	ttk::label $inner.tilde$row -text "~" -width 2
	ttk::entry $inner.emax$row -textvariable ed(trim,$colname,max) -width 12
	grid $inner.l$row $inner.emin$row $inner.tilde$row $inner.emax$row \
	    -padx 2 -pady 1 -sticky w
	incr row
    }

    # Update scroll region after layout
    bind $inner <Configure> [list $canvas_w configure -scrollregion \
				 [$canvas_w bbox all]]

    pack $canvas_w -side left -fill both -expand true
    pack $sf.vs -side right -fill y

    # Buttons
    set bf [ttk::frame $w.buttons]
    ttk::button $bf.ok -text {Apply} -command {set ed(ok) 1} -default active
    ttk::button $bf.cancel -text {Cancel} -command {set ed(ok) 0}
    ttk::button $bf.reset -text {Reset} -command {
	foreach col $ed(trim,cols) {
	    set ed(trim,$col,min) {}
	    set ed(trim,$col,max) {}
	}
    }
    ttk::button $bf.save -text {Save} -command {
	CatalogPanelTrimSaveFromEd
    }
    ttk::button $bf.load -text {Load} -command {
	CatalogPanelTrimLoadToEd
    }
    pack $bf.ok $bf.cancel $bf.reset $bf.save $bf.load \
	-side left -expand true -padx 2 -pady 4

    bind $w <Return> {set ed(ok) 1}

    # Fini
    ttk::separator $w.sep -orient horizontal
    pack $w.buttons $w.sep -side bottom -fill x
    pack $w.param -side top -fill both -expand true

    DialogWait $w ed(ok)
    destroy $w

    if {$ed(ok)} {
	# Copy trim values from ed to catpanel
	foreach colname $ed(trim,cols) {
	    ::ogf::cat::set trim,$colname,min $ed(trim,$colname,min)
	    ::ogf::cat::set trim,$colname,max $ed(trim,$colname,max)
	}
	CatalogPanelTrimApply
    }

    unset ed
}

proc CatalogPanelTrimSaveFromEd {} {
    global ed

    set prefdir [file join [file normalize ~] .ds9]
    if {![file isdirectory $prefdir]} {
	file mkdir $prefdir
    }
    set preffile [file join $prefdir sextract_trim.prf]
    if {[catch {set fd [open $preffile w]} err]} return

    foreach colname $ed(trim,cols) {
	puts $fd "$colname\t$ed(trim,$colname,min)\t$ed(trim,$colname,max)"
    }
    close $fd
}

proc CatalogPanelTrimLoadToEd {} {
    global ed

    set preffile [file join [file normalize ~] .ds9 sextract_trim.prf]
    if {![file exists $preffile]} return
    if {[catch {set fd [open $preffile r]} err]} return

    while {[gets $fd line] >= 0} {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set parts [split $line "\t"]
	if {[llength $parts] >= 3} {
	    set colname [lindex $parts 0]
	    set minval [lindex $parts 1]
	    set maxval [lindex $parts 2]
	    if {[lsearch -exact $ed(trim,cols) $colname] >= 0} {
		set ed(trim,$colname,min) $minval
		set ed(trim,$colname,max) $maxval
	    }
	}
    }
    close $fd
}

proc CatalogPanelTrimApply {} {
    global current

    if {![::ogf::cat::has]} return

    set lines [split [::ogf::cat::tsv] \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]
    set ncols [llength $headers]

    # Build list of active trim conditions
    set conditions {}
    for {set i 0} {$i < $ncols} {incr i} {
	set colname [string trim [lindex $headers $i]]
	set has_min 0
	set has_max 0
	set minval 0
	set maxval 0
	if {[::ogf::cat::exists trim,$colname,min] && [::ogf::cat::get trim,$colname,min] ne {}} {
	    if {[string is double -strict [::ogf::cat::get trim,$colname,min]]} {
		set has_min 1
		set minval [::ogf::cat::get trim,$colname,min]
	    }
	}
	if {[::ogf::cat::exists trim,$colname,max] && [::ogf::cat::get trim,$colname,max] ne {}} {
	    if {[string is double -strict [::ogf::cat::get trim,$colname,max]]} {
		set has_max 1
		set maxval [::ogf::cat::get trim,$colname,max]
	    }
	}
	if {$has_min || $has_max} {
	    lappend conditions [list $i $has_min $minval $has_max $maxval]
	}
    }

    # If no conditions, show all
    if {[llength $conditions] == 0} {
	::ogf::cat::set trim,active 0
	CatalogPanelLoadTSV [::ogf::cat::tsv] "all"
	::ogf::cat::set status "Trim cleared - showing all sources"
	return
    }

    # Filter rows
    set filtered $header
    set count 0
    set total 0
    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	incr total
	set fields [split $line "\t"]
	set pass 1

	foreach cond $conditions {
	    set cidx [lindex $cond 0]
	    set has_min [lindex $cond 1]
	    set minval [lindex $cond 2]
	    set has_max [lindex $cond 3]
	    set maxval [lindex $cond 4]

	    set val [string trim [lindex $fields $cidx]]
	    if {![string is double -strict $val]} {
		set pass 0
		break
	    }
	    if {$has_min && $val < $minval} {
		set pass 0
		break
	    }
	    if {$has_max && $val > $maxval} {
		set pass 0
		break
	    }
	}

	if {$pass} {
	    append filtered "\n$line"
	    incr count
	}
    }

    set jc {}
    foreach cond $conditions {
	lassign $cond cidx hmin minv hmax maxv
	set one "\"col\": [OGFJStr [string trim [lindex $headers $cidx]]]"
	if {$hmin} {append one ", \"min\": [OGFJStr $minv]"}
	if {$hmax} {append one ", \"max\": [OGFJStr $maxv]"}
	lappend jc "\{$one\}"
    }
    OGFSessLog catalog.trim auto {} -tool native -requires catalog \
	-title "Trim catalog: [llength $conditions] condition(s)" \
	-payload [dict create conditions_json "\[[join $jc {, }]\]"]
    ::ogf::cat::set trim,active 1
    CatalogPanelLoadTSV $filtered "trimmed"

    # Re-mark from authoritative data
    CatalogPanelCreateAllMarkers

    ::ogf::cat::set status "Trimmed: $count of $total sources match conditions"
}

proc CatalogPanelTrimSave {} {

    set prefdir [file join [file normalize ~] .ds9]
    if {![file isdirectory $prefdir]} {
	file mkdir $prefdir
    }
    set preffile [file join $prefdir sextract_trim.prf]
    if {[catch {set fd [open $preffile w]} err]} return

    set lines [split [::ogf::cat::tsv] \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]
    foreach h $headers {
	set colname [string trim $h]
	set minval {}
	set maxval {}
	if {[::ogf::cat::exists trim,$colname,min]} {
	    set minval [::ogf::cat::get trim,$colname,min]
	}
	if {[::ogf::cat::exists trim,$colname,max]} {
	    set maxval [::ogf::cat::get trim,$colname,max]
	}
	puts $fd "$colname\t$minval\t$maxval"
    }
    close $fd
}

proc CatalogPanelTrimLoad {} {

    set preffile [file join [file normalize ~] .ds9 sextract_trim.prf]
    if {![file exists $preffile]} return
    if {[catch {set fd [open $preffile r]} err]} return

    while {[gets $fd line] >= 0} {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set parts [split $line "\t"]
	if {[llength $parts] >= 3} {
	    set colname [lindex $parts 0]
	    ::ogf::cat::set trim,$colname,min [lindex $parts 1]
	    ::ogf::cat::set trim,$colname,max [lindex $parts 2]
	}
    }
    close $fd
}

