# star_psf/star_psf.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/star_psf/plugin.json.

proc CatalogPanelStarFinder {} {
    global current
    if {[info commands OGFAIBackendHook] ne {} && [OGFAIBackendHook star]} return  ;# ogf_ai.tcl: backend local|external

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return
    if {![::ogf::cat::has]} {
	::ogf::cat::set status "No sources — run Extract first"
	return
    }

    # Get FITS filename
    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	::ogf::cat::set status "No FITS image loaded"
	return
    }

    # Find script
    set script [CatalogPanelGetScript ds9_star_finder.py]
    if {![file exists $script]} {
	::ogf::cat::set status "ERROR: ds9_star_finder.py not found"
	return
    }

    # Save catalog to temp TSV
    set catfile [CatalogPanelSaveTempCatalog star]
    if {$catfile eq {}} {
	::ogf::cat::set status "Star Finder: cannot write catalog"
	return
    }

    # Build arguments
    set paramargs {}
    lappend paramargs "--catalog" $catfile

    # Pass PSF file if available
    if {[::ogf::cat::exists psf,file] && [::ogf::cat::get psf,file] ne {} &&
	[file exists [::ogf::cat::get psf,file]]} {
	lappend paramargs "--psf" [::ogf::cat::get psf,file]
    }

    ::ogf::cat::set status "AI Star Classification: classifying sources on [file tail $fn] ..."
    update idletasks

    # Run classification
    set errfile [file join [file normalize ~] .ds9 star_finder_stderr.txt]
    OGFSessLog stars.classify auto [list [OGFPython] $script $fn {*}$paramargs] -title {AI star classification} -post [dict create kind star]
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
	    ::ogf::cat::set status "Star Finder error: $last_err"
	    puts "Star Finder full stderr:\n$stderr_msg"
	} else {
	    ::ogf::cat::set status "Star Finder error: $err"
	}
	return
    }
    catch {file delete $errfile}
    catch {file delete $catfile}

    # Parse results and add columns
    CatalogPanelStarFinderParse $data
}

proc CatalogPanelStarFinderParse {data} {

    set lines [split $data \n]
    set n_classified 0

    # Parse header: #STAR_FINDER	N_CLASSIFIED=245	N_SOURCES=300
    foreach line $lines {
	if {[string match "#STAR_FINDER*" $line]} {
	    foreach field [split $line "\t"] {
		if {[string match "N_CLASSIFIED=*" $field]} {
		    set n_classified [string range $field 13 end]
		}
	    }
	    continue
	}
    }

    # Collect result lines (skip header/comments)
    set result_lines {}
    set header_line ""
    foreach line $lines {
	if {[string match "#*" $line]} continue
	if {[string match "NUMBER*" $line]} {
	    set header_line $line
	    continue
	}
	if {[string trim $line] eq {}} continue
	lappend result_lines $line
    }

    if {[llength $result_lines] == 0} {
	::ogf::cat::set status "Star Finder: no sources classified"
	return
    }

    # Rebuild TSV result for AddColumnsFromTSV
    set result_data $header_line
    foreach line $result_lines {
	append result_data "\n" $line
    }

    # Add AI_STAR and AI_STAR_CONF columns
    CatalogPanelAddColumnsFromTSV $result_data {AI_STAR AI_STAR_CONF}

    # Recolor markers
    CatalogPanelStarFinderColorMarkers $result_lines

    ::ogf::cat::set status "AI Star Classification: $n_classified sources classified"
}

proc CatalogPanelStarFinderColorMarkers {result_lines} {
    global current

    if {$current(frame) == {}} return
    set frame $current(frame)

    # Build color lookup: NUMBER -> color
    array set star_color {}
    foreach line $result_lines {
	set fields [split $line "\t"]
	if {[llength $fields] < 4} continue
	set src_num [lindex $fields 0]
	set color [lindex $fields 3]
	set star_color($src_num) $color
    }

    # Delete existing markers and recreate with star/galaxy colors
    OGFMarkDelete $frame

    if {![::ogf::cat::has]} return

    set lines [split [::ogf::cat::tsv] \n]
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

    # Build region strings with star/galaxy colors
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

	# Color: star=cyan, galaxy=yellow, unclassified=yellow
	set color yellow
	if {[info exists star_color($src_num)]} {
	    set color $star_color($src_num)
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

    ::ogf::cat::set markall,on 1
}

proc CatalogPanelPSFParamLoad {} {

    set preffile [file join [file normalize ~] .ds9 psf_deconv.prf]
    if {![file exists $preffile]} return
    if {[catch {set fd [open $preffile r]} err]} return
    while {[gets $fd line] >= 0} {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set parts [split $line]
	if {[llength $parts] >= 2} {
	    set key [lindex $parts 0]
	    set val [lindex $parts 1]
	    if {[::ogf::cat::exists psf,param,$key]} {
		::ogf::cat::set psf,param,$key $val
	    }
	}
    }
    close $fd
}

proc CatalogPanelPSFParamSave {} {

    set prefdir [file join [file normalize ~] .ds9]
    if {![file isdirectory $prefdir]} {
	file mkdir $prefdir
    }
    set preffile [file join $prefdir psf_deconv.prf]
    if {[catch {set fd [open $preffile w]} err]} return
    foreach pname {class-star-thresh max-ellipticity fwhm-sigma min-flux-snr \
		   psf-size rl-iterations wiener-nsr tikhonov-lambda tv-lambda \
		   clean-gain clean-niter clean-threshold mem-lambda mem-niter \
		   ext-core-mag-min ext-core-mag-max ext-wing-mag-max \
		   ext-core-size ext-wing-size ext-blend-inner ext-blend-outer \
		   ext-saturation-limit \
		   sim-telescope sim-instrument sim-filter sim-psf-size \
		   sim-oversample sim-jitter-sigma sim-focus-offset} {
	puts $fd "$pname [::ogf::cat::get psf,param,$pname]"
    }
    close $fd
}

proc CatalogPanelPSFGetScript {} {
    set bindir [file dirname [info nameofexecutable]]
    set script [file join $bindir ds9_psf_deconv.py]
    if {![file exists $script]} {
	set libdir [file join [file dirname $bindir] ds9 library]
	set script [file join $libdir ds9_psf_deconv.py]
    }
    return $script
}

proc CatalogPanelPSFGetFITS {} {
    global current

    set fn {}
    if {$current(frame) != {}} {
	catch {set fn [$current(frame) get fits file name full]}
    }
    set fn [string trim $fn "{}"]
    regsub {\[.*\]$} $fn {} fn
    return $fn
}

proc CatalogPanelFindStars {method} {
    global current

    if {![::ogf::cat::has]} {
	::ogf::cat::set status "Extract sources first before finding stars"
	return
    }

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	::ogf::cat::set status "No FITS image loaded"
	return
    }

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	::ogf::cat::set status "ERROR: ds9_psf_deconv.py not found"
	return
    }

    # Save catalog to temp TSV
    set catfile [file join [file normalize ~] .ds9 psf_catalog.tsv]
    catch {file mkdir [file dirname $catfile]}
    if {[catch {
	set fd [open $catfile w]
	puts $fd [::ogf::cat::tsv]
	close $fd
    } err]} {
	::ogf::cat::set status "Star finding error: cannot write catalog: $err"
	return
    }

    # Build arguments
    set paramargs {}
    lappend paramargs "--mode" "find_stars"
    lappend paramargs "--catalog" $catfile
    lappend paramargs "--method" $method
    lappend paramargs "--class-star-thresh" [::ogf::cat::get psf,param,class-star-thresh]
    lappend paramargs "--max-ellipticity" [::ogf::cat::get psf,param,max-ellipticity]
    lappend paramargs "--fwhm-sigma" [::ogf::cat::get psf,param,fwhm-sigma]
    lappend paramargs "--min-flux-snr" [::ogf::cat::get psf,param,min-flux-snr]

    ::ogf::cat::set status "Finding stars ($method) ..."
    update idletasks

    set errfile [file join [file normalize ~] .ds9 psf_stderr.txt]
    OGFSessLog psf.find_stars auto [list [OGFPython] $script $fn {*}$paramargs] -title {Find PSF stars ($method)} -post [dict create kind stars]
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
	    ::ogf::cat::set status "Star finding error: $last_err"
	    puts "Star finding stderr:\n$stderr_msg"
	} else {
	    ::ogf::cat::set status "Star finding error: $err"
	}
	return
    }
    catch {file delete $errfile}
    catch {file delete $catfile}

    # Parse output
    set lines [split $data \n]
    set star_indices {}
    set n_stars 0

    foreach line $lines {
	if {[string match "#PSF_STARS*" $line]} {
	    foreach field [split $line "\t"] {
		if {[string match "N_STARS=*" $field]} {
		    set n_stars [string range $field 8 end]
		}
	    }
	    continue
	}
	# Skip column header
	if {[string match "NUMBER*" $line]} continue
	if {$line eq {}} continue

	set fields [split $line "\t"]
	if {[llength $fields] >= 3} {
	    set num [lindex $fields 0]
	    lappend star_indices $num
	}
    }

    ::ogf::cat::set psf,star_indices $star_indices
    ::ogf::cat::set status "Found $n_stars stars ($method)"

    # Show star markers
    CatalogPanelShowStars
}

proc CatalogPanelShowStars {} {
    global current

    if {$current(frame) eq {}} return
    set frame $current(frame)

    # Clear existing star markers
    catch {$frame marker catalog tag psf_star delete}

    if {![::ogf::cat::exists psf,star_indices] || [::ogf::cat::get psf,star_indices] eq {}} {
	::ogf::cat::set status "No stars found — run Star Finding first"
	return
    }

    # Parse alldata to find star positions
    set lines [split [::ogf::cat::tsv] \n]
    if {[llength $lines] < 2} return

    set header [lindex $lines 0]
    set cols [split $header "\t"]
    set num_idx -1
    set x_idx -1
    set y_idx -1
    set a_idx -1
    set b_idx -1
    for {set i 0} {$i < [llength $cols]} {incr i} {
	set col [string trim [lindex $cols $i]]
	switch $col {
	    NUMBER {set num_idx $i}
	    X_IMAGE {set x_idx $i}
	    Y_IMAGE {set y_idx $i}
	    A_IMAGE {set a_idx $i}
	    B_IMAGE {set b_idx $i}
	}
    }
    if {$num_idx < 0 || $x_idx < 0 || $y_idx < 0} return

    set reg "image\n"
    set count 0
    foreach line [lrange $lines 1 end] {
	set fields [split $line "\t"]
	if {[llength $fields] <= $num_idx} continue
	set num [string trim [lindex $fields $num_idx]]
	if {[lsearch -exact [::ogf::cat::get psf,star_indices] $num] < 0} continue

	set x [string trim [lindex $fields $x_idx]]
	set y [string trim [lindex $fields $y_idx]]
	set a 5.0
	set b 5.0
	if {$a_idx >= 0} {set a [expr {max(3.0, [string trim [lindex $fields $a_idx]] * 2)}]}
	if {$b_idx >= 0} {set b [expr {max(3.0, [string trim [lindex $fields $b_idx]] * 2)}]}

	append reg "circle($x $y ${a}i) # color=purple width=2 dash=1 tag={psf_star} tag={psf_star.$num} select=0 edit=0 move=0 rotate=0 delete=1\n"
	incr count
    }

    if {$count > 0} {
	set psf_star_reg $reg
	catch {$frame marker catalog command ds9 var psf_star_reg}
	::ogf::cat::set status "Showing $count star markers"
    }
}

proc CatalogPanelClearStars {} {
    global current

    ::ogf::cat::set psf,star_indices {}
    if {$current(frame) ne {}} {
	catch {$current(frame) marker catalog tag psf_star delete}
    }
    ::ogf::cat::set status "Star markers cleared"
}

proc CatalogPanelBuildPSF {method} {
    global current

    if {![::ogf::cat::exists psf,star_indices] || [::ogf::cat::get psf,star_indices] eq {}} {
	::ogf::cat::set status "Find stars first before building PSF"
	return
    }

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	::ogf::cat::set status "No FITS image loaded"
	return
    }

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	::ogf::cat::set status "ERROR: ds9_psf_deconv.py not found"
	return
    }

    # Save catalog
    set catfile [file join [file normalize ~] .ds9 psf_catalog.tsv]
    catch {file mkdir [file dirname $catfile]}
    if {[catch {
	set fd [open $catfile w]
	puts $fd [::ogf::cat::tsv]
	close $fd
    } err]} {
	::ogf::cat::set status "PSF build error: cannot write catalog: $err"
	return
    }

    set star_list [join [::ogf::cat::get psf,star_indices] ","]

    set paramargs {}
    lappend paramargs "--mode" "build_psf"
    lappend paramargs "--catalog" $catfile
    lappend paramargs "--star-indices" $star_list
    lappend paramargs "--psf-method" $method
    lappend paramargs "--psf-size" [::ogf::cat::get psf,param,psf-size]
    lappend paramargs "--psf-output" [::ogf::cat::get psf,file]

    ::ogf::cat::set status "Building PSF ($method) from [llength [::ogf::cat::get psf,star_indices]] stars ..."
    update idletasks

    set errfile [file join [file normalize ~] .ds9 psf_stderr.txt]
    OGFSessLog psf.build auto [list [OGFPython] $script $fn {*}$paramargs] -title {Build PSF ($method)} -post [dict create kind psfbuilt] -requires [list stars]
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
	    ::ogf::cat::set status "PSF build error: $last_err"
	    puts "PSF build stderr:\n$stderr_msg"
	} else {
	    ::ogf::cat::set status "PSF build error: $err"
	}
	return
    }
    catch {file delete $errfile}
    catch {file delete $catfile}

    ::ogf::cat::set psf,has_psf 1

    # Parse info from output
    foreach line [split $data \n] {
	if {[string match "#PSF_BUILT*" $line]} {
	    set info_str [string range $line 10 end]
	    ::ogf::cat::set status "PSF built: $info_str"
	    break
	}
    }
}

proc CatalogPanelBuildExtendedPSF {} {
    global ed

    if {![::ogf::cat::has]} {
	::ogf::cat::set status "Extract sources first before building extended PSF"
	return
    }

    set w .extpsfdlg
    if {[winfo exists $w]} {
	raise $w
	return
    }

    toplevel $w
    wm title $w "Extended PSF"
    wm geometry $w 360x420

    # Copy current values
    foreach pname {ext-core-mag-min ext-core-mag-max ext-wing-mag-max \
		   ext-core-size ext-wing-size ext-blend-inner ext-blend-outer \
		   ext-saturation-limit} {
	set ed(psf,$pname) [::ogf::cat::get psf,param,$pname]
    }

    set f [ttk::frame $w.content]
    pack $f -fill both -expand true -padx 8 -pady 8

    set r 0

    # Core Stars section
    ttk::label $f.hcore -text "Core Stars" -font TkHeadingFont
    grid $f.hcore -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {4 2}
    incr r

    ttk::label $f.lmagmin -text "Magnitude range:"
    ttk::entry $f.emagmin -textvariable ed(psf,ext-core-mag-min) -width 8
    ttk::label $f.ltilde -text "~"
    ttk::entry $f.emagmax -textvariable ed(psf,ext-core-mag-max) -width 8
    grid $f.lmagmin -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.emagmin -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.ltilde  -row $r -column 2 -padx 2 -pady 2
    grid $f.emagmax -row $r -column 3 -sticky w -padx 2 -pady 2
    incr r

    ttk::label $f.lcsize -text "Cutout size:"
    ttk::entry $f.ecsize -textvariable ed(psf,ext-core-size) -width 8
    ttk::label $f.lcpx -text "px"
    grid $f.lcsize -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ecsize -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lcpx   -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    # Wing Stars section
    ttk::label $f.hwing -text "Wing Stars" -font TkHeadingFont
    grid $f.hwing -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}
    incr r

    ttk::label $f.lwmag -text "Max magnitude:"
    ttk::entry $f.ewmag -textvariable ed(psf,ext-wing-mag-max) -width 8
    grid $f.lwmag -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ewmag -row $r -column 1 -sticky w -padx 2 -pady 2
    incr r

    ttk::label $f.lwsize -text "Cutout size:"
    ttk::entry $f.ewsize -textvariable ed(psf,ext-wing-size) -width 8
    ttk::label $f.lwpx -text "px"
    grid $f.lwsize -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ewsize -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lwpx   -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    ttk::label $f.lsat -text "Saturation:"
    ttk::entry $f.esat -textvariable ed(psf,ext-saturation-limit) -width 8
    ttk::label $f.ladu -text "ADU"
    grid $f.lsat -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.esat -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.ladu -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    # Blending section
    ttk::label $f.hblend -text "Blending" -font TkHeadingFont
    grid $f.hblend -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}
    incr r

    ttk::label $f.lbin -text "Inner radius:"
    ttk::entry $f.ebin -textvariable ed(psf,ext-blend-inner) -width 8
    ttk::label $f.lbinpx -text "px"
    grid $f.lbin   -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ebin   -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lbinpx -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    ttk::label $f.lbout -text "Outer radius:"
    ttk::entry $f.ebout -textvariable ed(psf,ext-blend-outer) -width 8
    ttk::label $f.lbopx -text "px"
    grid $f.lbout  -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ebout  -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lbopx  -row $r -column 2 -sticky w -padx 2 -pady 2

    # Buttons
    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8

    ttk::button $bf.build -text "Build" \
	-command [list CatalogPanelBuildExtendedPSFExec $w]
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.build -side right -padx 4
}

proc CatalogPanelBuildExtendedPSFExec {w} {
    global ed

    # Apply params
    foreach pname {ext-core-mag-min ext-core-mag-max ext-wing-mag-max \
		   ext-core-size ext-wing-size ext-blend-inner ext-blend-outer \
		   ext-saturation-limit} {
	::ogf::cat::set psf,param,$pname $ed(psf,$pname)
    }
    CatalogPanelPSFParamSave

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	::ogf::cat::set status "No FITS image loaded"
	return
    }

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	::ogf::cat::set status "ERROR: ds9_psf_deconv.py not found"
	return
    }

    # Save catalog
    set catfile [file join [file normalize ~] .ds9 psf_catalog.tsv]
    catch {file mkdir [file dirname $catfile]}
    if {[catch {
	set fd [open $catfile w]
	puts $fd [::ogf::cat::tsv]
	close $fd
    } err]} {
	::ogf::cat::set status "Extended PSF error: cannot write catalog: $err"
	return
    }

    set paramargs {}
    lappend paramargs "--mode" "build_psf_extended"
    lappend paramargs "--catalog" $catfile
    lappend paramargs "--ext-core-mag-min" [::ogf::cat::get psf,param,ext-core-mag-min]
    lappend paramargs "--ext-core-mag-max" [::ogf::cat::get psf,param,ext-core-mag-max]
    lappend paramargs "--ext-wing-mag-max" [::ogf::cat::get psf,param,ext-wing-mag-max]
    lappend paramargs "--ext-core-size" [::ogf::cat::get psf,param,ext-core-size]
    lappend paramargs "--ext-wing-size" [::ogf::cat::get psf,param,ext-wing-size]
    lappend paramargs "--ext-blend-inner" [::ogf::cat::get psf,param,ext-blend-inner]
    lappend paramargs "--ext-blend-outer" [::ogf::cat::get psf,param,ext-blend-outer]
    lappend paramargs "--ext-saturation-limit" [::ogf::cat::get psf,param,ext-saturation-limit]
    lappend paramargs "--psf-output" [::ogf::cat::get psf,file]

    ::ogf::cat::set status "Building extended PSF ..."
    update idletasks

    set errfile [file join [file normalize ~] .ds9 psf_stderr.txt]
    OGFSessLog psf.build_extended auto [list [OGFPython] $script $fn {*}$paramargs] -title {Build extended PSF} -post [dict create kind psfbuilt]
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
	    ::ogf::cat::set status "Extended PSF error: $last_err"
	    puts "Extended PSF stderr:\n$stderr_msg"
	} else {
	    ::ogf::cat::set status "Extended PSF error: $err"
	}
	return
    }
    catch {file delete $errfile}
    catch {file delete $catfile}

    ::ogf::cat::set psf,has_psf 1

    foreach line [split $data \n] {
	if {[string match "#PSF_EXTENDED*" $line]} {
	    ::ogf::cat::set status "Extended PSF: [string range $line 15 end]"
	    break
	}
    }

    CatalogPanelViewPSF
}

proc CatalogPanelCheckSimAvail {} {

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	::ogf::cat::set psf,sim_webbpsf_ok 0
	::ogf::cat::set psf,sim_tinytim_ok 0
	return
    }

    # Use a dummy fits arg for check_sim mode
    if {[catch {set data [exec [OGFPython] $script dummy.fits --mode check_sim 2>/dev/null]} err]} {
	::ogf::cat::set psf,sim_webbpsf_ok 0
	::ogf::cat::set psf,sim_tinytim_ok 0
	return
    }

    foreach line [split $data \n] {
	if {[string match "#SIM_STATUS*" $line]} {
	    foreach part [split $line \t] {
		if {[string match "WEBBPSF=*" $part]} {
		    ::ogf::cat::set psf,sim_webbpsf_ok [string range $part 8 end]
		}
		if {[string match "TINYTIM=*" $part]} {
		    ::ogf::cat::set psf,sim_tinytim_ok [string range $part 8 end]
		}
	    }
	}
    }
}

proc CatalogPanelSimPSFWebbPSF {} {
    global ed

    set w .webbpsfdlg
    if {[winfo exists $w]} {
	raise $w
	return
    }

    # Check availability if not yet done
    if {[::ogf::cat::get psf,sim_webbpsf_ok] == -1} {
	CatalogPanelCheckSimAvail
    }

    toplevel $w
    wm title $w "WebbPSF (JWST)"
    wm geometry $w 360x380

    set ed(psf,sim-instrument) [::ogf::cat::get psf,param,sim-instrument]
    set ed(psf,sim-filter) [::ogf::cat::get psf,param,sim-filter]
    set ed(psf,sim-psf-size) [::ogf::cat::get psf,param,sim-psf-size]
    set ed(psf,sim-oversample) [::ogf::cat::get psf,param,sim-oversample]
    set ed(psf,sim-jitter-sigma) [::ogf::cat::get psf,param,sim-jitter-sigma]
    set ed(psf,sim-focus-offset) [::ogf::cat::get psf,param,sim-focus-offset]

    set f [ttk::frame $w.content]
    pack $f -fill both -expand true -padx 8 -pady 8

    set r 0

    # Availability status
    if {[::ogf::cat::get psf,sim_webbpsf_ok] == 1} {
	set statxt "WebbPSF: Available"
    } else {
	set statxt "WebbPSF: Not found (pip install webbpsf)"
    }
    ttk::label $f.status -text $statxt -foreground \
	[expr {[::ogf::cat::get psf,sim_webbpsf_ok] == 1 ? "green" : "red"}]
    grid $f.status -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {4 8}
    incr r

    # Instrument
    ttk::label $f.linst -text "Instrument:"
    ttk::combobox $f.cinst -textvariable ed(psf,sim-instrument) -width 14 \
	-values {NIRCAM MIRI NIRISS NIRSPEC FGS} -state readonly
    grid $f.linst -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.cinst -row $r -column 1 -sticky w -padx 2 -pady 2
    incr r

    # Filter
    ttk::label $f.lfilt -text "Filter:"
    ttk::combobox $f.cfilt -textvariable ed(psf,sim-filter) -width 14
    grid $f.lfilt -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.cfilt -row $r -column 1 -sticky w -padx 2 -pady 2
    incr r

    # Dynamic filter update
    bind $f.cinst <<ComboboxSelected>> [list CatalogPanelSimPSFUpdateFilters $f.cfilt jwst]
    # Initialize filter list
    CatalogPanelSimPSFUpdateFilters $f.cfilt jwst

    # PSF size
    ttk::label $f.lsz -text "PSF size:"
    ttk::entry $f.esz -textvariable ed(psf,sim-psf-size) -width 8
    ttk::label $f.lpx -text "px"
    grid $f.lsz -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.esz -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lpx -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    # Oversample
    ttk::label $f.lover -text "Oversample:"
    ttk::entry $f.eover -textvariable ed(psf,sim-oversample) -width 8
    grid $f.lover -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.eover -row $r -column 1 -sticky w -padx 2 -pady 2
    incr r

    # Jitter
    ttk::label $f.ljit -text "Jitter sigma:"
    ttk::entry $f.ejit -textvariable ed(psf,sim-jitter-sigma) -width 8
    ttk::label $f.ljitas -text "arcsec"
    grid $f.ljit   -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ejit   -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.ljitas -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    # Focus
    ttk::label $f.lfoc -text "Focus offset:"
    ttk::entry $f.efoc -textvariable ed(psf,sim-focus-offset) -width 8
    ttk::label $f.lfocw -text "waves"
    grid $f.lfoc  -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.efoc  -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lfocw -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    # Auto-detect button
    ttk::button $f.autodet -text "Auto-detect from FITS" \
	-command [list CatalogPanelSimPSFAutoDetect $f.cinst $f.cfilt jwst]
    grid $f.autodet -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}

    # Buttons
    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8

    ttk::button $bf.gen -text "Generate" \
	-command [list CatalogPanelSimPSFWebbPSFExec $w]
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.gen -side right -padx 4
}

proc CatalogPanelSimPSFWebbPSFExec {w} {
    global ed

    foreach pname {sim-instrument sim-filter sim-psf-size sim-oversample \
		   sim-jitter-sigma sim-focus-offset} {
	::ogf::cat::set psf,param,$pname $ed(psf,$pname)
    }
    ::ogf::cat::set psf,param,sim-telescope jwst
    CatalogPanelPSFParamSave

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	::ogf::cat::set status "No FITS image loaded"
	return
    }

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	::ogf::cat::set status "ERROR: ds9_psf_deconv.py not found"
	return
    }

    set paramargs {}
    lappend paramargs "--mode" "sim_psf"
    lappend paramargs "--sim-telescope" "jwst"
    lappend paramargs "--sim-instrument" [::ogf::cat::get psf,param,sim-instrument]
    lappend paramargs "--sim-filter" [::ogf::cat::get psf,param,sim-filter]
    lappend paramargs "--sim-psf-size" [::ogf::cat::get psf,param,sim-psf-size]
    lappend paramargs "--sim-oversample" [::ogf::cat::get psf,param,sim-oversample]
    lappend paramargs "--sim-jitter-sigma" [::ogf::cat::get psf,param,sim-jitter-sigma]
    lappend paramargs "--sim-focus-offset" [::ogf::cat::get psf,param,sim-focus-offset]
    lappend paramargs "--psf-output" [::ogf::cat::get psf,file]

    ::ogf::cat::set status "Generating WebbPSF ([::ogf::cat::get psf,param,sim-instrument] / [::ogf::cat::get psf,param,sim-filter]) ..."
    update idletasks

    set errfile [file join [file normalize ~] .ds9 psf_stderr.txt]
    OGFSessLog psf.webbpsf auto [list [OGFPython] $script $fn {*}$paramargs] -title {WebbPSF (JWST) PSF} -post [dict create kind psfbuilt]
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
	    ::ogf::cat::set status "WebbPSF error: $last_err"
	    puts "WebbPSF stderr:\n$stderr_msg"
	} else {
	    ::ogf::cat::set status "WebbPSF error: $err"
	}
	return
    }
    catch {file delete $errfile}

    ::ogf::cat::set psf,has_psf 1

    foreach line [split $data \n] {
	if {[string match "#PSF_SIM*" $line]} {
	    ::ogf::cat::set status "Sim PSF: [string range $line 9 end]"
	    break
	}
    }

    CatalogPanelViewPSF
}

proc CatalogPanelSimPSFTinyTim {} {
    global ed

    set w .tinytimdlg
    if {[winfo exists $w]} {
	raise $w
	return
    }

    if {[::ogf::cat::get psf,sim_tinytim_ok] == -1} {
	CatalogPanelCheckSimAvail
    }

    toplevel $w
    wm title $w "TinyTim (HST)"
    wm geometry $w 360x340

    set ed(psf,sim-instrument) [::ogf::cat::get psf,param,sim-instrument]
    set ed(psf,sim-filter) [::ogf::cat::get psf,param,sim-filter]
    set ed(psf,sim-psf-size) [::ogf::cat::get psf,param,sim-psf-size]
    set ed(psf,sim-oversample) [::ogf::cat::get psf,param,sim-oversample]
    set ed(psf,sim-focus-offset) [::ogf::cat::get psf,param,sim-focus-offset]

    set f [ttk::frame $w.content]
    pack $f -fill both -expand true -padx 8 -pady 8

    set r 0

    # Availability status
    if {[::ogf::cat::get psf,sim_tinytim_ok] == 1} {
	set statxt "TinyTim: Available"
    } else {
	set statxt "TinyTim: Not found (tiny1/tiny2/tiny3 not on PATH)"
    }
    ttk::label $f.status -text $statxt -foreground \
	[expr {[::ogf::cat::get psf,sim_tinytim_ok] == 1 ? "green" : "red"}]
    grid $f.status -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {4 8}
    incr r

    # Instrument
    ttk::label $f.linst -text "Instrument:"
    ttk::combobox $f.cinst -textvariable ed(psf,sim-instrument) -width 14 \
	-values {ACS_WFC ACS_HRC WFC3_UVIS WFC3_IR WFPC2} -state readonly
    grid $f.linst -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.cinst -row $r -column 1 -sticky w -padx 2 -pady 2
    incr r

    # Filter
    ttk::label $f.lfilt -text "Filter:"
    ttk::combobox $f.cfilt -textvariable ed(psf,sim-filter) -width 14
    grid $f.lfilt -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.cfilt -row $r -column 1 -sticky w -padx 2 -pady 2
    incr r

    bind $f.cinst <<ComboboxSelected>> [list CatalogPanelSimPSFUpdateFilters $f.cfilt hst]
    CatalogPanelSimPSFUpdateFilters $f.cfilt hst

    # PSF size
    ttk::label $f.lsz -text "PSF size:"
    ttk::entry $f.esz -textvariable ed(psf,sim-psf-size) -width 8
    ttk::label $f.lpx -text "px"
    grid $f.lsz -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.esz -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lpx -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    # Oversample
    ttk::label $f.lover -text "Oversample:"
    ttk::entry $f.eover -textvariable ed(psf,sim-oversample) -width 8
    grid $f.lover -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.eover -row $r -column 1 -sticky w -padx 2 -pady 2
    incr r

    # Focus
    ttk::label $f.lfoc -text "Focus offset:"
    ttk::entry $f.efoc -textvariable ed(psf,sim-focus-offset) -width 8
    ttk::label $f.lfocum -text "um"
    grid $f.lfoc   -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.efoc   -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lfocum -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    # Auto-detect
    ttk::button $f.autodet -text "Auto-detect from FITS" \
	-command [list CatalogPanelSimPSFAutoDetect $f.cinst $f.cfilt hst]
    grid $f.autodet -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}

    # Buttons
    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8

    ttk::button $bf.gen -text "Generate" \
	-command [list CatalogPanelSimPSFTinyTimExec $w]
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.gen -side right -padx 4
}

proc CatalogPanelSimPSFTinyTimExec {w} {
    global ed

    foreach pname {sim-instrument sim-filter sim-psf-size sim-oversample \
		   sim-focus-offset} {
	::ogf::cat::set psf,param,$pname $ed(psf,$pname)
    }
    ::ogf::cat::set psf,param,sim-telescope hst
    CatalogPanelPSFParamSave

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	::ogf::cat::set status "No FITS image loaded"
	return
    }

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	::ogf::cat::set status "ERROR: ds9_psf_deconv.py not found"
	return
    }

    set paramargs {}
    lappend paramargs "--mode" "sim_psf"
    lappend paramargs "--sim-telescope" "hst"
    lappend paramargs "--sim-instrument" [::ogf::cat::get psf,param,sim-instrument]
    lappend paramargs "--sim-filter" [::ogf::cat::get psf,param,sim-filter]
    lappend paramargs "--sim-psf-size" [::ogf::cat::get psf,param,sim-psf-size]
    lappend paramargs "--sim-oversample" [::ogf::cat::get psf,param,sim-oversample]
    lappend paramargs "--sim-focus-offset" [::ogf::cat::get psf,param,sim-focus-offset]
    lappend paramargs "--psf-output" [::ogf::cat::get psf,file]

    ::ogf::cat::set status "Generating TinyTim PSF ([::ogf::cat::get psf,param,sim-instrument] / [::ogf::cat::get psf,param,sim-filter]) ..."
    update idletasks

    set errfile [file join [file normalize ~] .ds9 psf_stderr.txt]
    OGFSessLog psf.tinytim auto [list [OGFPython] $script $fn {*}$paramargs] -title {TinyTim (HST) PSF} -post [dict create kind psfbuilt]
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
	    ::ogf::cat::set status "TinyTim error: $last_err"
	    puts "TinyTim stderr:\n$stderr_msg"
	} else {
	    ::ogf::cat::set status "TinyTim error: $err"
	}
	return
    }
    catch {file delete $errfile}

    ::ogf::cat::set psf,has_psf 1

    foreach line [split $data \n] {
	if {[string match "#PSF_SIM*" $line]} {
	    ::ogf::cat::set status "Sim PSF: [string range $line 9 end]"
	    break
	}
    }

    CatalogPanelViewPSF
}

proc CatalogPanelSimPSFUpdateFilters {cfilt telescope} {
    global ed

    # Filter lists per instrument
    array set jwst_filters {
	NIRCAM  {F070W F090W F115W F140M F150W F162M F200W F210M F250M F277W F300M F322W2 F335M F356W F360M F410M F430M F444W F460M F480M}
	MIRI    {F560W F770W F1000W F1065C F1130W F1140C F1280W F1500W F1550C F1800W F2100W F2300C F2550W}
	NIRISS  {F090W F115W F140M F150W F158M F200W F277W F356W F380M F430M F444W F480M}
	NIRSPEC {F070LP F100LP F170LP F290LP CLEAR}
	FGS     {FGS}
    }
    array set hst_filters {
	ACS_WFC   {F435W F475W F502N F550M F555W F606W F625W F658N F775W F814W F850LP}
	ACS_HRC   {F220W F250W F330W F435W F475W F555W F606W F625W F775W F814W F850LP}
	WFC3_UVIS {F218W F225W F275W F336W F390W F438W F475W F555W F606W F625W F775W F814W F850LP}
	WFC3_IR   {F098M F105W F110W F125W F140W F160W}
	WFPC2     {F300W F336W F439W F450W F555W F606W F675W F702W F791W F814W}
    }

    set inst $ed(psf,sim-instrument)
    set inst [string toupper $inst]

    if {$telescope eq "jwst" && [info exists jwst_filters($inst)]} {
	$cfilt configure -values $jwst_filters($inst)
    } elseif {$telescope eq "hst" && [info exists hst_filters($inst)]} {
	$cfilt configure -values $hst_filters($inst)
    } else {
	$cfilt configure -values {}
    }
}

proc CatalogPanelSimPSFAutoDetect {cinst cfilt telescope} {
    global ed

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	::ogf::cat::set status "No FITS image loaded"
	return
    }

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	::ogf::cat::set status "ERROR: ds9_psf_deconv.py not found"
	return
    }

    # Run sim_psf with all auto to just get detection output
    set paramargs {}
    lappend paramargs "--mode" "sim_psf"
    lappend paramargs "--sim-telescope" "auto"
    lappend paramargs "--sim-instrument" "auto"
    lappend paramargs "--sim-filter" "auto"
    lappend paramargs "--psf-output" "/dev/null"

    # We expect this to potentially fail (no webbpsf/tinytim), but stderr has detection info
    set errfile [file join [file normalize ~] .ds9 psf_stderr.txt]
    catch {exec python3 $script $fn {*}$paramargs 2>$errfile}

    set stderr_msg ""
    catch {
	set fd [open $errfile r]
	set stderr_msg [read $fd]
	close $fd
    }

    # Parse "Auto-detected: {'telescope': ..., 'instrument': ..., 'filter': ...}"
    if {[regexp {instrument.*?:\s*'([^']+)'} $stderr_msg -> inst_val]} {
	set ed(psf,sim-instrument) [string toupper $inst_val]
	$cinst set [string toupper $inst_val]
    }
    if {[regexp {filter.*?:\s*'([^']+)'} $stderr_msg -> filt_val]} {
	set ed(psf,sim-filter) [string toupper $filt_val]
	$cfilt set [string toupper $filt_val]
    }

    # Update filter list for the detected instrument
    CatalogPanelSimPSFUpdateFilters $cfilt $telescope

    ::ogf::cat::set status "Auto-detected: $ed(psf,sim-instrument) / $ed(psf,sim-filter)"
}

proc CatalogPanelViewPSF {} {

    if {![::ogf::cat::get psf,has_psf] || ![file exists [::ogf::cat::get psf,file]]} {
	::ogf::cat::set status "No PSF available — build PSF first"
	return
    }

    set w .psfviewer
    if {[winfo exists $w]} {
	raise $w
	CatalogPanelViewPSFRender $w
	return
    }

    toplevel $w
    wm title $w "PSF Viewer"

    # Image display
    ttk::label $w.img -anchor center
    pack $w.img -padx 8 -pady 8 -fill both -expand true

    # Info label
    ttk::label $w.info -text "" -anchor center
    pack $w.info -padx 8 -pady {0 4}

    # Buttons: Save / Load / Close
    ttk::frame $w.btn
    ttk::button $w.btn.save -text "Save PSF..." \
	-command CatalogPanelSavePSF
    ttk::button $w.btn.load -text "Load PSF..." \
	-command [list CatalogPanelViewPSFLoad $w]
    ttk::button $w.btn.close -text "Close" \
	-command [list destroy $w]
    pack $w.btn.save -side left -padx 4
    pack $w.btn.load -side left -padx 4
    pack $w.btn.close -side right -padx 4
    pack $w.btn -fill x -padx 8 -pady 8

    CatalogPanelViewPSFRender $w
}

proc CatalogPanelViewPSFRender {w} {

    set psffile [::ogf::cat::get psf,file]
    set tmpimg [file join [file normalize ~] .ds9 psf_view.ppm]

    # Write render script to temp file
    set tmpscript [file join [file normalize ~] .ds9 psf_render.py]
    set fd [open $tmpscript w]
    puts $fd {import numpy as np
from astropy.io import fits
from scipy.ndimage import zoom
import sys

psffile = sys.argv[1]
outfile = sys.argv[2]

with fits.open(psffile) as hdul:
    data = hdul[0].data.astype(np.float64)

h0, w0 = data.shape
vmin, vmax = float(data.min()), float(data.max())

# Asinh stretch
if vmax > vmin:
    norm = (data - vmin) / (vmax - vmin)
    stretched = np.arcsinh(norm * 10) / np.arcsinh(10)
    img = (stretched * 255).clip(0, 255).astype(np.uint8)
else:
    img = np.zeros_like(data, dtype=np.uint8)

# Scale up to at least 256x256
scale = max(1, 256 // max(img.shape))
if scale > 1:
    img = zoom(img, scale, order=0)

h, w = img.shape

# Write PPM (P6 RGB) — Tk reads this natively
with open(outfile, 'wb') as f:
    f.write(f'P6\n{w} {h}\n255\n'.encode())
    rgb = np.stack([img, img, img], axis=-1)
    f.write(rgb.tobytes())

print(f'{w0}x{h0}  peak={vmax:.4g}')
}
    close $fd

    if {[catch {set info [exec [OGFPython] $tmpscript $psffile $tmpimg]} err]} {
	catch {$w.info configure -text "Render error: $err"}
	catch {file delete $tmpscript}
	return
    }
    catch {file delete $tmpscript}

    # Load into Tk photo image
    catch {image delete psfviewimg}
    image create photo psfviewimg -file $tmpimg
    $w.img configure -image psfviewimg
    $w.info configure -text "PSF: $info"

    # Resize window to fit image + buttons
    set iw [image width psfviewimg]
    set ih [image height psfviewimg]
    set ww [expr {max($iw + 16, 280)}]
    set wh [expr {$ih + 90}]
    wm geometry $w ${ww}x${wh}
}

proc CatalogPanelViewPSFLoad {w} {

    set types {
	{{FITS Files} {.fits .fit .fts}}
	{{All Files} *}
    }
    set infile [tk_getOpenFile -filetypes $types \
		    -title "Load PSF FITS"]
    if {$infile eq {}} return

    file copy -force $infile [::ogf::cat::get psf,file]
    ::ogf::cat::set psf,has_psf 1
    ::ogf::cat::set status "PSF loaded from $infile"

    # Refresh the viewer
    CatalogPanelViewPSFRender $w
}

proc CatalogPanelSavePSF {} {

    if {![::ogf::cat::get psf,has_psf] || ![file exists [::ogf::cat::get psf,file]]} {
	::ogf::cat::set status "No PSF available — build PSF first"
	return
    }

    set types {
	{{FITS Files} {.fits .fit .fts}}
	{{All Files} *}
    }
    set outfile [tk_getSaveFile -filetypes $types \
		     -title "Save PSF FITS" \
		     -initialfile "psf.fits"]
    if {$outfile eq {}} return

    file copy -force [::ogf::cat::get psf,file] $outfile
    ::ogf::cat::set status "PSF saved to $outfile"
}

proc CatalogPanelLoadPSF {} {

    set types {
	{{FITS Files} {.fits .fit .fts}}
	{{All Files} *}
    }
    set infile [tk_getOpenFile -filetypes $types \
		    -title "Load PSF FITS"]
    if {$infile eq {}} return

    file copy -force $infile [::ogf::cat::get psf,file]
    ::ogf::cat::set psf,has_psf 1
    ::ogf::cat::set status "PSF loaded from $infile"

    # Refresh viewer if open
    if {[winfo exists .psfviewer]} {
	CatalogPanelViewPSFRender .psfviewer
    }
}

