# extract/extract.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/extract/plugin.json.

# Run source extraction on the currently loaded FITS image
proc CatalogPanelExtract {} {
    global ds9
    global current
    global loadParam

    # Find the ds9_sextract binary (platform-aware)
    set bindir [file dirname [info nameofexecutable]]
    set os $::tcl_platform(os)

    if {$os eq "Windows NT"} {
	set sextract [file join $bindir ds9_sextract.exe]
    } else {
	set sextract [file join $bindir ds9_sextract]
    }
    if {![file executable $sextract]} {
	::ogf::cat::set status "ERROR: ds9_sextract not found in $bindir"
	return
    }

    # Get current FITS filename
    set fn {}
    if {$current(frame) != {}} {
	catch {set fn [$current(frame) get fits file name full]}
    }
    if {$fn eq {}} {
	::ogf::cat::set status "No FITS image loaded"
	return
    }

    # Strip curly braces if present
    set fn [string trim $fn "{}"]

    # Strip FITS HDU extension specifier (e.g. [SCI], [1], [SCI,2])
    regsub {\[.*\]$} $fn {} fn

    if {![file exists $fn]} {
	::ogf::cat::set status "File not found: $fn"
	return
    }

    # Set log scale with optimized limits
    CatalogPanelSetLogScale

    ::ogf::cat::set status "Extracting sources from [file tail $fn] ..."
    update idletasks

    # Save extraction parameters so AI Merge uses the same values
    foreach pname {detect-thresh detect-minarea deblend-nthresh deblend-mincont \
		   mag-zeropoint back-size back-filtersize} {
	if {[::ogf::cat::exists param,$pname]} {
	    ::ogf::cat::set extract_param,$pname [::ogf::cat::get param,$pname]
	}
    }

    # Build parameter arguments list
    set paramargs {}
    foreach pname {detect-thresh detect-minarea deblend-nthresh deblend-mincont \
		   phot-aperture mag-zeropoint gain pixel-scale seeing-fwhm \
		   back-size back-filtersize \
		   phot-aperture-2 phot-aperture-3 phot-aperture-5 conv-filter} {
	if {[::ogf::cat::exists param,$pname]} {
	    lappend paramargs "--$pname" [::ogf::cat::get param,$pname]
	}
    }

    # Platform-specific library path setup and execution
    if {$os eq "Darwin"} {
	# macOS: set DYLD_LIBRARY_PATH
	set libpaths {}
	if {[info exists ::env(CONDA_PREFIX)]} {
	    lappend libpaths "$::env(CONDA_PREFIX)/lib"
	}
	set home_conda [file join [file normalize ~] miniconda3/lib]
	if {[file isdirectory $home_conda]} {
	    lappend libpaths $home_conda
	}
	if {[info exists ::env(DYLD_LIBRARY_PATH)]} {
	    lappend libpaths $::env(DYLD_LIBRARY_PATH)
	}
	if {[llength $libpaths] > 0} {
	    set ::env(DYLD_LIBRARY_PATH) [join $libpaths :]
	}
    } elseif {$os ne "Windows NT"} {
	# Linux/Unix: set LD_LIBRARY_PATH
	set libpaths {}
	if {[info exists ::env(CONDA_PREFIX)]} {
	    lappend libpaths "$::env(CONDA_PREFIX)/lib"
	}
	set home_conda [file join [file normalize ~] miniconda3/lib]
	if {[file isdirectory $home_conda]} {
	    lappend libpaths $home_conda
	}
	if {[info exists ::env(LD_LIBRARY_PATH)]} {
	    lappend libpaths $::env(LD_LIBRARY_PATH)
	}
	if {[llength $libpaths] > 0} {
	    set ::env(LD_LIBRARY_PATH) [join $libpaths :]
	}
    }
    # Windows: DLLs found via PATH automatically

    OGFSessLog extract auto [list $sextract $fn {*}$paramargs] -title "Extract sources (ds9_sextract)" \
	-tool sextract -post [dict create kind set]
    # Run extraction (cross-platform exec)
    if {[catch {set data [exec $sextract $fn {*}$paramargs 2>@stderr]} err]} {
	::ogf::cat::set status "Extraction error: $err"
	return
    }

    # Parse TSV output into table
    CatalogPanelLoadTSV $data [file tail $fn]
}

proc CatalogPanelParamDef {} {

    ::ogf::cat::set param,detect-thresh 1.5
    ::ogf::cat::set param,detect-minarea 5
    ::ogf::cat::set param,deblend-nthresh 32
    ::ogf::cat::set param,deblend-mincont 0.005
    ::ogf::cat::set param,phot-aperture 5.0
    ::ogf::cat::set param,mag-zeropoint 25.0
    ::ogf::cat::set param,gain 0.0
    ::ogf::cat::set param,pixel-scale 1.0
    ::ogf::cat::set param,seeing-fwhm 3.0
    ::ogf::cat::set param,back-size 64
    ::ogf::cat::set param,back-filtersize 3
    ::ogf::cat::set param,phot-aperture-2 4.0
    ::ogf::cat::set param,phot-aperture-3 6.0
    ::ogf::cat::set param,phot-aperture-5 10.0
    ::ogf::cat::set param,conv-filter default
    ::ogf::cat::set param,n-workers 0

    # Separate (deblend) parameters
    ::ogf::cat::set param,sep-deblend-nthresh 64
    ::ogf::cat::set param,sep-deblend-mincont 0.0001
    ::ogf::cat::set param,sep-detect-thresh 0.8
    ::ogf::cat::set param,sep-detect-minarea 3
    ::ogf::cat::set param,sep-radius-factor 3.0
    ::ogf::cat::set param,sep-back-size 32

    CatalogPanelParamLoad
}

proc CatalogPanelParamLoad {} {

    set preffile [file join [file normalize ~] .ds9 sextract.prf]
    if {![file exists $preffile]} return
    if {[catch {set fd [open $preffile r]} err]} return
    while {[gets $fd line] >= 0} {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set parts [split $line]
	if {[llength $parts] >= 2} {
	    set key [lindex $parts 0]
	    set val [lindex $parts 1]
	    if {[::ogf::cat::exists param,$key]} {
		::ogf::cat::set param,$key $val
	    }
	}
    }
    close $fd
}

proc CatalogPanelParamSave {} {

    set prefdir [file join [file normalize ~] .ds9]
    if {![file isdirectory $prefdir]} {
	file mkdir $prefdir
    }
    set preffile [file join $prefdir sextract.prf]
    if {[catch {set fd [open $preffile w]} err]} return
    foreach pname {detect-thresh detect-minarea deblend-nthresh deblend-mincont \
		   phot-aperture mag-zeropoint gain pixel-scale seeing-fwhm \
		   back-size back-filtersize \
		   phot-aperture-2 phot-aperture-3 phot-aperture-5 conv-filter \
		   n-workers \
		   sep-deblend-nthresh sep-deblend-mincont sep-detect-thresh \
		   sep-detect-minarea sep-radius-factor sep-back-size} {
	puts $fd "$pname [::ogf::cat::get param,$pname]"
    }
    close $fd
}

proc CatalogPanelParamDefaults {} {
    global ed

    set ed(detect-thresh) 1.5
    set ed(detect-minarea) 5
    set ed(deblend-nthresh) 32
    set ed(deblend-mincont) 0.005
    set ed(phot-aperture) 5.0
    set ed(mag-zeropoint) 25.0
    set ed(gain) 0.0
    set ed(pixel-scale) 1.0
    set ed(seeing-fwhm) 3.0
    set ed(back-size) 64
    set ed(back-filtersize) 3
    set ed(phot-aperture-2) 4.0
    set ed(phot-aperture-3) 6.0
    set ed(phot-aperture-5) 10.0
    set ed(conv-filter) default
    set ed(n-workers) 0
}

proc CatalogPanelDualExtract {} {
    global ds9 current ed

    set w {.dualextract}
    set ed(ok) 0

    DialogCreate $w {Dual-Image Extract} ed(ok)

    set f [ttk::frame $w.param]

    # Detection image
    ttk::label $f.ldet -text "Detection Image:" -anchor w
    ttk::entry $f.edet -textvariable ed(dual,detect) -width 40
    ttk::button $f.bdet -text "Browse..." -command {
	set ff [tk_getOpenFile -title "Detection Image" \
	    -filetypes {{{FITS Files} {.fits .fit .fts}} {{All Files} {*}}}]
	if {$ff ne {}} { set ed(dual,detect) $ff }
    }
    grid $f.ldet $f.edet $f.bdet -padx 4 -pady 2 -sticky w

    # Measurement image
    ttk::label $f.lmeas -text "Measurement Image:" -anchor w
    ttk::entry $f.emeas -textvariable ed(dual,measure) -width 40
    ttk::button $f.bmeas -text "Browse..." -command {
	set ff [tk_getOpenFile -title "Measurement Image" \
	    -filetypes {{{FITS Files} {.fits .fit .fts}} {{All Files} {*}}}]
	if {$ff ne {}} { set ed(dual,measure) $ff }
    }
    grid $f.lmeas $f.emeas $f.bmeas -padx 4 -pady 2 -sticky w

    # Set defaults from current frame
    set ed(dual,detect) [CatalogPanelGetFITS]
    set ed(dual,measure) {}

    # Buttons
    set bf [ttk::frame $w.buttons]
    ttk::button $bf.ok -text {Extract} -command {set ed(ok) 1} -default active
    ttk::button $bf.cancel -text {Cancel} -command {set ed(ok) 0}
    pack $bf.ok $bf.cancel -side left -expand true -padx 2 -pady 4

    bind $w <Return> {set ed(ok) 1}
    ttk::separator $w.sep -orient horizontal
    pack $w.buttons $w.sep -side bottom -fill x
    pack $w.param -side top -fill both -expand true

    DialogWait $w ed(ok) $f.edet
    set detect_img $ed(dual,detect)
    set measure_img $ed(dual,measure)
    destroy $w

    if {!$ed(ok)} { unset ed; return }
    unset ed

    if {$detect_img eq {} || ![file exists $detect_img]} {
	::ogf::cat::set status "Detection image not found"
	return
    }
    if {$measure_img eq {} || ![file exists $measure_img]} {
	::ogf::cat::set status "Measurement image not found"
	return
    }

    set script [CatalogPanelGetScript ds9_dual_extract.py]
    if {![file exists $script]} {
	::ogf::cat::set status "Script not found: ds9_dual_extract.py"
	return
    }

    ::ogf::cat::set status "Dual-image extraction..."
    update idletasks

    set args [list [OGFPython] $script \
	--detect-image $detect_img --measure-image $measure_img]
    if {[::ogf::cat::exists param,detect-thresh]} {
	lappend args --detect-thresh [::ogf::cat::get param,detect-thresh]
    }
    if {[::ogf::cat::exists param,detect-minarea]} {
	lappend args --detect-minarea [::ogf::cat::get param,detect-minarea]
    }
    if {[::ogf::cat::exists param,phot-aperture]} {
	lappend args --phot-aperture [::ogf::cat::get param,phot-aperture]
    }
    if {[::ogf::cat::exists param,mag-zeropoint]} {
	lappend args --mag-zeropoint [::ogf::cat::get param,mag-zeropoint]
    }

    OGFSessLog analysis.dual_extract auto $args -title {Dual-image extract} 
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	::ogf::cat::set status "Dual extract error: $err"
	return
    }

    CatalogPanelLoadTSV $data "dual-image"
    ::ogf::cat::set status "Dual-image extraction complete"
}

# Hook: automatically extract sources after FITS file is loaded
proc CatalogPanelAutoExtract {} {
    if {[::ogf::cat::exists tbl]} {
	after 500 CatalogPanelExtract
    }
}

