# Multi-component fit plugin (plugins/multifit): Tcl side = outputs of the CLI steps (frames, montage window, tables, selection -> object list).
# Catalog keys multifit,* (registered in ogf_core.tcl); no direct catpanel access.

proc OGFMultifitFile {name} {return [file join [OGFSessWorkDir] multifit multifit_$name]}

proc OGFMultifitAfter {} {
    foreach {k f} {results_file results.tsv model_file model.fits residual_file residual.fits montage_file montage.png} {
	::ogf::cat::set multifit,$k [OGFMultifitFile $f]
    }
    if {[::ogf::params::get multifit show-frames]} {catch {OGFMultifitOpenFrames {residual.fits model.fits}}}
}

proc OGFMultifitOpenFrames {names} {
    global current scale
    set orig $current(frame)
    set n 0
    foreach nm $names {
	set f [OGFMultifitFile $nm]
	if {![file exists $f]} continue
	CreateFrame
	if {[catch {LoadFitsFile $f {} {}} err]} {::ogf::log ERROR "multifit: cannot load $f: $err"; continue}
	set scale(mode) zscale
	ChangeScaleMode
	incr n
    }
    catch {GotoFrame $orig}
    return $n
}

proc OGFMultifitFrames {} {
    set n [OGFMultifitOpenFrames {residual.fits model.fits}]
    ::ogf::status "Multi-fit: $n frame(s) opened (residual, model)"
}

# fit only the selected catalog rows: the NUMBERs go into the (recorded) parameter `objects`
proc OGFMultifitSelected {} {
    set sel [::ogf::cat::selection]
    if {$sel eq {}} {::ogf::status "Multi-fit: select catalog rows first"; return}
    ::ogf::params::put multifit objects [join $sel \;]
    catch {::ogf::params::save multifit}
    ::ogf::step::run multifit fit
}

proc OGFMultifitTable {} {
    set f [OGFMultifitFile results.tsv]
    if {![file exists $f]} {::ogf::status "Multi-fit: no results yet - run the fit first"; return}
    set fd [open $f r]; set txt [read $fd]; close $fd
    OGFTextWindow "Multi-fit components ($f)" $txt
}

proc OGFMultifitMontage {} {
    set png [OGFMultifitFile montage.png]
    if {![file exists $png]} {::ogf::status "Multi-fit: no montage yet - run the fit first"; return}
    set w .ogfmultifitmontage
    if {[winfo exists $w]} {destroy $w}
    toplevel $w
    wm title $w "Multi-fit: data / model / residual"
    catch {image delete ogfmultifitimg}
    image create photo ogfmultifitimg -file $png
    label $w.l -image ogfmultifitimg
    pack $w.l -fill both -expand 1
    ttk::button $w.close -text Close -command [list destroy $w]
    pack $w.close -pady 3
    return $w
}

# ---- automatic structure decomposition (step autodecomp: CLI multifit.py --model decomp --columns ad) -----------------------------------
proc OGFAutoDecompAfter {} {
    foreach {k f} {results_file results.tsv model_file model.fits residual_file residual.fits montage_file montage.png decomp_file decomp.tsv decomp_plot decomp_plot.png} {
	if {[file exists [OGFMultifitFile $f]]} {::ogf::cat::set multifit,$k [OGFMultifitFile $f]}
    }
    set f [OGFMultifitFile decomp.tsv]
    if {[file exists $f]} {
	set fd [open $f r]; set L [split [string trim [read $fd]] "\n"]; close $fd
	set h [split [lindex $L 0] "\t"]
	set it [lsearch $h AD_TYPE]
	array set n {1 0 2 0 3 0}
	foreach l [lrange $L 1 end] {
	    set t [lindex [split $l "\t"] $it]
	    if {[info exists n($t)]} {incr n($t)}
	}
	::ogf::cat::set multifit,decomp_counts [list $n(1) $n(2) $n(3)]
	::ogf::status "Auto decomposition: $n(1) single Sersic, $n(2) bulge+disc, $n(3) nucleus+bulge+disc"
    }
    if {[::ogf::params::get multifit show-frames]} {catch {OGFMultifitOpenFrames {residual.fits model.fits}}}
}

proc OGFAutoDecompTable {} {
    set f [OGFMultifitFile decomp.tsv]
    if {![file exists $f]} {::ogf::status "Auto decomposition: nothing yet - run the step first"; return}
    set fd [open $f r]; set txt [read $fd]; close $fd
    OGFTextWindow "Auto decomposition ($f)" $txt
}

proc OGFAutoDecompPlot {} {
    set png [OGFMultifitFile decomp_plot.png]
    if {![file exists $png]} {::ogf::status "Auto decomposition: no plot yet - run the step first"; return}
    set w .ogfdecompplot
    if {[winfo exists $w]} {destroy $w}
    toplevel $w
    wm title $w "Auto decomposition: isophote profiles and fitted components"
    catch {image delete ogfdecompimg}
    image create photo ogfdecompimg -file $png
    label $w.l -image ogfdecompimg
    pack $w.l -fill both -expand 1
    ttk::button $w.close -text Close -command [list destroy $w]
    pack $w.close -pady 3
    return $w
}
