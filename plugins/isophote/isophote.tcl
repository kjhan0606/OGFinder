# Isophote plugin (plugins/isophote): Tcl side = apply the outputs of the CLI step (frames, plot window, table window).
# State is kept in catalog keys isophote,* through ::ogf::cat (registered in ogf_core.tcl); no direct catpanel access.

proc OGFIsophoteFile {kind} {
    return [file join [OGFSessWorkDir] isophote_$kind]
}

# job done: remember the files, optionally show model and residual
proc OGFIsophoteAfter {} {
    foreach {k f} {model_file model.fits resid_file resid.fits table_file profiles.tsv plot_file plot.png json_file profiles.json} {
	::ogf::cat::set isophote,$k [OGFIsophoteFile $f]
    }
    if {[::ogf::params::get isophote show-frames]} {catch {OGFIsophoteShowFrames}}
}

proc OGFIsophoteLoadFrame {fn} {
    CreateFrame
    if {[catch {LoadFitsFile $fn {} {}} err]} {::ogf::log ERROR "isophote: cannot load $fn: $err"; return 0}
    global scale
    set scale(mode) zscale
    ChangeScaleMode
    return 1
}

# model (cmodel) and residual as two new frames; the original frame is made current again
proc OGFIsophoteShowFrames {} {
    global current
    set mf [OGFIsophoteFile model.fits]
    set rf [OGFIsophoteFile resid.fits]
    if {![file exists $mf]} {::ogf::status "Isophote: no model yet - run Fit Isophotes first"; return}
    set orig $current(frame)
    set n 0
    foreach f [list $mf $rf] {if {[file exists $f] && [OGFIsophoteLoadFrame $f]} {lappend ::ogf_iso_frames $current(frame); incr n}}
    catch {GotoFrame $orig}
    ::ogf::status "Isophote: $n frame(s) opened (model, residual)"
}

# plot window: the PNG written by the CLI (matplotlib Agg) shown in a Tk label
proc OGFIsophotePlot {} {
    set pf [OGFIsophoteFile plot.png]
    if {![file exists $pf]} {::ogf::status "Isophote: no plot yet - run Fit Isophotes first"; return}
    set w .ogfisoplot
    if {[winfo exists $w]} {destroy $w}
    toplevel $w
    wm title $w "Isophote profiles"
    catch {image delete ogfisoimg}
    image create photo ogfisoimg -file $pf
    label $w.l -image ogfisoimg
    pack $w.l -fill both -expand 1
    ttk::button $w.close -text Close -command [list destroy $w]
    pack $w.close -pady 3
    return $w
}

proc OGFIsophoteTable {} {
    set tf [OGFIsophoteFile profiles.tsv]
    if {![file exists $tf]} {::ogf::status "Isophote: no table yet - run Fit Isophotes first"; return}
    set fd [open $tf r]; set txt [read $fd]; close $fd
    OGFTextWindow "Isophote profiles ($tf)" $txt
}
