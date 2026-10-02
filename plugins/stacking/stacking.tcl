# Stacking plugin (plugins/stacking): Tcl side = keys, plot window, profile table window, stack image in a frame.
proc OGFStackFile {name} {return [file join [OGFSessWorkDir] stacking stacking_$name]}

proc OGFStackAfter {} {
    foreach {k f} {stack_file stack.fits err_file err.fits cube_file cutouts.fits profile_file profile.tsv json_file summary.json plot_file plot.png} {
	set p [OGFStackFile $f]
	if {[file exists $p]} {::ogf::cat::set stacking,$k $p}
    }
    set j [OGFStackFile summary.json]
    if {[file exists $j]} {
	set fd [open $j r]; set js [read $fd]; close $fd
	foreach k {n_stacked aperture_sum aperture_err} {
	    if {[regexp "\"$k\": (\[-0-9.eE+\]+)" $js -> v]} {::ogf::cat::set stacking,$k $v}
	}
	::ogf::status "Stack: [::ogf::cat::get stacking,n_stacked ?] sources, aperture sum [::ogf::cat::get stacking,aperture_sum ?] +- [::ogf::cat::get stacking,aperture_err ?]"
    }
    if {[::ogf::params::get stacking show-stack]} {catch {OGFStackShow}}
}

proc OGFStackShow {} {
    global current scale
    set f [OGFStackFile stack.fits]
    if {![file exists $f]} {::ogf::status "Stacking: nothing yet - run Stack Sources first"; return 0}
    set orig $current(frame)
    CreateFrame
    if {[catch {LoadFitsFile $f {} {}} err]} {::ogf::log ERROR "stacking: cannot load $f: $err"; catch {GotoFrame $orig}; return 0}
    set scale(mode) zscale
    ChangeScaleMode
    catch {GotoFrame $orig}
    ::ogf::status "Stacking: stack image opened in a new frame"
    return 1
}

proc OGFStackTable {} {
    set f [OGFStackFile profile.tsv]
    if {![file exists $f]} {::ogf::status "Stacking: no profile yet"; return}
    set fd [open $f r]; set txt [read $fd]; close $fd
    OGFTextWindow "Stack radial profile" $txt
}

proc OGFStackPlot {} {
    set png [OGFStackFile plot.png]
    if {![file exists $png]} {::ogf::status "Stacking: no plot yet - run Stack Sources first"; return}
    set w .ogfstackplot
    if {[winfo exists $w]} {destroy $w}
    toplevel $w
    wm title $w "Stack and mean radial profile"
    catch {image delete ogfstackimg}
    image create photo ogfstackimg -file $png
    label $w.l -image ogfstackimg
    pack $w.l -fill both -expand 1
    ttk::button $w.close -text Close -command [list destroy $w]
    pack $w.close -pady 3
    return $w
}
