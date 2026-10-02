# Completeness plugin (plugins/completeness): Tcl side = keys, plot window, table window.  Catalog state only through ::ogf::cat.
proc OGFCompletenessFile {name} {return [file join [OGFSessWorkDir] $name]}

proc OGFCompletenessAfter {} {
    foreach {k f} {json_file completeness.json curve_file completeness_curve.tsv plot_file completeness_plot.png} {
	::ogf::cat::set completeness,$k [OGFCompletenessFile $f]
    }
    set jf [OGFCompletenessFile completeness.json]
    if {![file exists $jf]} return
    set fd [open $jf r]; set txt [read $fd]; close $fd
    # the two limits are the first "lim50" / "lim90" keys of the result (top level)
    foreach k {lim50 lim90} {
	if {[regexp "\"$k\":\\s*(\[-0-9.eE+\]+|null)" $txt -> v]} {::ogf::cat::set completeness,$k $v}
    }
    foreach k {lim50 lim90} {if {[::ogf::cat::get completeness,$k ?] eq "null"} {::ogf::cat::set completeness,$k none}}
    ::ogf::status "Completeness: 50% limit [::ogf::cat::get completeness,lim50 ?], 90% limit [::ogf::cat::get completeness,lim90 ?] mag (none = not reached in the magnitude range)"
}

proc OGFCompletenessPlot {} {
    set pf [OGFCompletenessFile completeness_plot.png]
    if {![file exists $pf]} {::ogf::status "Completeness: no plot yet - run the simulation first"; return}
    set w .ogfcompplot
    if {[winfo exists $w]} {destroy $w}
    toplevel $w
    wm title $w "Completeness curve"
    catch {image delete ogfcompimg}
    image create photo ogfcompimg -file $pf
    label $w.l -image ogfcompimg
    pack $w.l -fill both -expand 1
    ttk::button $w.close -text Close -command [list destroy $w]
    pack $w.close -pady 3
    return $w
}

proc OGFCompletenessTable {} {
    set tf [OGFCompletenessFile completeness_curve.tsv]
    if {![file exists $tf]} {::ogf::status "Completeness: no table yet"; return}
    set fd [open $tf r]; set txt [read $fd]; close $fd
    OGFTextWindow "Completeness curve" $txt
}
