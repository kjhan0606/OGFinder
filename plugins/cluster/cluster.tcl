# Cluster and lensing (plugins/cluster): outputs of cluster.py -> catalog keys cluster,*, CMD plot, summary, density-map frame.

proc OGFClusterFile {name} {return [file join [OGFSessWorkDir] cluster $name]}

proc OGFClusterAfter {} {
    foreach {k f} {rs_png cluster_rs.png density_file cluster_density.fits peaks_file cluster_peaks.tsv json_file cluster_summary.json} {
	set p [OGFClusterFile $f]
	if {[file exists $p]} {::ogf::cat::set cluster,$k $p}
    }
    if {[::ogf::params::get cluster show-frames]} {catch {OGFClusterFrame}}
}

proc OGFClusterPlot {} {
    set png [OGFClusterFile cluster_rs.png]
    if {![file exists $png]} {::ogf::status "Cluster: no plot yet - run Red Sequence + Members first"; return}
    set w .ogfclusterplot
    if {[winfo exists $w]} {destroy $w}
    toplevel $w
    wm title $w "Colour-magnitude diagram and red sequence"
    catch {image delete ogfclusterimg}
    image create photo ogfclusterimg -file $png
    label $w.l -image ogfclusterimg
    pack $w.l -fill both -expand 1
    ttk::button $w.close -text Close -command [list destroy $w]
    pack $w.close -pady 3
    return $w
}

proc OGFClusterSummary {} {
    set j [OGFClusterFile cluster_summary.json]
    if {![file exists $j]} {::ogf::status "Cluster: nothing yet - run a cluster step first"; return}
    set fd [open $j r]; set js [read $fd]; close $fd
    OGFTextWindow "Cluster summary" $js
}

proc OGFClusterFrame {} {
    global current scale
    set f [OGFClusterFile cluster_density.fits]
    if {![file exists $f]} {::ogf::status "Cluster: no density map yet - run Overdensity Map first"; return 0}
    set orig $current(frame)
    CreateFrame
    if {[catch {LoadFitsFile $f {} {}} err]} {::ogf::log ERROR "cluster: cannot load $f: $err"; catch {GotoFrame $orig}; return 0}
    set scale(mode) minmax
    ChangeScaleMode
    catch {GotoFrame $orig}
    ::ogf::status "Cluster: overdensity significance map opened in a new frame"
    return 1
}
