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

# ---- depth / completeness maps (depthmap.py): files in <work>/depth
proc OGFDepthFile {name} {return [file join [OGFSessWorkDir] depth $name]}

proc OGFDepthAfter {} {
    foreach {k f} {mag_map depth_mag_map.fits sb_map depth_sb_map.fits rms_map depth_rms_map.fits tiles_file depth_tiles.tsv area_file depth_area.tsv
		   lim50_map compmap_lim50_map.fits lim90_map compmap_lim90_map.fits regions_file compmap_regions.tsv regions_map compmap_regions.fits
		   json_file depth_summary.json comp_json compmap_summary.json plot_file depth_plot.png comp_plot compmap_plot.png} {
	set p [OGFDepthFile $f]
	if {[file exists $p]} {::ogf::cat::set depth,$k $p}
    }
    set jf [OGFDepthFile depth_summary.json]
    if {[file exists $jf]} {
	set fd [open $jf r]; set txt [read $fd]; close $fd
	if {[regexp {"mag_limit": \{\s*"median": ([-0-9.eE+]+)} $txt -> v]} {::ogf::cat::set depth,mag_limit_median $v}
	if {[regexp {"sb_limit": \{[^\}]*?"median": ([-0-9.eE+]+)} $txt -> v]} {::ogf::cat::set depth,sb_limit_median $v}
	::ogf::status "Depth: median limiting magnitude [::ogf::cat::get depth,mag_limit_median ?], SB limit [::ogf::cat::get depth,sb_limit_median ?]"
    }
}

proc OGFDepthShow {} {
    global current scale
    set orig $current(frame)
    set n 0
    foreach nm {depth_mag_map.fits depth_sb_map.fits compmap_lim50_map.fits} {
	set f [OGFDepthFile $nm]
	if {![file exists $f]} continue
	CreateFrame
	if {[catch {LoadFitsFile $f {} {}} err]} {::ogf::log ERROR "depth: cannot load $f: $err"; continue}
	set scale(mode) minmax
	ChangeScaleMode
	incr n
    }
    catch {GotoFrame $orig}
    ::ogf::status "Depth maps: $n frame(s) opened"
    return $n
}

proc OGFDepthTable {} {
    set txt {}
    foreach nm {compmap_regions.tsv depth_tiles.tsv} {
	set f [OGFDepthFile $nm]
	if {![file exists $f]} continue
	set fd [open $f r]; append txt "== $nm ==\n[read $fd]\n"; close $fd
    }
    if {$txt eq {}} {::ogf::status "Depth: nothing yet - run Depth Maps first"; return}
    OGFTextWindow "Depth tiles / completeness regions" $txt
}

proc OGFDepthPlot {} {
    set png [OGFDepthFile compmap_plot.png]
    if {![file exists $png]} {set png [OGFDepthFile depth_plot.png]}
    if {![file exists $png]} {::ogf::status "Depth: no plot yet - run Depth Maps first"; return}
    set w .ogfdepthplot
    if {[winfo exists $w]} {destroy $w}
    toplevel $w
    wm title $w "Depth map and completeness by region"
    catch {image delete ogfdepthimg}
    image create photo ogfdepthimg -file $png
    label $w.l -image ogfdepthimg
    pack $w.l -fill both -expand 1
    ttk::button $w.close -text Close -command [list destroy $w]
    pack $w.close -pady 3
    return $w
}
