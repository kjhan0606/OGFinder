# Transient light curves (plugins/lightcurves): result table and light-curve viewer.

proc OGFLcFile {name} {return [file join [OGFSessWorkDir] lightcurves $name]}

proc OGFLcAfter {} {
    foreach {k f} {results_file lc_results.tsv json_file lc_results.json} {
	set p [OGFLcFile $f]
	if {[file exists $p]} {::ogf::cat::set lightcurves,$k $p}
    }
}

proc OGFLcIds {} {
    set f [OGFLcFile lc_results.tsv]
    set ids {}
    if {![file exists $f]} {return $ids}
    set fd [open $f r]; set lines [split [read $fd] \n]; close $fd
    foreach l [lrange $lines 1 end] {if {$l ne {}} {lappend ids [lindex [split $l \t] 0]}}
    return $ids
}

proc OGFLcRender {id} {
    set cat [::ogf::cat::temp_file .tsv]
    set out [OGFLcFile view_$id.png]
    file mkdir [file dirname $out]
    set script [file join [::ogf::step::plugin_dir lightcurves] lightcurves.py]
    set a [list --task plot --id $id --out $out --work [file join [OGFSessWorkDir] lightcurves]]
    foreach {p opt} {source --source lc-file --lc-file moving-dir --moving-dir flux-unit --flux-unit zp --zp min-prob --min-prob min-points --min-points min-snr --min-snr} {
	lappend a $opt [::ogf::params::get lightcurves $p]
    }
    if {[::ogf::params::get lightcurves no-host]} {lappend a --no-host}
    if {[catch {exec [OGFPython] $script --catalog $cat {*}$a} err]} {
	::ogf::status "Light curves: $id: [lindex [split $err \n] end]"
	return {}
    }
    return $out
}

proc OGFLcViewer {{id {}}} {
    set ids [OGFLcIds]
    if {$ids eq {}} {::ogf::status "Light curves: nothing yet - run Classify Light Curves first"; return {}}
    if {$id eq {}} {set id [lindex $ids 0]}
    set w .ogflcview
    if {![winfo exists $w]} {
	toplevel $w
	wm title $w "Light-curve viewer"
	label $w.l
	pack $w.l -fill both -expand 1
	frame $w.b
	ttk::button $w.b.prev -text "< Prev" -command [list OGFLcStep -1]
	ttk::combobox $w.b.id -width 12 -values $ids -textvariable ::ogf_lc_id
	ttk::button $w.b.go -text Show -command {OGFLcShow $::ogf_lc_id}
	ttk::button $w.b.next -text "Next >" -command [list OGFLcStep 1]
	ttk::button $w.b.close -text Close -command [list destroy $w]
	pack $w.b.prev $w.b.id $w.b.go $w.b.next $w.b.close -side left -padx 3
	pack $w.b -pady 3
	bind $w.b.id <<ComboboxSelected>> {OGFLcShow $::ogf_lc_id}
	bind $w.b.id <Return> {OGFLcShow $::ogf_lc_id}
    } else {
	$w.b.id configure -values $ids
    }
    OGFLcShow $id
    return $w
}

proc OGFLcShow {id} {
    set ::ogf_lc_id $id
    set png [OGFLcRender $id]
    if {$png eq {} || ![file exists $png]} return
    catch {image delete ogflcimg}
    image create photo ogflcimg -file $png
    .ogflcview.l configure -image ogflcimg
    wm title .ogflcview "Light-curve viewer - $id"
}

proc OGFLcStep {d} {
    set ids [OGFLcIds]
    set i [lsearch -exact $ids $::ogf_lc_id]
    set i [expr {$i < 0 ? 0 : ($i + $d) % [llength $ids]}]
    OGFLcShow [lindex $ids $i]
}

proc OGFLcTable {} {
    set f [OGFLcFile lc_results.tsv]
    if {![file exists $f]} {::ogf::status "Light curves: nothing yet - run Classify Light Curves first"; return}
    set fd [open $f r]; set txt [read $fd]; close $fd
    OGFTextWindow "Light curves: results" $txt
}
