# Reproducibility bundle (plugins/repro): gathers the GUI state and calls plugins/repro/repro.py.

proc OGFReproScript {} {return [file join [::ogf::step::plugin_dir repro] repro.py]}

proc OGFReproFile {} {
    set f [::ogf::params::get repro bundle-file]
    if {$f eq {}} {set f [file join [OGFSessWorkDir] ogf_bundle.zip]}
    return $f
}

proc OGFReproParamsJSON {} {
    set parts {}
    foreach id [lsort [::ogf::reg::ids]] {
	set m [::ogf::reg::get $id]
	set kv {}
	foreach p [::ogf::json::get $m params] {
	    set n [dict get $p name]
	    if {[catch {::ogf::params::get $id $n} v]} continue
	    lappend kv "[OGFJStr $n]: [OGFJStr $v]"
	}
	if {[llength $kv]} {lappend parts "[OGFJStr $id]: \{[join $kv {, }]\}"}
    }
    return "\{[join $parts ",\n"]\}"
}

proc OGFReproBundle {} {
    set out [OGFReproFile]
    file mkdir [file dirname $out]
    set d [file join [OGFSessWorkDir] repro_tmp]
    file mkdir $d
    set cat [CatalogPanelSaveTempCatalog .tsv]
    set cf [file join $d catalog_final.tsv]
    file copy -force $cat $cf
    set fh [open [file join $d params.json] w]; puts $fh [OGFReproParamsJSON]; close $fh
    set recs {}
    set images {}
    foreach rec [::ogf::session::steps] {
	lappend recs [OGFSessRecordJSON $rec]
	catch {foreach im [dict get $rec images] {if {$im ni $images && [file isfile $im]} {lappend images $im}}}
    }
    set fh [open [file join $d steps.json] w]; puts $fh "\[[join $recs ",\n"]\]"; close $fh
    set a [list [OGFPython] [OGFReproScript] bundle --root [OGFSessRoot] --out $out --kind catalog --catalog $cf --params [file join $d params.json] --steps [file join $d steps.json] \
	--note [::ogf::params::get repro note]]
    if {[llength $images]} {lappend a --inputs [join $images ,]}
    if {[::ogf::params::get repro include-inputs]} {lappend a --include-inputs}
    set sess [file join $d session.py]
    if {[llength $recs] && ![catch {OGFSessExport $sess}] && [file exists $sess]} {
	set a [lreplace $a [lsearch $a catalog] [lsearch $a catalog] session]
	lappend a --session-script $sess
    }
    if {[catch {exec {*}$a} msg]} {
	::ogf::status "Bundle: failed: [lindex [split $msg \n] end]"
	return {}
    }
    ::ogf::status "Bundle written: $out"
    return $out
}

proc OGFReproVerify {} {
    set f [OGFReproFile]
    if {![file exists $f]} {::ogf::status "Verify: no bundle at $f"; return}
    set a [list [OGFPython] [OGFReproScript] verify $f --root [OGFSessRoot] --python [OGFPython] --rtol [::ogf::params::get repro rtol]]
    if {![::ogf::params::get repro rerun]} {lappend a --no-rerun}
    if {[::ogf::params::get repro inputs-dir] ne {}} {lappend a --inputs-dir [::ogf::params::get repro inputs-dir]}
    set rc [catch {exec {*}$a} out]
    set ::ogf::repro::last_rc [expr {$rc ? 1 : 0}]
    set ::ogf::repro::last_out $out
    OGFTextWindow "Bundle verification" $out
    ::ogf::status "Verify: [lindex [split [string trim $out] \n] end]"
    return [expr {!$rc}]
}

proc OGFReproInfo {} {
    set f [OGFReproFile]
    if {![file exists $f]} {::ogf::status "Bundle info: no bundle at $f"; return}
    OGFTextWindow "Bundle info" [exec [OGFPython] [OGFReproScript] info $f --root [OGFSessRoot]]
}

namespace eval ::ogf::repro {variable last_rc {} last_out {}}
