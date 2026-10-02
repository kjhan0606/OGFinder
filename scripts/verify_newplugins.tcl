# GUI test of the analysis plugins added after the plugin restructuring (isophote, completeness, ...).  Real ds9, real m51 extraction,
# every step run through ::ogf::step::run (real job runner + recorder), outputs checked, session exported and replayed headless.
#   rm -rf ~/ds9.auto ~/ds9.auto.dir; HOME=/tmp/np_home OGF_NP_OUT=/tmp/np.txt OGF_NP_DIR=/tmp/np_work DISPLAY=:77 \
#       bin/ds9 /workspace/fits/m51.fits -geometry 1300x950 -source scripts/verify_newplugins.tcl
# OGF_NP_ONLY=isophote,completeness limits the sections.  Output: PASS/FAIL lines + SUMMARY failures=N.
global catpanel current ds9
set ::fh [open $::env(OGF_NP_OUT) w]; set ::nf 0
set ::dir $::env(OGF_NP_DIR); file mkdir $::dir
proc R {tag ok {d {}}} {puts $::fh "[expr {$ok?{PASS}:{FAIL}}] $tag $d"; flush $::fh; if {!$ok} {incr ::nf}}
proc bgerror {m} {puts $::fh "BGERROR $m"; flush $::fh; incr ::nf}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}
proc geom {} {
    update idletasks; update
    set tf [winfo parent $::catpanel(tbl)]
    return "[winfo rooty $tf] [winfo height $tf] [winfo height $::catpanel(infoarea)] [winfo width .]x[winfo height .]"
}
proc wait_job {{ms 300000}} {set t [clock milliseconds]; while {[::ogf::job::busy] && [clock milliseconds]-$t < $ms} {update; after 50}; wait_idle 300}
proc steps_since {n} {lmap r [lrange [::ogf::session::steps] $n end] {list [dict get $r step] [dict get $r class]}}
proc col_values {name} {return [::ogf::cat::values $name]}
proc nonempty {name} {set n 0; foreach v [col_values $name] {if {$v ne {}} {incr n}}; return $n}
proc run_step {plugin step} {
    set s0 [llength [::ogf::session::steps]]
    set ok [::ogf::step::run $plugin $step]
    wait_job
    return [list $ok [steps_since $s0]]
}
proc want {only name} {return [expr {$only eq {} || $name in $only}]}

proc sec_isophote {} {
    set f0 [llength $::ds9(frames)]
    ::ogf::params::put isophote max-objects 3
    ::ogf::params::put isophote step 0.15
    lassign [run_step isophote fit] ok recs
    R isophote_ran $ok $recs
    R isophote_recorded [expr {[lindex $recs 0 0] eq "analysis.isophote"}] $recs
    R isophote_columns [expr {[::ogf::cat::columns] ne {} && "ISO_EPS_HL" in [::ogf::cat::columns] && "ISO_SHAPE" in [::ogf::cat::columns]}]
    set n [nonempty ISO_NISO]
    R isophote_rows_filled [expr {$n >= 1 && $n <= 3}] "rows=$n"
    set w [OGFSessWorkDir]
    foreach f {isophote_model.fits isophote_resid.fits isophote_profiles.tsv isophote_plot.png isophote_profiles.json} {
	R isophote_file_$f [expr {[file exists [file join $w $f]] && [file size [file join $w $f]] > 100}]
    }
    R isophote_cat_key [expr {[file exists [::ogf::cat::get isophote,model_file {}]]}]
    R isophote_frames [expr {[llength $::ds9(frames)] == $f0 + 2}] "frames [llength $::ds9(frames)] (was $f0)"
    R isophote_frame_restored [expr {$::current(frame) eq [lindex $::ds9(frames) 0]}] $::current(frame)
    R isophote_argv_templated [expr {[string match {*@{WORK}/isophote_model.fits*} [dict get [lindex [::ogf::session::steps] end] argv_t]]}]
    set pw [OGFIsophotePlot]
    R isophote_plot_window [expr {[winfo exists $pw] && [image width ogfisoimg] > 300 && [image height ogfisoimg] > 200}] "[image width ogfisoimg]x[image height ogfisoimg]"
    destroy $pw
    OGFIsophoteTable; update
    R isophote_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc run {} {
    global catpanel
    set only [expr {[info exists ::env(OGF_NP_ONLY)] ? [split $::env(OGF_NP_ONLY) ,] : {}}]
    R geometry_before [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    CatalogPanelExtract
    set t0 [clock milliseconds]; while {![::ogf::cat::has] && [clock milliseconds]-$t0 < 90000} {update; after 100}
    wait_idle 500
    R extracted [expr {[::ogf::cat::nrows] > 50}] "rows=[::ogf::cat::nrows]"
    foreach sec {isophote} {
	if {[want $only $sec]} {
	    if {[catch {sec_$sec} err]} {R ${sec}_error 0 "$err [string range $::errorInfo 0 300]"}
	}
    }
    # export the recorded session and replay it headless: the new columns must be byte-identical to the GUI ones
    set sp [file join $::dir ogfinder_session.py]
    CatalogPanelSessionSave $sp
    set fd [open [file join $::dir gui_catalog.tsv] w]; puts -nonewline $fd [::ogf::cat::tsv]; close $fd
    R session_exported [file exists $sp]
    R final_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    puts $::fh "SUMMARY failures=$::nf"; close $::fh
}
after 3000 {
    if {[catch run err]} {puts $::fh "FAIL run_error $err $::errorInfo"; incr ::nf; puts $::fh "SUMMARY failures=$::nf"; close $::fh}
    exit 0
}
