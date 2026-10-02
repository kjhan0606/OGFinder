# Equivalence test for the plugin steps that were converted from Tcl procs to `cli` argv templates (item 3).
# For each converted step it runs (A) the legacy proc, which is kept as the oracle, with `exec` stubbed, and (B) the manifest step through
# ::ogf::step::run (real job runner, real recorder) with OGFINDER_PYTHON = scripts/fake_python_for_templates.sh, and compares
#   the argv that reaches the interpreter, the recorder record (step, class, title, templated argv, post, requires), the catalog after the step,
#   and the final status text.   HOME must be a scratch directory.
#   HOME=/tmp/ogf_clitpl_home OGF_CLITPL_OUT=/tmp/clitpl.txt FAKE_LOG=/tmp/clitpl_fake.log OGFINDER_PYTHON=$PWD/scripts/fake_python_for_templates.sh \
#     DISPLAY=:77 bin/ds9 /workspace/fits/m51.fits -geometry 1300x950 -source scripts/verify_cli_templates.tcl
global catpanel current
set ::fh [open $::env(OGF_CLITPL_OUT) w]; set ::nf 0
set ::HOME [file normalize ~]
proc R {tag ok {d {}}} {puts $::fh "[expr {$ok?{PASS}:{FAIL}}] $tag $d"; flush $::fh; if {!$ok} {incr ::nf}}
proc bgerror {m} {puts $::fh "BGERROR $m"; flush $::fh}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}
proc geom {} {
    update idletasks; update
    set tf [winfo parent $::catpanel(tbl)]
    return "[winfo rooty $tf] [winfo height $tf] [winfo height $::catpanel(infoarea)] [winfo width .]x[winfo height .]"
}
rename exec ::real_exec
set ::EXEC {}; set ::STUB 0
proc exec {args} {
    set a {}
    foreach x $args {if {[string match 2>* $x]} continue; lappend a $x}
    if {$::STUB && [llength $a] > 1 && [string match *.py [lindex $a 1]]} {
	lappend ::EXEC $a
	return [::real_exec {*}$a]          ;# the fake interpreter prints the canned TSV
    }
    return [::real_exec {*}$args]
}
proc norm {s} {string map [list $::HOME <HOME> /workspace/fits <FITS> [file normalize bin] <BIN>] $s}
proc rec_sig {from} {
    set out {}
    foreach r [lrange [::ogf::session::steps] $from end] {
	lappend out "[dict get $r step] [dict get $r class] | [norm [dict get $r argv_t]] | post=[dict get $r post] title=[dict get $r title] requires=[dict get $r requires] inputs=[norm [dict get $r inputs]] outputs=[norm [dict get $r outputs]]"
    }
    return $out
}
proc wait_job {} {set t [clock milliseconds]; while {[::ogf::job::busy] && [clock milliseconds]-$t < 60000} {update; after 50}; wait_idle 200}

proc compare {name legacy plugin step} {
    global catpanel
    set snap $::BASE
    ::ogf::cat::load_tsv $snap restored; update; wait_idle 100
    # ---- A: legacy proc
    set s0 [llength [::ogf::session::steps]]; set e0 [llength $::EXEC]; set ::STUB 1
    set rc [catch {uplevel #0 $legacy} err]
    update; wait_idle 200; set ::STUB 0
    set A_argv [lrange $::EXEC $e0 end]; set A_rec [rec_sig $s0]; set A_tsv $catpanel(alldata); set A_status $catpanel(status)
    R ${name}_legacy_ran [expr {!$rc && [llength $A_argv] == 1}] "$err [llength $A_argv]"
    # ---- B: manifest step
    ::ogf::cat::load_tsv $snap restored
    update; wait_idle 100
    set s1 [llength [::ogf::session::steps]]
    set fl $::env(FAKE_LOG); set n0 0
    catch {set f [open $fl r]; set n0 [llength [split [string trim [read $f]] \n]]; close $f}
    set ok [::ogf::step::run $plugin $step]
    wait_job
    set B_argv {}
    catch {set f [open $fl r]; set B_argv [lrange [split [string trim [read $f]] \n] $n0 end]; close $f}
    set B_rec [rec_sig $s1]; set B_tsv $catpanel(alldata); set B_status $catpanel(status)
    R ${name}_step_started [expr {$ok == 1 && [llength $B_argv] == 1}] "ok=$ok n=[llength $B_argv]"
    set a_line [join [lrange [lindex $A_argv 0] 1 end] { }]; set b_line [lindex $B_argv 0]      ;# the fake interpreter logs its arguments without argv[0]
    R ${name}_argv_identical [expr {$a_line eq $b_line}] "\n   A: [norm $a_line]\n   B: [norm $b_line]"
    R ${name}_record_identical [expr {$A_rec eq $B_rec}] "\n   A: [norm $A_rec]\n   B: [norm $B_rec]"
    R ${name}_catalog_identical [expr {$A_tsv eq $B_tsv && $A_tsv ne $snap}] "cols=[llength [split [lindex [split $B_tsv \n] 0] \t]]"
    R ${name}_status_identical [expr {$A_status eq $B_status}] "A='$A_status' B='$B_status'"
}

proc run {} {
    global catpanel current
    update; wait_idle 500
    if {![string match /tmp/* $::HOME]} {puts $::fh "SUMMARY-ABORT HOME not scratch"; close $::fh; exit 2}
    R geometry_before [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    CatalogPanelExtract
    set t0 [clock milliseconds]; while {$catpanel(alldata) eq {} && [clock milliseconds]-$t0 < 60000} {update; after 100}
    wait_idle 300
    R extracted [expr {[::ogf::cat::nrows] > 50}] [::ogf::cat::nrows]
    set ::BASE $catpanel(alldata)
    set catpanel(param,n-workers) 2; set catpanel(param,mag-zeropoint) 24.5
    set catpanel(bd,param,max-sources) 50; set catpanel(bd,param,free-bulge-n) 0; set catpanel(bd,param,mag-zeropoint) 25.0; set catpanel(bd,param,pixel-scale) 0.05
    R fake_python [expr {[info exists ::env(OGFINDER_PYTHON)] && [string match *fake_python_for_templates.sh [OGFPython]]}] [OGFPython]
    compare sersic CatalogPanelSersicFit morphology sersic
    compare morphometry CatalogPanelMorphometry morphology morphometry
    compare bulge_disk CatalogPanelBulgeDisk morphology bulge_disk
    set psf [file join $::HOME fake_psf.fits]; set fd [open $psf w]; puts $fd x; close $fd
    set catpanel(psf,file) $psf; set catpanel(psf,has_psf) 1
    compare bulge_disk_psf CatalogPanelBulgeDisk morphology bulge_disk
    set catpanel(bd,param,free-bulge-n) 1; set catpanel(param,mag-zeropoint) 21.5; set catpanel(param,n-workers) 5
    compare bulge_disk_free_n CatalogPanelBulgeDisk morphology bulge_disk
    compare sersic_changed_params CatalogPanelSersicFit morphology sersic
    compare psf_phot CatalogPanelPSFPhotometry photometry psf_phot
    compare crowded CatalogPanelCrowdedPhot photometry crowded
    # preconditions: the step refuses (no job started) without a PSF, like the legacy proc
    set catpanel(psf,has_psf) 0; set catpanel(psf,file) {}
    set n0 [llength [::ogf::session::steps]]
    R psf_step_refuses_without_psf [expr {[::ogf::step::run photometry psf_phot] == 0 && [llength [::ogf::session::steps]] == $n0}] $catpanel(status)
    R geometry_after [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    puts $::fh "SUMMARY failures=$::nf"; close $::fh; exit
}
after 3000 run
