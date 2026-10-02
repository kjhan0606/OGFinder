# Drives the Moving Objects GUI steps headlessly (align -> difference -> link -> identify -> orbit -> transients -> light curve ->
# export) in a real ds9 on the cached BB89 exposures and saves the recorded session script.  Used by scripts/verify_moving_session_auto.sh
# to make the reference outputs for the replay comparison (run_all_checks: moving_session).  Needs network (SkyBoT, Horizons, Gaia).
#   MOV_OUT=/tmp/mov HOME=/tmp/movhome DISPLAY=:77 bin/ds9 -source scripts/verify_moving_record.tcl
global ogfmov catpanel ogfsess
set out $::env(MOV_OUT)
file mkdir $out
set lg [open [file join $out gui.log] w]
proc L {m} {puts $::lg $m; flush $::lg; puts "MOV: $m"}
proc waitidle {{ms 1800000}} {
    global ogfmov
    set t0 [clock milliseconds]
    after 300
    while {$ogfmov(busy)} {
	update; after 50
	if {[clock milliseconds]-$t0 > $ms} {error timeout}
    }
    update
}
proc go {} {
    global ogfmov ogfsess catpanel
    OGFSessInit
    set ogfmov(workdir) [file join [file normalize ~] .ds9 moving_work]
    file delete -force $ogfmov(workdir)
    set ogfmov(files) [lsort [glob ~/.ds9/mast_cache/j8pu38*_flc.fits]]
    set ogfmov(ra) 150.1375; set ogfmov(dec) 2.3610
    set seq {OGFMovAlign OGFMovDifference OGFMovLink OGFMovIdentify}
    foreach p $seq {
	set t0 [clock milliseconds]
	$p
	waitidle
	L "$p done [expr {[clock milliseconds]-$t0}] ms status=$ogfmov(status)"
    }
    # choose the tracklet that SkyBoT identified as 2015 BB89, else the best-ranked
    set wd [OGFMovWork]
    set tr 0
    set fd [open [file join $wd identified.tsv]]; set txt [read $fd]; close $fd
    foreach l [lrange [split $txt \n] 1 end] {
	set f [split $l \t]
	if {[string match *BB89* [lindex $f 2]]} {set tr [lindex $f 0]; break}
    }
    L "orbit tracklet $tr"
    set ogfmov(tracklet) $tr; set ogfmov(designation) {}
    OGFMovOrbitRun; waitidle; L "orbit done status=$ogfmov(status)"
    OGFMovTransients_run
    OGFMovLightCurve; waitidle; L "lightcurve done"
    OGFMovExport; waitidle; L "export done"
    set n [CatalogPanelSessionSave [file join $::env(MOV_OUT) ogfinder_session.py]]
    L "session saved [llength $ogfsess(steps)] steps"
    set fd [open [file join $::env(MOV_OUT) session_log.txt] w]; puts $fd [OGFSessLogText]; close $fd
    exit
}
proc OGFMovTransients_run {} {
    global ogfmov
    # headless equivalent of the dialog: same argv as OGFMovTransients builds
    set wd [OGFMovWork]
    set extra [list --mode transients --workdir $wd --min-epochs $ogfmov(minep) --snr $ogfmov(snr)]
    OGFMovRun moving.transients "Transient candidates" $extra OGFMovTransientsDone 0 [list [file join $wd transients.tsv] [file join $wd transients.reg]]
    waitidle; L "transients done"
}
after 3000 {if {[catch go err]} {L "ERROR $err $::errorInfo"; exit 1}}
