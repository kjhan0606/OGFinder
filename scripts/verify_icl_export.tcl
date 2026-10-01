# ICL bash-export smoke test: ds9 m51.fits -source verify_icl_export.tcl   (env OGF_VERIFY_OUT)
global current catpanel ds9
set out $::env(OGF_VERIFY_OUT)
file mkdir $out
set ::vlog [open [file join $out icl_export_steps.log] w]
proc vlog {m} {puts $::vlog $m; flush $::vlog; puts "VERIFY: $m"}
proc wait_until {cond {ms 600000}} {
    set t0 [clock milliseconds]
    while {![uplevel #0 $cond]} {update; after 50; if {[clock milliseconds]-$t0 > $ms} {error "timeout"}}
}
proc tk_getSaveFile {args} {return $::env(OGF_VERIFY_OUT)/icl_pipeline.sh}
proc run {} {
    global catpanel ds9
    set catpanel(param,detect-thresh) 2.0
    CatalogPanelExtract
    CatalogPanelICLBackground polynomial
    set ::icl_bcg_entry_id 1
    global catpanel_fdata
    set ff [lindex $ds9(frames) 0]
    set cd $catpanel_fdata($ff,alldata)
    set hdr [split [lindex [split $cd \n] 0] \t]
    CatalogPanelICLSetCenterByID .none $cd [lsearch $hdr NUMBER] [lsearch $hdr X_IMAGE] [lsearch $hdr Y_IMAGE]
    set catpanel(icl,param,rmax) 150
    CatalogPanelICLProfile
    CatalogPanelICLMeasure
    vlog "cmdlog entries: [llength $catpanel(icl,cmdlog)]"
    CatalogPanelICLExportScript
    vlog "export status: $catpanel(status)"
    CatalogPanelSaveCatalogTo [file join $::env(OGF_VERIFY_OUT) icl_measure_export.tsv]
    # session recorder in parallel
    vlog "session steps: [llength $::ogfsess(steps)]"
}
after 3000 {
    if {[catch run err]} {vlog "ERROR: $err\n$::errorInfo"; exit 3}
    vlog DONE
    exit 0
}
