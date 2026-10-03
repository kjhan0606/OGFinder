# GUI test of the trails plugin ("Remove Trails" button): detect, mask (shared mask manager), overlay, catalogue flags, result panel, Undo, preview,
# interpolate option, chip button geometry.  Run through scripts/verify_trails.sh (synthetic field with two injected trails).
global catpanel current ds9 ogfmask ogfui
set ::fh [open $::env(OGF_TR_OUT) w]; set ::nf 0
set ::dir $::env(OGF_TR_DIR); file mkdir $::dir
set ::sim $::env(OGF_TR_SIM)
proc R {tag ok {d {}}} {puts $::fh "[expr {$ok?{PASS}:{FAIL}}] $tag $d"; flush $::fh; if {!$ok} {incr ::nf}}
proc bgerror {m} {puts $::fh "BGERROR $m"; flush $::fh; incr ::nf}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}
proc geom {} {
    update idletasks; update
    set tf [winfo parent $::catpanel(tbl)]
    return "[winfo rooty $tf] [winfo height $tf] [winfo height $::catpanel(infoarea)] [winfo width .]x[winfo height .]"
}
proc wait_job {{ms 300000}} {set t [clock milliseconds]; while {[::ogf::job::busy] && [clock milliseconds]-$t < $ms} {update; after 50}; wait_idle 400}
proc steps_since {n} {lmap r [lrange [::ogf::session::steps] $n end] {list [dict get $r step] [dict get $r class]}}
proc run_step {plugin step} {
    set s0 [llength [::ogf::session::steps]]
    set ok [::ogf::step::run $plugin $step]
    wait_job
    return [list $ok [steps_since $s0]]
}
proc readjson {f} {set fd [open $f r]; set t [read $fd]; close $fd; return [::ogf::json::parse $t]}
proc maskstats {} {
    set p [OGFMaskPaths]
    set out [exec $::env(OGFINDER_PYTHON) [CatalogPanelGetScript ds9_mask.py] [dict get $p fits] --mode stats --mask [dict get $p mask]]
    return $out
}
proc ntrail {} {if {[regexp {N_TRAIL=(\d+)} [maskstats] -> n]} {return $n}; return -1}

proc run {} {
    global catpanel current ds9 ogfmask ogfui
    set truth [readjson [file join $::sim truth.json]]
    set g0 [geom]
    CatalogPanelExtract
    set t0 [clock milliseconds]; while {![::ogf::cat::has] && [clock milliseconds]-$t0 < 90000} {update; after 100}
    wait_idle 500
    R extracted [expr {[::ogf::cat::nrows] >= 20}] "rows=[::ogf::cat::nrows]"
    # ---- the button: chip in the Detect tab with a run button
    R chip_registered [expr {[::ogf::reg::get trails] ne {}}]
    R chip_run_button [expr {[info exists ogfui(run,trails)] && [winfo exists $ogfui(run,trails)]}]
    R chip_menu_has_remove [expr {[info exists ogfui(menu,trails)] && [$ogfui(menu,trails) index "Remove Trails"] ne "none"}]
    R chip_menu_has_undo [expr {[$ogfui(menu,trails) index "Undo Trail Removal"] ne "none"}]
    R chip_label_text [expr {[string match "Remove Trails*" [[winfo parent $ogfui(run,trails)].m cget -text]]}] [[winfo parent $ogfui(run,trails)].m cget -text]
    OGFUIShowTab Detect
    wait_idle 300
    R chip_mapped_in_detect_tab [winfo ismapped $ogfui(run,trails)]
    R geometry_unchanged_by_chip [expr {[geom] eq $g0}] "$g0 -> [geom]"
    # ---- one click: invoke the run button
    set s0 [llength [::ogf::session::steps]]
    $ogfui(run,trails) invoke
    wait_job
    set recs [steps_since $s0]
    R remove_recorded [expr {[lindex $recs 0 0] eq "trails.remove"}] $recs
    set w [file join [OGFSessWorkDir] trails]
    foreach f {trails.json trails_mask.fits trails_overlay.reg trails_flags.tsv} {R file_$f [expr {[file exists [file join $w $f]] && [file size [file join $w $f]] > 20}]}
    set d [readjson [file join $w trails.json]]
    R two_trails_found [expr {[dict get $d n_trails] == 2}] "n=[dict get $d n_trails] best=[dict get $d best_zscore]"
    set tr [dict get $d trails]
    set widths {}
    foreach t $tr {lappend widths [format %.1f [dict get $t fwhm]]}
    R widths_reasonable [expr {[lindex $widths 0] > 2 && [lindex $widths 0] < 9}] $widths
    R trail_mask_bit [expr {[ntrail] > 1000}] [maskstats]
    R overlay_on [expr {$ogfmask(overlay) == 1}]
    set mp [dict get [OGFMaskPaths] bool]
    R icl_lsbg_mask_synced [expr {[::ogf::cat::get icl,mask_file {}] eq $mp && [::ogf::cat::get lsbg,mask_file {}] eq $mp}]
    set cols [::ogf::cat::columns]
    R catalog_columns [expr {"TRAIL_FLAG" in $cols && "TRAIL_ID" in $cols && "TRAIL_DIST" in $cols && "TRAIL_FLUX" in $cols && "TRAIL_FRAC" in $cols}] $cols
    set fl [::ogf::cat::values TRAIL_FLAG]
    set nflag 0; foreach v $fl {if {$v ne {} && $v > 0} {incr nflag}}
    R some_objects_flagged [expr {$nflag >= 3 && $nflag < [llength $fl]}] "flagged $nflag of [llength $fl]"
    R result_panel [winfo exists .ogftrails]
    # ---- Undo from the panel
    .ogftrails.b.undo invoke
    wait_idle 600
    R undo_restores_mask [expr {[ntrail] == 0}] [maskstats]
    set fl [::ogf::cat::values TRAIL_FLAG]
    set nflag 0; foreach v $fl {if {$v ne {} && $v > 0} {incr nflag}}
    R undo_resets_columns [expr {$nflag == 0}] "flagged $nflag"
    R panel_closed_after_undo [expr {![winfo exists .ogftrails]}]
    # ---- preview does not touch the mask
    lassign [run_step trails preview] ok recs
    R preview_ran [expr {$ok && [lindex $recs 0 0] eq "trails.preview"}] $recs
    R preview_leaves_mask [expr {[ntrail] == 0}]
    # ---- interpolate option
    ::ogf::params::put trails fill interpolate
    ::ogf::params::put trails fill-noise 1
    ::ogf::params::put trails confirm-panel 0
    lassign [run_step trails remove] ok recs
    R interpolate_ran [expr {$ok}] $recs
    set d [readjson [file join $w trails.json]]
    set fi [file join $w [file rootname [file tail [dict get [OGFMaskPaths] fits]]]_trailfree.fits]
    R filled_image_written [expr {[file exists $fi] && [file size $fi] > 1000}] $fi
    R mask_set_again [expr {[ntrail] > 1000}]
    R no_panel_when_disabled [expr {![winfo exists .ogftrails]}]
    # stack option writes a MEF
    ::ogf::params::put trails fill stack
    lassign [run_step trails remove] ok recs
    R mef_written [expr {[llength [glob -nocomplain [file join $w *_trails_mef.fits]]] == 1}]
    # undo twice then redo through the mask manager
    OGFMaskRun undo {} {}; OGFMaskRun undo {} {}
    R mask_undo_chain [expr {[ntrail] == 0}] [maskstats]
    # Auto Mask keeps the trail bit
    OGFTrailsClear
    R clear_trails [expr {[ntrail] == 0}]
    set sp [file join $::dir ogfinder_session.py]
    CatalogPanelSessionSave $sp
    R session_exported [file exists $sp]
    R final_geometry [expr {[geom] eq $g0}] "$g0 -> [geom]"
    puts $::fh "SUMMARY failures=$::nf"; close $::fh
}
after 3000 {
    if {[catch run err]} {puts $::fh "FAIL run_error $err $::errorInfo"; incr ::nf; puts $::fh "SUMMARY failures=$::nf"; close $::fh}
    exit 0
}
