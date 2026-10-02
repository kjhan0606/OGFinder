# Layout-invariant check (run_all_checks: geometry).  At -geometry 1300x950 the catalog table must stay at y=181, h=769 with an
# info area of 154 px in EVERY state of the panel, and detach / reattach must give main-window widths 1300 -> 736 -> 1300.
#   rm -rf ~/ds9.auto ~/ds9.auto.dir; DISPLAY=:77 OGF_GEO_OUT=/tmp/geo.txt bin/ds9 /workspace/fits/m51.fits -geometry 1300x950 -source scripts/verify_geometry.tcl
# The numbers are the ones measured on this box (X11 fonts, no window manager); other platforms/fonts will differ (docs/windows_macos_build.md).
set ::fh [open $::env(OGF_GEO_OUT) w]; set ::nf 0; set ::nstate 0
proc R {tag ok {d {}}} {puts $::fh "[expr {$ok?{PASS}:{FAIL}}] $tag $d"; flush $::fh; if {!$ok} {incr ::nf}}
proc bgerror {m} {puts $::fh "BGERROR $m"; flush $::fh}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}
proc geom {} {
    update idletasks; update
    set tf [winfo parent $::catpanel(tbl)]
    return "[winfo rooty $tf] [winfo height $tf] [winfo height $::catpanel(infoarea)]"
}
proc main_w {} {update idletasks; update; return [winfo width .]}
proc check {tag} {
    incr ::nstate
    set g [geom]
    R "geo_$tag" [expr {$g eq "181 769 154"}] "table_y/table_h/info_h = $g (want 181 769 154)"
}
proc run {} {
    global catpanel ds9 current ogfui
    wait_idle 500
    check start
    R main_width_start [expr {[main_w] == 1300}] [main_w]
    # the six tabs
    foreach t $::ogf::tabs {OGFUIShowTab $t; wait_idle 200; check tab_[OGFUITabId $t]}
    OGFUIShowTab Detect
    # a catalog, a text search and a column filter
    CatalogPanelExtract
    set t0 [clock milliseconds]; while {![::ogf::cat::has] && [clock milliseconds]-$t0 < 90000} {update; after 100}
    wait_idle 400
    check after_extract
    set catpanel(search_var) 1; CatalogPanelFilter; wait_idle 200; check search
    set catpanel(search_var) {}; CatalogPanelFilter
    # the time-domain views
    set rows [::ogf::cat::rows]
    set mrows {}; set k 0
    foreach i {3 4 5} {set r [lindex $rows $i]; incr k; lappend mrows [list $k [dict get $r ALPHA_J2000] [dict get $r DELTA_J2000] 22 3 1 1 1 1 ok]}
    ::ogf::td::set_rows moving {id ra dec mag n rate_arcsec_h pa_deg rms_arcsec score status} $mrows
    set r [lindex $rows 7]
    ::ogf::td::set_rows transient {id ra dec mag n_det files snr_max host_id host_sep host_z offset_re heuristic} [list [list 1 [dict get $r ALPHA_J2000] [dict get $r DELTA_J2000] 24 2 x 9 {} {} {} {} {}]]
    foreach k {moving transient all galaxies} {::ogf::td::show $k; wait_idle 300; check td_$k}
    # review columns on a td view and the filter
    ::ogf::td::show moving; ::ogf::cat::select M1; wait_idle 200; OGFReviewSet accept; wait_idle 200; check td_review
    ::ogf::review::show accept; wait_idle 200; check td_review_filter
    ::ogf::review::show all; ::ogf::td::show galaxies; wait_idle 300
    # a selection and a multi-selection (summary text changes the info area content, not its height)
    set nums [::ogf::cat::values NUMBER]
    ::ogf::cat::select [lindex $nums 0]; wait_idle 200; check select_one
    ::ogf::cat::select [lrange $nums 1 4] add; wait_idle 200; check select_multi
    CatalogPanelClearSelection
    # tile / single
    CreateFrame; LoadFitsFile /workspace/fits/m51.fits {} {}; wait_idle 800
    OGFUIDisplay tile; wait_idle 1000
    R tile_on [OGFTileIsOn] "frames=$ds9(frames)"
    check tile
    OGFUIDisplay single; wait_idle 800; check single_again
    # detach / reattach (the Tools menu entry's command)
    set w0 [main_w]
    set catpanel(detached) 1; CatalogPanelToggleDetach; wait_idle 600
    set w1 [main_w]
    check detached
    R width_detached [expr {$w0 == 1300 && $w1 == 736}] "main width $w0 -> $w1 (want 1300 -> 736)"
    set catpanel(detached) 0; CatalogPanelToggleDetach; wait_idle 600
    set w2 [main_w]
    check reattached
    R width_reattached [expr {$w2 == 1300}] "main width after reattach $w2 (want 1300)"
    R states_checked [expr {$::nstate >= 20}] "$::nstate states"
    puts $::fh "SUMMARY failures=$::nf"; close $::fh; exit
}
after 3000 {if {[catch run err]} {puts $::fh "ERROR $err $::errorInfo"; close $::fh; exit 1}}
