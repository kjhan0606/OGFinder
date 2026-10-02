# GUI test of the review feature on moving / transient / detection rows (plugins/report + ::ogf::td::review_*).
# Real ds9, real extraction of m51; the moving / transient rows are synthetic (placed on real catalog positions).
#   rm -rf ~/ds9.auto ~/ds9.auto.dir; DISPLAY=:77 OGF_RVTD_OUT=/tmp/rvtd.txt OGF_RVTD_DIR=/tmp/rvtddir bin/ds9 /workspace/fits/m51.fits -geometry 1300x950 -source scripts/verify_review_td.tcl
set ::fh [open $::env(OGF_RVTD_OUT) w]; set ::nf 0
set ::dir $::env(OGF_RVTD_DIR); file mkdir $::dir
proc R {tag ok {d {}}} {puts $::fh "[expr {$ok?{PASS}:{FAIL}}] $tag $d"; flush $::fh; if {!$ok} {incr ::nf}}
proc bgerror {m} {puts $::fh "BGERROR $m"; flush $::fh}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}
proc geom {} {
    update idletasks; update
    set tf [winfo parent $::catpanel(tbl)]
    return "[winfo rooty $tf] [winfo height $tf] [winfo height $::catpanel(infoarea)] [winfo width .]x[winfo height .]"
}
proc col {name} {return [::ogf::cat::table_col $name]}
# values of a table column, in table order
proc colvals {name} {
    set c [col $name]; if {$c < 0} {return {}}
    set out {}; for {set r 1} {$r <= [::ogf::cat::table_nrows]} {incr r} {lappend out [::ogf::cat::cell $r $c]}
    return $out
}
proc steps_since {n} {lmap r [lrange [::ogf::session::steps] $n end] {list [dict get $r step] [dict get $r class]}}
rename tk_messageBox ::rv_real_mb
# (the rerun test places one synthetic moving row far off the image: ds9 then warns "Bad Coordinate mapping"; the stub records it)
proc tk_messageBox {args} {puts $::fh "MESSAGEBOX $args"; flush $::fh; return yes}
proc rv_cb {ok out} {set ::rv_done [list $ok $out]}
proc run {} {
    global catpanel
    R geometry_before [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    CatalogPanelExtract
    set t0 [clock milliseconds]; while {![::ogf::cat::has] && [clock milliseconds]-$t0 < 60000} {update; after 100}
    wait_idle 300
    set ::rows [::ogf::cat::rows]; set rows $::rows
    R extracted [expr {[llength $rows] > 50}] [llength $rows]
    # synthetic rows on real galaxy positions (ra dec of rows 5..9), so the host / overlap decorators and the pixel grid work
    proc pos {i} {set r [lindex $::rows $i]; return [list [dict get $r ALPHA_J2000] [dict get $r DELTA_J2000]]}
    set mrows {}; set k 0
    foreach i {5 6 7} {
        lassign [pos $i] ra de; incr k
        lappend mrows [list $k $ra $de 22.$k 3 12.5 90 0.1 0.8 ok]
    }
    ::ogf::td::set_rows moving {id ra dec mag n rate_arcsec_h pa_deg rms_arcsec score status} $mrows
    lassign [pos 8] ra de
    ::ogf::td::set_rows transient {id ra dec mag n_det files snr_max host_id host_sep host_z offset_re heuristic} [list [list 1 $ra $de 24.5 2 x 9.0 {} {} {} {} {}]]
    lassign [pos 9] ra de
    ::ogf::td::set_rows detection {id ra dec mag ex chip cls snr sign flux_e_s elong a_pix} [list [list 1 $ra $de 23 1 1 star 9 1 5 1.1 2]]
    set crc0 [::ogf::cat::tsv]
    set ng0 [::ogf::cat::nrows]
    # ---- Moving view, nothing reviewed yet
    ::ogf::td::show moving
    wait_idle 200
    R moving_rows [expr {[::ogf::cat::table_nrows] == 3}] [::ogf::cat::table_nrows]
    R no_review_columns_yet [expr {[col REVIEW] < 0}]
    R ready_in_td_view [::ogf::review::_ready x]
    ::ogf::cat::select M1
    wait_idle 200
    R select_td_row [expr {[::ogf::cat::selection] eq {M1}}] [::ogf::cat::selection]
    set s0 [llength [::ogf::session::steps]]
    OGFReviewSet accept
    R accept_stored [expr {[lindex [::ogf::td::review_get M1] 0] eq {accept}}] [::ogf::td::review_get M1]
    set cols [lmap c [lrange [colvals NUMBER] 0 -1] {set c}]
    set hdr {}; for {set c 1} {$c <= [::ogf::cat::table_ncols]} {incr c} {lappend hdr [::ogf::cat::header $c]}
    R columns_appended [expr {[lrange $hdr end-2 end] eq {REVIEW REVIEW_NOTE REVIEW_TIME}}] [lrange $hdr end-2 end]
    R accept_in_table [expr {[lindex [colvals REVIEW] 0] eq {accept} && [lindex [colvals REVIEW] 1] eq {}}] [colvals REVIEW]
    R time_format [regexp {^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d$} [lindex [colvals REVIEW_TIME] 0]] [colvals REVIEW_TIME]
    R selection_kept [expr {[::ogf::cat::selection] eq {M1}}] [::ogf::cat::selection]
    ::ogf::review::set_status reject {M2}
    ::ogf::review::set_status uncertain {M3}
    ::ogf::review::set_note "a note\twith tab" {M3}
    set c [::ogf::review::counts]
    R counts_td [expr {[dict get $c accept] == 1 && [dict get $c reject] == 1 && [dict get $c uncertain] == 1 && [dict get $c none] == 0}] $c
    R note_clean [expr {[lindex [colvals REVIEW_NOTE] 2] eq "a note with tab"}] [colvals REVIEW_NOTE]
    R still_moving_view [expr {[::ogf::td::kind] eq {moving}}] [::ogf::td::kind]
    # ---- session steps are informational, same step names as for galaxies
    set st [steps_since $s0]
    R steps_recorded [expr {$st eq {{review.set note} {review.set note} {review.set note} {review.note note}}}] $st
    # ---- the galaxy catalog is untouched
    R galaxy_catalog_untouched [expr {[::ogf::cat::tsv] eq $crc0 && [::ogf::cat::nrows] == $ng0 && "REVIEW" ni [::ogf::cat::columns]}]
    # ---- filters work on the time-domain rows
    ::ogf::review::show accept
    R filter_accept [expr {[colvals NUMBER] eq {M1}}] [colvals NUMBER]
    ::ogf::review::show reject
    R filter_reject [expr {[colvals NUMBER] eq {M2}}] [colvals NUMBER]
    ::ogf::review::show notrej
    R filter_notrej [expr {[colvals NUMBER] eq {M1 M3}}] [colvals NUMBER]
    set catpanel(search_var) M3
    CatalogPanelFilter
    R filter_and_search [expr {[colvals NUMBER] eq {M3}}] [colvals NUMBER]
    set catpanel(search_var) {}
    ::ogf::review::show all
    R filter_all [expr {[colvals NUMBER] eq {M1 M2 M3}}] [colvals NUMBER]
    # ---- tint
    set tags [$catpanel(tbl) tag names]
    R tint_tags [expr {"rv_accept" in $tags && "rv_reject" in $tags}] $tags
    set rc [col REVIEW]
    R tint_cell [expr {[$catpanel(tbl) tag includes rv_accept 1,$rc] && [$catpanel(tbl) tag includes rv_reject 2,$rc] && [$catpanel(tbl) tag includes rv_uncertain 3,$rc]}]
    # ---- review all shown (after a search for two rows)
    rename tk_messageBox ::rv_orig_mb; proc tk_messageBox {args} {return yes}
    set catpanel(search_var) {M}
    CatalogPanelFilter
    OGFReviewSetShown uncertain
    rename tk_messageBox {}; rename ::rv_orig_mb tk_messageBox
    set catpanel(search_var) {}
    ::ogf::review::show all
    R review_all_shown [expr {[colvals REVIEW] eq {uncertain uncertain uncertain}}] [colvals REVIEW]
    # ---- sort keeps the decisions with their rows
    ::ogf::review::set_status accept M2
    ::ogf::td::sort rate_arcsec_h descending
    set nu [colvals NUMBER]; set rv [colvals REVIEW]
    R sort_keeps_decisions [expr {[lindex $rv [lsearch -exact $nu M2]] eq {accept} && [lindex $rv [lsearch -exact $nu M1]] eq {uncertain}}] "$nu / $rv"
    ::ogf::td::sort id ascending
    # ---- other kinds and the All view
    ::ogf::td::show transient
    wait_idle 200
    R transient_no_review_yet [expr {[colvals REVIEW] eq {{}}}] [colvals REVIEW]
    ::ogf::cat::select T1; wait_idle 100
    OGFReviewSet reject
    R transient_reject [expr {[colvals REVIEW] eq {reject}}] [colvals REVIEW]
    ::ogf::td::show detection
    ::ogf::cat::select D1; wait_idle 100
    OGFReviewSet accept
    R detection_accept [expr {[colvals REVIEW] eq {accept}}] [colvals REVIEW]
    ::ogf::td::show all
    wait_idle 200
    set nu [colvals NUMBER]; set rv [colvals REVIEW]
    R all_view_has_review [expr {[col REVIEW] > 0 && [lindex $rv [lsearch -exact $nu M2]] eq {accept} && [lindex $rv [lsearch -exact $nu T1]] eq {reject}}] "$nu / $rv"
    R all_view_galaxies_unreviewed [expr {[llength [lsearch -all -exact $rv {}]] >= $ng0}] [llength $nu]
    # a galaxy decision made from the All view goes to the catalog and the kind filter stays on All
    set g1 [lindex [::ogf::cat::values NUMBER] 0]
    ::ogf::review::set_status accept $g1
    R galaxy_from_all_view [expr {[::ogf::td::kind] eq {all} && [dict get [lindex [::ogf::cat::rows] 0] REVIEW] eq {accept}}] "kind=[::ogf::td::kind]"
    set nu [colvals NUMBER]; set rv [colvals REVIEW]
    R all_view_shows_galaxy_decision [expr {[lindex $rv [lsearch -exact $nu $g1]] eq {accept}}]
    # mixed selection: one galaxy and one moving row in one call
    ::ogf::review::set_status reject [list $g1 M1]
    R mixed_galaxy [expr {[dict get [lindex [::ogf::cat::rows] 0] REVIEW] eq {reject}}]
    R mixed_td [expr {[lindex [::ogf::td::review_get M1] 0] eq {reject}}]
    # ---- clear
    ::ogf::review::set_status clear {M1 M2 M3 T1 D1 }
    R clear_td [expr {[::ogf::td::review_count] == 0}] [::ogf::td::review_count]
    # ---- save / load of the time-domain decisions
    ::ogf::review::set_status accept {M1}; ::ogf::review::set_status reject {M2}; ::ogf::review::set_status uncertain {T1}
    ::ogf::review::set_note "keep this" {M1}
    set f [file join $::dir td_review.tsv]
    R td_save [expr {[OGFReviewTDSave $f] == 3 && [file exists $f]}]
    set fd [open $f r]; set txt [read $fd]; close $fd
    R td_save_header [string match "key\tkind\tra\tdec\tREVIEW\tREVIEW_NOTE\tREVIEW_TIME*" $txt]
    ::ogf::review::set_status clear {M1 M2 T1}
    R td_wiped [expr {[::ogf::td::review_count] == 0}]
    R td_load [expr {[OGFReviewTDLoad $f] == 3}]
    R td_restored [expr {[::ogf::td::review_get M1] ne {} && [lindex [::ogf::td::review_get M1] 0] eq {accept} && [lindex [::ogf::td::review_get M1] 1] eq {keep this} && [lindex [::ogf::td::review_get M2] 0] eq {reject}}] [::ogf::td::review_get M1]
    # a decision whose row has a different position is not applied
    set bad [file join $::dir bad.tsv]
    set fd [open $bad w]; puts $fd "key\tkind\tra\tdec\tREVIEW\tREVIEW_NOTE\tREVIEW_TIME\nM1\tmoving\t1.0\t2.0\taccept\t\t\nM9\tmoving\t1.0\t2.0\taccept\t\t"; close $fd
    ::ogf::review::set_status clear {M1 M2 T1}
    lassign [::ogf::td::review_import $bad] ok skip
    R td_import_position_checked [expr {$ok == 0 && $skip == 2}] "ok=$ok skip=$skip"
    # ---- a new run replaces the rows: decisions whose row changed are dropped, unchanged ones stay
    ::ogf::review::set_status accept {M1}; ::ogf::review::set_status reject {M2}
    lassign [pos 5] ra5 de5
    ::ogf::td::set_rows moving {id ra dec mag n rate_arcsec_h pa_deg rms_arcsec score status} [list [list 1 $ra5 $de5 22.1 3 12.5 90 0.1 0.8 ok] [list 2 1.5 2.5 22.2 3 12.5 90 0.1 0.8 ok]]
    R rerun_keeps_same_row [expr {[lindex [::ogf::td::review_get M1] 0] eq {accept}}] [::ogf::td::review_get M1]
    R rerun_drops_changed_row [expr {[::ogf::td::review_get M2] eq {}}] [::ogf::td::review_get M2]
    ::ogf::td::clear moving
    R clear_kind_drops_decisions [expr {[::ogf::td::review_get M1] eq {}}]
    # ---- report export from a time-domain view
    ::ogf::td::set_rows moving {id ra dec mag n rate_arcsec_h pa_deg rms_arcsec score status} $mrows
    ::ogf::td::show moving
    ::ogf::review::set_status accept {M1 M2}
    ::ogf::review::set_status reject {M3}
    ::ogf::review::show accept
    set out [file join $::dir report_td.html]; file delete $out
    set ::rv_done {}
    set s1 [llength [::ogf::session::steps]]
    set started [OGFReportExport $out rv_cb]
    R export_started $started
    set tw [clock milliseconds]
    while {$::rv_done eq {} && [clock milliseconds]-$tw < 60000} {update; after 50}
    if {$::rv_done eq {}} {set ::rv_done [list 0 timeout]}
    R export_ok [expr {[lindex $::rv_done 0] == 1}] [string range $::rv_done 0 160]
    R export_file [expr {[file exists $out] && [file size $out] > 3000}] [expr {[file exists $out] ? [file size $out] : 0}]
    R export_stdout [regexp {^REPORT\t\S+\t2\t2} [lindex $::rv_done 1]] [lindex $::rv_done 1]
    if {[file exists $out]} {
        set fd [open $out r]; set h [read $fd]; close $fd
        R export_html_rows [expr {[llength [regexp -all -inline {class="cand st-accept"} $h]] == 2 && ![string match {*class="cand st-reject"*} $h]}]
        set ids {}
        foreach m [lrange [split [string map [list {<tr class="cand} "\x01"] $h] "\x01"] 1 end] {if {[regexp {<td>(M\d)</td>} $m -> id]} {lappend ids $id}}
        R export_html_ids [expr {[lsort $ids] eq {M1 M2}}] $ids
    }
    R export_step [expr {[lrange [steps_since $s1] end end] eq {{report.export note}}}] [steps_since $s1]
    ::ogf::review::show all
    # ---- back to the galaxy view: galaxy decisions made in the All view survived
    ::ogf::td::show galaxies
    wait_idle 200
    R galaxies_view_back [expr {![::ogf::td::active]}]
    R galaxy_decision_kept [expr {[dict get [lindex [::ogf::cat::rows] 0] REVIEW] eq {reject}}]
    # ---- session export lists the informational steps
    set py [file join $::dir session.py]
    CatalogPanelSessionSave $py
    set fd [open $py r]; set pyt [read $fd]; close $fd
    R session_lists_review [expr {[string match "*review.set*" $pyt] && [string match "*review.td_save*" $pyt]}]
    R geometry_after [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    puts $::fh "SUMMARY failures=$::nf"; close $::fh; exit
}
after 3000 {if {[catch run err]} {puts $::fh "ERROR $err $::errorInfo"; close $::fh; exit 1}}
