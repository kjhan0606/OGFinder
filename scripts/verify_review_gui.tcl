# GUI test of the review feature and the report export (plugins/report).  Real ds9, real extraction of m51, real report.py job.
#   rm -rf ~/ds9.auto ~/ds9.auto.dir; DISPLAY=:77 OGF_REVIEW_OUT=/tmp/rv.txt OGF_REVIEW_DIR=/tmp/rvdir bin/ds9 /workspace/fits/m51.fits -geometry 1300x950 -source scripts/verify_review_gui.tcl
# Output lines: PASS/FAIL tag detail, SUMMARY failures=N.  HTML contents are checked by plugins/report/tests (Python); here only the
# GUI side: columns, table filter, tint, session steps, persistence in save/load, layout invariants, the export job.
set ::fh [open $::env(OGF_REVIEW_OUT) w]; set ::nf 0
set ::dir $::env(OGF_REVIEW_DIR); file mkdir $::dir
proc R {tag ok {d {}}} {puts $::fh "[expr {$ok?{PASS}:{FAIL}}] $tag $d"; flush $::fh; if {!$ok} {incr ::nf}}
proc bgerror {m} {puts $::fh "BGERROR $m"; flush $::fh}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}
proc geom {} {
    update idletasks; update
    set tf [winfo parent $::catpanel(tbl)]
    return "[winfo rooty $tf] [winfo height $tf] [winfo height $::catpanel(infoarea)] [winfo width .]x[winfo height .]"
}
proc table_numbers {} {
    global catpanel; set db $catpanel(tbldb); global $db
    set nc [OGFTableCol NUMBER]; set out {}
    for {set r 1} {$r < [$catpanel(tbl) cget -rows]} {incr r} {if {[info exists ${db}($r,$nc)]} {lappend out [set ${db}($r,$nc)]}}
    return $out
}
proc rv_cb {ok out} {set ::rv_done [list $ok $out]}
proc steps_since {n} {lmap r [lrange [::ogf::session::steps] $n end] {list [dict get $r step] [dict get $r class]}}
proc run {} {
    global catpanel
    R geometry_before [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    CatalogPanelExtract
    set t0 [clock milliseconds]; while {![::ogf::cat::has] && [clock milliseconds]-$t0 < 60000} {update; after 100}
    wait_idle 300
    set n [::ogf::cat::nrows]
    R extracted [expr {$n > 50}] "rows=$n"
    R plugin_enabled [expr {"report" in [::ogf::reg::ids]}]
    R plugin_nothing_missing [expr {[::ogf::reg::missing report 1] eq {}}] [::ogf::reg::missing report 1]
    set nums [::ogf::cat::values NUMBER]
    set s0 [llength [::ogf::session::steps]]
    set c0 [::ogf::cat::columns]
    R no_review_columns_yet [expr {"REVIEW" ni $c0}]
    # ---- set by selection through the real select path
    ::ogf::cat::select [lindex $nums 0]
    wait_idle 200
    R select_one [expr {[::ogf::cat::selection] eq [list [lindex $nums 0]]}]
    OGFReviewSet accept
    set cols [::ogf::cat::columns]
    R columns_appended [expr {[lrange $cols end-2 end] eq {REVIEW REVIEW_NOTE REVIEW_TIME}}] [lrange $cols end-2 end]
    R columns_original_untouched [expr {[lrange $cols 0 [expr {[llength $c0]-1}]] eq $c0}]
    set row [lindex [::ogf::cat::rows] 0]
    R accept_value [expr {[dict get $row REVIEW] eq "accept"}]
    R time_format [regexp {^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d$} [dict get $row REVIEW_TIME]] [dict get $row REVIEW_TIME]
    R selection_kept [expr {[::ogf::cat::selection] eq [list [lindex $nums 0]]}] [::ogf::cat::selection]
    ::ogf::review::set_status reject [list [lindex $nums 1] [lindex $nums 2]]
    ::ogf::review::set_status uncertain [lindex $nums 3]
    ::ogf::review::set_note "a note\twith tab and \"quotes\"" [lindex $nums 3]
    set c [::ogf::review::counts]
    R counts [expr {[dict get $c accept] == 1 && [dict get $c reject] == 2 && [dict get $c uncertain] == 1 && [dict get $c none] == $n-4}] $c
    set r3 [lindex [::ogf::cat::rows] 3]
    R note_stored_clean [expr {[dict get $r3 REVIEW_NOTE] eq "a note with tab and \"quotes\""}] [dict get $r3 REVIEW_NOTE]
    R rows_unchanged [expr {[::ogf::cat::nrows] == $n}]
    R bad_status_errors [catch {::ogf::review::set_status maybe $nums}]
    # ---- session steps (informational)
    set st [steps_since $s0]
    R steps_recorded [expr {$st eq {{review.set note} {review.set note} {review.set note} {review.note note}}}] $st
    set last [lindex [::ogf::session::steps] end]
    R step_payload [expr {[dict get $last payload] eq {numbers_list 4 note {a note with tab and "quotes"}}}] [dict get $last payload]
    # ---- the catalog fingerprint ignores review columns (replay check)
    set plain [OGFSessPlainTSV $catpanel(alldata)]
    R plain_tsv_drops_review [expr {[llength [split [lindex [split $plain \n] 0] \t]] == [llength $c0]}]
    R plain_tsv_is_original [expr {[OGFSessPlainTSV [join [lrange [split $catpanel(alldata) \n] 0 3] \n]] ne {}}]
    # ---- filter
    ::ogf::review::show accept
    wait_idle 200
    R filter_accept [expr {[table_numbers] eq [list [lindex $nums 0]]}] [table_numbers]
    R filter_status [string match "Filtered: 1 sources*REVIEW*" $catpanel(status)] $catpanel(status)
    R shown_numbers_api [expr {[::ogf::cat::shown_numbers] eq [list [lindex $nums 0]]}]
    ::ogf::review::show reject
    R filter_reject [expr {[lsort [table_numbers]] eq [lsort [lrange $nums 1 2]]}]
    ::ogf::review::show notrej
    R filter_not_rejected [expr {[llength [table_numbers]] == $n-2}] [llength [table_numbers]]
    ::ogf::review::show none
    R filter_none [expr {[llength [table_numbers]] == $n-4}] [llength [table_numbers]]
    # text search and review filter are ANDed
    ::ogf::review::show reject
    set xs [dict get [lindex [::ogf::cat::rows] 1] X_IMAGE]
    set catpanel(search_var) $xs
    CatalogPanelFilter
    R filter_and_search [expr {[table_numbers] eq [list [lindex $nums 1]]}] "search $xs -> [table_numbers]"
    set catpanel(search_var) {}
    # filter survives a sort and a reload
    CatalogPanelSort MAG_AUTO ascending
    R filter_survives_sort [expr {[lsort [table_numbers]] eq [lsort [lrange $nums 1 2]]}] [table_numbers]
    # tint
    set tbl $catpanel(tbl)
    set cells [$tbl tag cell rv_reject]
    R tint_reject_cells [expr {[llength $cells] == 4}] "[llength $cells] cells"
    ::ogf::review::show all
    R filter_all [expr {[llength [table_numbers]] == $n}]
    R tint_after_all [expr {[llength [$tbl tag cell rv_accept]] == 2 && [llength [$tbl tag cell rv_uncertain]] == 2}] \
        "acc=[llength [$tbl tag cell rv_accept]] unc=[llength [$tbl tag cell rv_uncertain]]"
    # review all shown: filter to a search result and mark
    ::ogf::review::show none
    rename tk_messageBox ::rv_orig_mb; proc tk_messageBox {args} {return yes}
    set catpanel(search_var) {}
    OGFReviewSetShown uncertain
    rename tk_messageBox {}; rename ::rv_orig_mb tk_messageBox
    ::ogf::review::show all
    set c [::ogf::review::counts]
    R review_all_shown [expr {[dict get $c uncertain] == $n-3 && [dict get $c none] == 0}] $c
    # clear
    ::ogf::review::set_status clear [lrange $nums 4 end]
    set c [::ogf::review::counts]
    R clear [expr {[dict get $c none] == $n-4}] $c
    # ---- persistence: save TSV / CSV, wipe, load
    set tsv [file join $::dir saved.tsv]
    proc tk_getSaveFile {args} {return $::tsvf}
    proc tk_getOpenFile {args} {return $::tsvf}
    set ::tsvf $tsv; CatalogPanelSaveCatalog
    set fd [open $tsv r]; set d [read $fd]; close $fd
    R saved_has_columns [string match "*\tREVIEW\tREVIEW_NOTE\tREVIEW_TIME" [lindex [split $d \n] 0]]
    CatalogPanelClear
    R cleared [expr {![::ogf::cat::has]}]
    CatalogPanelLoadCatalog
    R loaded_has_review [expr {[dict get [lindex [::ogf::cat::rows] 0] REVIEW] eq "accept" && [dict get [lindex [::ogf::cat::rows] 3] REVIEW] eq "uncertain"}]
    R loaded_note [expr {[dict get [lindex [::ogf::cat::rows] 3] REVIEW_NOTE] eq "a note with tab and \"quotes\""}]
    set csv [file join $::dir saved.csv]; set ::tsvf $csv; CatalogPanelSaveCatalog
    CatalogPanelClear; CatalogPanelLoadCatalog
    R csv_roundtrip [expr {[dict get [lindex [::ogf::cat::rows] 1] REVIEW] eq "reject"}]
    rename tk_getSaveFile {}; rename tk_getOpenFile {}
    # a new extraction (no REVIEW column) drops the review filter
    ::ogf::review::show accept
    ::ogf::cat::load_tsv "NUMBER\tX_IMAGE\tY_IMAGE\n1\t1\t1\n2\t2\t2" "other"
    R filter_dropped_on_new_catalog [expr {![::ogf::cat::filter_active]}] [::ogf::cat::filters]
    R new_catalog_all_rows [expr {[llength [table_numbers]] == 2}]
    # ---- export job
    CatalogPanelExtract
    set t0 [clock milliseconds]; while {[::ogf::cat::nrows] < 50 && [clock milliseconds]-$t0 < 60000} {update; after 100}
    wait_idle 300
    set nums [::ogf::cat::values NUMBER]
    ::ogf::review::set_status accept [lrange $nums 0 2]
    ::ogf::review::set_status reject [lrange $nums 3 4]
    ::ogf::review::show accept
    set out [file join $::dir report.html]; file delete $out
    set ::rv_done {}
    set s1 [llength [::ogf::session::steps]]
    set started [OGFReportExport $out rv_cb]
    R export_started $started
    after 60000 {if {$::rv_done eq {}} {set ::rv_done [list 0 timeout]}}
    vwait ::rv_done
    R export_ok [expr {[lindex $::rv_done 0] == 1}] [string range $::rv_done 0 120]
    R export_file [expr {[file exists $out] && [file size $out] > 5000}] [expr {[file exists $out] ? [file size $out] : 0}]
    R export_stdout [regexp {^REPORT\t\S+\t3\t3\t3} [lindex $::rv_done 1]] [lindex $::rv_done 1]
    R export_argv_table [expr {"--only-numbers" in $::ogf::report::last(argv) && "--include" in $::ogf::report::last(argv)}]
    R export_step [expr {[lrange [steps_since $s1] end end] eq {{report.export note}}}] [steps_since $s1]
    R export_status [string match "Report written: report.html (3 candidates, 3 cut-outs*" $catpanel(status)] $catpanel(status)
    set pj [file join [OGFSessWorkDir] report_provenance.json]
    R provenance_file [file exists $pj]
    if {[file exists $pj]} {
        set fd [open $pj r]; set pjt [read $fd]; close $fd
        R provenance_has_steps [expr {[string match "*review.set*" $pjt] && [string match "*extract*" $pjt]}]
        R provenance_filter [string match "*REVIEW in {accept}*" $pjt]
    }
    # ---- session export still works and lists the informational steps
    set py [file join $::dir session.py]
    CatalogPanelSessionSave $py
    set fd [open $py r]; set pyt [read $fd]; close $fd
    R session_lists_review [expr {[string match "*review.set*" $pyt] && [string match "*\"class\": \"note\"*" $pyt]}]
    # ---- layout invariants
    R geometry_after [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    puts $::fh "SUMMARY failures=$::nf"; close $::fh; exit
}
after 3000 {if {[catch run err]} {puts $::fh "ERROR $err $::errorInfo"; close $::fh; exit 1}}
