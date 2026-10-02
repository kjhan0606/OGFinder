#  OGFinder plugin "report": per-candidate review state + HTML/PDF candidate report.
#
#  Review state is stored as three CATALOG COLUMNS, so it is saved/loaded with the catalog (TSV/CSV/FITS export), follows the catalog
#  through the frame switch and shows up in the table:
#       REVIEW        accept | reject | uncertain | (empty)
#       REVIEW_NOTE   free text (tabs/newlines are folded to spaces)
#       REVIEW_TIME   YYYY-MM-DD HH:MM:SS of the last decision
#  Every change is recorded in the session log as an INFORMATIONAL step (class "note": the exported script lists it and skips it):
#       review.set   payload  numbers_list / status / note        (accept, reject, uncertain, clear)
#       review.note  payload  numbers_list / note
#       report.export payload file / rows / ...                   (a report was written)
#  The table can be filtered by review state (::ogf::cat::filter_set REVIEW ...) and the review cells are tinted.
#  "Export Report..." runs report.py (this directory) in the job runner and writes one self-contained HTML (optionally a PDF).
#  Keys used: review,show (see the registry in ogf_core.tcl).  All catalog access goes through ::ogf::cat.
#  See docs/reports_review.md.

package provide DS9 1.0

namespace eval ::ogf::report {variable last; array set last {out {} argv {}}}
namespace eval ::ogf::review {
    variable states {accept reject uncertain}
    variable colors {accept #cfeccf reject #f3c9c9 uncertain #fbe7b0}
    variable cols {REVIEW REVIEW_NOTE REVIEW_TIME}
}

proc OGFReviewInit {} {
    ::ogf::cat::set review,show all
    ::ogf::cat::on_table_filled ::ogf::review::colorize
}

# ------------------------------------------------------------------ helpers
proc ::ogf::review::_clean {s} {
    return [string trim [regsub -all {[\t\r\n]+} $s { }]]
}
# is the galaxy table the one on screen?  (the review columns live in the galaxy catalog)
proc ::ogf::review::_ready {what} {
    if {![::ogf::cat::has]} {::ogf::status "$what: extract or load a catalog first"; return 0}
    if {![catch {::ogf::td::active} a] && $a} {
        ::ogf::status "$what: switch the kind filter to Galaxies first (review applies to the galaxy catalog)"
        return 0
    }
    if {[lsearch -exact [::ogf::cat::columns] NUMBER] < 0} {::ogf::status "$what: the catalog has no NUMBER column"; return 0}
    return 1
}
# the review decision of a NUMBER ("" when none)
proc ::ogf::review::get {num} {
    foreach d [::ogf::cat::rows] {
        if {[dict get $d NUMBER] eq $num} {return [expr {[dict exists $d REVIEW] ? [dict get $d REVIEW] : {}}]}
    }
    return {}
}
# counts {accept N reject N uncertain N none N}
proc ::ogf::review::counts {} {
    set c [dict create accept 0 reject 0 uncertain 0 none 0]
    foreach v [::ogf::cat::values REVIEW] {
        if {$v in {accept reject uncertain}} {dict incr c $v} else {dict incr c none}
    }
    if {![llength [::ogf::cat::values REVIEW]]} {dict set c none [::ogf::cat::nrows]}
    return $c
}

# ------------------------------------------------------------------ set / clear / note
# nums: NUMBERs (default: the selection); status: accept|reject|uncertain|clear.  note: if given, replaces the note ("" removes it).
proc ::ogf::review::set_status {status {nums {}} {note {}} {setnote 0}} {
    variable states
    if {$status ni [concat $states clear]} {error "review status must be accept, reject, uncertain or clear (got \"$status\")"}
    if {![_ready "Review"]} {return 0}
    if {![llength $nums]} {set nums [::ogf::cat::selection]}
    if {![llength $nums]} {::ogf::status "Review: select a candidate in the table (or on the image) first"; return 0}
    set now [clock format [clock seconds] -format "%Y-%m-%d %H:%M:%S"]
    set note [_clean $note]
    set changes [dict create]
    # column order REVIEW, REVIEW_NOTE, REVIEW_TIME (a missing column is appended in this order)
    set havenote [expr {[lsearch -exact [::ogf::cat::columns] REVIEW_NOTE] >= 0}]
    foreach n $nums {
        set d [dict create]
        if {$status eq "clear"} {
            dict set d REVIEW {}; dict set d REVIEW_NOTE {}; dict set d REVIEW_TIME {}
        } else {
            dict set d REVIEW $status
            if {$setnote} {dict set d REVIEW_NOTE $note} elseif {!$havenote} {dict set d REVIEW_NOTE {}}
            dict set d REVIEW_TIME $now
        }
        dict set changes $n $d
    }
    set k [::ogf::cat::set_cells $changes]
    if {$k == 0} {::ogf::status "Review: none of [llength $nums] selected NUMBER(s) is in the catalog"; return 0}
    set verb [expr {$status eq "clear" ? "cleared" : $status}]
    OGFSessLog review.set note {} -tool internal \
        -title "Review: $verb [llength $nums] candidate[expr {[llength $nums] == 1 ? {} : {s}}]" \
        -payload [dict create numbers_list $nums status $status note $note] \
        -note "manual review decision stored in the catalog columns REVIEW / REVIEW_NOTE / REVIEW_TIME; informational, not replayed"
    refresh_filter
    ::ogf::status "Review: $verb $k candidate[expr {$k == 1 ? {} : {s}}] ([format_counts])"
    return $k
}
proc ::ogf::review::set_note {note {nums {}}} {
    if {![_ready "Review note"]} {return 0}
    if {![llength $nums]} {set nums [::ogf::cat::selection]}
    if {![llength $nums]} {::ogf::status "Review note: select a candidate first"; return 0}
    set note [_clean $note]
    set changes [dict create]
    foreach n $nums {dict set changes $n [dict create REVIEW_NOTE $note]}
    if {[lsearch -exact [::ogf::cat::columns] REVIEW] < 0} {foreach n $nums {dict set changes $n REVIEW {}; dict set changes $n REVIEW_TIME {}}}
    set k [::ogf::cat::set_cells $changes]
    OGFSessLog review.note note {} -tool internal -title "Review: note on [llength $nums] candidate[expr {[llength $nums] == 1 ? {} : {s}}]" \
        -payload [dict create numbers_list $nums note $note] -note "free-text note stored in the catalog column REVIEW_NOTE; informational"
    refresh_filter
    ::ogf::status "Review: note set on $k candidate[expr {$k == 1 ? {} : {s}}]"
    return $k
}
proc ::ogf::review::format_counts {} {
    set c [counts]
    return "[dict get $c accept] accepted, [dict get $c reject] rejected, [dict get $c uncertain] uncertain, [dict get $c none] unreviewed"
}

# ------------------------------------------------------------------ filter
# what: all | accept | reject | uncertain | none | notrej | accept+uncertain
proc ::ogf::review::show {what} {
    if {![::ogf::cat::has]} {::ogf::status "Review filter: no catalog"; return}
    switch -- $what {
        all {::ogf::cat::filter_clear REVIEW}
        accept - reject - uncertain {::ogf::cat::filter_set REVIEW [list $what]}
        none {::ogf::cat::filter_set REVIEW {{}}}
        notrej {::ogf::cat::filter_set REVIEW {accept uncertain {}}}
        accept+uncertain {::ogf::cat::filter_set REVIEW {accept uncertain}}
        default {error "unknown review filter \"$what\""}
    }
    ::ogf::cat::set review,show $what
    CatalogPanelFilter
}
proc ::ogf::review::refresh_filter {} {
    # set_cells reloads the table through LoadTSV, which re-applies the filters; keep the key honest
    if {![dict exists [::ogf::cat::filters] REVIEW]} {::ogf::cat::set review,show all}
}

# ------------------------------------------------------------------ table tint
proc ::ogf::review::colorize {} {
    variable colors
    if {![::ogf::cat::exists tbl] || [catch {winfo exists [::ogf::cat::get tbl]} ok] || !$ok} return
    set tbl [::ogf::cat::get tbl]
    foreach {s c} $colors {
        catch {$tbl tag delete rv_$s}
        $tbl tag configure rv_$s -background $c -foreground black
    }
    set rc [::ogf::cat::table_col REVIEW]
    if {$rc < 0} return
    set nc [::ogf::cat::table_col NUMBER]
    set cells [dict create accept {} reject {} uncertain {}]
    set nr [expr {[::ogf::cat::table_nrows] + 1}]
    for {set r 1} {$r < $nr} {incr r} {
        if {![::ogf::cat::cell_exists $r $rc]} continue
        set v [::ogf::cat::cell $r $rc]
        if {[dict exists $cells $v]} {
            dict lappend cells $v $r,$rc
            if {$nc > 0} {dict lappend cells $v $r,$nc}
        }
    }
    dict for {s l} $cells {if {[llength $l]} {$tbl tag cell rv_$s {*}$l}}
}

# ------------------------------------------------------------------ menu procs
proc OGFReviewSet {status} {::ogf::review::set_status $status}
proc OGFReviewShow {what} {::ogf::review::show $what}
proc OGFReviewSetShown {status} {
    if {![::ogf::review::_ready "Review"]} return
    set nums [::ogf::cat::shown_numbers]
    if {![llength $nums]} {::ogf::status "Review: the table shows no rows"; return}
    if {[llength $nums] > 1 && [tk_messageBox -type yesno -icon question -title "Review all shown rows" \
            -message "Mark all [llength $nums] rows shown in the table as $status?"] ne "yes"} return
    ::ogf::review::set_status $status $nums
}
proc OGFReviewNoteDialog {} {
    if {![::ogf::review::_ready "Review note"]} return
    set nums [::ogf::cat::selection]
    if {![llength $nums]} {::ogf::status "Review note: select a candidate first"; return}
    set cur {}
    if {[llength $nums] == 1} {
        foreach d [::ogf::cat::rows] {if {[dict get $d NUMBER] eq [lindex $nums 0] && [dict exists $d REVIEW_NOTE]} {set cur [dict get $d REVIEW_NOTE]}}
    }
    set r [OGFForm "Review note" [list [list note "Note ([llength $nums] candidate[expr {[llength $nums] == 1 ? {} : {s}}])" $cur]]]
    if {$r eq {}} return
    ::ogf::review::set_note [dict get $r note]
}

# ------------------------------------------------------------------ provenance + export
# JSON with what the report lists as provenance: session steps, versions, parameters
proc ::ogf::report::provenance_json {} {
    global ogfsess ds9
    set root [OGFSessRoot]
    set git {}; set dirty {}
    catch {set git [string trim [OGFsess_exec git -C $root rev-parse HEAD]]}
    catch {set dirty [expr {[string trim [OGFsess_exec git -C $root status --porcelain -uno -- ds9/library plugins moving ai_bridge]] ne {}}]}
    set parts {}
    lappend parts "  [OGFJStr git_head]: [OGFJStr $git]"
    lappend parts "  [OGFJStr git_dirty_sources]: [OGFJStr $dirty]"
    lappend parts "  [OGFJStr ds9_version]: [OGFJStr [expr {[info exists ds9(version)] ? $ds9(version) : {}}]]"
    lappend parts "  [OGFJStr ds9_exe]: [OGFJStr [info nameofexecutable]]"
    lappend parts "  [OGFJStr python]: [OGFJStr [OGFPython]]"
    lappend parts "  [OGFJStr session_started]: [OGFJStr [expr {[info exists ogfsess(t0)] ? [clock format $ogfsess(t0) -format %Y-%m-%dT%H:%M:%S%z] : {}}]]"
    lappend parts "  [OGFJStr table_filter]: [OGFJStr [::ogf::report::table_filter_text]]"
    # steps
    set js {}
    foreach rec [::ogf::session::steps] {
        set l {}
        foreach k {seq step title class wall ms failed note} {lappend l "[OGFJStr $k]: [OGFJStr [dict get $rec $k]]"}
        lappend l "[OGFJStr argv]: [OGFJList [dict get $rec argv]]"
        lappend l "[OGFJStr payload]: [OGFJDict [dict get $rec payload]]"
        lappend js "    \{[join $l {, }]\}"
    }
    lappend parts "  [OGFJStr steps]: \[\n[join $js ",\n"]\n  \]"
    # plugin parameters
    set pj {}
    foreach id [::ogf::reg::ids] {
        set specs [::ogf::json::get [::ogf::reg::get $id] params]
        if {$specs eq {} || $id eq "report"} continue
        set l {}
        foreach s $specs {
            set n [dict get $s name]
            if {[catch {::ogf::params::get $id $n} v]} continue
            lappend l "[OGFJStr $n]: [OGFJStr $v]"
        }
        if {[llength $l]} {lappend pj "    [OGFJStr $id]: \{[join $l {, }]\}"}
    }
    lappend parts "  [OGFJStr plugin_params]: \{\n[join $pj ",\n"]\n  \}"
    set ex {}
    foreach k [concat [::ogf::cat::keys param,*] [::ogf::cat::keys extract_param,*]] {
        lappend ex "[OGFJStr $k]: [OGFJStr [::ogf::cat::get $k]]"
    }
    lappend parts "  [OGFJStr extraction_params]: \{[join $ex {, }]\}"
    return "\{\n[join $parts ",\n"]\n\}\n"
}
proc ::ogf::report::table_filter_text {} {
    set t {}
    set pat [::ogf::cat::get search_var {}]
    if {$pat ne {}} {lappend t "text search \"$pat\""}
    if {[::ogf::cat::filter_active]} {lappend t [::ogf::cat::filter_text]}
    return [join $t {; }]
}

# the argv of report.py for the stored parameters.  Pure function (tests call it): out = HTML path, files = paths of the inputs
proc ::ogf::report::build_argv {out catfile provfile image {numfile {}}} {
    set P {::ogf::params::get report}
    set argv [list [OGFPython] [file join [::ogf::step::plugin_dir report] report.py] \
        --catalog $catfile --image $image --output $out --provenance $provfile \
        --include [{*}$P include] \
        --title [{*}$P title] --author [{*}$P author] --max-rows [{*}$P max_rows] \
        --cutout-size [{*}$P cutout_size] --thumb-zoom [{*}$P thumb_zoom] --stretch [{*}$P stretch] \
        --table-filter [table_filter_text]]
    if {[{*}$P include] eq "table" && $numfile ne {}} {lappend argv --only-numbers $numfile}
    if {[{*}$P sort_by] ne {}} {lappend argv --sort-by [{*}$P sort_by]}
    if {![{*}$P thumbnails]} {lappend argv --no-thumbs}
    if {![{*}$P marker]} {lappend argv --no-marker}
    if {![{*}$P all_columns]} {lappend argv --no-all-columns}
    if {[{*}$P columns] ne {}} {lappend argv --columns [string map {{ } {}} [{*}$P columns]]}
    if {[{*}$P pdf]} {lappend argv --pdf}
    return $argv
}

# OUT: HTML path (a file dialog when empty).  Returns 1 when the job was started.
proc OGFReportExport {{out {}} {done {}}} {
    if {![::ogf::cat::has]} {::ogf::status "Export Report: extract or load a catalog first"; return 0}
    if {![catch {::ogf::td::active} a] && $a} {::ogf::status "Export Report: switch the kind filter to Galaxies first"; return 0}
    set image [::ogf::cat::image_file]
    if {$out eq {}} {
        set base [expr {$image ne {} ? [file rootname [file tail $image]] : {catalog}}]
        set out [tk_getSaveFile -title "Export Report" -defaultextension .html -initialfile "${base}_report.html" \
            -filetypes {{{HTML report} {.html}} {{All files} *}}]
        if {$out eq {}} return
    }
    set work [OGFSessWorkDir]
    file mkdir $work
    set provf [file join $work report_provenance.json]
    set numf [file join $work report_numbers.txt]
    if {[catch {
        set fd [open $provf w]; fconfigure $fd -encoding utf-8; puts -nonewline $fd [::ogf::report::provenance_json]; close $fd
        set fd [open $numf w]; puts $fd [join [::ogf::cat::shown_numbers] \n]; close $fd
    } err]} {::ogf::status "Export Report: $err"; return 0}
    set catf [::ogf::cat::temp_file report]
    set argv [::ogf::report::build_argv $out $catf $provf $image $numf]
    set ::ogf::report::last(out) $out
    set ::ogf::report::last(argv) $argv
    return [::ogf::job::run $argv -title "Export Report" -plugin report \
        -done [list ::ogf::report::finished $out $done]]
}
proc ::ogf::report::finished {out done ok output ms} {
    if {!$ok} {
        ::ogf::status "Export Report failed: [string range [lindex [split [string trim $output] \n] end] 0 140]"
        if {$done ne {}} {catch {uplevel #0 [list {*}$done 0 $output]}}
        return
    }
    set f [split [string trim $output] \t]
    OGFSessLog report.export note {} -tool internal -title "Export report [file tail $out]" \
        -payload [dict create file [file tail $out] rows [lindex $f 2] thumbnails [lindex $f 4]] \
        -note "report written (HTML, self-contained); informational, nothing to replay"
    ::ogf::status "Report written: [file tail $out] ([lindex $f 2] candidates, [lindex $f 4] cut-outs[expr {[lindex $f 5] ne {-} && [lindex $f 5] ne {} ? ", PDF" : {}}])"
    if {$done ne {}} {catch {uplevel #0 [list {*}$done 1 $output]}}
}
