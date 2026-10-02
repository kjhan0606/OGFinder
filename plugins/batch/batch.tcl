# Batch processing (plugins/batch): runs plugins/batch/batch.py in the background, shows progress, opens the summary.

namespace eval ::ogf::batch {
    variable fd {} pid {} running 0 log {} outdir {} rc {}
}

proc OGFBatchScript {} {return [file join [::ogf::step::plugin_dir batch] batch.py]}

proc OGFBatchArgs {} {
    set a [list --fields [::ogf::params::get batch fields-file] --recipe [::ogf::params::get batch recipe-file] --outdir [::ogf::params::get batch outdir] \
	--jobs [::ogf::params::get batch jobs] --retries [::ogf::params::get batch retries] --timeout [::ogf::params::get batch step-timeout] \
	--root [OGFSessRoot] --python [OGFPython]]
    if {[::ogf::params::get batch resume]} {lappend a --resume}
    if {[::ogf::params::get batch only] ne {}} {lappend a --only [::ogf::params::get batch only]}
    return $a
}

proc OGFBatchExample {} {
    set d [::ogf::params::get batch outdir]
    if {$d eq {}} {::ogf::status "Batch: set the output directory first"; return}
    file mkdir $d
    set f [file join $d recipe_example.json]
    exec [OGFPython] [OGFBatchScript] --example-recipe $f
    ::ogf::params::put batch recipe-file $f
    ::ogf::status "Batch: example recipe written to $f"
    return $f
}

proc OGFBatchRun {} {
    if {$::ogf::batch::running} {::ogf::status "Batch: already running"; return}
    foreach {p what} {fields-file "fields list" recipe-file "recipe" outdir "output directory"} {
	if {[::ogf::params::get batch $p] eq {}} {::ogf::status "Batch: set the $what first"; return}
    }
    set ::ogf::batch::outdir [::ogf::params::get batch outdir]
    file mkdir $::ogf::batch::outdir
    set ::ogf::batch::log {}
    set cmd [list [OGFPython] [OGFBatchScript] {*}[OGFBatchArgs]]
    set fd [open [concat | $cmd [list 2>@1]] r]
    set ::ogf::batch::fd $fd
    set ::ogf::batch::running 1
    set ::ogf::batch::pid [pid $fd]
    fconfigure $fd -blocking 0 -buffering line
    fileevent $fd readable ::ogf::batch::readable
    OGFBatchProgress
    ::ogf::status "Batch: started (pid $::ogf::batch::pid)"
    after 500 ::ogf::batch::poll
}

proc ::ogf::batch::readable {} {
    variable fd; variable running; variable log
    if {[catch {gets $fd line} n]} {set n -1}
    if {$n >= 0} {
	append log $line\n
	if {[winfo exists .ogfbatch.log]} {.ogfbatch.log insert end $line\n; .ogfbatch.log see end}
    }
    if {[eof $fd]} {
	fconfigure $fd -blocking 1
	set rc 0
	if {[catch {close $fd} err opts]} {
	    set ec [dict get $opts -errorcode]
	    set rc [expr {[lindex $ec 0] eq "CHILDSTATUS" ? [lindex $ec 2] : 1}]
	}
	set fd {}
	set running 0
	set ::ogf::batch::rc $rc
	::ogf::batch::update_view
	::ogf::status "Batch: finished (exit $rc) - $::ogf::batch::outdir/summary.tsv"
    }
}

proc ::ogf::batch::poll {} {
    if {$::ogf::batch::running} {
	::ogf::batch::update_view
	after 700 ::ogf::batch::poll
    }
}

# status.json -> tree rows; returns {done total}
proc ::ogf::batch::update_view {} {
    variable outdir
    set f [file join $outdir status.json]
    if {![file exists $f] || ![winfo exists .ogfbatch]} {return {0 0}}
    if {[catch {set fh [open $f r]; set txt [read $fh]; close $fh; set st [::ogf::json::parse $txt]}]} {return {0 0}}
    set t .ogfbatch.tree
    set done 0; set total 0
    foreach name [lsort [dict keys $st]] {
	set d [dict get $st $name]
	set state [::ogf::json::get $d state]; set step [::ogf::json::get $d step]
	set i [::ogf::json::get $d index 0]; set n [::ogf::json::get $d total 0]
	incr total
	if {$state in {ok failed}} {incr done}
	set vals [list $name $state $step "$i/$n"]
	if {[$t exists $name]} {$t item $name -values $vals} else {$t insert {} end -id $name -values $vals}
	$t tag add $state $name
    }
    .ogfbatch.pb configure -maximum [expr {max($total,1)}] -value $done
    .ogfbatch.status configure -text "$done of $total fields finished[expr {$::ogf::batch::running ? {} : {  (finished)}}]"
    return [list $done $total]
}

proc OGFBatchProgress {} {
    set w .ogfbatch
    if {![winfo exists $w]} {
	toplevel $w
	wm title $w "Batch progress"
	ttk::treeview $w.tree -columns {field state step progress} -show headings -height 10
	foreach {c t wd} {field Field 140 state State 80 step Step 180 progress Steps 70} {$w.tree heading $c -text $t; $w.tree column $c -width $wd}
	$w.tree tag configure failed -foreground red
	$w.tree tag configure ok -foreground darkgreen
	ttk::progressbar $w.pb -mode determinate
	ttk::label $w.status -text ""
	text $w.log -height 8 -width 80
	frame $w.b
	ttk::button $w.b.cancel -text Cancel -command ::ogf::batch::cancel
	ttk::button $w.b.resume -text "Resume" -command {if {!$::ogf::batch::running} {::ogf::params::put batch resume 1; OGFBatchRun}}
	ttk::button $w.b.sum -text "Summary table" -command OGFBatchSummary
	ttk::button $w.b.close -text Close -command [list destroy $w]
	pack $w.b.cancel $w.b.resume $w.b.sum $w.b.close -side left -padx 3
	pack $w.tree -fill both -expand 1
	pack $w.pb -fill x -padx 4 -pady 2
	pack $w.status
	pack $w.log -fill both -expand 1
	pack $w.b -pady 3
    }
    if {$::ogf::batch::outdir eq {}} {set ::ogf::batch::outdir [::ogf::params::get batch outdir]}
    ::ogf::batch::update_view
    return $w
}

proc ::ogf::batch::cancel {} {
    variable pid; variable running
    if {$running && $pid ne {}} {
	catch {exec kill $pid}
	::ogf::status "Batch: cancelled - run again with Resume to continue"
    }
}

proc OGFBatchSummary {} {
    set d [expr {$::ogf::batch::outdir ne {} ? $::ogf::batch::outdir : [::ogf::params::get batch outdir]}]
    set f [file join $d summary.tsv]
    if {![file exists $f]} {::ogf::status "Batch: no summary yet"; return}
    set fh [open $f r]; set txt [read $fh]; close $fh
    OGFTextWindow "Batch: summary" $txt
}
