# Cross-match (plugins/xmatch): summary and pair-table windows.

proc OGFXmatchFile {name} {return [file join [OGFSessWorkDir] xmatch $name]}

proc OGFXmatchAfter {} {
    foreach {k f} {summary_file xmatch_summary.json pairs_file xmatch_pairs.tsv} {
	set p [OGFXmatchFile $f]
	if {[file exists $p]} {::ogf::cat::set xmatch,$k $p}
    }
}

proc OGFXmatchSummary {} {
    set f [OGFXmatchFile xmatch_summary.json]
    if {![file exists $f]} {::ogf::status "Cross-match: nothing yet - run Cross-match first"; return}
    set fd [open $f r]; set txt [read $fd]; close $fd
    OGFTextWindow "Cross-match: summary" $txt
}

proc OGFXmatchPairs {} {
    set f [OGFXmatchFile xmatch_pairs.tsv]
    if {![file exists $f]} {::ogf::status "Cross-match: nothing yet - run Cross-match first"; return}
    set fd [open $f r]; set txt [read $fd]; close $fd
    OGFTextWindow "Cross-match: pairs" $txt
}
