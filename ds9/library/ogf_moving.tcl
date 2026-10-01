#  OGFinder: "Moving Objects" menu (asteroids / moving objects / static transients in multi-epoch imaging).
#  All computation is done by ds9_moving.py (package <root>/moving); this file only provides the menu, dialogs,
#  the result table, region markers and the orbit-fit window.  Every computational step is ONE argv list with explicit
#  inputs/outputs, recorded through OGFSessLog when the session recorder exists.
#  See docs/moving_objects.md

package provide DS9 1.0

proc OGFMovingInit {} {
    global ogfmov
    array set ogfmov {
	workdir {} files {} ra {} dec {} radius 1.0 provider mast collections HST,JWST,PS1
	download 0 filter {} localdir {} snr 8.0 tol 1.0 halfpix 700 status {} busy 0 catalog {}
	tracklet 0 designation {} mjdmin {} mjdmax {} minep 2
    }
    set ogfmov(workdir) [file join [file normalize ~] .ds9 moving_work]
}

# ---------------------------------------------------------------- menu
# Single top-level menu "Moving Objects", entries in workflow order.
proc OGFMovingMenu {bar} {
    set mb $bar.moving
    ttk::menubutton $mb -text "Moving Objects" -menu $mb.m -style CatMenu.TMenubutton
    menu $mb.m -tearoff 0
    $mb.m add command -label "Fetch Reference..." -command OGFMovFetchDialog
    $mb.m add command -label "Align" -command OGFMovAlign
    $mb.m add command -label "Difference" -command OGFMovDifference
    $mb.m add command -label "Detect" -command OGFMovDetect
    $mb.m add command -label "Link Tracklets" -command OGFMovLink
    $mb.m add command -label "Identify Known Objects" -command OGFMovIdentify
    $mb.m add command -label "Orbit Fit..." -command OGFMovOrbitDialog
    $mb.m add command -label "Transient Candidates..." -command OGFMovTransients
    $mb.m add command -label "Light Curve" -command OGFMovLightCurve
    $mb.m add command -label "Export..." -command OGFMovExport
    $mb.m add separator
    $mb.m add command -label "Select Exposures..." -command OGFMovSelectFiles
    $mb.m add command -label "Show Results Table" -command OGFMovTable
    $mb.m add command -label "Work Directory..." -command OGFMovWorkdir
    return $mb
}

proc OGFMovScript {} {
    return [CatalogPanelGetScript ds9_moving.py]
}

proc OGFMovStatus {msg} {
    global ogfmov catpanel
    set ogfmov(status) $msg
    catch {set catpanel(status) "Moving: $msg"}
    catch {.ogfmovtbl.st configure -text $msg}
    update idletasks
}

proc OGFMovWork {} {
    global ogfmov
    file mkdir $ogfmov(workdir)
    return $ogfmov(workdir)
}

proc OGFMovWorkdir {} {
    global ogfmov
    set d [tk_chooseDirectory -initialdir $ogfmov(workdir) -title "Moving Objects work directory"]
    if {$d ne {}} {set ogfmov(workdir) $d}
}

proc OGFMovSelectFiles {} {
    global ogfmov
    set init [file join [file normalize ~] .ds9 mast_cache]
    set fl [tk_getOpenFile -multiple 1 -initialdir $init -title "Select calibrated exposures (flc/flt/cal/FITS)" \
	-filetypes {{{FITS} {.fits .fit .fits.gz}} {{All files} *}}]
    if {$fl ne {}} {set ogfmov(files) $fl; OGFMovStatus "[llength $fl] exposure file(s) selected"}
}

# current image centre (fk5 degrees) as default position
proc OGFMovDefaultPos {} {
    global current ogfmov
    if {$ogfmov(ra) ne {}} return
    if {$current(frame) eq {}} return
    if {[catch {set p [$current(frame) get pan wcs fk5 degrees]}]} return
    if {[llength $p] == 2} {
	set ogfmov(ra) [format %.5f [lindex $p 0]]; set ogfmov(dec) [format %.5f [lindex $p 1]]
    }
}

# ---------------------------------------------------------------- generic runner (async, logged)
# argv: full command (python, script, args...).  done: script called with {ok output} when finished.
proc OGFMovRun {step title extra done {network 0} {outputs {}}} {
    global ogfmov
    if {$ogfmov(busy)} {OGFMovStatus "another Moving Objects step is still running"; return}
    set script [OGFMovScript]
    if {![file exists $script]} {OGFMovStatus "ds9_moving.py not found"; return}
    set argv [list [OGFPython] $script {*}$extra]
    set seq 0
    if {[info commands OGFSessLog] ne {}} {
	catch {set seq [OGFSessLog $step auto $argv -title $title -network $network -outputs $outputs -tool python]}
    }
    set ogfmov(busy) 1
    set ogfmov(out) {}
    set ogfmov(seq) $seq
    set ogfmov(t0) [clock milliseconds]
    OGFMovStatus "$title ..."
    OGFPrepareLibPath
    if {[catch {set fd [open |[concat $argv {2>@1}] r]} err]} {
	set ogfmov(busy) 0
	OGFMovStatus "cannot start: [string range $err 0 100]"
	return
    }
    fconfigure $fd -blocking 0 -buffering line
    fileevent $fd readable [list OGFMovReadable $fd $done $title]
}

proc OGFMovReadable {fd done title} {
    global ogfmov
    if {[catch {set chunk [read $fd]}]} {set chunk {}}
    append ogfmov(out) $chunk
    foreach line [split $chunk \n] {
	if {[string match "#MOVING*" $line]} {
	    if {[string match "*provider_unavailable*" $line]} {
		regexp {"reason": "([^"]*)"} $line -> why
		OGFMovStatus "LSST provider unavailable: [string range $why 0 140]"
	    }
	} elseif {$line ne {} && ![string match "*Warning*" $line] && ![string match "  *" $line]} {
	    OGFMovStatus "$title: [string range $line 0 110]"
	}
    }
    if {[eof $fd]} {
	fconfigure $fd -blocking 1
	set rc [catch {close $fd} err]
	set ogfmov(busy) 0
	set ms [expr {[clock milliseconds] - $ogfmov(t0)}]
	if {$ogfmov(seq) > 0} {
	    catch {OGFSessSet $ogfmov(seq) ms $ms}
	    if {$rc} {catch {OGFSessSet $ogfmov(seq) failed 1}}
	}
	if {$rc} {OGFMovStatus "$title failed: [string range $err 0 140]"} else {OGFMovStatus "$title done ([format %.1f [expr {$ms/1000.0}]] s)"}
	catch {{*}$done [expr {!$rc}] $ogfmov(out)}
    }
}

# ---------------------------------------------------------------- Fetch Reference
proc OGFMovFetchDialog {} {
    global ogfmov
    OGFMovDefaultPos
    set w .ogfmovfetch
    catch {destroy $w}
    toplevel $w
    wm title $w "Moving Objects: Fetch Reference"
    wm transient $w .
    set r 0
    foreach {k lab} {ra "RA (deg)" dec "Dec (deg)" radius "Search radius (arcmin)" collections "MAST collections"
	filter "Filter contains" localdir "Rubin DP1 local dir (optional)"} {
	ttk::label $w.l$r -text $lab
	ttk::entry $w.e$r -textvariable ogfmov($k) -width 22
	grid $w.l$r -row $r -column 0 -sticky w -padx 8 -pady 2
	grid $w.e$r -row $r -column 1 -sticky we -padx 8 -pady 2
	incr r
    }
    ttk::label $w.lp -text "Reference provider"
    ttk::frame $w.pf
    foreach {v t} {mast MAST lsst LSST both both} {
	ttk::radiobutton $w.pf.$v -text $t -variable ogfmov(provider) -value $v
	pack $w.pf.$v -side left -padx 4
    }
    grid $w.lp -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $w.pf -row $r -column 1 -sticky w -padx 8; incr r
    ttk::checkbutton $w.dl -text "Download calibrated exposures (MAST, up to 6; large files)" -variable ogfmov(download)
    grid $w.dl -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady 2; incr r
    ttk::label $w.note -wraplength 380 -justify left -text \
	"LSST/Rubin: public no-credential access only. DP1 images need a Rubin Science Platform login, so the LSST provider is reported unavailable unless you supply downloaded DP1 FITS above; MAST/Pan-STARRS is then used. Nothing is stored or sent."
    grid $w.note -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady 4; incr r
    ttk::frame $w.bb
    ttk::button $w.bb.ok -text Fetch -command {OGFMovFetchRun}
    ttk::button $w.bb.c -text Close -command {destroy .ogfmovfetch}
    pack $w.bb.ok $w.bb.c -side left -padx 4
    grid $w.bb -row $r -column 0 -columnspan 2 -pady 6
}

proc OGFMovFetchRun {} {
    global ogfmov
    if {![string is double -strict $ogfmov(ra)] || ![string is double -strict $ogfmov(dec)]} {OGFMovStatus "RA/Dec required"; return}
    set wd [OGFMovWork]
    set extra [list --mode fetch --workdir $wd --ra $ogfmov(ra) --dec $ogfmov(dec) --radius $ogfmov(radius) \
	--provider $ogfmov(provider) --collections $ogfmov(collections)]
    if {$ogfmov(filter) ne {}} {lappend extra --filter $ogfmov(filter)}
    if {$ogfmov(localdir) ne {}} {lappend extra --lsst-local-dir $ogfmov(localdir)}
    if {$ogfmov(download)} {lappend extra --download}
    OGFMovRun moving.fetch "Fetch reference ($ogfmov(provider))" $extra OGFMovFetchDone 1 [list [file join $wd candidates.tsv]]
}

proc OGFMovFetchDone {ok out} {
    global ogfmov
    if {!$ok} return
    set wd [OGFMovWork]
    set fj [file join $wd files.json]
    if {$ogfmov(download) && [file exists $fj]} {
	set fd [open $fj r]; set txt [read $fd]; close $fd
	set fl [regexp -all -inline {"path": "([^"]+)"} $txt]
	set files {}
	foreach {m p} $fl {lappend files $p}
	if {[llength $files]} {set ogfmov(files) $files}
    }
    OGFMovShowCandidates [file join $wd candidates.tsv]
}

proc OGFMovReadTSV {fn} {
    if {![file exists $fn]} {return {{} {}}}
    set fd [open $fn r]; set txt [read $fd]; close $fd
    set lines [split [string trimright $txt \n] \n]
    set cols [split [lindex $lines 0] \t]
    set rows {}
    foreach l [lrange $lines 1 end] {lappend rows [split $l \t]}
    return [list $cols $rows]
}

# ---------------------------------------------------------------- results table (own window)
proc OGFMovTable {} {
    set w .ogfmovtbl
    if {![winfo exists $w]} {
	toplevel $w
	wm title $w "Moving Objects: results"
	wm geometry $w 760x340
	ttk::label $w.st -text {} -anchor w
	ttk::treeview $w.t -show headings -selectmode browse -yscrollcommand [list $w.sb set]
	ttk::scrollbar $w.sb -command [list $w.t yview]
	pack $w.st -side bottom -fill x
	pack $w.sb -side right -fill y
	pack $w.t -side left -fill both -expand 1
	bind $w.t <<TreeviewSelect>> OGFMovRowSelected
	bind $w.t <Double-Button-1> OGFMovRowSelected
    }
    raise $w
    return $w
}

proc OGFMovFillTable {kind cols rows} {
    global ogfmov
    set w [OGFMovTable]
    set ogfmov(tblkind) $kind
    $w.t delete [$w.t children {}]
    $w.t configure -columns $cols
    foreach c $cols {$w.t heading $c -text $c; $w.t column $c -width 80 -anchor w}
    foreach r $rows {$w.t insert {} end -values $r}
    wm title $w "Moving Objects: $kind ([llength $rows])"
}

proc OGFMovShowCandidates {fn} {
    lassign [OGFMovReadTSV $fn] cols rows
    set keep {provider obs_id instrument filter texp mjd epoch pixscale depth proposal kind}
    set idx {}
    foreach k $keep {set i [lsearch $cols $k]; if {$i >= 0} {lappend idx $i}}
    set rr {}
    foreach r $rows {set x {}; foreach i $idx {lappend x [lindex $r $i]}; lappend rr $x}
    set cc {}; foreach i $idx {lappend cc [lindex $cols $i]}
    OGFMovFillTable "reference candidates" $cc $rr
}

proc OGFMovShowMovers {} {
    set fn [file join [OGFMovWork] movers.tsv]
    lassign [OGFMovReadTSV $fn] cols rows
    set keep {id status n rate_arcsec_h pa_deg rms_arcsec score ra dec}
    set idx {}
    foreach k $keep {set i [lsearch $cols $k]; if {$i >= 0} {lappend idx $i}}
    set rr {}
    foreach r $rows {
	set x {}
	foreach i $idx {
	    set v [lindex $r $i]
	    if {[string is double -strict $v] && [lindex $cols $i] ni {id n}} {set v [format %.4g $v]}
	    lappend x $v
	}
	lappend rr $x
    }
    set cc {}; foreach i $idx {lappend cc [lindex $cols $i]}
    OGFMovFillTable movers $cc $rr
    OGFMovLoadRegions [file join [OGFMovWork] movers.reg]
}

proc OGFMovShowTransients {} {
    set fn [file join [OGFMovWork] transients.tsv]
    lassign [OGFMovReadTSV $fn] cols rows
    set keep {id ra dec n_det files snr_max host_sep_arcsec offset_re host_z heuristic}
    set idx {}
    foreach k $keep {set i [lsearch $cols $k]; if {$i >= 0} {lappend idx $i}}
    set rr {}
    foreach r $rows {
	set x {}
	foreach i $idx {set v [lindex $r $i]; if {[string is double -strict $v] && [lindex $cols $i] ni {id n_det files}} {set v [format %.5g $v]}; lappend x $v}
	lappend rr $x
    }
    set cc {}; foreach i $idx {lappend cc [lindex $cols $i]}
    OGFMovFillTable transients $cc $rr
    OGFMovLoadRegions [file join [OGFMovWork] transients.reg]
}

proc OGFMovShowDetections {} {
    global ogfmov
    set fn [file join [OGFMovWork] detections.tsv]
    lassign [OGFMovReadTSV $fn] cols rows
    set keep {id ex chip cls snr sign flux_e_s elong a_pix ra dec}
    set idx {}
    foreach k $keep {set i [lsearch $cols $k]; if {$i >= 0} {lappend idx $i}}
    set ci [lsearch $cols snr]; set cl [lsearch $cols cls]; set sg [lsearch $cols sign]
    set rr {}
    foreach r $rows {
	if {[lindex $r $sg] < 0 && !$ogfmov(shownegative)} continue
	if {[lindex $r $ci] < $ogfmov(snr)} continue
	set x {}
	foreach i $idx {set v [lindex $r $i]; if {[string is double -strict $v] && [lindex $cols $i] ni {id ex sign}} {set v [format %.4g $v]}; lappend x $v}
	lappend rr $x
    }
    set cc {}; foreach i $idx {lappend cc [lindex $cols $i]}
    OGFMovFillTable detections $cc $rr
    return [llength $rr]
}

proc OGFMovRowSelected {} {
    global ogfmov current
    set w .ogfmovtbl
    set sel [$w.t selection]
    if {$sel eq {}} return
    set cols [$w.t cget -columns]
    set vals [$w.t item [lindex $sel 0] -values]
    set ra {}; set dec {}
    foreach c $cols v $vals {
	if {$c eq "ra"} {set ra $v}; if {$c eq "dec"} {set dec $v}
	if {$c eq "id"} {set id $v}
    }
    set info [join [lmap c $cols v $vals {string cat $c = $v}] "   "]
    OGFMovStatus [string range $info 0 220]
    set ogfmov(selected) [expr {[info exists id] ? $id : {}}]
    if {$ogfmov(tblkind) eq "movers" && [info exists id]} {set ogfmov(tracklet) $id}
    if {[string is double -strict $ra] && [string is double -strict $dec] && $current(frame) ne {}} {
	catch {$current(frame) pan to wcs fk5 degrees $ra $dec}
    }
}

proc OGFMovLoadRegions {fn} {
    global current
    if {![file exists $fn] || $current(frame) eq {}} return
    if {[catch {MarkerLoadFile $fn $current(frame) ds9 wcs fk5} err]} {OGFMovStatus "region load failed: $err"}
}

# ---------------------------------------------------------------- pipeline steps
proc OGFMovNeedFiles {} {
    global ogfmov
    if {[llength $ogfmov(files)] == 0} {
	OGFMovSelectFiles
    }
    if {[llength $ogfmov(files)] == 0} {OGFMovStatus "no exposures selected"; return 0}
    return 1
}

proc OGFMovAlign {} {
    global ogfmov
    if {![OGFMovNeedFiles]} return
    set wd [OGFMovWork]
    OGFMovRun moving.align "Align exposures (Gaia DR3 / relative)" \
	[concat [list --mode align --workdir $wd --files] $ogfmov(files)] {OGFMovDone} 1 [list [file join $wd align_report.json]]
}

proc OGFMovDifference {} {
    global ogfmov
    if {![OGFMovNeedFiles]} return
    OGFMovDefaultPos
    set wd [OGFMovWork]
    set extra [concat [list --mode difference --workdir $wd --half-pix $ogfmov(halfpix) --snr 5 --files] $ogfmov(files)]
    if {[string is double -strict $ogfmov(ra)] && [string is double -strict $ogfmov(dec)]} {
	lappend extra --ra $ogfmov(ra) --dec $ogfmov(dec)
    }
    OGFMovRun moving.difference "Difference (template + ZOGY) and detection" $extra OGFMovDifferenceDone 0 \
	[list [file join $wd detections.tsv] [file join $wd diff]]
}

proc OGFMovDifferenceDone {ok out} {
    if {$ok} {OGFMovDetect}
}

# "Detect" shows the detections of the last Difference run (S/N and class filter), no recomputation
proc OGFMovDetect {} {
    global ogfmov
    set fn [file join [OGFMovWork] detections.tsv]
    if {![file exists $fn]} {OGFMovStatus "run Difference first"; return}
    if {![info exists ogfmov(shownegative)]} {set ogfmov(shownegative) 0}
    set n [OGFMovShowDetections]
    OGFMovStatus "$n detections with S/N >= $ogfmov(snr) (classes: point, trail, artefact_*; see 'cls')"
}

proc OGFMovLink {} {
    global ogfmov
    if {![OGFMovNeedFiles]} return
    set wd [OGFMovWork]
    OGFMovRun moving.link "Link tracklets" \
	[concat [list --mode link --workdir $wd --snr $ogfmov(snr) --tol $ogfmov(tol) --max-per-exposure 900 --max-tracklets 400 --files] $ogfmov(files)] \
	OGFMovLinkDone 0 [list [file join $wd tracklets.json] [file join $wd movers.tsv] [file join $wd movers.reg]]
}

proc OGFMovLinkDone {ok out} {
    if {$ok} {OGFMovShowMovers}
}

proc OGFMovIdentify {} {
    set wd [OGFMovWork]
    if {![file exists [file join $wd tracklets.json]]} {OGFMovStatus "run Link Tracklets first"; return}
    OGFMovRun moving.identify "Identify known objects (SkyBoT + Horizons)" \
	[list --mode identify --workdir $wd --radius 6 --tol 1.5 --max-tracklets 15] OGFMovIdentifyDone 1 [list [file join $wd identified.tsv]]
}

proc OGFMovIdentifyDone {ok out} {
    if {$ok} {OGFMovShowMovers}
}

proc OGFMovTransients {} {
    global ogfmov
    set w .ogfmovtr
    catch {destroy $w}
    set fl [list [list catalog "Host catalogue TSV (optional)" $ogfmov(catalog)] [list minep "Min. exposure files" $ogfmov(minep)]]
    set r [OGFForm "Moving Objects: transient candidates" $fl {} \
	"Static, point-like, positive detections seen in >= N different exposure files, not part of a mover tracklet. Host = nearest catalogue galaxy within 5 arcsec."]
    if {$r eq {}} return
    set wd [OGFMovWork]
    set ogfmov(catalog) [dict get $r catalog]
    set extra [list --mode transients --workdir $wd --min-epochs [dict get $r minep] --snr $ogfmov(snr)]
    if {$ogfmov(catalog) ne {}} {lappend extra --catalog $ogfmov(catalog)}
    OGFMovRun moving.transients "Transient candidates" $extra OGFMovTransientsDone 0 [list [file join $wd transients.tsv] [file join $wd transients.reg]]
}

proc OGFMovTransientsDone {ok out} {
    if {$ok} {OGFMovShowTransients}
}

proc OGFMovLightCurve {} {
    global ogfmov
    if {![OGFMovNeedFiles]} return
    set wd [OGFMovWork]
    OGFMovRun moving.lightcurve "Light curves (forced photometry on difference images)" \
	[concat [list --mode lightcurve --workdir $wd --files] $ogfmov(files)] OGFMovLCDone 0 [list [file join $wd lightcurves.json]]
}

proc OGFMovLCDone {ok out} {
    if {!$ok} return
    set w .ogfmovlc
    catch {destroy $w}
    toplevel $w
    wm title $w "Moving Objects: light curves"
    text $w.t -width 90 -height 18 -font TkFixedFont
    pack $w.t -fill both -expand 1
    foreach l [split $out \n] {if {[string match "T*:*" $l]} {$w.t insert end "$l\n"}}
}

proc OGFMovExport {} {
    global ogfmov
    set wd [OGFMovWork]
    set oj [lindex [lsort [glob -nocomplain -directory $wd orbit_*.json]] end]
    if {$oj eq {}} {OGFMovStatus "run Orbit Fit first"; return}
    set extra [list --mode export --workdir $wd --orbit-json [file tail $oj]]
    if {[string match "orbit_tracklet*" [file tail $oj]]} {
	regexp {tracklet(\d+)} $oj -> n
	lappend extra --tracklet $n
    }
    OGFMovRun moving.export "Export orbit (MPC 80-col + JSON, local files only)" $extra OGFMovExportDone 0 [list [file join $wd elements_export.json]]
}

proc OGFMovExportDone {ok out} {
    if {$ok} {OGFMovStatus "exported to [OGFMovWork] (nothing submitted to the MPC)"}
}

proc OGFMovDone {ok out} {}

# ---------------------------------------------------------------- Orbit fit dialog
proc OGFMovOrbitDialog {} {
    global ogfmov
    set w .ogfmovorb
    catch {destroy $w}
    toplevel $w
    wm title $w "Moving Objects: Orbit Fit"
    ttk::frame $w.top
    ttk::label $w.top.l1 -text "Tracklet id"
    ttk::entry $w.top.e1 -textvariable ogfmov(tracklet) -width 6
    ttk::label $w.top.l2 -text "or MPC designation / number"
    ttk::entry $w.top.e2 -textvariable ogfmov(designation) -width 12
    ttk::label $w.top.l3 -text "MJD range"
    ttk::entry $w.top.e3 -textvariable ogfmov(mjdmin) -width 9
    ttk::entry $w.top.e4 -textvariable ogfmov(mjdmax) -width 9
    ttk::button $w.top.go -text "Fit" -command OGFMovOrbitRun
    pack $w.top.l1 $w.top.e1 $w.top.l2 $w.top.e2 $w.top.l3 $w.top.e3 $w.top.e4 $w.top.go -side left -padx 3
    pack $w.top -side top -fill x -pady 4
    ttk::label $w.note -wraplength 640 -justify left -text \
	"Tracklet fit: ASSIST (DE440 + 16 massive asteroids + GR) differential correction seeded by statistical ranging; HST parallax from Horizons. A designation downloads the public MPC observations (debiased, default station weights) and fits them."
    pack $w.note -side top -fill x -padx 6
    text $w.res -width 84 -height 11 -font TkFixedFont
    pack $w.res -side top -fill x -padx 6 -pady 4
    canvas $w.plot -width 640 -height 220 -background white
    pack $w.plot -side top -padx 6 -pady 4
    ttk::frame $w.bb
    ttk::button $w.bb.ex -text "Export (MPC 80-col + JSON)" -command OGFMovExport
    ttk::button $w.bb.cl -text Close -command [list destroy $w]
    pack $w.bb.ex $w.bb.cl -side left -padx 4
    pack $w.bb -side top -pady 4
}

proc OGFMovOrbitRun {} {
    global ogfmov
    set wd [OGFMovWork]
    set extra [list --mode orbit --workdir $wd]
    if {$ogfmov(designation) ne {}} {
	lappend extra --designation $ogfmov(designation)
	if {$ogfmov(mjdmin) ne {}} {lappend extra --mjd-min $ogfmov(mjdmin)}
	if {$ogfmov(mjdmax) ne {}} {lappend extra --mjd-max $ogfmov(mjdmax)}
	set label $ogfmov(designation)
    } else {
	lappend extra --tracklet $ogfmov(tracklet)
	set label tracklet$ogfmov(tracklet)
    }
    set ogfmov(orbitlabel) $label
    OGFMovRun moving.orbit "Orbit fit ($label)" $extra OGFMovOrbitDone [expr {$ogfmov(designation) ne {}}] \
	[list [file join $wd orbit_$label.json]]
}

proc OGFMovOrbitDone {ok out} {
    global ogfmov
    if {!$ok} return
    set w .ogfmovorb
    if {![winfo exists $w]} {OGFMovOrbitDialog}
    set fn [file join [OGFMovWork] orbit_$ogfmov(orbitlabel).kv]
    if {![file exists $fn]} return
    set fd [open $fn r]; set txt [read $fd]; close $fd
    $w.res delete 1.0 end
    set res {}
    $w.res insert end [format "%-6s %16s %12s\n" element value "+/- 1-sigma"]
    foreach l [split $txt \n] {
	set f [split $l \t]
	switch [lindex $f 0] {
	    elem {
		if {[string is double -strict [lindex $f 2]] && [lindex $f 2] ne "nan"} {
		    $w.res insert end [format "%-6s %16.8g %12.3g\n" [lindex $f 1] [lindex $f 2] [lindex $f 3]]
		} else {
		    $w.res insert end [format "%-6s %16s %12s\n" [lindex $f 1] "undetermined" "-"]
		}
	    }
	    info {$w.res insert end "[lindex $f 1] = [lindex $f 2]\n"}
	    rq {$w.res insert end [format "ranging %-2s 16/50/84%%: %.4g / %.4g / %.4g\n" [lindex $f 1] [lindex $f 2] [lindex $f 3] [lindex $f 4]]}
	    class {$w.res insert end [format "class probability %-28s %.3f\n" [lindex $f 1] [lindex $f 2]]}
	    res {lappend res [lrange $f 1 end]}
	}
    }
    OGFMovPlotResiduals $w.plot $res
}

proc OGFMovPlotResiduals {c res} {
    $c delete all
    if {[llength $res] == 0} return
    set W [$c cget -width]; set H [$c cget -height]
    set l 50; set r 10; set t 14; set b 28
    set tmin 1e99; set tmax -1e99; set ymax 0.0
    foreach r_ $res {
	lassign $r_ mjd dra dde act stn
	set tmin [expr {min($tmin,$mjd)}]; set tmax [expr {max($tmax,$mjd)}]
	set ymax [expr {max($ymax,abs($dra),abs($dde))}]
    }
    if {$tmax <= $tmin} {set tmax [expr {$tmin+1e-3}]}
    set ymax [expr {$ymax*1.15}]
    if {$ymax <= 0} {set ymax 1.0}
    set x0 $l; set x1 [expr {$W-$r}]; set y0 [expr {$H-$b}]; set y1 $t
    $c create rectangle $x0 $y1 $x1 $y0 -outline gray50
    set ym [expr {($y0+$y1)/2.0}]
    $c create line $x0 $ym $x1 $ym -fill gray70 -dash {2 2}
    $c create text [expr {$l-4}] $y1 -anchor e -text [format %.2g $ymax] -font TkSmallCaptionFont
    $c create text [expr {$l-4}] $y0 -anchor e -text [format %.2g -$ymax] -font TkSmallCaptionFont
    $c create text [expr {$l-4}] $ym -anchor e -text 0 -font TkSmallCaptionFont
    $c create text [expr {($x0+$x1)/2}] [expr {$H-4}] -anchor s -font TkSmallCaptionFont \
	-text [format "UTC MJD %.4f .. %.4f   (O-C, arcsec;  o = RA*cos(dec),  + = Dec,  grey = rejected)" $tmin $tmax]
    foreach r_ $res {
	lassign $r_ mjd dra dde act stn
	set x [expr {$x0 + ($mjd-$tmin)/($tmax-$tmin)*($x1-$x0)}]
	set col [expr {$act ? "#1f4e9c" : "#aaaaaa"}]
	set col2 [expr {$act ? "#c0392b" : "#aaaaaa"}]
	set ya [expr {$ym - $dra/$ymax*($ym-$y1)}]; set yb [expr {$ym - $dde/$ymax*($ym-$y1)}]
	$c create oval [expr {$x-3}] [expr {$ya-3}] [expr {$x+3}] [expr {$ya+3}] -outline $col
	$c create line [expr {$x-3}] $yb [expr {$x+3}] $yb -fill $col2
	$c create line $x [expr {$yb-3}] $x [expr {$yb+3}] -fill $col2
    }
}
