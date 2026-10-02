#  OGFinder plugin "moving": asteroids / moving objects / static transients in multi-epoch imaging.
#  All computation is done by ds9_moving.py (package <root>/moving); this file provides the step procs, the
#  Fetch Reference dialog, the feed of the shared time-domain table (kinds moving / transient / detection, see
#  ogf_td.tcl) and the single "Time-domain details" window (Orbit and Light curve tabs).  Every computational step
#  is ONE argv list with explicit inputs/outputs, recorded through OGFSessLog when the session recorder exists
#  (unchanged by the move into the plugin system).
#  See docs/moving_objects.md and docs/plugins.md

package provide DS9 1.0

proc OGFMovingInit {} {
    global ogfmov
    array set ogfmov {
	workdir {} files {} ra {} dec {} radius 1.0 provider mast collections HST,JWST,PS1
	download 0 filter {} localdir {} snr 8.0 tol 1.0 halfpix 700 status {} busy 0 catalog {}
	tracklet 0 designation {} mjdmin {} mjdmax {} minep 2 shownegative 0 lcmode mag lc {} orbitlabel {} candcols {} candrows {} detwin .ogftd selkind {} selected {}
    }
    set ogfmov(workdir) [file join [file normalize ~] .ds9 moving_work]
    # the three row kinds of the time-domain table
    ::ogf::td::register_kind moving -label Moving -prefix M -color magenta -order 10 \
	-columns {n rate_arcsec_h pa_deg rms_arcsec score status overlap overlap_id} -decorate OGFMovDecorateMovers \
	-point {point=circle 7}
    ::ogf::td::register_kind transient -label Transients -prefix T -color red -order 20 \
	-columns {n_det files snr_max host_id host_sep host_z offset_re heuristic} -decorate OGFMovDecorateTransients \
	-point {point=cross 12}
    ::ogf::td::register_kind detection -label Detections -prefix D -color white -order 30 -show_in_all 0 \
	-columns {ex chip cls snr sign flux_e_s elong a_pix} -point {point=boxcircle 8}
    ::ogf::td::on_select OGFMovOnSelect
    # results that exist in the work directory from an earlier session are shown on request (no auto-load)
}

proc OGFMovScript {} {
    return [CatalogPanelGetScript ds9_moving.py]
}

proc OGFMovStatus {msg} {
    global ogfmov
    set ogfmov(status) $msg
    catch {::ogf::cat::set status "Moving: $msg"}
    catch {$::ogfmov(detwin).st configure -text $msg}
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
    if {[winfo exists $w]} {raise $w; return}
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
    grid $w.bb -row $r -column 0 -columnspan 2 -pady 6; incr r
    # reference candidates found by the fetch are listed in this dialog (they are observations, not sky objects)
    ttk::treeview $w.cand -show headings -height 6 -selectmode browse
    grid $w.cand -row $r -column 0 -columnspan 2 -sticky we -padx 8 -pady {0 8}
    grid columnconfigure $w 1 -weight 1
    if {$ogfmov(candrows) ne {}} {OGFMovFillCandidates}
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

proc OGFMovShowCandidates {fn} {
    global ogfmov
    lassign [OGFMovReadTSV $fn] cols rows
    set keep {provider obs_id instrument filter texp mjd epoch pixscale depth proposal kind}
    set idx {}
    foreach k $keep {set i [lsearch $cols $k]; if {$i >= 0} {lappend idx $i}}
    set rr {}
    foreach r $rows {set x {}; foreach i $idx {lappend x [lindex $r $i]}; lappend rr $x}
    set cc {}; foreach i $idx {lappend cc [lindex $cols $i]}
    set ogfmov(candcols) $cc; set ogfmov(candrows) $rr
    OGFMovFillCandidates
    OGFMovStatus "[llength $rr] reference candidates"
}

proc OGFMovFillCandidates {} {
    global ogfmov
    set w .ogfmovfetch
    if {![winfo exists $w]} {OGFMovFetchDialog; return}
    set t $w.cand
    $t delete [$t children {}]
    $t configure -columns $ogfmov(candcols)
    foreach c $ogfmov(candcols) {$t heading $c -text $c; $t column $c -width 70 -anchor w}
    foreach r $ogfmov(candrows) {$t insert {} end -values $r}
}

# ---------------------------------------------------------------- results -> the shared time-domain table
proc OGFMovFmt {cols i v} {
    if {[string is double -strict $v] && [lindex $cols $i] ni {id n ex sign n_det files}} {return [format %.5g $v]}
    return $v
}

# {ra dec ra dec ...} per mover id from movers.reg ("# vector(...) text={M<id> ...}" followed by the circles)
proc OGFMovTracks {fn} {
    set tr [dict create]
    if {![file exists $fn]} {return $tr}
    set fd [open $fn r]; set txt [read $fd]; close $fd
    set id {}
    foreach l [split $txt \n] {
	if {[regexp {text=\{M(\d+) } $l -> id]} continue
	if {$id ne {} && [regexp {^circle\(([-0-9.]+),([-0-9.]+),} $l -> ra de]} {dict lappend tr $id $ra $de}
    }
    return $tr
}

proc OGFMovShowMovers {} {
    set wd [OGFMovWork]
    lassign [OGFMovReadTSV [file join $wd movers.tsv]] cols rows
    set tr [OGFMovTracks [file join $wd movers.reg]]
    set rr {}
    foreach r $rows {
	set x {}
	set i 0
	foreach v $r {lappend x [OGFMovFmt $cols $i $v]; incr i}
	set id [lindex $r [lsearch $cols id]]
	lappend x [expr {[dict exists $tr $id] ? [dict get $tr $id] : {}}]
	lappend rr $x
    }
    ::ogf::td::set_rows moving [concat $cols _pos] $rr 1
    OGFMovStatus "[llength $rr] moving-object tracklets"
}

proc OGFMovShowTransients {} {
    lassign [OGFMovReadTSV [file join [OGFMovWork] transients.tsv]] cols rows
    set rr {}
    foreach r $rows {
	set x {}
	set i 0
	foreach v $r {lappend x [OGFMovFmt $cols $i $v]; incr i}
	lappend rr $x
    }
    ::ogf::td::set_rows transient $cols $rr 1
    OGFMovStatus "[llength $rr] transient candidates"
}

proc OGFMovShowDetections {} {
    global ogfmov
    lassign [OGFMovReadTSV [file join [OGFMovWork] detections.tsv]] cols rows
    set ci [lsearch $cols snr]; set sg [lsearch $cols sign]; set fi [lsearch $cols flux_e_s]; set zi [lsearch $cols zp_ab]
    set keep {id ex chip cls snr sign flux_e_s elong a_pix ra dec}
    set idx {}
    foreach k $keep {set i [lsearch $cols $k]; if {$i >= 0} {lappend idx $i}}
    set rr {}
    foreach r $rows {
	if {[lindex $r $sg] < 0 && !$ogfmov(shownegative)} continue
	if {[lindex $r $ci] < $ogfmov(snr)} continue
	set x {}
	foreach i $idx {lappend x [OGFMovFmt $cols $i [lindex $r $i]]}
	set mag {}
	set f [lindex $r $fi]; set z [lindex $r $zi]
	if {[string is double -strict $f] && $f > 0 && [string is double -strict $z]} {set mag [format %.3f [expr {$z - 2.5*log10($f)}]]}
	lappend x $mag
	lappend rr $x
    }
    set cc {}; foreach i $idx {lappend cc [lindex $cols $i]}
    lappend cc mag
    ::ogf::td::set_rows detection $cc $rr 1
    return [llength $rr]
}

# computed columns -----------------------------------------------------------
# transients: nearest galaxy of the CURRENT galaxy catalog within 5 arcsec (empty without a catalog)
proc OGFMovDecorateTransients {cols rows} {
    global ogftd
    set ir [lsearch $cols ra]; set id [lsearch $cols dec]
    set drop {host_id host_sep host_sep_arcsec}
    set keepi {}
    set ncols {}
    foreach c $cols i [lsearch -all -not -exact $cols __none] {
	if {$c in $drop} continue
	lappend keepi $i; lappend ncols $c
    }
    lappend ncols host_id host_sep
    set have [::ogf::td::have_galaxies]
    set out {}
    foreach r $rows {
	set x {}
	foreach i $keepi {lappend x [lindex $r $i]}
	set hid {}; set hs {}
	if {$have} {
	    set near [::ogf::td::galaxy_near [lindex $r $ir] [lindex $r $id] $ogftd(hostradius)]
	    if {[llength $near]} {lassign [lindex $near 0] hid hs; set hs [format %.2f $hs]}
	}
	lappend x $hid $hs
	lappend out $x
    }
    return [list $ncols $out]
}

# moving objects: overlap flag when any epoch position lies inside a galaxy of the catalog (ISO radius, >= 1 arcsec)
proc OGFMovDecorateMovers {cols rows} {
    set ip [lsearch $cols _pos]
    set ir [lsearch $cols ra]; set id [lsearch $cols dec]
    set have [::ogf::td::have_galaxies]
    set mx [expr {$have ? [::ogf::td::galaxy_max_extent] : 0.0}]
    set ncols [concat $cols overlap overlap_id]
    set out {}
    foreach r $rows {
	set flag {}; set gid {}
	if {$have} {
	    set flag no
	    set pos [expr {$ip >= 0 && [lindex $r $ip] ne {} ? [lindex $r $ip] : [list [lindex $r $ir] [lindex $r $id]]}]
	    foreach {ra de} $pos {
		foreach g [::ogf::td::galaxy_near $ra $de [expr {max($mx,1.0)}]] {
		    lassign $g gn gs ge
		    if {$gs <= max($ge,1.0)} {set flag yes; set gid $gn; break}
		}
		if {$flag eq "yes"} break
	    }
	}
	lappend out [concat $r [list $flag $gid]]
    }
    return [list $ncols $out]
}

# a row of a Moving / Transients / Detections view was selected (called by the table layer)
proc OGFMovOnSelect {key raw} {
    global ogfmov
    set kind [dict get $raw kind]
    set info [join [lmap {c v} [dict remove $raw kind _pos host_sep_arcsec NUMBER] {string cat $c = $v}] "   "]
    OGFMovStatus [string range $info 0 220]
    set ogfmov(selected) [dict get $raw id]
    set ogfmov(selkind) $kind
    if {$kind eq "moving"} {set ogfmov(tracklet) [dict get $raw id]}
    if {[winfo exists $ogfmov(detwin)]} {OGFMovDetailsUpdate}
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

# ZOGY options from the plugin parameters (Tools > Plugin settings > Moving); only NON-default values are added, so the default argv is unchanged
proc OGFMovZogyArgs {} {
    set a {}
    if {![catch {::ogf::params::get moving source-noise} v] && $v} {lappend a --source-noise}
    if {![catch {::ogf::params::get moving astrom} v] && $v ni {{} off}} {lappend a --astrom $v}
    if {![catch {::ogf::params::get moving psf-tile} v] && [string is integer -strict $v] && $v > 0} {lappend a --psf-tile $v}
    if {![catch {::ogf::params::get moving psf-source} v] && $v ni {{} tiles}} {lappend a --psf-source $v}
    if {![catch {::ogf::params::get moving template-psf} v] && $v ni {{} target}} {lappend a --template-psf $v}
    return $a
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
    lappend extra {*}[OGFMovZogyArgs]
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
    global ogfmov
    if {!$ok} return
    set ogfmov(lctext) {}
    foreach l [split $out \n] {if {[string match "T*:*" $l]} {append ogfmov(lctext) "$l\n"}}
    set ogfmov(lc) [dict create]
    set fn [file join [OGFMovWork] lightcurves.json]
    if {[file exists $fn]} {
	if {![catch {set fd [open $fn r]; set txt [read $fd]; close $fd; set lcs [::ogf::json::parse $txt]}]} {
	    foreach e $lcs {dict set ogfmov(lc) [dict get $e id] $e}
	}
    }
    OGFMovDetails lightcurve
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

# ---------------------------------------------------------------- Time-domain details (ONE window, two tabs)
# Orbit fit for the selected moving object / designation, light curve for the selected transient.  The window is
# created once and reused (raised / updated); it is an ordinary toplevel, i.e. it can be moved off the main window.
proc OGFMovDetails {{tab {}}} {
    global ogfmov
    set w $ogfmov(detwin)
    if {![winfo exists $w]} {
	toplevel $w
	wm title $w "Time-domain details"
	wm geometry $w 700x620
	wm protocol $w WM_DELETE_WINDOW [list wm withdraw $w]
	ttk::notebook $w.nb
	ttk::frame $w.nb.orbit
	ttk::frame $w.nb.lc
	$w.nb add $w.nb.orbit -text Orbit
	$w.nb add $w.nb.lc -text "Light curve"
	OGFMovBuildOrbitTab $w.nb.orbit
	OGFMovBuildLCTab $w.nb.lc
	ttk::label $w.st -text {} -anchor w
	pack $w.st -side bottom -fill x
	pack $w.nb -side top -fill both -expand 1
    }
    wm deiconify $w
    raise $w
    if {$tab eq "orbit"} {$w.nb select $w.nb.orbit} elseif {$tab eq "lightcurve"} {$w.nb select $w.nb.lc}
    OGFMovDetailsUpdate
    return $w
}

# follow the selected row: moving -> Orbit tab fields, transient -> Light curve plot
proc OGFMovDetailsUpdate {} {
    global ogfmov
    set w $ogfmov(detwin)
    if {![winfo exists $w]} return
    switch -- $ogfmov(selkind) {
	moving {$w.nb select $w.nb.orbit}
	transient {$w.nb select $w.nb.lc}
    }
    OGFMovPlotLC
}

proc OGFMovBuildOrbitTab {f} {
    global ogfmov
    ttk::frame $f.top
    ttk::label $f.top.l1 -text "Tracklet id"
    ttk::entry $f.top.e1 -textvariable ogfmov(tracklet) -width 6
    ttk::label $f.top.l2 -text "or MPC designation / number"
    ttk::entry $f.top.e2 -textvariable ogfmov(designation) -width 12
    ttk::label $f.top.l3 -text "MJD range"
    ttk::entry $f.top.e3 -textvariable ogfmov(mjdmin) -width 9
    ttk::entry $f.top.e4 -textvariable ogfmov(mjdmax) -width 9
    ttk::button $f.top.go -text "Fit" -command OGFMovOrbitRun
    pack $f.top.l1 $f.top.e1 $f.top.l2 $f.top.e2 $f.top.l3 $f.top.e3 $f.top.e4 $f.top.go -side left -padx 3
    pack $f.top -side top -fill x -pady 4
    ttk::label $f.note -wraplength 640 -justify left -text \
	"Tracklet fit: ASSIST (DE440 + 16 massive asteroids + GR) differential correction seeded by statistical ranging; HST parallax from Horizons. A designation downloads the public MPC observations (debiased, default station weights) and fits them."
    pack $f.note -side top -fill x -padx 6
    text $f.res -width 84 -height 11 -font TkFixedFont
    pack $f.res -side top -fill x -padx 6 -pady 4
    canvas $f.plot -width 640 -height 220 -background white
    pack $f.plot -side top -padx 6 -pady 4
    ttk::frame $f.bb
    ttk::button $f.bb.ex -text "Export (MPC 80-col + JSON)" -command OGFMovExport
    pack $f.bb.ex -side left -padx 4
    pack $f.bb -side top -pady 4
}

proc OGFMovBuildLCTab {f} {
    global ogfmov
    ttk::label $f.hdr -textvariable ogfmov(lchdr) -anchor w
    pack $f.hdr -side top -fill x -padx 6 -pady 4
    canvas $f.plot -width 640 -height 260 -background white
    pack $f.plot -side top -padx 6 -pady 4
    text $f.t -width 90 -height 10 -font TkFixedFont
    pack $f.t -side top -fill both -expand 1 -padx 6 -pady 4
    ttk::label $f.hint -text "Select a transient row in the table (Time-domain > Transients), then run Light Curve." -anchor w
    pack $f.hint -side bottom -fill x -padx 6
}

# light curve of the selected transient (mag vs MJD; epochs without a detection are upper limits -> flux plot)
proc OGFMovPlotLC {} {
    global ogfmov
    set f $ogfmov(detwin).nb.lc
    if {![winfo exists $f]} return
    $f.plot delete all
    $f.t delete 1.0 end
    if {[info exists ogfmov(lctext)]} {$f.t insert end $ogfmov(lctext)}
    set id $ogfmov(selected)
    set ogfmov(lchdr) "No transient selected"
    if {$ogfmov(selkind) ne "transient" || ![info exists ogfmov(lc)] || ![dict exists $ogfmov(lc) $id]} {
	if {$ogfmov(selkind) eq "transient"} {set ogfmov(lchdr) "T$id: no light curve yet (run Light Curve)"}
	return
    }
    set e [dict get $ogfmov(lc) $id]
    set t [dict get $e t]; set fl [dict get $e flux]; set er [dict get $e err]
    set ogfmov(lchdr) "T$id: [llength $t] epochs, forced photometry on the difference images (flux, e-/s)"
    set c $f.plot
    set W [$c cget -width]; set H [$c cget -height]
    set l 56; set r 12; set tp 14; set b 30
    set tmin [tcl::mathfunc::min {*}$t]; set tmax [tcl::mathfunc::max {*}$t]
    if {$tmax <= $tmin} {set tmax [expr {$tmin + 1e-3}]}
    set ymin 1e99; set ymax -1e99
    foreach y $fl s $er {
	if {![string is double -strict $y] || $y in {nan inf -inf}} continue
	set ymin [expr {min($ymin, $y-abs($s))}]; set ymax [expr {max($ymax, $y+abs($s))}]
    }
    if {$ymax <= $ymin} {set ymax [expr {$ymin + 1}]}
    set pad [expr {0.1*($ymax-$ymin)}]; set ymin [expr {$ymin-$pad}]; set ymax [expr {$ymax+$pad}]
    set x0 $l; set x1 [expr {$W-$r}]; set y0 [expr {$H-$b}]; set y1 $tp
    $c create rectangle $x0 $y1 $x1 $y0 -outline gray50
    $c create text [expr {$l-4}] $y1 -anchor e -text [format %.3g $ymax] -font TkSmallCaptionFont
    $c create text [expr {$l-4}] $y0 -anchor e -text [format %.3g $ymin] -font TkSmallCaptionFont
    $c create text [expr {($x0+$x1)/2}] [expr {$H-4}] -anchor s -font TkSmallCaptionFont -text [format "MJD %.4f .. %.4f" $tmin $tmax]
    if {$ymin < 0 && $ymax > 0} {
	set yz [expr {$y0 - (0-$ymin)/($ymax-$ymin)*($y0-$y1)}]
	$c create line $x0 $yz $x1 $yz -fill gray70 -dash {2 2}
    }
    foreach tt $t y $fl s $er {
	if {![string is double -strict $y] || $y in {nan inf -inf}} continue
	set x [expr {$x0 + ($tt-$tmin)/($tmax-$tmin)*($x1-$x0)}]
	set ya [expr {$y0 - ($y-$ymin)/($ymax-$ymin)*($y0-$y1)}]
	set yl [expr {$y0 - ($y-abs($s)-$ymin)/($ymax-$ymin)*($y0-$y1)}]
	set yh [expr {$y0 - ($y+abs($s)-$ymin)/($ymax-$ymin)*($y0-$y1)}]
	$c create line $x $yl $x $yh -fill "#c0392b"
	$c create oval [expr {$x-3}] [expr {$ya-3}] [expr {$x+3}] [expr {$ya+3}] -outline "#1f4e9c" -fill white
    }
}

# Orbit Fit... step: open the details window on the Orbit tab (tracklet of the selected row, or a designation)
proc OGFMovOrbitDialog {} {
    OGFMovDetails orbit
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
    OGFMovDetails orbit
    set w $ogfmov(detwin).nb.orbit
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
