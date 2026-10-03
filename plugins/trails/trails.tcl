# Satellite / aircraft trail removal (plugins/trails): Tcl side = hooks of the "Remove Trails" button (chip in the Detect tab).
# The work is done by plugins/trails/trails.py (CLI step, recorded by the session recorder); this file only prepares the mask path, shows the
# result (mask overlay on, catalogue columns, result panel with Undo) and implements Undo / Clear / the stacking dialog.
# The trail pixels live in the shared mask manager (bit 32 of ~/.ds9/mask_<base>.fits) so ICL / LSBG / stacking see them and Mask > Edit Mask > Undo works.

proc OGFTrailsWork {} {return [file join [OGFSessWorkDir] trails]}

proc OGFTrailsResult {} {
    set f [file join [OGFTrailsWork] trails.json]
    if {![file exists $f]} {return {}}
    set fd [open $f r]; set js [read $fd]; close $fd
    if {[catch {::ogf::json::parse $js} d]} {return {}}
    return $d
}

# step "before": tell the CLI which flag mask to update
proc OGFTrailsBefore {} {
    set p [OGFMaskPaths]
    if {$p eq {}} {error "no FITS image loaded"}
    ::ogf::cat::set trails,mask_path [dict get $p mask]
    ::ogf::cat::set trails,image [dict get $p fits]
    catch {file delete [file join [OGFTrailsWork] trails_flags.tsv]}
}

# step "after" of Remove Trails
proc OGFTrailsAfter {} {
    global ogfmask
    set d [OGFTrailsResult]
    if {$d eq {}} {::ogf::status "Remove Trails: no result file"; return}
    set n [::ogf::json::get $d n_trails 0]
    ::ogf::cat::set trails,n_trails $n
    ::ogf::cat::set trails,masked_fraction [::ogf::json::get $d masked_fraction 0]
    ::ogf::cat::set trails,result_file [file join [OGFTrailsWork] trails.json]
    set nfl [::ogf::json::get $d n_flagged 0]
    if {$n == 0} {
	::ogf::status "Remove Trails: no trail found (threshold [::ogf::json::get $d threshold ?], best score [format %.1f [::ogf::json::get $d best_zscore 0]]) - mask unchanged"
	OGFTrailsPanel 0 0 0
	return
    }
    # mask: sync ICL / LSBG, overlay on
    catch {OGFMaskAfterEdit}
    if {[::ogf::params::get trails show-overlay]} {
	set ogfmask(overlay) 1
	catch {OGFMaskRefreshOverlay}
    }
    # catalogue columns
    set ff [file join [OGFTrailsWork] trails_flags.tsv]
    if {[file exists $ff] && [::ogf::cat::has]} {
	set fd [open $ff r]; set txt [read $fd]; close $fd
	catch {::ogf::cat::add_columns $txt {TRAIL_FLAG TRAIL_ID TRAIL_DIST TRAIL_FLUX TRAIL_FRAC}}
	::ogf::cat::set trails,flags_file $ff
    }
    set fi [::ogf::json::get $d filled_image]
    if {$fi ne {} && [file exists $fi]} {
	::ogf::cat::set trails,filled_file $fi
	if {[::ogf::params::get trails show-filled]} {catch {OGFTrailsShowFilled}}
    }
    ::ogf::status "Remove Trails: $n trail(s) masked ([format %.2f [expr {100.0*[::ogf::json::get $d masked_fraction 0]}]] % of the image), $nfl catalogue object(s) flagged - Undo: Remove Trails menu or Mask > Edit Mask > Undo"
    if {[::ogf::params::get trails confirm-panel]} {OGFTrailsPanel $n [::ogf::json::get $d masked_fraction 0] $nfl}
}

# preview: draw the trails as line regions, no mask change
proc OGFTrailsPreviewAfter {} {
    global current
    set d [OGFTrailsResult]
    set n [::ogf::json::get $d n_trails 0]
    set reg [file join [OGFTrailsWork] trails_overlay.reg]
    if {$n > 0 && [file exists $reg]} {
	catch {$current(frame) marker load ds9 $reg}
    }
    ::ogf::status "Detect Trails: $n trail(s) found (preview only, mask unchanged)"
    set txt "Trails found: $n\n"
    foreach t [::ogf::json::get $d trails] {
	append txt [format "#%d  (%.1f,%.1f) -> (%.1f,%.1f)  angle %.2f deg  length %.0f px  FWHM %.1f px  mask half-width %.1f px  score %.1f\n" \
	    [::ogf::json::get $t id] [::ogf::json::get $t x1] [::ogf::json::get $t y1] [::ogf::json::get $t x2] [::ogf::json::get $t y2] \
	    [::ogf::json::get $t theta_deg] [::ogf::json::get $t length] [::ogf::json::get $t fwhm] [::ogf::json::get $t halfwidth] [::ogf::json::get $t zscore]]
    }
    OGFTextWindow "Detected trails (preview)" $txt
}

proc OGFTrailsPanel {n frac nfl} {
    set w .ogftrails
    catch {destroy $w}
    toplevel $w
    wm title $w "Remove Trails"
    wm transient $w .
    if {$n == 0} {
	set msg "No trail found.\nThe mask was not changed."
    } else {
	set msg [format "%d trail(s) masked (%.2f %% of the image).\n%d catalogue object(s) flagged (TRAIL_FLAG)." $n [expr {100.0*$frac}] $nfl]
    }
    ttk::label $w.l -text $msg -justify left
    pack $w.l -padx 12 -pady 8
    ttk::frame $w.b
    pack $w.b -pady {0 8}
    ttk::button $w.b.undo -text "Undo" -command OGFTrailsUndo
    ttk::button $w.b.set -text "Settings..." -command [list OGFParamDialog trails]
    ttk::button $w.b.close -text "Close" -command [list destroy $w]
    pack $w.b.undo $w.b.set $w.b.close -side left -padx 4
    if {$n == 0} {$w.b.undo state disabled}
}

# undo = the mask manager's undo (restores the mask before Remove Trails) + neutral catalogue columns
proc OGFTrailsUndo {} {
    if {[OGFMaskRun undo {} {}] eq {}} {return}
    catch {OGFMaskAfterEdit}
    if {[::ogf::cat::has]} {
	set nums [::ogf::cat::values NUMBER]
	set txt "NUMBER\tTRAIL_FLAG\tTRAIL_ID\tTRAIL_DIST\tTRAIL_FLUX\tTRAIL_FRAC\n"
	foreach nn $nums {append txt "$nn\t0\t0\t-1\t0\t-1\n"}
	catch {::ogf::cat::add_columns $txt {TRAIL_FLAG TRAIL_ID TRAIL_DIST TRAIL_FLUX TRAIL_FRAC}}
    }
    catch {destroy .ogftrails}
    ::ogf::status "Remove Trails undone (mask restored, TRAIL_* columns reset)"
}

proc OGFTrailsClear {} {
    if {[OGFMaskRun trails-clear {} {}] eq {}} {return}
    catch {OGFMaskAfterEdit}
    ::ogf::status "Trail bits cleared from the mask"
}

proc OGFTrailsShowFilled {} {
    global current
    set f [::ogf::cat::get trails,filled_file {}]
    if {$f eq {} || ![file exists $f]} {::ogf::status "Trails: no interpolated image yet (set Fill = interpolate)"; return 0}
    set orig $current(frame)
    CreateFrame
    if {[catch {LoadFitsFile $f {} {}} err]} {::ogf::log ERROR "trails: cannot load $f: $err"; catch {GotoFrame $orig}; return 0}
    catch {GotoFrame $orig}
    ::ogf::status "Trails: interpolated image opened in a new frame"
    return 1
}

# stack registered frames with their trails excluded (median / sigma-clipped mean)
proc OGFTrailsStackDialog {} {
    set files [tk_getOpenFile -title "Registered frames to stack (select several)" -multiple 1 -filetypes {{FITS {.fits .fit .fts .gz}} {All *}}]
    if {[llength $files] < 2} {::ogf::status "Trails stack: select at least two frames"; return}
    set r [OGFForm "Stack excluding trails" {{method "Method (sigclip|median|mean)" sigclip}}]
    if {$r eq {}} return
    set m [dict get $r method]
    if {$m ni {sigclip median mean}} {set m sigclip}
    set argv [list [OGFPython] [file join [::ogf::step::plugin_dir trails] trails.py] [lindex $files 0] --task stack --frames {*}$files \
	--work [file join [OGFSessWorkDir] trails_stack] --stack-method $m \
	--threshold [::ogf::params::get trails threshold] --max-trails [::ogf::params::get trails max-trails]]
    ::ogf::job::run $argv -step trails.stack -class manual -title "Stack excluding trails" -plugin trails \
	-done [list OGFTrailsStackDone]
}

proc OGFTrailsStackDone {ok output ms} {
    if {!$ok} return
    set p [file join [OGFSessWorkDir] trails_stack trails_stack.fits]
    if {[file exists $p]} {
	::ogf::cat::set trails,stack_file $p
	::ogf::status "Stack excluding trails written: $p"
	catch {CreateFrame; LoadFitsFile $p {} {}}
    }
}
