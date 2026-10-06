# GUI test of plugins/ds10core (stand-alone shell of the shared Astrafex Core).  Run through scripts/verify_ds10core.sh.
global catpanel current ds9 ogfui
set ::fh [open $::env(OGF_DS_OUT) w]; set ::nf 0
set ::dir $::env(OGF_DS_DIR)
set ::shots $::env(OGF_DS_SHOTS)
proc R {tag ok {d {}}} {puts $::fh "[expr {$ok?{PASS}:{FAIL}}] $tag $d"; flush $::fh; if {!$ok} {incr ::nf}}
proc bgerror {m} {puts $::fh "BGERROR $m"; flush $::fh; incr ::nf}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}
proc wait_job {{ms 300000}} {set t [clock milliseconds]; while {[::ogf::job::busy] && [clock milliseconds]-$t < $ms} {update; after 50}; wait_idle 500}
proc steps_since {n} {lmap r [lrange [::ogf::session::steps] $n end] {list [dict get $r step] [dict get $r class]}}
proc run_step {plugin step} {
    set s0 [llength [::ogf::session::steps]]
    set ok [::ogf::step::run $plugin $step]
    wait_job
    return [list $ok [steps_since $s0]]
}
proc text_window {} {if {[winfo exists .ogftext.t]} {return [.ogftext.t get 1.0 end]}; return {}}
proc shot {name} {
    if {$::shots eq {}} return
    update idletasks; update; wait_idle 400
    file mkdir $::shots
    catch {exec $::env(OGFINDER_PYTHON) scripts/xshot.py $::env(OGF_DS_DISPLAY) [file join $::shots $name.png]} m
    puts $::fh "SHOT $name $m"; flush $::fh
}

proc run {} {
    global catpanel current ds9 ogfui
    wait_idle 800
    # ---- registration (plugin manifest, chip, menu, licence tags are in the manifest)
    R registered [expr {[::ogf::reg::get ds10core] ne {}}]
    set m [::ogf::reg::get ds10core]
    R licence_tag [expr {[dict exists $m license] || [string match *license* $m]}]
    R chip_run_button [expr {[info exists ogfui(run,ds10core)] && [winfo exists $ogfui(run,ds10core)]}]
    OGFUIShowTab Measure
    wait_idle 400
    # the Measure tab has more chips than fit the panel: this one may sit in the "More" menu (ogfui(overflow,Measure)); it is reachable either way
    set vis [expr {[info exists ogfui(run,ds10core)] && [winfo ismapped $ogfui(run,ds10core)]}]
    set inmore [expr {[info exists ogfui(overflow,Measure)] && "ds10core" in $ogfui(overflow,Measure)}]
    R chip_reachable_in_measure_tab [expr {$vis || $inmore}] "visible=$vis in_more_menu=$inmore"
    foreach lbl {"Forced photometry (bands -> colours, optional PSF match)" "Region statistics (exact pixel membership)" "Pixel table at (X, Y)" "Run astrafex-script / replay.py / bundle"} {
	R menu_has_[string map {{ } _ ( {} ) {} , {} / _ > _ - _} $lbl] [expr {[info exists ogfui(menu,ds10core)] && [$ogfui(menu,ds10core) index $lbl] ne "none"}]
    }
    R param_defaults [expr {[::ogf::params::get ds10core aperture-radii] eq "0.18,0.3,0.48" && [::ogf::params::get ds10core cog-radius] == 0.72}]
    R core_found [expr {![catch {exec $::env(OGFINDER_PYTHON) plugins/ds10core/ds10.py --where} where]}] $where
    shot astrafex_01_chip_menu_measure_tab
    # the plugin's own menu (all steps) posted next to the panel
    set mm [expr {$inmore ? "$ogfui(more,Measure).m.p_ds10core" : $ogfui(menu,ds10core)}]
    catch {$mm post [expr {[winfo rootx .] + 760}] [expr {[winfo rooty .] + 120}]}
    shot astrafex_01b_plugin_menu
    catch {$mm unpost}
    # ---- header calibration of the open image (primary action of the chip)
    lassign [run_step ds10core calib] ok recs
    set txt [text_window]
    R calib_step_ok [expr {$ok}] $recs
    R calib_text [expr {[string match "*AB zero point: 25.9463*" $txt] && [string match "*F160W*" $txt] && [string match "*origin: header*" $txt]}] [string range $txt 0 200]
    R calib_recorded [expr {[lindex $recs 0 0] eq "ds10core.calib"}] $recs
    shot astrafex_02_calibration_result
    # ---- maps (nothing to find next to the synthetic image: the step must run and say so)
    lassign [run_step ds10core maps] ok recs
    set txt [text_window]
    R maps_step [expr {$ok && [string match "*maps*" $txt]}] [string range $txt 0 120]
    # ---- regions
    ::ogf::params::put ds10core regions-file [file join $::dir r.reg]
    lassign [run_step ds10core regions-stats] ok recs
    set txt [text_window]
    R regions_stats_text [expr {$ok && [string match "*circle*" $txt] && [string match "*box*" $txt] && [string match "*annulus*" $txt] && [string match "*n_pix*" $txt]}] [string range $txt 0 300]
    shot astrafex_03_region_statistics
    lassign [run_step ds10core regions-mask] ok recs
    set mf [file join [OGFSessWorkDir] ds10core regions_mask.fits]
    R regions_mask_file [expr {$ok && [file exists $mf] && [file size $mf] > 1000}] $mf
    # ---- pixel table
    ::ogf::params::put ds10core px 150
    ::ogf::params::put ds10core py 148
    ::ogf::params::put ds10core size 9
    lassign [run_step ds10core pixtab] ok recs
    set txt [text_window]
    R pixtab_text [expr {$ok && [llength [split [string trim $txt] "\n"]] >= 9}] [string range $txt 0 120]
    shot astrafex_04_pixel_table
    # ---- forced photometry on the three synthetic bands
    ::ogf::params::put ds10core band-images "[file join $::dir syn_f105w.fits] [file join $::dir syn_f125w.fits] [file join $::dir syn_f160w.fits]"
    ::ogf::params::put ds10core thresh 3.0
    ::ogf::params::put ds10core smooth-fwhm 2.0
    ::ogf::params::put ds10core minarea 4
    lassign [run_step ds10core forced] ok recs
    set txt [text_window]
    set wd [file join [OGFSessWorkDir] ds10core forced_out]
    R forced_ok [expr {$ok && [file exists [file join $wd forced_catalog.tsv]] && [file exists [file join $wd script.json]]}] $recs
    R forced_text [expr {[string match "*forced photometry:*" $txt] && [string match "*F105W*" $txt] && [string match "*header*" $txt]}] [string range $txt 0 300]
    set fd [open [file join $wd forced_catalog.tsv] r]; set hdr [gets $fd]; set nrow 0; while {[gets $fd line] >= 0} {incr nrow}; close $fd
    R forced_catalog_columns [expr {[string match "*COLOR_F105W_F160W*" $hdr] && [string match "*MAG_AUTO_F125W*" $hdr] && $nrow >= 5}] "rows=$nrow"
    shot astrafex_05_forced_photometry_result
    # ---- forced photometry again with PSF matching for colours
    ::ogf::params::put ds10core psf-match 1
    ::ogf::params::put ds10core psf-method fourier
    ::ogf::params::put ds10core psf-target broadest
    lassign [run_step ds10core forced] ok recs
    set txt [text_window]
    set wd2 [file join [OGFSessWorkDir] ds10core forced_out]
    R forced_psf_ok [expr {$ok && [file exists [file join $wd2 forced_summary.json]]}] $recs
    set summ {}
    if {[file exists [file join $wd2 forced_summary.json]]} {
	set fd [open [file join $wd2 forced_summary.json] r]; set summ [read $fd]; close $fd
    }
    R forced_psf_summary [expr {[regexp {"enabled": true} $summ] || [regexp {"enabled":true} $summ]}] [string range $summ 0 200]
    shot psf_forced_match_result
    ::ogf::params::put ds10core psf-match 0
    # ---- shared-core LSBG finder with the SVM classification stage (keep-all: every candidate gets P_LSBG / LSBG_CLASS)
    foreach {k v} {lsbg-mu-eff-min 10 lsbg-r-eff-min 0.05 lsbg-detect-minarea 10 lsbg-classify 0 lsbg-classify-keep-all 1 lsbg-classify-threshold 0} {::ogf::params::put ds10core $k $v}
    lassign [run_step ds10core lsbg] ok recs
    set lcat [lindex [glob -nocomplain [file join [OGFSessWorkDir] lsbg_catalog.tsv] [file join [OGFSessWorkDir] ds10core lsbg_catalog.tsv] [file join [OGFSessWorkDir] ds10core * lsbg_catalog.tsv]] 0]
    set lsum [file join [file dirname $lcat] lsbg_summary.json]
    set lhdr {}; set lsumtxt {}
    if {$lcat ne {} && [file exists $lcat]} {set fd [open $lcat r]; set lhdr [gets $fd]; close $fd}
    if {[file exists $lsum]} {set fd [open $lsum r]; set lsumtxt [read $fd]; close $fd}
    R lsbg_classify_ok [expr {$ok && [string match "*P_LSBG*LSBG_CLASS*" $lhdr]}] "$recs [string range $lhdr end-60 end]"
    R lsbg_classify_summary [expr {[string match "*\"classifier\"*" $lsumtxt] && [string match "*lsbg*" [lindex $recs 0 0]]}] [string range $lsumtxt 0 300]
    catch {destroy .ogftext}; wait_idle 1500
    R lsbg_classify_table [expr {[::ogf::cat::nrows] > 0 && [lsearch [::ogf::cat::columns] P_LSBG] >= 0}] "rows=[::ogf::cat::nrows]"
    catch {[::ogf::cat::get tbl] xview moveto 1.0}; wait_idle 300
    shot lsbg_classify_desktop_1440x900
    set dlg [OGFParamDialog ds10core -group LSBG -title {LSBG finder (shared core) - Settings}]
    R lsbg_classify_dialog [expr {$dlg ne {} && [winfo exists $dlg]}] $dlg
    if {$dlg ne {}} {catch {foreach t [$dlg.body.nb tabs] {if {[$dlg.body.nb tab $t -text] eq {LSBG}} {$dlg.body.nb select $t}}}; catch {wm geometry $dlg 640x700+790+60}; wait_idle 800; shot lsbg_classify_desktop_params_1440x900; catch {destroy $dlg}}
    ::ogf::params::put ds10core lsbg-classify-keep-all 0
    # ---- run the generated step list (web-compatible astrafex-script/1) as a script
    ::ogf::params::put ds10core script [file join $wd script.json]
    ::ogf::params::put ds10core files-dir $::dir
    lassign [run_step ds10core run-script] ok recs
    set rr [file join [OGFSessWorkDir] ds10core run_out run_report.json]
    R run_script_ok [expr {$ok && [file exists $rr]}] $recs
    if {[file exists $rr]} {
	set fd [open $rr r]; set rep [read $fd]; close $fd
	R run_report_ok [expr {[regexp {"ok": true} $rep]}]
    }
    # ---- the session recorded the steps with their exact argument vectors
    set kinds [lmap r [::ogf::session::steps] {dict get $r step}]
    R session_records [expr {"ds10core.calib" in $kinds && "ds10core.forced" in $kinds && "ds10core.run_script" in $kinds}] $kinds
    # ---- open an offline bundle of the web app: restore + show the first image and its regions in ds9
    set bundle [file join [pwd] plugins ds10core tests data offline_bundle_small.zip]
    ::ogf::params::put ds10core bundle $bundle
    ::ogf::params::put ds10core open-dest [file join $::dir offline ws]
    lassign [run_step ds10core open-bundle] ok recs
    set txt [text_window]
    set wsd [file join $::dir offline ws]
    R open_bundle_step [expr {$ok && [string match "*opened offline bundle*" $txt] && [string match "*field96.fits*" $txt] && [string match "*steps with stored results: s1*" $txt]}] [string range $txt 0 200]
    R open_bundle_workspace [expr {[file exists [file join $wsd workspace.json]] && [file exists [file join $wsd files field96.fits]] && [file exists [file join $wsd regions regions.reg]]}]
    wait_idle 800
    set shown [catch {$current(frame) get fits file name root base} fnm]
    R open_bundle_image_in_ds9 [expr {$shown == 0 && [string match "*field96*" $fnm]}] "frame file: $fnm"
    set nmark -1
    catch {$current(frame) marker select all; set nmark [$current(frame) get marker select number]; $current(frame) marker unselect all}
    R open_bundle_regions_in_ds9 [expr {$nmark == 1}] "markers=$nmark"
    R open_bundle_workspace_recorded [expr {[::ogf::cat::get ds10core,workspace ?] eq $wsd}]
    shot astrafex_offline_open_bundle_ds9
    # the restored script runs in the workspace and reuses the stored result
    ::ogf::params::put ds10core script [file join $wsd session_script.json]
    ::ogf::params::put ds10core files-dir [file join $wsd files]
    lassign [run_step ds10core run-script] ok recs
    set rr [file join [OGFSessWorkDir] ds10core run_out run_report.json]
    set rep {}; if {[file exists $rr]} {set fd [open $rr r]; set rep [read $fd]; close $fd}
    R open_bundle_script_restored [expr {$ok && [regexp {"status": "restored"} $rep]}] [string range $rep 0 200]
    shot astrafex_offline_open_bundle_script_run
    # the way back: pack the workspace into a return bundle for the web import (nothing was run locally here, so only the web step is listed)
    ::ogf::params::put ds10core pack-workspace $wsd
    ::ogf::params::put ds10core pack-out [file join $::dir offline back.zip]
    lassign [run_step ds10core pack-return] ok recs
    set txt [text_window]
    R pack_return_step [expr {$ok && [file exists [file join $::dir offline back.zip]] && [string match "*s1: web*" $txt]}] [string range $txt 0 200]
    puts $::fh "SUMMARY failures=$::nf"; close $::fh
    exit
}
after 3000 {if {[catch run err]} {puts $::fh "BGERROR $err"; puts $::fh "SUMMARY failures=99"; close $::fh; exit}}
