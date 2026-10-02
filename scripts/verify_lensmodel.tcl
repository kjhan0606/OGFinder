# GUI test of the lens-model plugin on a synthetic lens with known parameters (plugins/lensmodel/lensmodel.py --task simulate).
#   rm -rf ~/ds9.auto ~/ds9.auto.dir; HOME=/tmp/lm_home OGF_LM_OUT=/tmp/lm.txt OGF_LM_DIR=/tmp/lm_work OGF_LM_SIM=/tmp/lm_sim DISPLAY=:77 \
#       bin/ds9 /tmp/lm_sim/sim_lens.fits -geometry 1300x950 -source scripts/verify_lensmodel.tcl
# Part 1 (recorded, replayable): positions mode - fit, curves, predict, magnification map, source plane.  Session exported.
# Part 2 (not recorded for replay): the simulated catalog is loaded and the NUMBER mode is used (columns LENS_MU ..., selection helpers, dialog).
global catpanel current ds9
set ::fh [open $::env(OGF_LM_OUT) w]; set ::nf 0
set ::dir $::env(OGF_LM_DIR); file mkdir $::dir
set ::sim $::env(OGF_LM_SIM)
proc R {tag ok {d {}}} {puts $::fh "[expr {$ok?{PASS}:{FAIL}}] $tag $d"; flush $::fh; if {!$ok} {incr ::nf}}
proc bgerror {m} {puts $::fh "BGERROR $m"; flush $::fh; incr ::nf}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}
proc geom {} {
    update idletasks; update
    set tf [winfo parent $::catpanel(tbl)]
    return "[winfo rooty $tf] [winfo height $tf] [winfo height $::catpanel(infoarea)] [winfo width .]x[winfo height .]"
}
proc wait_job {{ms 300000}} {set t [clock milliseconds]; while {[::ogf::job::busy] && [clock milliseconds]-$t < $ms} {update; after 50}; wait_idle 300}
proc steps_since {n} {lmap r [lrange [::ogf::session::steps] $n end] {list [dict get $r step] [dict get $r class]}}
proc run_step {plugin step} {
    set s0 [llength [::ogf::session::steps]]
    set ok [::ogf::step::run $plugin $step]
    wait_job
    return [list $ok [steps_since $s0]]
}
proc readjson {f} {set fd [open $f r]; set t [read $fd]; close $fd; return [::ogf::json::parse $t]}
proc nregions {} {
    return [expr {[llength [split [$::current(frame) marker list ds9 image fk5 degrees yes] ";"]] - 2}]
}

proc run {} {
    global catpanel current ds9
    set truth [readjson [file join $::sim sim_truth.json]]
    set tp [dict get $truth params]
    R geometry_before [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    CatalogPanelExtract
    set t0 [clock milliseconds]; while {![::ogf::cat::has] && [clock milliseconds]-$t0 < 90000} {update; after 100}
    wait_idle 500
    R extracted [expr {[::ogf::cat::nrows] >= 1}] "rows=[::ogf::cat::nrows]"
    # ---- part 1: positions mode (replayable)
    set pos {}
    foreach im [dict get $truth images] {}
    set cat [open [file join $::sim sim_catalog.tsv] r]; set lines [split [string trim [read $cat]] \n]; close $cat
    set pl {}
    foreach l [lrange $lines 2 end] {set f [split $l \t]; lappend pl "[lindex $f 1],[lindex $f 2]"}
    ::ogf::params::put lensmodel positions [join $pl ";"]
    set lc [split [lindex $lines 1] \t]
    ::ogf::params::put lensmodel lens-x [lindex $lc 1]
    ::ogf::params::put lensmodel lens-y [lindex $lc 2]
    ::ogf::params::put lensmodel sigma-pos 0.005
    ::ogf::params::put lensmodel z-lens 0.5
    ::ogf::params::put lensmodel z-source 2.0
    ::ogf::params::put lensmodel show-overlay 1
    set f0 [llength $ds9(frames)]
    lassign [run_step lensmodel fit] ok recs
    R lens_fit_ran $ok $recs
    R lens_fit_recorded [expr {[lindex $recs 0 0] eq "analysis.lens_fit"}] $recs
    set w [file join [OGFSessWorkDir] lensmodel]
    foreach f {lens_model.json lens_images.tsv lens_curves.json lens_overlay.reg} {
	R lens_file_$f [expr {[file exists [file join $w $f]] && [file size [file join $w $f]] > 50}]
    }
    set res [readjson [file join $w lens_model.json]]
    set p [dict get $res params]
    set dth [expr {abs([dict get $p theta_E] - [dict get $tp theta_E])/[dict get $tp theta_E]}]
    R lens_theta_E_recovered [expr {$dth < 2e-3}] "theta_E [dict get $p theta_E] vs [dict get $tp theta_E] (rel $dth)"
    R lens_q_recovered [expr {abs([dict get $p q] - [dict get $tp q]) < 0.03}] "q [dict get $p q] vs [dict get $tp q]"
    R lens_chi2_ok [expr {[dict get $res dof] == 1 && [dict get $res rms_arcsec] < 0.02}] "dof [dict get $res dof] rms [dict get $res rms_arcsec]"
    R lens_four_images [expr {[dict get $res n_images_predicted] == 4}] [dict get $res n_images_predicted]
    R lens_overlay_loaded [expr {[nregions] == 7}] "regions [nregions]"
    OGFLensOverlay
    R lens_overlay_replaced_not_duplicated [expr {[nregions] == 7}] "regions [nregions]"
    lassign [run_step lensmodel curves] ok recs
    R lens_curves_ran [expr {$ok && [lindex $recs 0 0] eq "analysis.lens_curves"}] $recs
    ::ogf::params::put lensmodel z-source 3.0
    lassign [run_step lensmodel predict] ok recs
    R lens_predict_ran [expr {$ok && [lindex $recs 0 0] eq "analysis.lens_predict"}] $recs
    set res2 [readjson [file join $w lens_model.json]]
    R lens_predict_new_z [expr {[dict get $res2 z_source] == 3.0 && [dict get [dict get $res2 physical] mass_Einstein_Msun] != [dict get [dict get $res physical] mass_Einstein_Msun]}]
    ::ogf::params::put lensmodel z-source 2.0
    lassign [run_step lensmodel predict] ok recs
    set nfr [llength $ds9(frames)]
    lassign [run_step lensmodel magmap] ok recs
    R lens_magmap_ran [expr {$ok && [lindex $recs 0 0] eq "analysis.lens_magmap"}] $recs
    R lens_magmap_file [file exists [file join $w lens_magnification.fits]]
    R lens_magmap_frame [expr {[llength $ds9(frames)] == $nfr + 1}] "frames [llength $ds9(frames)]"
    GotoFrame [lindex $ds9(frames) 0]
    ::ogf::params::put lensmodel lens-mask 0.0
    lassign [run_step lensmodel source] ok recs
    R lens_source_ran [expr {$ok && [lindex $recs 0 0] eq "analysis.lens_source"}] $recs
    R lens_source_file [file exists [file join $w lens_source.fits]]
    R lens_source_frame [expr {[llength $ds9(frames)] == $nfr + 2}] "frames [llength $ds9(frames)]"
    set res3 [readjson [file join $w lens_model.json]]
    set sr [dict get $res3 source_reconstruction]
    set sx [lindex [dict get $sr source_centroid_arcsec] 0]; set sy [lindex [dict get $sr source_centroid_arcsec] 1]
    set bx [lindex [dict get $truth source_arcsec] 0]; set by [lindex [dict get $truth source_arcsec] 1]
    # the four point-source images (flux 100 |mu|) dominate the back-projected flux: the centroid must be the true source position
    R lens_source_centroid [expr {hypot($sx-$bx, $sy-$by) < 0.03}] "centroid ($sx,$sy) truth ($bx,$by)"
    GotoFrame [lindex $ds9(frames) 0]
    set sp [file join $::dir ogfinder_session.py]
    CatalogPanelSessionSave $sp
    R session_exported [file exists $sp]
    # ---- part 2: NUMBER mode with the simulated catalog (not replayable)
    set fd [open [file join $::sim sim_catalog.tsv] r]; set data [read $fd]; close $fd
    ::ogf::cat::load_tsv $data sim
    wait_idle 800
    R sim_catalog_loaded [expr {[::ogf::cat::nrows] == 5}] "rows=[::ogf::cat::nrows]"
    ::ogf::cat::select {2 3 4 5}
    wait_idle 300
    R pick_images [expr {[OGFLensPickImages] && [::ogf::params::get lensmodel image-numbers] eq "2,3,4,5"}] [::ogf::params::get lensmodel image-numbers]
    ::ogf::cat::select {1}
    wait_idle 300
    R pick_lens [expr {[OGFLensPickLens] && [::ogf::params::get lensmodel lens-number] == 1}] [::ogf::params::get lensmodel lens-number]
    ::ogf::params::put lensmodel positions {}
    ::ogf::params::put lensmodel flux-column FLUX
    ::ogf::params::put lensmodel flux-sigma 0.05
    lassign [run_step lensmodel fit] ok recs
    R number_fit_ran $ok $recs
    set cols [::ogf::cat::columns]
    R number_columns [expr {"LENS_MU" in $cols && "LENS_RES" in $cols && "LENS_DT" in $cols && "LENS_PARITY" in $cols}]
    set mus [lsearch -all -inline -not [::ogf::cat::values LENS_MU] {}]
    R number_rows_filled [expr {[llength $mus] == 4}] "mu values [llength $mus]"
    set truemu {}
    foreach im [dict get $truth images] {lappend truemu [dict get $im mu]}
    set okmu 1
    foreach m $mus {set best 1e9; foreach t $truemu {set best [expr {min($best, abs($m-$t))}]}; if {$best > 0.05} {set okmu 0}}
    R number_mu_match_truth $okmu "mu $mus truth $truemu"
    set w2 [OGFLensDialog]
    wait_idle 300
    R dialog_open [winfo exists .ogflens]
    R dialog_images_table [expr {[llength [.ogflens.nb.img.t children {}]] == 4}] [llength [.ogflens.nb.img.t children {}]]
    R dialog_curves_table [expr {[llength [.ogflens.nb.dly.t children {}]] >= 1}]
    R dialog_summary_text [expr {[string match "*chi2*" [.ogflens.nb.sum.t get 1.0 end]]}]
    destroy .ogflens
    R final_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    puts $::fh "SUMMARY failures=$::nf"; close $::fh
}
after 3000 {
    if {[catch run err]} {puts $::fh "FAIL run_error $err $::errorInfo"; incr ::nf; puts $::fh "SUMMARY failures=$::nf"; close $::fh}
    exit 0
}
