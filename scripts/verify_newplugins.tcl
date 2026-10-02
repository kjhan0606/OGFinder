# GUI test of the analysis plugins added after the plugin restructuring (isophote, completeness, ...).  Real ds9, real m51 extraction,
# every step run through ::ogf::step::run (real job runner + recorder), outputs checked, session exported and replayed headless.
#   rm -rf ~/ds9.auto ~/ds9.auto.dir; HOME=/tmp/np_home OGF_NP_OUT=/tmp/np.txt OGF_NP_DIR=/tmp/np_work DISPLAY=:77 \
#       bin/ds9 /workspace/fits/m51.fits -geometry 1300x950 -source scripts/verify_newplugins.tcl
# OGF_NP_ONLY=isophote,completeness limits the sections.  Output: PASS/FAIL lines + SUMMARY failures=N.
global catpanel current ds9
set ::fh [open $::env(OGF_NP_OUT) w]; set ::nf 0
set ::dir $::env(OGF_NP_DIR); file mkdir $::dir
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
proc col_values {name} {return [::ogf::cat::values $name]}
proc nonempty {name} {set n 0; foreach v [col_values $name] {if {$v ne {}} {incr n}}; return $n}
proc run_step {plugin step} {
    set s0 [llength [::ogf::session::steps]]
    set ok [::ogf::step::run $plugin $step]
    wait_job
    return [list $ok [steps_since $s0]]
}
proc want {only name} {return [expr {$only eq {} || $name in $only}]}

proc sec_isophote {} {
    set f0 [llength $::ds9(frames)]
    ::ogf::params::put isophote max-objects 3
    ::ogf::params::put isophote step 0.15
    lassign [run_step isophote fit] ok recs
    R isophote_ran $ok $recs
    R isophote_recorded [expr {[lindex $recs 0 0] eq "analysis.isophote"}] $recs
    R isophote_columns [expr {[::ogf::cat::columns] ne {} && "ISO_EPS_HL" in [::ogf::cat::columns] && "ISO_SHAPE" in [::ogf::cat::columns]}]
    set n [nonempty ISO_NISO]
    R isophote_rows_filled [expr {$n >= 1 && $n <= 3}] "rows=$n"
    set w [OGFSessWorkDir]
    foreach f {isophote_model.fits isophote_resid.fits isophote_profiles.tsv isophote_plot.png isophote_profiles.json} {
	R isophote_file_$f [expr {[file exists [file join $w $f]] && [file size [file join $w $f]] > 100}]
    }
    R isophote_cat_key [expr {[file exists [::ogf::cat::get isophote,model_file {}]]}]
    R isophote_frames [expr {[llength $::ds9(frames)] == $f0 + 2}] "frames [llength $::ds9(frames)] (was $f0)"
    R isophote_frame_restored [expr {$::current(frame) eq [lindex $::ds9(frames) 0]}] $::current(frame)
    R isophote_argv_templated [expr {[string match {*@{WORK}/isophote_model.fits*} [dict get [lindex [::ogf::session::steps] end] argv_t]]}]
    set pw [OGFIsophotePlot]
    R isophote_plot_window [expr {[winfo exists $pw] && [image width ogfisoimg] > 300 && [image height ogfisoimg] > 200}] "[image width ogfisoimg]x[image height ogfisoimg]"
    destroy $pw
    OGFIsophoteTable; update
    R isophote_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc sec_completeness {} {
    ::ogf::params::put completeness n-bins 5
    ::ogf::params::put completeness per-bin 12
    ::ogf::params::put completeness per-image 6
    ::ogf::params::put completeness mag-min 12.5
    ::ogf::params::put completeness mag-max 17.5
    ::ogf::params::put completeness crop-size 400
    lassign [run_step completeness measure] ok recs
    R completeness_ran $ok $recs
    R completeness_recorded [expr {[lindex $recs 0 0] eq "analysis.completeness_sim"}] $recs
    set cols [::ogf::cat::columns]
    R completeness_columns [expr {"COMPL_FRAC" in $cols && "COMPL_LIM50" in $cols && "COMPL_LIM90" in $cols}]
    R completeness_rows_filled [expr {[nonempty COMPL_FRAC] > 50}] "rows=[nonempty COMPL_FRAC]"
    set fr [lsearch -all -inline -not [col_values COMPL_FRAC] {}]
    R completeness_fraction_range [expr {[tcl::mathfunc::min {*}$fr] >= 0 && [tcl::mathfunc::max {*}$fr] <= 1}] "[tcl::mathfunc::min {*}$fr] .. [tcl::mathfunc::max {*}$fr]"
    set w [OGFSessWorkDir]
    foreach f {completeness.json completeness_curve.tsv completeness_plot.png catalog_meta.json} {
	R completeness_file_$f [expr {[file exists [file join $w $f]] && [file size [file join $w $f]] > 50}]
    }
    R completeness_keys [expr {[::ogf::cat::exists completeness,lim50] && [string is double -strict [::ogf::cat::get completeness,lim50]]}] [::ogf::cat::get completeness,lim50 ?]
    set pw [OGFCompletenessPlot]
    R completeness_plot_window [expr {[winfo exists $pw] && [image width ogfcompimg] > 300}] "[image width ogfcompimg]x[image height ogfcompimg]"
    destroy $pw
    # catalog metadata travels with a saved catalog
    set sv [file join $::dir saved.tsv]
    CatalogPanelSaveCatalogTo $sv
    R completeness_meta_sidecar [file exists $sv.meta.json] $sv.meta.json
    R completeness_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc sec_psfex {} {
    ::ogf::params::put psfex snr-min 8
    ::ogf::params::put psfex fwhm 2.2
    ::ogf::params::put psfex order 1
    lassign [run_step psfex build] ok recs
    R psfex_ran $ok $recs
    R psfex_recorded [expr {[lindex $recs 0 0] eq "analysis.psfex"}] $recs
    set cols [::ogf::cat::columns]
    R psfex_columns [expr {"PSFM_FWHM" in $cols && "PSFM_E" in $cols && "PSFM_PA" in $cols && "PSFM_NSTAR" in $cols}]
    R psfex_rows_filled [expr {[nonempty PSFM_FWHM] > 100}] "rows=[nonempty PSFM_FWHM]"
    set fw [lsearch -all -inline -not [col_values PSFM_FWHM] {}]
    R psfex_fwhm_range [expr {[tcl::mathfunc::min {*}$fw] > 1.0 && [tcl::mathfunc::max {*}$fw] < 8.0}] "[tcl::mathfunc::min {*}$fw] .. [tcl::mathfunc::max {*}$fw]"
    set w [file join [OGFSessWorkDir] psfex]
    foreach f {model.json model.fits center.fits stars.tsv maps.tsv info.json diag.png} {
	R psfex_file_$f [expr {[file exists [file join $w psfex_$f]] && [file size [file join $w psfex_$f]] > 50}]
    }
    R psfex_keys [expr {[::ogf::cat::exists psfex,model_file] && [file exists [::ogf::cat::get psfex,model_file]]}] [::ogf::cat::get psfex,model_file ?]
    set pw [OGFPsfexDiag]
    R psfex_plot_window [expr {[winfo exists $pw] && [image width ogfpsfeximg] > 300}] "[image width ogfpsfeximg]x[image height ogfpsfeximg]"
    destroy $pw
    run_step psfex use
    R psfex_use [expr {[::ogf::cat::get psf,model ?] eq [::ogf::cat::get psfex,model_file ?] && [file exists [::ogf::cat::get psf,file ?]]}] [::ogf::cat::get psf,model ?]
    # PSF photometry now carries --psf-model (argv unchanged when unset)
    set argv [::ogf::step::build_argv photometry [::ogf::reg::step photometry psf_phot] [::ogf::step::context photometry psfphot]]
    R psfex_photometry_argv [expr {[lsearch $argv --psf-model] >= 0}] [lrange $argv end-1 end]
    ::ogf::cat::set psf,model {}
    set argv [::ogf::step::build_argv photometry [::ogf::reg::step photometry psf_phot] [::ogf::step::context photometry psfphot]]
    R psfex_photometry_argv_unset [expr {[lsearch $argv --psf-model] < 0}]
    R psfex_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc sec_multifit {} {
    set f0 [llength $::ds9(frames)]
    ::ogf::params::put multifit max-objects 10
    ::ogf::params::put multifit model sersic
    ::ogf::params::put multifit neighbours fit
    ::ogf::params::put multifit psf-fwhm 3.0
    lassign [run_step multifit fit] ok recs
    R multifit_ran $ok $recs
    R multifit_recorded [expr {[lindex $recs 0 0] eq "analysis.multifit"}] $recs
    set cols [::ogf::cat::columns]
    R multifit_columns [expr {"GF_MAG" in $cols && "GF_RE" in $cols && "GF_N" in $cols && "GF_CHI2" in $cols && "GF_FLAG" in $cols && "GF_BT" in $cols}]
    R multifit_rows_filled [expr {[nonempty GF_MAG] >= 8 && [nonempty GF_MAG] <= 10}] "rows=[nonempty GF_MAG]"
    set w [file join [OGFSessWorkDir] multifit]
    foreach f {results.tsv model.fits residual.fits montage.png psf.json} {
	R multifit_file_$f [expr {[file exists [file join $w multifit_$f]] && [file size [file join $w multifit_$f]] > 50}]
    }
    R multifit_frames [expr {[llength $::ds9(frames)] == $f0 + 2}] "frames [llength $::ds9(frames)] (was $f0)"
    R multifit_frame_restored [expr {$::current(frame) eq [lindex $::ds9(frames) 0]}] $::current(frame)
    R multifit_keys [expr {[::ogf::cat::exists multifit,residual_file] && [file exists [::ogf::cat::get multifit,residual_file]]}]
    set pw [OGFMultifitMontage]
    R multifit_montage_window [expr {[winfo exists $pw] && [image width ogfmultifitimg] > 300}] "[image width ogfmultifitimg]x[image height ogfmultifitimg]"
    destroy $pw
    # selected rows only
    set nums [lrange [::ogf::cat::values NUMBER] 20 22]
    ::ogf::cat::select $nums replace 0
    update
    set s0 [llength [::ogf::session::steps]]
    OGFMultifitSelected
    wait_job
    R multifit_selected_only [expr {[nonempty GF_MAG] == 3}] "rows=[nonempty GF_MAG]"
    R multifit_selected_param [expr {[::ogf::params::get multifit objects] eq [join $nums \;]}] [::ogf::params::get multifit objects]
    ::ogf::params::put multifit objects {}
    R multifit_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc sec_morphext {} {
    ::ogf::params::put morph_ext mx-max-sources 60
    ::ogf::params::put morph_ext mx-curves 3
    lassign [run_step morph_ext morph_ext] ok recs
    R morphext_ran $ok $recs
    R morphext_recorded [expr {[lindex $recs 0 0] eq "analysis.morph_ext"}] $recs
    set cols [::ogf::cat::columns]
    R morphext_columns [expr {"MX_RP" in $cols && "MX_KRON_MAG" in $cols && "MX_SMOOTH" in $cols && "MX_GINI_P" in $cols && "MX_CONC" in $cols}]
    R morphext_rows_filled [expr {[nonempty MX_RP] >= 10 && [nonempty MX_KRON_R] >= 30}] "rp=[nonempty MX_RP] kron=[nonempty MX_KRON_R]"
    set rp [lsearch -all -inline -not [col_values MX_RP] {}]
    R morphext_rp_sane [expr {[tcl::mathfunc::min {*}$rp] > 1.0 && [tcl::mathfunc::max {*}$rp] < 200}] "[tcl::mathfunc::min {*}$rp] .. [tcl::mathfunc::max {*}$rp]"
    set w [file join [OGFSessWorkDir] morph_ext]
    foreach f {morph_ext_growth.tsv morph_ext_curves.png} {
	R morphext_file_$f [expr {[file exists [file join $w $f]] && [file size [file join $w $f]] > 200}]
    }
    R morphext_keys [expr {[::ogf::cat::exists morphext,growth_file]}]
    set pw [OGFMorphExtCurves]
    R morphext_plot_window [expr {[winfo exists $pw] && [image width ogfmorphextimg] > 300}] "[image width ogfmorphextimg]x[image height ogfmorphextimg]"
    destroy $pw
    R morphext_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc sec_noisemodel {} {
    ::ogf::params::put noisemodel n-aper 300
    ::ogf::params::put noisemodel radii 2,4,6,8
    ::ogf::params::put noisemodel aperture fixed
    ::ogf::params::put noisemodel aper-radius 4
    ::ogf::params::put noisemodel correct 1
    lassign [run_step noisemodel model] ok recs
    R noisemodel_ran $ok $recs
    R noisemodel_recorded [expr {[lindex $recs 0 0] eq "analysis.noisemodel"}] $recs
    set cols [::ogf::cat::columns]
    R noisemodel_columns [expr {"NM_RMS" in $cols && "NM_FLUXERR" in $cols && "NM_CORR" in $cols && "NM_MAGERR_AP" in $cols}]
    R noisemodel_rows_filled [expr {[nonempty NM_FLUXERR] >= 20 && [nonempty NM_RMS] >= 20}] "err=[nonempty NM_FLUXERR] rms=[nonempty NM_RMS]"
    set corr [lsearch -all -inline -not [col_values NM_CORR] {}]
    R noisemodel_corr_sane [expr {[tcl::mathfunc::min {*}$corr] > 0.5 && [tcl::mathfunc::max {*}$corr] < 20}] "[tcl::mathfunc::min {*}$corr] .. [tcl::mathfunc::max {*}$corr]"
    set w [file join [OGFSessWorkDir] noisemodel]
    foreach f {noisemodel_bkg.fits noisemodel_rms.fits noisemodel_sub.fits noisemodel_summary.json noisemodel_curve.tsv noisemodel_noise_curve.png} {
	R noisemodel_file_$f [expr {[file exists [file join $w $f]] && [file size [file join $w $f]] > 100}]
    }
    R noisemodel_keys [expr {[::ogf::cat::exists noisemodel,rms_file] && [::ogf::cat::exists noisemodel,plot_file]}]
    set pw [OGFNoisePlot]
    R noisemodel_plot_window [expr {[winfo exists $pw] && [image width ogfnoiseimg] > 300}] "[image width ogfnoiseimg]x[image height ogfnoiseimg]"
    destroy $pw
    set f0 [llength $::ds9(frames)]
    set n [OGFNoiseFrames]
    R noisemodel_frames [expr {$n == 2 && [llength $::ds9(frames)] == $f0 + 2}] "$n frames, $f0 -> [llength $::ds9(frames)]"
    R noisemodel_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc sec_sedcodes {} {
    set root [file normalize [file join [::ogf::step::plugin_dir sedcodes] .. ..]]
    set mock [file join $root sed_adapters mocks]
    ::ogf::params::put sedcodes mag-columns {MAG_AUTO:F606W,MAG_ISOCOR:F775W,MAG_APER:F850LP,MAG_APER_2:F105W,MAG_APER_3:F125W,MAG_APER_5:F160W}
    ::ogf::params::put sedcodes max-objects 40
    ::ogf::params::put sedcodes engine external
    ::ogf::params::put sedcodes command "[OGFPython] [file join $mock mock_eazy.py]"
    ::ogf::params::put sedcodes filters-res [file join $mock FILTER.RES.toy]
    ::ogf::params::put sedcodes z-step 0.05
    ::ogf::params::put sedcodes z-max 4.0
    lassign [run_step sedcodes photoz] ok recs
    R sedcodes_photoz_ran $ok $recs
    R sedcodes_photoz_recorded [expr {[lindex $recs 0 0] eq "analysis.sedcodes_photoz"}] $recs
    set cols [::ogf::cat::columns]
    R sedcodes_photoz_columns [expr {"EZ_Z" in $cols && "EZ_Z16" in $cols && "EZ_Z84" in $cols && "EZ_CHI2" in $cols}]
    R sedcodes_photoz_rows [expr {[nonempty EZ_Z] >= 30}] "filled=[nonempty EZ_Z]"
    set zs [lsearch -all -inline -not [col_values EZ_Z] {}]
    R sedcodes_photoz_range [expr {[tcl::mathfunc::min {*}$zs] >= 0.0 && [tcl::mathfunc::max {*}$zs] <= 4.0}] "[tcl::mathfunc::min {*}$zs] .. [tcl::mathfunc::max {*}$zs]"
    ::ogf::params::put sedcodes sed-code cigale
    ::ogf::params::put sedcodes z-column EZ_Z
    ::ogf::params::put sedcodes command "[OGFPython] [file join $mock mock_cigale.py]"
    lassign [run_step sedcodes sedfit] ok recs
    R sedcodes_sedfit_ran $ok $recs
    R sedcodes_sedfit_recorded [expr {[lindex $recs 0 0] eq "analysis.sedcodes_sedfit"}] $recs
    set cols [::ogf::cat::columns]
    R sedcodes_sedfit_columns [expr {"SC_LOGM" in $cols && "SC_AV" in $cols && "SC_CHI2" in $cols}]
    R sedcodes_sedfit_rows [expr {[nonempty SC_LOGM] >= 30}] "filled=[nonempty SC_LOGM]"
    set lm [lsearch -all -inline -not [col_values SC_LOGM] {}]
    R sedcodes_logm_finite [expr {[tcl::mathfunc::min {*}$lm] > 0 && [tcl::mathfunc::max {*}$lm] < 20}] "[tcl::mathfunc::min {*}$lm] .. [tcl::mathfunc::max {*}$lm]"
    R sedcodes_file [expr {[file exists [file join [OGFSessWorkDir] sedcodes sedcodes_results.json]] && [::ogf::cat::exists sedcodes,results_file]}]
    R sedcodes_work_files [expr {[file exists [file join [OGFSessWorkDir] sedcodes cigale pcigale.ini]] && [file exists [file join [OGFSessWorkDir] sedcodes eazy catalog.cat]]}]
    OGFSedcodesCheck; update
    R sedcodes_check_window [expr {[winfo exists .ogftext] && [string match "*== eazy ==*== prospector ==*" [.ogftext.t get 1.0 end]]}]
    catch {destroy .ogftext}
    OGFSedcodesProfiles; update
    R sedcodes_profiles [expr {[file exists [::ogf::cat::get sedcodes,profile_file {}]]}]
    R sedcodes_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc sec_daophot {} {
    set f0 [llength $::ds9(frames)]
    ::ogf::params::put daophot fwhm 3.5
    ::ogf::params::put daophot n-iter 2
    ::ogf::params::put daophot neighbour-iter 1
    ::ogf::params::put daophot psf-order 1
    # manual PSF-star list editing from the catalog selection (recorded parameters)
    set nums [lrange [::ogf::cat::values NUMBER] 0 2]
    ::ogf::cat::select $nums replace 0
    update
    OGFDaophotAddPsf
    R daophot_psf_add_param [expr {[llength [split [::ogf::params::get daophot psf-add] \;]] == 3}] [::ogf::params::get daophot psf-add]
    OGFDaophotRemPsf
    R daophot_psf_swap [expr {[::ogf::params::get daophot psf-add] eq {} && [llength [split [::ogf::params::get daophot psf-remove] \;]] == 3}] "[::ogf::params::get daophot psf-add] | [::ogf::params::get daophot psf-remove]"
    ::ogf::params::put daophot psf-remove {}
    foreach st {find phot pickpsf psf} {
	lassign [run_step daophot $st] ok recs
	R daophot_${st}_ran $ok $recs
	R daophot_${st}_recorded [expr {[lindex $recs 0 0] eq "analysis.daophot_$st"}] $recs
    }
    set w [OGFSessWorkDir]
    foreach f {daophot_find.tsv daophot_phot.tsv daophot_psfstars.tsv daophot_psf.json daophot_psf.fits} {
	R daophot_file_$f [expr {[file exists [file join $w daophot $f]] && [file size [file join $w daophot $f]] > 50}]
    }
    lassign [run_step daophot fit] ok recs
    R daophot_fit_ran $ok $recs
    R daophot_fit_recorded [expr {[lindex $recs 0 0] eq "analysis.daophot_fit"}] $recs
    set cols [::ogf::cat::columns]
    R daophot_columns [expr {"DAO_MAG" in $cols && "DAO_CHI" in $cols && "DAO_FLAG" in $cols && "DAO_MAG_APC" in $cols}]
    R daophot_rows_filled [expr {[nonempty DAO_MAG] > 30}] "rows=[nonempty DAO_MAG]"
    foreach f {daophot_stars.tsv daophot_resid.fits daophot_diag.png daophot_result.json} {
	R daophot_file_$f [expr {[file exists [file join $w daophot $f]] && [file size [file join $w daophot $f]] > 100}]
    }
    R daophot_frames [expr {[llength $::ds9(frames)] == $f0 + 1}] "frames [llength $::ds9(frames)] (was $f0)"
    R daophot_frame_restored [expr {$::current(frame) eq [lindex $::ds9(frames) 0]}] $::current(frame)
    foreach st {substar apcorr} {
	lassign [run_step daophot $st] ok recs
	R daophot_${st}_ran $ok $recs
    }
    R daophot_sub_file [file exists [file join $w daophot daophot_sub.fits]]
    OGFDaophotDiag; update
    R daophot_diag_window [expr {[winfo exists .ogfdaodiag] && [image width .ogfdaodiag_img] > 300 && [image height .ogfdaodiag_img] > 200}] "[image width .ogfdaodiag_img]x[image height .ogfdaodiag_img]"
    destroy .ogfdaodiag
    OGFDaophotTable; update
    R daophot_argv_templated [expr {[string match {*@{WORK}*} [dict get [lindex [::ogf::session::steps] end] argv_t]]}]
    R daophot_cat_key [file exists [::ogf::cat::get daophot,stars_file {}]]
    R daophot_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

# every tab's chip row must fit into the panel width (chips that do not fit live in the "More" menu)
proc chips_fit {} {
    global ogfui
    set cur $ogfui(tab)
    set cw [winfo width $ogfui(content)]
    set bad {}
    set more {}
    foreach t $::ogf::tabs {
	OGFUIShowTab $t; update
	set right 0
	foreach c [winfo children $ogfui(tabframe,$t)] {
	    if {[winfo ismapped $c]} {set r [expr {[winfo x $c] + [winfo reqwidth $c]}]; if {$r > $right} {set right $r}}
	}
	if {$right > $cw} {lappend bad "$t:$right>$cw"}
	if {[info exists ogfui(overflow,$t)]} {lappend more "$t=[join $ogfui(overflow,$t) /]"}
    }
    OGFUIShowTab $cur; update
    R chips_fit [expr {$bad eq {}}] "content=$cw overflow-menus: $more $bad"
    # the plugin that moved into a More menu is still reachable through the cascade
    if {[info exists ogfui(more,Measure)]} {
	set m [$ogfui(more,Measure) cget -menu]
	# every plugin that moved into the More menu stays reachable: each entry is a cascade whose sub-menu exists
	set okm [expr {[$m index end] ne "none" && [$m index end] >= 0}]
	for {set i 0} {$i <= [$m index end]} {incr i} {
	    set sub [lindex [$m entryconfigure $i -menu] end]
	    if {$sub eq {} || ![winfo exists $sub]} {set okm 0}
	}
	R chips_more_menu $okm "[$m index end] entries; first [lindex [$m entryconfigure 0 -menu] end]"

    }
}

proc run {} {
    global catpanel
    set only [expr {[info exists ::env(OGF_NP_ONLY)] ? [split $::env(OGF_NP_ONLY) ,] : {}}]
    R geometry_before [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    CatalogPanelExtract
    set t0 [clock milliseconds]; while {![::ogf::cat::has] && [clock milliseconds]-$t0 < 90000} {update; after 100}
    wait_idle 500
    R extracted [expr {[::ogf::cat::nrows] > 50}] "rows=[::ogf::cat::nrows]"
    foreach sec {isophote completeness daophot psfex multifit morphext noisemodel sedcodes} {
	if {[want $only $sec]} {
	    if {[catch {sec_$sec} err]} {R ${sec}_error 0 "$err [string range $::errorInfo 0 300]"}
	}
    }
    # export the recorded session and replay it headless: the new columns must be byte-identical to the GUI ones
    set sp [file join $::dir ogfinder_session.py]
    CatalogPanelSessionSave $sp
    set fd [open [file join $::dir gui_catalog.tsv] w]; puts -nonewline $fd [::ogf::cat::tsv]; close $fd
    R session_exported [file exists $sp]
    chips_fit
    R final_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    puts $::fh "SUMMARY failures=$::nf"; close $::fh
}
after 3000 {
    if {[catch run err]} {puts $::fh "FAIL run_error $err $::errorInfo"; incr ::nf; puts $::fh "SUMMARY failures=$::nf"; close $::fh}
    exit 0
}
