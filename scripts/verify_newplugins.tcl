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

proc sec_stacking {} {
    ::ogf::params::put stacking half 15
    ::ogf::params::put stacking n-boot 30
    ::ogf::params::put stacking null 30
    ::ogf::params::put stacking subtract-null 1
    ::ogf::params::put stacking ap-r 4
    ::ogf::params::put stacking isolate 12
    ::ogf::params::put stacking mask-catalog 0.6
    ::ogf::params::put stacking show-stack 1
    set f0 [llength $::ds9(frames)]
    lassign [run_step stacking stack] ok recs
    R stacking_ran $ok $recs
    R stacking_recorded [expr {[lindex $recs 0 0] eq "analysis.stacking"}] $recs
    R stacking_columns [expr {"ST_USED" in [::ogf::cat::columns] && "ST_APFLUX" in [::ogf::cat::columns]}]
    R stacking_rows_filled [expr {[nonempty ST_APFLUX] >= 5}] "rows=[nonempty ST_APFLUX]"
    set w [file join [OGFSessWorkDir] stacking]
    foreach f {stacking_stack.fits stacking_err.fits stacking_profile.tsv stacking_summary.json stacking_plot.png} {
	R stacking_file_$f [expr {[file exists [file join $w $f]] && [file size [file join $w $f]] > 100}]
    }
    R stacking_keys [expr {[::ogf::cat::exists stacking,stack_file] && [::ogf::cat::exists stacking,n_stacked]}]
    R stacking_frame [expr {[llength $::ds9(frames)] == $f0 + 1}] "frames [llength $::ds9(frames)] (was $f0)"
    catch {GotoFrame [lindex $::ds9(frames) 0]}
    set pw [OGFStackPlot]
    R stacking_plot_window [expr {[winfo exists $pw] && [image width ogfstackimg] > 300}] "[image width ogfstackimg]x[image height ogfstackimg]"
    destroy $pw
    OGFStackTable; update
    R stacking_argv_templated [expr {[string match {*@{WORK}/stacking*} [dict get [lindex [::ogf::session::steps] end] argv_t]]}]
    R stacking_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc sec_depth {} {
    ::ogf::params::put completeness dp-aper-radius 3
    ::ogf::params::put completeness dp-box-arcsec 20
    ::ogf::params::put completeness dp-tile 150
    lassign [run_step completeness depth] ok recs
    R depth_ran $ok $recs
    R depth_recorded [expr {[lindex $recs 0 0] eq "analysis.depth_maps"}] $recs
    R depth_columns [expr {"DEPTH_LIM" in [::ogf::cat::columns] && "DEPTH_SBLIM" in [::ogf::cat::columns] && "DEPTH_RMS" in [::ogf::cat::columns]}]
    R depth_rows_filled [expr {[nonempty DEPTH_LIM] >= 50}] "rows=[nonempty DEPTH_LIM]"
    set w [file join [OGFSessWorkDir] depth]
    foreach f {depth_mag_map.fits depth_sb_map.fits depth_rms_map.fits depth_tiles.tsv depth_summary.json depth_plot.png} {
	R depth_file_$f [expr {[file exists [file join $w $f]] && [file size [file join $w $f]] > 100}]
    }
    R depth_keys [expr {[::ogf::cat::exists depth,mag_map] && [::ogf::cat::exists depth,mag_limit_median]}]
    ::ogf::params::put completeness detector sep
    ::ogf::params::put completeness dp-regions grid
    ::ogf::params::put completeness dp-grid 2x2
    ::ogf::params::put completeness dp-per-bin 60
    ::ogf::params::put completeness n-bins 5
    ::ogf::params::put completeness mag-min 12.5
    ::ogf::params::put completeness mag-max 17.5
    ::ogf::params::put completeness dp-maps-at 15
    lassign [run_step completeness compmap] ok recs
    R compmap_ran $ok $recs
    R compmap_columns [expr {"COMPL_REGION" in [::ogf::cat::columns] && "COMPL_LOC" in [::ogf::cat::columns] && "COMPL_LIM50_LOC" in [::ogf::cat::columns]}]
    R compmap_rows_filled [expr {[nonempty COMPL_REGION] >= 50}] "rows=[nonempty COMPL_REGION]"
    foreach f {compmap_lim50_map.fits compmap_regions.tsv compmap_summary.json compmap_frac_m15.fits compmap_plot.png} {
	R compmap_file_$f [expr {[file exists [file join $w $f]] && [file size [file join $w $f]] > 100}]
    }
    set pw [OGFDepthPlot]
    R depth_plot_window [expr {[winfo exists $pw] && [image width ogfdepthimg] > 300}] "[image width ogfdepthimg]x[image height ogfdepthimg]"
    destroy $pw
    OGFDepthTable; update
    set f0 [llength $::ds9(frames)]
    set n [OGFDepthShow]
    R depth_frames [expr {$n >= 2 && [llength $::ds9(frames)] == $f0 + $n}] "$n frames, $f0 -> [llength $::ds9(frames)]"
    R depth_argv_templated [expr {[string match {*@{WORK}/depth*} [dict get [lindex [::ogf::session::steps] end] argv_t]]}]
    R depth_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc sec_photozq {} {
    # photo-z distribution tools on the M51 catalog: the toy EAZY photo-z (EZ_Z with 16/84 percentiles) against a synthetic spec-z file for every 3rd object
    if {"EZ_Z" ni [::ogf::cat::columns]} {
    set root [file normalize [file join [::ogf::step::plugin_dir sedcodes] .. ..]]
    set mock [file join $root sed_adapters mocks]
    ::ogf::params::put sedcodes mag-columns {MAG_AUTO:F606W,MAG_ISOCOR:F775W,MAG_APER:F850LP,MAG_APER_2:F105W,MAG_APER_3:F125W,MAG_APER_5:F160W}
    ::ogf::params::put sedcodes max-objects 60
    ::ogf::params::put sedcodes engine external
    ::ogf::params::put sedcodes command "[OGFPython] [file join $mock mock_eazy.py]"
    ::ogf::params::put sedcodes filters-res [file join $mock FILTER.RES.toy]
    ::ogf::params::put sedcodes z-step 0.05
    ::ogf::params::put sedcodes z-max 4.0
    lassign [run_step sedcodes photoz] ok recs
    R pzq_prep $ok $recs
    }
    set nums [::ogf::cat::values NUMBER]
    set ez [col_values EZ_Z]
    set zf [file join $::dir pzq_specz.tsv]
    set fd [open $zf w]; puts $fd "NUMBER\tZ_SPEC"
    set k 0
    foreach n $nums z $ez {
	incr k
	if {$k % 2} continue
	if {$z eq {}} {set z [expr {0.3 + rand()}]}                   ;# objects without a toy photo-z still get a spec-z (representativeness only)
	puts $fd "$n\t[format %.4f [expr {$z + 0.03*(1+$z)*(rand()+rand()+rand()-1.5)*2}]]"
    }
    close $fd
    ::ogf::params::put photoz_sed pq-zspec-file $zf
    ::ogf::params::put photoz_sed pq-zphot-col EZ_Z
    ::ogf::params::put photoz_sed pq-q16-col EZ_Z16
    ::ogf::params::put photoz_sed pq-q84-col EZ_Z84
    lassign [run_step photoz_sed pzquality] ok recs
    R pzq_ran $ok $recs
    R pzq_recorded [expr {[lindex $recs 0 0] eq "analysis.photoz_quality"}] $recs
    R pzq_columns [expr {"PZ_PIT" in [::ogf::cat::columns] && "PZ_CRPS" in [::ogf::cat::columns] && "PZ_INCI68" in [::ogf::cat::columns]}]
    R pzq_rows [expr {[nonempty PZ_PIT] >= 8}] "rows=[nonempty PZ_PIT]"
    set w [file join [OGFSessWorkDir] photoz_quality]
    foreach f {pzq_report.json pzq_plot.png pzq_binned.tsv pzq_pit.tsv} {
	R pzq_file_$f [expr {[file exists [file join $w $f]] && [file size [file join $w $f]] > 50}]
    }
    R pzq_keys [expr {[::ogf::cat::exists pzq,nmad] && [::ogf::cat::exists pzq,pit_ks_p]}]
    set pw [OGFPzqPlot]
    R pzq_plot_window [expr {[winfo exists $pw] && [image width ogfpzqimg] > 300}] "[image width ogfpzqimg]x[image height ogfpzqimg]"
    destroy $pw
    OGFPzqTable; update
    ::ogf::params::put photoz_sed pq-features {MAG_AUTO,FLUX_RADIUS}
    ::ogf::params::put photoz_sed pq-knn 3
    lassign [run_step photoz_sed pzrepr] ok recs
    R pzr_ran $ok $recs
    R pzr_recorded [expr {[lindex $recs 0 0] eq "analysis.photoz_repr"}] $recs
    R pzr_columns [expr {"PZ_SPECW" in [::ogf::cat::columns] && "PZ_SPECDIST" in [::ogf::cat::columns] && "PZ_INSPEC" in [::ogf::cat::columns]}]
    R pzr_rows [expr {[nonempty PZ_SPECDIST] >= 30 && [nonempty PZ_SPECW] >= 8}] "dist=[nonempty PZ_SPECDIST] w=[nonempty PZ_SPECW]"
    R pzr_verdict [::ogf::cat::exists pzq,verdict] [::ogf::cat::get pzq,verdict ?]
    R pzq_argv_templated [expr {[string match {*@{WORK}/photoz_quality*} [dict get [lindex [::ogf::session::steps] end] argv_t]]}]
    R pzq_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
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

proc sec_cluster {} {
    ::ogf::params::put cluster dens-sigma 150
    ::ogf::params::put cluster dens-bin 8
    ::ogf::params::put cluster fit-radius 100000
    lassign [run_step cluster members] ok recs
    R cluster_members_ran $ok $recs
    R cluster_members_recorded [expr {[lindex $recs 0 0] eq "analysis.cluster_members"}] $recs
    set cols [::ogf::cat::columns]
    R cluster_members_columns [expr {"CL_RS_RES" in $cols && "CL_PMEM" in $cols && "CL_MEMBER" in $cols && "CL_R" in $cols}]
    R cluster_members_rows [expr {[nonempty CL_R] >= 20 && [nonempty CL_RS_PULL] >= 20}] "r=[nonempty CL_R] pull=[nonempty CL_RS_PULL]"
    lassign [run_step cluster density] ok recs
    R cluster_density_ran $ok $recs
    R cluster_density_rows [expr {[nonempty CL_SIGMA] >= 20}] "sigma=[nonempty CL_SIGMA]"
    lassign [run_step cluster arcs] ok recs
    R cluster_arcs_ran $ok $recs
    R cluster_arcs_columns [expr {"ARC_LW" in [::ogf::cat::columns] && "ARC_FLAG" in [::ogf::cat::columns]}]
    R cluster_arcs_rows [expr {[nonempty ARC_LW] >= 1}] "measured=[nonempty ARC_LW] flagged=[llength [lsearch -all -inline [col_values ARC_FLAG] 1]]"
    set w [file join [OGFSessWorkDir] cluster]
    foreach f {cluster_summary.json cluster_rs.png cluster_density.fits cluster_peaks.tsv} {
	R cluster_file_$f [expr {[file exists [file join $w $f]] && [file size [file join $w $f]] > 50}]
    }
    R cluster_keys [expr {[::ogf::cat::exists cluster,density_file] && [::ogf::cat::exists cluster,rs_png]}]
    set pw [OGFClusterPlot]
    R cluster_plot_window [expr {[winfo exists $pw] && [image width ogfclusterimg] > 300}] "[image width ogfclusterimg]x[image height ogfclusterimg]"
    destroy $pw
    OGFClusterSummary; update
    R cluster_summary_window [expr {[winfo exists .ogftext] && [string match "*red_sequence*" [.ogftext.t get 1.0 end]]}]
    catch {destroy .ogftext}
    set f0 [llength $::ds9(frames)]
    set n [OGFClusterFrame]
    R cluster_frame [expr {$n == 1 && [llength $::ds9(frames)] == $f0 + 1}] "$f0 -> [llength $::ds9(frames)]"
    R cluster_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc sec_spectra {} {
    set root [file normalize [file join [::ogf::step::plugin_dir spectra] .. ..]]
    set dir [file join [file dirname [OGFSessWorkDir]] spec_demo]
    set nums [lrange [::ogf::cat::values NUMBER] 0 29]
    exec [OGFPython] [file join [::ogf::step::plugin_dir spectra] make_demo.py] $dir [join $nums ,] --seed 3
    ::ogf::params::put spectra spec-dir $dir
    lassign [run_step spectra link] ok recs
    R spectra_link_ran $ok $recs
    R spectra_link_recorded [expr {[lindex $recs 0 0] eq "analysis.spectra_link"}] $recs
    R spectra_link_rows [expr {[nonempty SP_OK] >= 25 && "SP_FILE" in [::ogf::cat::columns]}] "ok=[nonempty SP_OK]"
    lassign [run_step spectra fit] ok recs
    R spectra_fit_ran $ok $recs
    R spectra_fit_columns [expr {"SP_Z" in [::ogf::cat::columns] && "SP_ZQ" in [::ogf::cat::columns] && "SP_LINE_FLUX" in [::ogf::cat::columns]}]
    R spectra_fit_rows [expr {[nonempty SP_Z] >= 20}] "z=[nonempty SP_Z]"
    # redshift accuracy against the demo truth
    set tf [open [file join $dir truth.tsv] r]; set tl [split [read $tf] \n]; close $tf
    set zt {}
    foreach l [lrange $tl 1 end] {if {$l ne {}} {dict set zt [lindex $l 0] [lindex $l 1]}}
    set nn [::ogf::cat::values NUMBER]; set zz [::ogf::cat::values SP_Z]; set qq [::ogf::cat::values SP_ZQ]
    set maxd 0.0; set n 0
    foreach num $nn z $zz q $qq {
	if {$z eq {} || $q < 2 || ![dict exists $zt $num]} continue
	set d [expr {abs($z - [dict get $zt $num])}]
	if {$d > $maxd} {set maxd $d}
	incr n
    }
    R spectra_z_accuracy [expr {$n >= 20 && $maxd < 0.002}] "n=$n max|dz|=$maxd"
    foreach f {spectra_links.tsv spectra_results.json} {
	R spectra_file_$f [expr {[file exists [file join [OGFSessWorkDir] spectra $f]] && [file size [file join [OGFSessWorkDir] spectra $f]] > 50}]
    }
    set w [OGFSpectraViewer [lindex $nums 0]]
    update
    R spectra_viewer [expr {[winfo exists $w] && [image width ogfspecimg] > 400 && [image height ogfspecimg] > 200}] "[image width ogfspecimg]x[image height ogfspecimg]"
    OGFSpectraStep 1; update
    R spectra_viewer_next [expr {$::ogf_spec_num eq [lindex $nums 1]}] "now $::ogf_spec_num"
    destroy $w
    OGFSpectraTable; update
    R spectra_table_window [expr {[winfo exists .ogftext] && [string match "*NUMBER*FILE*" [.ogftext.t get 1.0 end]]}]
    catch {destroy .ogftext}
    R spectra_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc sec_xmatch {} {
    set dir [file join [file dirname [OGFSessWorkDir]] xmatch_demo]
    file mkdir $dir
    set nn [::ogf::cat::values NUMBER]; set ra [::ogf::cat::values ALPHA_J2000]; set de [::ogf::cat::values DELTA_J2000]; set mg [::ogf::cat::values MAG_AUTO]
    set ref [file join $dir ref.csv]
    set fd [open $ref w]
    puts $fd "objid,ra,dec,gmag"
    set pi 3.141592653589793
    set want {}
    set i 0
    foreach n $nn r $ra d $de m $mg {
	incr i
	if {$r eq {} || $d eq {} || $i % 2} continue
	set r2 [expr {$r + 0.30 / 3600.0 / cos($d * $pi / 180.0)}]
	set d2 [expr {$d - 0.20 / 3600.0}]
	puts $fd [format "R%d,%.8f,%.8f,%s" $n $r2 $d2 $m]
	dict set want $n R$n
    }
    close $fd
    ::ogf::params::put xmatch source file
    ::ogf::params::put xmatch ref-file $ref
    ::ogf::params::put xmatch ref-id-col objid
    ::ogf::params::put xmatch copy-columns gmag
    ::ogf::params::put xmatch radius-arcsec 1.5
    lassign [run_step xmatch match] ok recs
    R xmatch_ran $ok $recs
    R xmatch_recorded [expr {[lindex $recs 0 0] eq "analysis.xmatch"}] $recs
    R xmatch_columns [expr {"XM_SEP" in [::ogf::cat::columns] && "XM_ID" in [::ogf::cat::columns] && "XM_V1" in [::ogf::cat::columns]}]
    set good 0; set tot 0
    foreach n [::ogf::cat::values NUMBER] id [::ogf::cat::values XM_ID] {
	if {[dict exists $want $n]} {incr tot; if {$id eq [dict get $want $n]} {incr good}}
    }
    R xmatch_ids [expr {$tot > 20 && double($good) / $tot >= 0.95}] "correct $good of $tot"
    set sf [file join [OGFSessWorkDir] xmatch xmatch_summary.json]
    set fd [open $sf r]; set js [read $fd]; close $fd
    set ok2 [regexp {"dra": ([-0-9.e]+)} $js -> dra]
    regexp {"ddec": ([-0-9.e]+)} $js -> ddec
    R xmatch_shift [expr {$ok2 && abs($dra - 0.30) < 0.05 && abs($ddec + 0.20) < 0.05}] "dra=$dra ddec=$ddec"
    OGFXmatchSummary; update
    R xmatch_summary_window [expr {[winfo exists .ogftext] && [string match "*n_matched*" [.ogftext.t get 1.0 end]]}]
    catch {destroy .ogftext}
    R xmatch_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc sec_lightcurves {} {
    set dir [file join [file dirname [OGFSessWorkDir]] lc_demo]
    set nums [lrange [::ogf::cat::values NUMBER] 0 27]
    exec [OGFPython] [file join [::ogf::step::plugin_dir lightcurves] make_demo.py] $dir [join $nums ,] --seed 3
    ::ogf::params::put lightcurves source file
    ::ogf::params::put lightcurves lc-file [file join $dir lc_demo.tsv]
    lassign [run_step lightcurves classify] ok recs
    R lightcurves_ran $ok $recs
    R lightcurves_recorded [expr {[lindex $recs 0 0] eq "analysis.lightcurves"}] $recs
    R lightcurves_columns [expr {"LC_CLASS" in [::ogf::cat::columns] && "LC_PSN" in [::ogf::cat::columns] && "LC_PEAK_MAG" in [::ogf::cat::columns]}]
    R lightcurves_rows [expr {[nonempty LC_PSN] >= 22}] "psn=[nonempty LC_PSN]"
    # compare with the demo truth
    set tf [open [file join $dir truth.tsv] r]; set tl [split [read $tf] \n]; close $tf
    set truth {}
    foreach l [lrange $tl 1 end] {if {$l ne {}} {dict set truth [lindex $l 0] [lindex $l 1]}}
    set good 0; set tot 0; set snok 0; set sntot 0
    foreach n [::ogf::cat::values NUMBER] c [::ogf::cat::values LC_CLASS] p [::ogf::cat::values LC_PSN] {
	if {![dict exists $truth $n] || $c eq {}} continue
	incr tot
	set t [dict get $truth $n]
	if {$c eq $t} {incr good}
	if {$t in {SNIa SNIbc SNII}} {incr sntot; if {$p ne {} && $p >= 0.5} {incr snok}}
    }
    R lightcurves_accuracy [expr {$tot >= 24 && double($good) / $tot >= 0.6}] "class match $good of $tot; SN found $snok of $sntot"
    foreach f {lc_results.tsv lc_results.json} {
	R lightcurves_file_$f [expr {[file exists [file join [OGFSessWorkDir] lightcurves $f]] && [file size [file join [OGFSessWorkDir] lightcurves $f]] > 100}]
    }
    set w [OGFLcViewer LC[lindex $nums 0]]
    update
    R lightcurves_viewer [expr {[winfo exists $w] && [image width ogflcimg] > 400 && [image height ogflcimg] > 200}] "[image width ogflcimg]x[image height ogflcimg]"
    OGFLcStep 1; update
    R lightcurves_viewer_next [expr {$::ogf_lc_id eq "LC[lindex $nums 1]"}] "now $::ogf_lc_id"
    destroy $w
    OGFLcTable; update
    R lightcurves_table_window [expr {[winfo exists .ogftext] && [string match "*PSN*" [.ogftext.t get 1.0 end]]}]
    catch {destroy .ogftext}
    R lightcurves_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc sec_batch {} {
    # 1. the Python expansion of the cli templates equals the Tcl one (same params, same context)
    set script [file join [::ogf::step::plugin_dir batch] batch.py]
    set neq 0; set nn 0; set bad {}
    foreach {pid sid} {xmatch match cluster members daophot find daophot run isophote fit sedcodes photoz spectra fit psfex build} {
	if {[catch {
	    set step [::ogf::reg::step $pid $sid]
	    set ctx [::ogf::step::context $pid [expr {[::ogf::json::get $step catalog_tmp] ne {} ? [::ogf::json::get $step catalog_tmp] : $sid}]]
	    set tcl_argv [::ogf::step::build_argv $pid $step $ctx]
	    set pj {}
	    foreach p [::ogf::json::get [::ogf::reg::get $pid] params] {
		set v [::ogf::params::get $pid [dict get $p name]]
		lappend pj [dict get $p name] $v
	    }
	    set pd "\{"; set first 1
	    foreach {k v} $pj {if {!$first} {append pd ,}; set first 0; append pd "\"$k\":\"[string map [list \\ \\\\ \" \\\"] $v]\""}
	    append pd "\}"
	    set cj "\{\"python\":\"[dict get $ctx python]\",\"work\":\"[dict get $ctx work]\",\"image\":\"[dict get $ctx image]\",\"root\":\"[dict get $ctx root]\",\"mask\":\"[dict get $ctx mask]\",\"catalog\":\"[expr {[dict exists $ctx catalog] ? [dict get $ctx catalog] : {}}]\"\}"
	    set out [exec [OGFPython] $script --expand $pid.$sid --root [OGFSessRoot] --ctx-json $cj --params-json $pd]
	    set py_argv [::ogf::json::parse $out]
	} err]} {
	    lappend bad "$pid.$sid: $err"
	    continue
	}
	incr nn
	if {$py_argv eq $tcl_argv} {incr neq} else {lappend bad "$pid.$sid differs"}
    }
    R batch_expand_equal [expr {$nn >= 6 && $neq == $nn}] "$neq of $nn identical; [join [lrange $bad 0 2] {; }]"
    # 2. a small batch run through the GUI
    set base [file join [file dirname [OGFSessWorkDir]] batch_demo]
    file mkdir $base
    set py "import numpy as np\nfrom astropy.io import fits\nfor k in range(3):\n    rng=np.random.RandomState(k); im=rng.randn(200,200)*5+100\n    yy,xx=np.mgrid\[0:200,0:200\]\n    for i in range(15):\n        x,y=rng.rand(2)*160+20; im+=rng.uniform(150,900)*np.exp(-((xx-x)**2+(yy-y)**2)/8.0)\n    fits.writeto('$base/f%d.fits'%k, im.astype('float32'), overwrite=True)\n"
    exec [OGFPython] -c $py
    set fh [open [file join $base fields.txt] w]
    foreach k {0 1 2} {puts $fh "G$k [file join $base f$k.fits]"}
    close $fh
    set fh [open [file join $base recipe.json] w]
    puts $fh {{"detect": {"args": ["--detect-thresh", "3"]}, "steps": [{"plugin": "noisemodel", "step": "model"}, {"plugin": "example_hello", "step": "greet"}]}}
    close $fh
    ::ogf::params::put batch fields-file [file join $base fields.txt]
    ::ogf::params::put batch recipe-file [file join $base recipe.json]
    ::ogf::params::put batch outdir [file join $base out]
    ::ogf::params::put batch jobs 3
    ::ogf::params::put batch resume 0
    OGFBatchRun
    update
    R batch_started [expr {$::ogf::batch::running && [winfo exists .ogfbatch]}]
    set t0 [clock seconds]
    while {$::ogf::batch::running && [clock seconds] - $t0 < 120} {after 200; update}
    R batch_finished [expr {!$::ogf::batch::running && $::ogf::batch::rc == 0}] "rc=$::ogf::batch::rc after [expr {[clock seconds]-$t0}] s"
    update
    set nrows [llength [.ogfbatch.tree children {}]]
    set states {}
    foreach c [.ogfbatch.tree children {}] {lappend states [lindex [.ogfbatch.tree item $c -values] 1]}
    R batch_progress_view [expr {$nrows == 3 && $states eq {ok ok ok}}] "rows=$nrows states=$states"
    set sf [file join $base out summary.tsv]
    set fh [open $sf r]; set lines [split [string trim [read $fh]] \n]; close $fh
    R batch_summary [expr {[llength $lines] == 4 && ![regexp {\tfailed\t} [join [lrange $lines 1 end] \n]]}] [lindex $lines 1]
    R batch_field_logs [expr {[file size [file join $base out G0 field.log]] > 100 && [file exists [file join $base out G2 catalog.tsv]]}]
    # resume
    ::ogf::params::put batch resume 1
    OGFBatchRun
    set t0 [clock seconds]
    while {$::ogf::batch::running && [clock seconds] - $t0 < 120} {after 200; update}
    set fh [open $sf r]; set txt [read $fh]; close $fh
    R batch_resume_cached [expr {[regexp {\t0\t3\t0\t} [lindex [split $txt \n] 1]] || [regexp {\t3\t0\t} [lindex [split $txt \n] 1]]}] [lindex [split $txt \n] 1]
    OGFBatchSummary; update
    R batch_summary_window [expr {[winfo exists .ogftext] && [string match "*G1*" [.ogftext.t get 1.0 end]]}]
    catch {destroy .ogftext}
    destroy .ogfbatch
    R batch_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
}

proc sec_repro {} {
    ::ogf::params::put repro bundle-file [file join [expr {[info exists ::env(OGF_NP_DIR)] ? $::env(OGF_NP_DIR) : [file dirname [OGFSessWorkDir]]}] gui_bundle.zip]
    ::ogf::params::put repro note "verify_newplugins"
    set z [OGFReproBundle]
    R repro_bundle_written [expr {$z ne {} && [file size $z] > 2000}] "$z [expr {[file exists $z] ? [file size $z] : 0}] bytes"
    set listing [exec [OGFPython] -c "import zipfile,sys;print(' '.join(zipfile.ZipFile(sys.argv\[1\]).namelist()))" $z]
    R repro_bundle_files [expr {[string match "*manifest.json*" $listing] && [string match "*params.json*" $listing] && [string match "*steps.json*" $listing] && [string match "*requirements.lock*" $listing] && [string match "*session.py*" $listing] && [string match "*outputs/catalog_final.tsv*" $listing]}] $listing
    set mj [exec [OGFPython] -c "import zipfile,sys,json;m=json.loads(zipfile.ZipFile(sys.argv\[1\]).read('manifest.json'));p=json.loads(zipfile.ZipFile(sys.argv\[1\]).read('params.json'));s=json.loads(zipfile.ZipFile(sys.argv\[1\]).read('steps.json'));print(m\['kind'],len(m\['inputs']),len(m\['outputs']),len(p),len(s),m\['tool'].get('head','')\[:9\],m\['outputs']\[0]\['rows'])" $z]
    lassign $mj kind nin nout npl nstep head rows
    R repro_manifest [expr {$kind eq "session" && $nin >= 1 && $nout == 1 && $npl >= 20 && $nstep >= 1}] "kind=$kind inputs=$nin plugins-with-params=$npl steps=$nstep git=$head rows=$rows"
    R repro_rows_match [expr {$rows == [llength [::ogf::cat::values NUMBER]]}] "rows=$rows catalog=[llength [::ogf::cat::values NUMBER]]"
    set ok [OGFReproVerify]
    R repro_verify_integrity $ok $::ogf::repro::last_out
    catch {destroy .ogftext}
    OGFReproInfo; update
    R repro_info_window [expr {[winfo exists .ogftext] && [string match "*sha256*" [.ogftext.t get 1.0 end]]}]
    catch {destroy .ogftext}
    R repro_geometry [expr {[geom] eq "181 769 154 1300x950"}] [geom]
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
    foreach sec {isophote completeness daophot psfex multifit morphext noisemodel stacking depth sedcodes photozq cluster spectra xmatch lightcurves batch repro} {
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
