# Equivalence test of the "headless" cli blocks (item 1 of the second work package): for each legacy dialog / proc step that now carries a
# `headless` template, run (A) the legacy proc (dialogs are filled by a timer; exec stubbed to the fake interpreter) and (B) the manifest
# step through ::ogf::step::run ID STEP headless (real job runner + recorder, fake interpreter) and compare the argv that reaches the interpreter,
# the recorder record (step, class, templated argv, post, requires, network; the title only where the legacy title is dynamic), the catalog and the status.
#   HOME=/tmp/ogf_clihl_home OGF_CLIHL_OUT=/tmp/clihl.txt FAKE_LOG=/tmp/clihl_fake.log OGFINDER_PYTHON=$PWD/scripts/fake_python_for_templates.sh \
#     DISPLAY=:77 bin/ds9 /workspace/fits/m51.fits -geometry 1300x950 -source scripts/verify_cli_headless.tcl
global catpanel current
set ::fh [open $::env(OGF_CLIHL_OUT) w]; set ::nf 0
set ::HOME [file normalize ~]
proc R {tag ok {d {}}} {puts $::fh "[expr {$ok?{PASS}:{FAIL}}] $tag $d"; flush $::fh; if {!$ok} {incr ::nf}}
proc bgerror {m} {puts $::fh "BGERROR $m"; flush $::fh}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}
proc geom {} {
    update idletasks; update
    set tf [winfo parent $::catpanel(tbl)]
    return "[winfo rooty $tf] [winfo height $tf] [winfo height $::catpanel(infoarea)] [winfo width .]x[winfo height .]"
}
rename exec ::real_exec
set ::EXEC {}; set ::STUB 0
proc exec {args} {
    set a {}
    foreach x $args {if {[string match 2>* $x]} continue; lappend a $x}
    if {$::STUB && [llength $a] > 1 && [string match *.py [lindex $a 1]]} {
	lappend ::EXEC $a
	return [::real_exec {*}$a]
    }
    return [::real_exec {*}$args]
}
proc norm {s} {string map [list $::HOME <HOME> /workspace/fits <FITS> [file normalize bin] <BIN>] $s}
proc rec_sig {from {title 1}} {
    set out {}
    foreach r [lrange [::ogf::session::steps] $from end] {
	set t [expr {$title ? "title=[dict get $r title]" : "title=-"}]
	lappend out "[dict get $r step] [dict get $r class] tool=[dict get $r tool] | [norm [dict get $r argv_t]] | post=[dict get $r post] $t requires=[dict get $r requires] network=[dict get $r network]"
    }
    return $out
}
proc back_to_image {} {
    for {set i 0} {$i < 6 && [file tail [CatalogPanelGetFITS]] ne "m51.fits"} {incr i} {catch {DeleteCurrentFrame}; wait_idle 150}
}
proc wait_job {} {set t [clock milliseconds]; while {[::ogf::job::busy] && [clock milliseconds]-$t < 60000} {update; after 50}; wait_idle 200}

# name legacy-script filler-script plugin step ?title-compare? ?kind?   kind: cat (add columns / set catalog) | none (no catalog change expected)
proc compare {name legacy filler plugin step {title 1} {kind cat} {cmpstatus 1}} {
    global catpanel
    set snap $::BASE
    ::ogf::cat::load_tsv $snap restored; update; wait_idle 100
    set ::catpanel(status) {}
    set s0 [llength [::ogf::session::steps]]; set e0 [llength $::EXEC]; set ::STUB 1
    if {$filler ne {}} {after 500 $filler}
    set rc [catch {uplevel #0 $legacy} err]
    update; wait_idle 300; set ::STUB 0
    set A_argv [lrange $::EXEC $e0 end]; set A_rec [rec_sig $s0 $title]; set A_tsv $catpanel(alldata); set A_status $catpanel(status)
    back_to_image
    R ${name}_legacy_ran [expr {!$rc && [llength [::ogf::session::steps]] == $s0 + 1}] "rc=$rc $err steps+=[expr {[llength [::ogf::session::steps]]-$s0}]"
    # ---- B: headless block
    ::ogf::cat::load_tsv $snap restored
    update; wait_idle 100
    set ::catpanel(status) {}
    set s1 [llength [::ogf::session::steps]]
    set fl $::env(FAKE_LOG); set n0 0
    catch {set f [open $fl r]; set n0 [llength [split [string trim [read $f]] \n]]; close $f}
    set ok [::ogf::step::run $plugin $step headless]
    wait_job
    set B_argv {}
    catch {set f [open $fl r]; set B_argv [lrange [split [string trim [read $f]] \n] $n0 end]; close $f}
    set B_rec [rec_sig $s1 $title]; set B_tsv $catpanel(alldata); set B_status $catpanel(status)
    back_to_image
    R ${name}_step_started [expr {$ok == 1 && [llength [::ogf::session::steps]] == $s1 + 1}] "ok=$ok"
    if {[llength $A_argv] == 1} {
	set a_line [join [lrange [lindex $A_argv 0] 1 end] { }]; set b_line [lindex $B_argv 0]
	R ${name}_argv_identical [expr {$a_line eq $b_line}] "\n   A: [norm $a_line]\n   B: [norm $b_line]"
    }
    R ${name}_record_identical [expr {$A_rec eq $B_rec}] "\n   A: [norm $A_rec]\n   B: [norm $B_rec]"
    if {$kind eq "cat"} {R ${name}_catalog_identical [expr {$A_tsv eq $B_tsv && $A_tsv ne $snap}] "cols=[llength [split [lindex [split $B_tsv \n] 0] \t]]"}
    if {$cmpstatus} {R ${name}_status_identical [expr {$A_status eq $B_status}] "A='$A_status' B='$B_status'"}
}

proc run {} {
    global catpanel current
    update; wait_idle 500
    if {![string match /tmp/* $::HOME]} {puts $::fh "SUMMARY-ABORT HOME not scratch"; close $::fh; exit 2}
    R geometry_before [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    R fake_python [string match *fake_python_for_templates.sh [OGFPython]] [OGFPython]
    # ---- extract: both paths run the real ds9_sextract on the image; the catalogs must be identical
    set catpanel(param,detect-thresh) 2.5; set catpanel(param,back-size) 48
    set ::BASE {}
    set ::catpanel(alldata) {}
    set s0 [llength [::ogf::session::steps]]
    CatalogPanelExtract
    set t0 [clock milliseconds]; while {$catpanel(alldata) eq {} && [clock milliseconds]-$t0 < 60000} {update; after 100}
    wait_idle 300
    set A_tsv $catpanel(alldata); set A_rec [rec_sig $s0]; set A_status $catpanel(status)
    R extract_legacy_rows [expr {[::ogf::cat::nrows] > 50}] [::ogf::cat::nrows]
    ::ogf::cat::load_tsv {} cleared; set catpanel(alldata) {}; update
    set s1 [llength [::ogf::session::steps]]
    set ok [::ogf::step::run extract extract headless]
    set t0 [clock milliseconds]; while {([::ogf::job::busy] || $catpanel(alldata) eq {}) && [clock milliseconds]-$t0 < 60000} {update; after 100}
    wait_idle 300
    set B_rec [rec_sig $s1]
    R extract_step_started [expr {$ok == 1}] ok=$ok
    R extract_record_identical [expr {$A_rec eq $B_rec}] "\n   A: [norm $A_rec]\n   B: [norm $B_rec]"
    R extract_catalog_identical [expr {$A_tsv eq $catpanel(alldata) && $A_tsv ne {}}] "rows=[::ogf::cat::nrows]"
    R extract_param_copy [expr {[::ogf::cat::get extract_param,detect-thresh] == 2.5}] [::ogf::cat::get extract_param,detect-thresh]
    set ::BASE $catpanel(alldata)
    set catpanel(param,n-workers) 2; set catpanel(param,mag-zeropoint) 24.5

    # ---- dual extract (dialog filled by a timer)
    set meas /workspace/fits/m51.fits
    set mp [OGFParamsDualSet $meas]
    compare dual_extract CatalogPanelDualExtract [list apply {{m} {set ::ed(dual,measure) $m; set ::ed(ok) 1}} $meas] extract dual
    # ---- segmap, completeness, multiband, crossmatch
    compare segmap CatalogPanelSegmentationMap {} photometry segmap 1 none
    ::ogf::params::put photometry comp-ninject 120; ::ogf::params::put photometry comp-magmin 21.0; ::ogf::params::put photometry comp-magmax 27.5; ::ogf::params::put photometry comp-nbins 7
    compare completeness CatalogPanelCompleteness [list apply {{} {set ::ed(comp,ninject) 120; set ::ed(comp,magmin) 21.0; set ::ed(comp,magmax) 27.5; set ::ed(comp,nbins) 7; set ::ed(ok) 1}}] photometry completeness
    ::ogf::params::put photometry mb-bands "g:/workspace/fits/g.fits,r:/workspace/fits/r.fits"
    compare multiband CatalogPanelMultiBand [list apply {{} {.multibandphot.param.tbands insert 1.0 "g:/workspace/fits/g.fits\nr:/workspace/fits/r.fits\n"; set ::ed(ok) 1}}] photometry multiband 1
    ::ogf::params::put photometry xm-catalog 2MASS; ::ogf::params::put photometry xm-radius 3.5
    compare crossmatch CatalogPanelCrossMatch [list apply {{} {set ::ed(xm,catalog) 2MASS; set ::ed(xm,radius) 3.5; set ::ed(ok) 1}}] photometry crossmatch 0
    # ---- photo-z / SED: the legacy dialogs are opened and their fields filled, then their Run proc is called
    ::ogf::params::put photoz_sed photoz-bands F435W,F606W,F814W; ::ogf::params::put photoz_sed photoz-mag-columns MAG_A,MAG_B,MAG_C
    set ::BANDS_PZ {F435W,F606W,F814W}
    compare photoz {CatalogPanelPhotoZ; set d .catphotoz; $d.bands delete 0 end; $d.bands insert 0 F435W,F606W,F814W; $d.mags delete 0 end; $d.mags insert 0 MAG_A,MAG_B,MAG_C; CatalogPanelPhotoZRun $d} {} photoz_sed photoz
    ::ogf::params::put photoz_sed sed-backend bagpipes; ::ogf::params::put photoz_sed sed-bands F435W,F606W; ::ogf::params::put photoz_sed sed-mag-columns MAG_A,MAG_B; ::ogf::params::put photoz_sed sed-photoz-column PHOTO_Z
    compare sed {CatalogPanelSEDFit; set d .catsedfit; $d.backend set bagpipes; $d.bands delete 0 end; $d.bands insert 0 F435W,F606W; $d.mags delete 0 end; $d.mags insert 0 MAG_A,MAG_B; $d.pzcol delete 0 end; $d.pzcol insert 0 PHOTO_Z; CatalogPanelSEDFitRun $d} {} photoz_sed sed
    # ---- deconvolution with a PSF file
    set psf [file join $::HOME fake_psf.fits]; set fd [open $psf w]; puts $fd x; close $fd
    set catpanel(psf,file) $psf; set catpanel(psf,has_psf) 1
    foreach alg {rl rl_tv wiener tikhonov clean mem rl_accelerated} {
	::ogf::params::put deconv algorithm $alg
	compare deconv_$alg [list CatalogPanelDeconvolve $alg] {} deconv deconvolve 0 none
    }
    # ---- ICL: background (no shared mask yet), profile on the raw image (no bgsub file), measurements (profile file prepared)
    ::ogf::cat::set icl,center_x 150.5; ::ogf::cat::set icl,center_y 120.0; ::ogf::cat::set icl,param,bkg-method chebyshev; ::ogf::cat::set icl,param,bkg-order 4
    compare icl_background {CatalogPanelICLBackground chebyshev} {} icl background 1 none
    compare icl_background_iter {::ogf::cat::set icl,param,bkg-iterative 1; CatalogPanelICLBackground chebyshev} {} icl background 1 none
    ::ogf::cat::set icl,param,bkg-iterative 0
    compare icl_profile CatalogPanelICLProfile {} icl profile 1 cat
    set d [file join $::HOME .ds9]; set fd [open [file join $d icl_profile_m51.tsv] w]; puts $fd "R\tSB"; close $fd; ::ogf::cat::set icl,has_profile 1
    compare icl_measure CatalogPanelICLMeasure {} icl measure 1 cat
    # ---- shared Auto Mask (legacy: CatalogPanelMaskAuto 0 with the ogfmask(p,*) values) and the LSBG full pipeline
    foreach {k v} {p,detect-thresh 4.0 p,minarea 7 p,expand-factor 1.8 p,max-dilate-radius 25 p,bright-star-mag-limit 17.5 p,bright-star-radius-scale 11.0 p,mag-threshold 23.5 p,lsb-protect 1} {set ::ogfmask($k) $v}
    foreach {n v} {auto-detect-thresh 4.0 auto-minarea 7 auto-expand-factor 1.8 auto-max-dilate-radius 25 auto-bright-star-mag-limit 17.5 auto-bright-star-radius-scale 11.0 auto-mag-threshold 23.5 auto-lsb-protect 1} {::ogf::params::put mask $n $v}
    R mask_param_keys [expr {$::ogfmask(p,detect-thresh) == 4.0 && [::ogf::params::get mask auto-lsb-protect] == 1}]
    compare mask_auto {CatalogPanelMaskAuto 0} {} mask auto 1 none 0
    set d [file join $::HOME .ds9]; file mkdir $d
    foreach f {mask_m51.fits mask_m51_bool.fits} {set fd [open [file join $d $f] w]; puts $fd x; close $fd}
    file copy -force /workspace/fits/m51.fits [file join $d lsbg_cleaned_m51.fits]; file copy -force /workspace/fits/m51.fits [file join $d lsbg_cleaned.fits]     ;# so that LoadFitsFile of the 'cleaned image' does not raise a dialog
    ::ogf::cat::set lsbg,param,svm-classify 1; ::ogf::cat::set lsbg,param,svm-checkpoint /nonexistent_ckpt.pt; ::ogf::cat::set lsbg,param,lsb-protect 0; ::ogf::cat::set lsbg,param,sersic-fit 0
    compare lsbg_run CatalogPanelLSBGRunAll {} lsbg run_all 1
    ::ogf::cat::set lsbg,param,svm-classify 0; ::ogf::cat::set lsbg,param,lsb-protect 1; ::ogf::cat::set lsbg,param,sersic-fit 1; ::ogf::cat::set lsbg,param,multiscale 0
    compare lsbg_run_defaultsvm CatalogPanelLSBGRunAll {} lsbg run_all 1
    # preconditions
    ::ogf::params::put photometry mb-bands {}
    R geometry_after [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    puts $::fh "SUMMARY failures=$::nf"; close $::fh; exit
}
proc OGFParamsDualSet {m} {::ogf::params::put extract dual-measure-image $m; return $m}
after 3000 run
