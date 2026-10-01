# Characterisation test for the catalog-panel store (item 5 of the decoupling work).  It records, for every migrated plugin feature,
#   * the argv that reaches exec (CLI args), with the python scripts replaced by canned TSV output,
#   * the catpanel keys afterwards (values reaching the store) and the status text,
#   * the contents of the .prf files written by the parameter save procs (HOME is a scratch dir) and what the load procs restore,
#   * the session-recorder steps (step, class, templated argv, post) the feature logged,
# and writes them to OGF_CAT_OUT.  scripts/verify_cat_behavior.sh diffs that file against scripts/golden/cat_behavior.golden, which
# was produced BEFORE any plugin was migrated to ::ogf::cat.   Usage (HOME must be a scratch directory!):
#   HOME=/tmp/ogf_cat_home OGF_CAT_OUT=/tmp/cat.txt DISPLAY=:77 bin/ds9 /workspace/fits/m51.fits -geometry 1300x950 -source scripts/verify_cat_behavior.tcl
global catpanel current ds9 ogfsess
set ::fh [open $::env(OGF_CAT_OUT) w]
set ::HOME [file normalize ~]
proc P {args} {puts $::fh [join $args { }]; flush $::fh}
proc bgerror {m} {P "BGERROR $m"}
rename tk_messageBox ::orig_mb
proc tk_messageBox {args} {P "MESSAGEBOX [dict get $args -message]"; return ok}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}

# ---- exec stub: record argv; python scripts return canned output, everything else runs for real
set ::EXEC {}
rename exec ::real_exec
proc exec {args} {
    set a {}
    foreach x $args {if {[string match 2>* $x]} continue; lappend a $x}
    lappend ::EXEC $a
    set py [expr {[llength $a] > 1 && [string match *.py [lindex $a 1]] ? [file tail [lindex $a 1]] : {}}]
    if {$py eq {}} {return [::real_exec {*}$args]}
    return [canned $py $a]
}
proc cat_numbers {a} {
    set i [lsearch -exact $a --catalog]
    if {$i < 0} {return {1 2 3}}
    set fd [open [lindex $a [expr {$i+1}]] r]; set d [read $fd]; close $fd
    set out {}
    foreach l [lrange [split $d \n] 1 end] {if {[string trim $l] ne {}} {lappend out [lindex [split $l \t] 0]}}
    return $out
}
proc tsv_cols {a cols} {
    set lines [list [join [linsert $cols 0 NUMBER] \t]]
    foreach n [cat_numbers $a] {
	set row [list $n]; set k 0
	foreach c $cols {lappend row [format %.3f [expr {($n % 7) + 0.25*[incr k]}]]}
	lappend lines [join $row \t]
    }
    return [join $lines \n]
}
proc canned {py a} {
    switch -- $py {
	ds9_photo_z.py {return [tsv_cols $a {PHOTO_Z PHOTO_Z_ERR PHOTO_Z_Q68 PHOTO_Z_OUTLIER}]}
	ds9_sed_fit.py {return [tsv_cols $a {LOG_MASS LOG_MASS_ERR LOG_AGE LOG_AGE_ERR LOG_Z AV SFR}]}
	ds9_bulge_disk.py {return [tsv_cols $a {BT_RATIO BULGE_RE BULGE_MAG DISK_RS DISK_MAG BD_CHI2 BD_FLAG}]}
	ds9_morphometry.py {return [tsv_cols $a {CONC ASYM GINI M20 R_PETRO}]}
	ds9_sersic.py {return [tsv_cols $a {SERSIC_N SERSIC_RE SERSIC_IE SERSIC_ELLIP SERSIC_THETA SERSIC_CHI2}]}
	ds9_psf_phot.py {return [tsv_cols $a {FLUX_PSF FLUXERR_PSF MAG_PSF MAGERR_PSF CHI2_PSF X_PSF Y_PSF}]}
	ds9_crowded_phot.py {return [tsv_cols $a {FLUX_CROWD FLUXERR_CROWD MAG_CROWD X_CROWD Y_CROWD N_NEIGHBORS}]}
	ds9_crossmatch.py {return [tsv_cols $a {MATCH_DIST MATCH_ID}]}
	ds9_galaxy_morph.py {
	    set lines [list "#GALAXY_MORPH\tN_CLASSIFIED=3\tN_SOURCES=3" "NUMBER\tT\tD\tC\ta\tb\tc\td\te\tf\tCOLOR"]
	    foreach n [lrange [cat_numbers $a] 0 2] {lappend lines [join [list $n E elliptical 0.9 x 1 x 1 x 1 green] \t]}
	    return [join $lines \n]
	}
	ds9_icl.py {
	    set m [lindex $a [expr {[lsearch -exact $a --mode]+1}]]
	    if {$m eq "background"} {return "OK"}
	    if {$m eq "profile"} {return "R_PIX\tSB\tSB_ERR\n1.0\t22.5\t0.1\n2.0\t23.5\t0.1"}
	    return "QUANTITY\tVALUE\nR_ISO\t12.5"
	}
	ds9_lsbg.py {
	    set m [lindex $a [expr {[lsearch -exact $a --mode]+1}]]
	    return "NUMBER\tX_IMAGE\tY_IMAGE\tMAG\tMODE\n1\t100\t100\t24.5\t$m\n2\t200\t150\t25.5\t$m\n3\t300\t90\t26.1\t$m"
	}
	ds9_separate.py {
	    return "#SEPARATE\tN_SUB=2\nNUMBER\tX_IMAGE\tY_IMAGE\tA_IMAGE\tB_IMAGE\tTHETA_IMAGE\tISO_RADIUS\tFLUX_AUTO\tMAG_AUTO\tFLAGS\n2.1\t41.5\t42.5\t3\t2\t10\t5\t900\t18.1\t0\n2.2\t45.5\t44.5\t3\t2\t10\t5\t800\t18.2\t0"
	}
	ds9_ai_merge.py {
	    return "#AI_MERGE\tN_GROUPS=2\tTHRESHOLD=0.70\nGROUP\tN_MEMBERS\tCONFIDENCE\tMEMBERS_X\tMEMBERS_Y\tMEMBERS_NUM\n0\t2\t0.91\t10.0,70.0\t10.0,70.0\t1,7\n1\t2\t0.80\t80.0,77.0\t80.0,88.0\t8,10"
	}
	ds9_add_source.py {
	    return "#ADD_SOURCE\tFOUND=1\nNUMBER\tX_IMAGE\tY_IMAGE\tA_IMAGE\tB_IMAGE\tTHETA_IMAGE\tELLIPTICITY\tKRON_RADIUS\tISO_RADIUS\tFLUX_AUTO\tFLUX_APER\tMAG_AUTO\tMAG_APER\tPEAK\tCLASS_STAR\tFLAGS\tFWHM_IMAGE\n9\t77.0\t88.0\t2\t1.5\t5\t0.25\t3.5\t4\t500\t450\t19.0\t19.1\t30\t0.9\t0\t3.0"
	}
	ds9_completeness.py {return "MAG\tFRAC\n20.0\t0.99\n21.0\t0.9"}
	ds9_multiband.py {return "NUMBER\tMAG_G\tMAG_R\n1\t20.1\t19.8"}
	default {return ""}
    }
}

# ---- recorder helpers
proc steps_sig {from} {
    global ogfsess
    set out {}
    foreach r [lrange $ogfsess(steps) $from end] {
	lappend out "STEP [dict get $r step] [dict get $r class] | [norm [dict get $r argv_t]] | post=[dict get $r post] title=[dict get $r title] requires=[dict get $r requires]"
    }
    return $out
}
proc norm {s} {
    set s [string map [list $::HOME <HOME> /workspace/fits <FITS> [file normalize bin] <BIN>] $s]
    regsub -all {/tmp/[A-Za-z0-9_.-]+} $s <TMP> s
    return $s
}
proc snap {globs} {
    global catpanel
    set out {}
    foreach g $globs {
	foreach k [lsort [array names catpanel $g]] {
	    set v $catpanel($k)
	    if {$k eq "alldata"} {set v "<tsv rows=[expr {[llength [split $v \n]]-1}] cols=[llength [split [lindex [split $v \n] 0] \t]]>"}
	    lappend out "  catpanel($k) = [norm $v]"
	}
    }
    return [join $out \n]
}
proc cols {} {global catpanel; return [lrange [split [lindex [split $catpanel(alldata) \n] 0] \t] 0 end]}
proc feature {name body {globs {status}}} {
    global ogfsess
    set s0 [llength $ogfsess(steps)]; set e0 [llength $::EXEC]
    P "=== $name"
    if {[catch {uplevel #0 $body} err]} {P "  ERROR: [norm $err]"}
    update; wait_idle 100
    foreach a [lrange $::EXEC $e0 end] {P "  EXEC [norm $a]"}
    foreach s [steps_sig $s0] {P "  $s"}
    set sn [snap $globs]; if {$sn ne {}} {P $sn}
}
proc prf {name} {
    set f [file join $::HOME .ds9 $name]
    if {![file exists $f]} {P "  PRF $name: (absent)"; return}
    set fd [open $f r]; set d [read $fd]; close $fd
    P "  PRF $name:"; foreach l [split [string trim $d] \n] {P "    $l"}
}

proc run {} {
    global catpanel current ds9 ogfsess
    update; wait_idle 500
    P "HOME-scratch [expr {[string match /tmp/* $::HOME] ? {yes} : {NO - refusing}}]"
    if {![string match /tmp/* $::HOME]} {P SUMMARY-ABORT; close $::fh; exit 2}
    # ---------------------------------------------------------------- extraction (real ds9_sextract)
    feature extract {set catpanel(param,detect-thresh) 3.0; CatalogPanelExtract; \
	set t0 [clock milliseconds]; while {$catpanel(alldata) eq {} && [clock milliseconds]-$t0 < 60000} {update; after 100}} \
	{alldata extract_param,* param,detect-thresh param,mag-zeropoint status}
    P "  columns: [llength [cols]] first=[lindex [cols] 0] n_rows=[expr {[llength [split $catpanel(alldata) \n]]-2}]"
    # ---------------------------------------------------------------- parameter stores: legacy write -> prf -> wipe -> load
    feature params_extract {
	set catpanel(param,detect-thresh) 2.5; set catpanel(param,back-size) 48; set catpanel(param,n-workers) 3
	CatalogPanelParamSave
    } {}
    prf sextract.prf
    feature params_extract_load {
	set catpanel(param,detect-thresh) 9; set catpanel(param,back-size) 9; set catpanel(param,n-workers) 9
	CatalogPanelParamLoad
    } {param,detect-thresh param,back-size param,n-workers}
    feature params_psf {
	set catpanel(psf,param,psf-size) 31; set catpanel(psf,param,rl-iterations) 25; set catpanel(psf,param,clean-gain) 0.2
	CatalogPanelPSFParamSave
    } {}
    prf psf_deconv.prf
    feature params_psf_load {
	set catpanel(psf,param,psf-size) 1; set catpanel(psf,param,rl-iterations) 1; CatalogPanelPSFParamLoad
    } {psf,param,psf-size psf,param,rl-iterations psf,param,clean-gain}
    feature params_bd {
	set catpanel(bd,param,max-sources) 17; set catpanel(bd,param,free-bulge-n) 1; set catpanel(bd,param,mag-zeropoint) 26.5; set catpanel(bd,param,pixel-scale) 0.05
	CatalogPanelBDParamSave
    } {}
    prf bulge_disk.prf
    feature params_bd_load {
	set catpanel(bd,param,max-sources) 1; CatalogPanelBDParamLoad
    } {bd,param,*}
    feature params_photoz {
	set catpanel(photoz,param,bands) g,r,i; set catpanel(photoz,param,mag-columns) MAG_APER_2; set catpanel(photoz,param,checkpoint) /nonexistent.pt
	CatalogPanelPhotoZParamSave
    } {}
    prf photo_z.prf
    feature params_photoz_load {
	set catpanel(photoz,param,bands) x; CatalogPanelPhotoZParamLoad
    } {photoz,param,*}
    feature params_sed {
	set catpanel(sed,param,backend) prospector; set catpanel(sed,param,bands) g,r; set catpanel(sed,param,photoz-column) PHOTO_Z
	CatalogPanelSEDParamSave
    } {}
    prf sed_fit.prf
    feature params_sed_load {
	set catpanel(sed,param,backend) x; CatalogPanelSEDParamLoad
    } {sed,param,*}
    # declarative-dialog path: ::ogf::params put/get/save/restore
    feature params_api {
	::ogf::params::put morphology max-sources 23
	::ogf::params::put morphology pixel-scale 0.11
	set ::pv [list [::ogf::params::get morphology max-sources] [::ogf::params::get morphology pixel-scale]]
	::ogf::params::save morphology
	::ogf::params::put deconv rl-iterations 44
	::ogf::params::save deconv
	::ogf::params::put extract detect-thresh 4.5
	::ogf::params::save extract
	P "  api values: $::pv"
    } {bd,param,max-sources bd,param,pixel-scale psf,param,rl-iterations param,detect-thresh}
    prf bulge_disk.prf; prf psf_deconv.prf; prf sextract.prf
    # restore sane values for the runs below
    set catpanel(param,detect-thresh) 3.0; set catpanel(param,n-workers) 2; set catpanel(param,mag-zeropoint) 24.5
    set catpanel(bd,param,max-sources) 50; set catpanel(bd,param,free-bulge-n) 0; set catpanel(bd,param,mag-zeropoint) 25.0; set catpanel(bd,param,pixel-scale) 0.05
    # ---------------------------------------------------------------- analysis steps (python canned)
    feature sersic {CatalogPanelSersicFit} {status param,mag-zeropoint}
    P "  cols+: [lrange [cols] end-5 end]"
    feature morphometry {CatalogPanelMorphometry}
    P "  cols+: [lrange [cols] end-4 end]"
    feature bulge_disk {CatalogPanelBulgeDisk}
    P "  cols+: [lrange [cols] end-6 end]"
    feature photoz_dialog_run {
	set catpanel(photoz,param,bands) g,r,i; set catpanel(photoz,param,mag-columns) {}
	CatalogPanelPhotoZ
	.catphotoz.bands delete 0 end; .catphotoz.bands insert 0 "u,g"
	CatalogPanelPhotoZRun .catphotoz
    } {status photoz,param,*}
    P "  cols+: [lrange [cols] end-3 end]"
    P "  photoz dialog destroyed: [expr {![winfo exists .catphotoz]}]"
    prf photo_z.prf
    feature sed_dialog_run {
	CatalogPanelSEDFit
	set w .catsedfit
	if {[winfo exists $w]} {
	    # the run proc is the "Run" button's command; find it
	    set cmd [$w.btns.run cget -command]
	    P "  run command: [norm $cmd]"
	    uplevel #0 $cmd
	}
    } {status sed,param,*}
    P "  cols+: [lrange [cols] end-6 end]"
    # PSF-dependent steps need a PSF file
    set psf [file join $::HOME fake_psf.fits]; set fd [open $psf w]; puts $fd x; close $fd
    set catpanel(psf,file) $psf; set catpanel(psf,has_psf) 1
    feature psf_photometry {CatalogPanelPSFPhotometry}
    P "  cols+: [lrange [cols] end-6 end]"
    feature crowded_phot {CatalogPanelCrowdedPhot}
    P "  cols+: [lrange [cols] end-5 end]"
    foreach alg {rl wiener clean} {
	feature deconvolve_$alg [list CatalogPanelDeconvolve $alg] {status psf,param,rl-iterations psf,param,wiener-nsr psf,param,clean-gain}
    }
    feature crossmatch_dialog {after 400 {set ed(ok) 1}; CatalogPanelCrossMatch}
    P "  cols+: [lrange [cols] end-1 end]"
    feature completeness_dialog {after 400 {set ed(ok) 1}; CatalogPanelCompleteness}
    # galaxy morphology (CNN) result parsing -> morph keys, new columns, marker colours
    feature galaxy_morphology {CatalogPanelGalaxyMorphology} {status morph,* markall,on}
    P "  cols+: [lrange [cols] end-1 end]"
    # mask: presets read the icl / lsbg parameters of the store
    feature mask_presets {
	set catpanel(icl,param,detect-thresh) 2.2; set catpanel(icl,param,expand-factor) 1.7; set catpanel(icl,param,max-dilate-radius) 40
	set catpanel(icl,param,bright-star-mag-limit) 17; set catpanel(icl,param,bright-star-radius-scale) 2.5
	OGFMaskInit; OGFMaskPresetICL
	set s1 "[array get ::ogfmask p,*]"
	set catpanel(lsbg,param,mask-detect-thresh) 1.1; set catpanel(lsbg,param,mask-detect-minarea) 7; set catpanel(lsbg,param,mask-expand-factor) 1.3
	set catpanel(lsbg,param,max-dilate-radius) 33; set catpanel(lsbg,param,bright-star-mag-limit) 16; set catpanel(lsbg,param,bright-star-radius-scale) 2.1
	set catpanel(lsbg,param,mask-mag-threshold) 21; set catpanel(lsbg,param,lsb-protect) 1
	OGFMaskPresetLSBG
	set s2 "[array get ::ogfmask p,*]"
	P "  mask preset icl : [lsort -stride 2 $s1]"
	P "  mask preset lsbg: [lsort -stride 2 $s2]"
    } {}
    feature mask_auto_args {
	set catpanel(lsbg,param,pixel-scale) 0.2; set catpanel(lsbg,param,lsb-mu-threshold) 26.5
	set a [OGFMaskAutoArgs]
	P "  auto args: [norm $a]"
	set catpanel(icl,cmdlog) {}
	OGFMaskSyncPipelines [dict create bool /tmp/zz_mask.fits base /tmp/zz]
	P "  sync keys: [lsort [array names catpanel *mask*]]"
    } {icl,mask_file icl,has_mask lsbg,mask_file lsbg,has_mask icl,fits_base_mask lsbg,fits_base_mask}
    # bands (the dialog-free procs)
    feature bands_register {
	OGFBandsInit
	CatalogPanelBandsRegisterFrame $current(frame) gband 25.5 {} [CatalogPanelGetFITS]
	P "  bands: $::ogfband(names) detect=$::ogfband(detect)"
    } {status param,mag-zeropoint param,pixel-scale}
    feature bands_set_detect {CatalogPanelBandsSetDetect gband} {status param,mag-zeropoint param,pixel-scale}
    # segmentation map (output loaded into a new frame)
    feature segmap {
	set ::segout [file join $::HOME segout.fits]; file copy -force /workspace/fits/m51.fits $::segout
	rename canned ::canned0
	proc canned {py a} {if {$py eq "ds9_segmap.py"} {return "OK 5 $::segout"}; return [::canned0 $py $a]}
	CatalogPanelSegmentationMap
	rename canned {}; rename ::canned0 canned
    } {status param,detect-thresh}
    # selection-dependent reads
    feature selection_reads {
	set catpanel(sel,nums) {3 5 8}
	P "  AI row count=[OGFAIRowCount] selection=$catpanel(sel,nums)"
	set catpanel(sel,nums) {}
    } {}
    # ---------------------------------------------------------------- ICL / LSBG (parameter files, derived file keys, exec args, state keys)
    set catpanel(icl,param,rmin) 3; set catpanel(icl,param,rmax) 150
    feature icl_params_save {CatalogPanelICLParamSave} {}
    prf icl.prf
    feature icl_params_load {set catpanel(icl,param,rmax) 1; set catpanel(icl,param,bkg-order) 9; CatalogPanelICLParamLoad} {icl,param,rmax icl,param,bkg-order icl,param,rmin}
    feature icl_update_files {set catpanel(icl,fits_base) {}; CatalogPanelICLUpdateFiles [CatalogPanelGetFITS]} {icl,fits_base icl,mask_file icl,masked_file icl,bkg_file icl,bgsub_file icl,profile_file icl,has_mask icl,has_bkg icl,has_profile}
    feature icl_background {CatalogPanelICLBackground median} {status icl,has_bkg icl,cmdlog icl,bkg_file icl,bgsub_file}
    feature icl_profile {set catpanel(icl,center_x) 120.5; set catpanel(icl,center_y) 130.5; CatalogPanelICLProfile} {status icl,has_profile icl,center_x icl,center_y icl,cmdlog}
    set fd [open $catpanel(icl,profile_file) w]; puts $fd x; close $fd; set catpanel(icl,has_profile) 1
    feature icl_measure {CatalogPanelICLMeasure} {status icl,cmdlog}
    feature icl_annuli {CatalogPanelICLDrawAnnuli} {}
    set catpanel(lsbg,param,detect-thresh) 2.1; set catpanel(lsbg,param,pixel-scale) 0.2
    feature lsbg_params_save {CatalogPanelLSBGParamSave} {}
    prf lsbg.prf
    feature lsbg_params_load {set catpanel(lsbg,param,detect-thresh) 1; CatalogPanelLSBGParamLoad} {lsbg,param,detect-thresh lsbg,param,pixel-scale}
    feature lsbg_update_files {set catpanel(lsbg,fits_base) {}; CatalogPanelLSBGUpdateFiles [CatalogPanelGetFITS]} {lsbg,fits_base lsbg,mask_file lsbg,masked_file lsbg,bkg_file lsbg,cleaned_file lsbg,segmap_file lsbg,catalog_file lsbg,has_*}
    # create the intermediate files the later LSBG steps require
    foreach f [list $catpanel(lsbg,cleaned_file) $catpanel(lsbg,segmap_file) $catpanel(lsbg,mask_file)] {set fd [open $f w]; puts $fd x; close $fd}
    feature lsbg_detect {CatalogPanelLSBGDetect} {status lsbg,has_detect lsbg,detect_data lsbg,cmdlog alldata}
    # Detect opened a new frame (CreateFrame resets the per-frame panel state): go back to the image frame and rebuild the file keys
    GotoFrame [lindex $ds9(frames) 0]; wait_idle 300
    set catpanel(lsbg,fits_base) {}; CatalogPanelLSBGUpdateFiles [CatalogPanelGetFITS]
    foreach f [list $catpanel(lsbg,cleaned_file) $catpanel(lsbg,segmap_file)] {set fd [open $f w]; puts $fd x; close $fd}
    feature lsbg_photometry {CatalogPanelLSBGPhotometry} {status lsbg,has_catalog lsbg,cmdlog alldata}
    feature lsbg_filter {CatalogPanelLSBGFilter} {status lsbg,cmdlog alldata}
    # ---------------------------------------------------------------- CLI script export / import (cmdlog -> .sh -> replay)
    proc script_roundtrip {pipeline} {
	set f [file join $::HOME ${pipeline}_pipeline.sh]
	file delete $f
	set ::cli_file $f
	foreach c {tk_getSaveFile tk_getOpenFile tk_messageBox} {rename $c ::cli_orig_$c}
	proc tk_getSaveFile {args} {return $::cli_file}
	proc tk_getOpenFile {args} {return $::cli_file}
	proc tk_messageBox {args} {return no}
	catch {CatalogPanelExportCLIScript $pipeline} e1
	if {[file exists $f]} {
	    set fd [open $f r]; set d [read $fd]; close $fd
	    P "  script lines=[llength [split $d \n]] bytes=[string length $d]"
	    foreach l [split $d \n] {if {[string match "python3 *" $l] || [regexp {^[A-Z_]+=} $l]} {P "    [norm $l]"}}
	}
	catch {CatalogPanelImportCLIScript $pipeline} e2
	foreach c {tk_getSaveFile tk_getOpenFile tk_messageBox} {rename $c {}; rename ::cli_orig_$c $c}
	if {$e1 ne {}} {P "  export error: [norm $e1]"}
	if {$e2 ne {}} {P "  import error: [norm $e2]"}
    }
    feature cli_export_import_icl {script_roundtrip icl} {status icl,cmdlog alldata}
    feature cli_export_import_lsbg {script_roundtrip lsbg} {status lsbg,cmdlog alldata}
    # ---------------------------------------------------------------- objects: merge / separate / delete / add / AI-merge / save+load
    set synth "NUMBER\tX_IMAGE\tY_IMAGE\tA_IMAGE\tB_IMAGE\tTHETA_IMAGE\tISO_RADIUS\tFLUX_AUTO\tMAG_AUTO\tNPIX_ISO\tFLAGS"
    for {set i 1} {$i <= 8} {incr i} {append synth "\n$i\t[expr {$i*10.0}]\t[expr {$i*10.0}]\t3.0\t2.0\t[expr {$i*5.0}]\t6.0\t[expr {1000.0*$i}]\t[format %.2f [expr {20.0-$i*0.1}]]\t[expr {20+$i}]\t0"}
    CatalogPanelLoadTSV $synth synth
    update; wait_idle 200
    set catpanel(markall,on) 0
    proc select_row {num} {
	global catpanel
	set db $catpanel(tbldb); global $db
	for {set r 1} {$r < [$catpanel(tbl) cget -rows]} {incr r} {
	    if {[info exists ${db}($r,1)] && [set ${db}($r,1)] eq $num} {$catpanel(tbl) selection clear all; $catpanel(tbl) selection set $r,1; return $r}
	}
	return -1
    }
    feature objects_selected_source {
	P "  row=[select_row 4]"
	P "  src=[lsort -stride 2 [CatalogPanelGetSelectedSource]]"
    } {}
    feature objects_merge {
	set catpanel(merge,list) {2 3 5}; set catpanel(merge,active) 1
	CatalogPanelMergeSources
	P "  merged rows: [lrange [split [::ogf::cat::tsv] \n] end end]"
    } {status merge,active merge,list alldata}
    feature objects_merge_cancel {
	set catpanel(merge,list) {1 2}; set catpanel(merge,active) 1
	CatalogPanelMergeCancel
    } {status merge,active merge,list}
    feature objects_escape_key {
        set catpanel(merge,active) 1; set catpanel(merge,list) {1}
	CatalogPanelEscapeKey
	set catpanel(ai,active) 0
	set catpanel(sel,nums) {1}
	CatalogPanelEscapeKey
    } {status merge,active merge,list sel,nums ai,active}
    feature objects_delete {
	P "  row=[select_row 6]"
	CatalogPanelDeleteSelected
	P "  numbers: [join [lmap l [lrange [split [::ogf::cat::tsv] \n] 1 end] {lindex [split $l \t] 0}] ,]"
    } {status alldata}
    feature objects_separate {
	P "  row=[select_row 4]"
	CatalogPanelSeparateSelected
	P "  numbers: [join [lmap l [lrange [split [::ogf::cat::tsv] \n] 1 end] {lindex [split $l \t] 0}] ,]"
    } {status alldata}
    set catpanel(add_objects_mode) 0
    feature objects_add_disabled {CatalogPanelAddObjectAtPosition $current(frame) 77 88} {status add_objects_mode}
    set catpanel(add_objects_mode) 1
    feature objects_add {
	CatalogPanelAddObjectAtPosition $current(frame) 77 88
	P "  last row: [lindex [split [::ogf::cat::tsv] \n] end]"
    } {status alldata add_objects_mode param,detect-thresh}
    set ::tsvf [file join $::HOME saved_cat.tsv]
    feature objects_save_load {
	rename tk_getSaveFile ::orig_gsf; proc tk_getSaveFile {args} {return $::tsvf}
	CatalogPanelSeparateSave
	rename tk_getSaveFile {}; rename ::orig_gsf tk_getSaveFile
	set fd [open $::tsvf r]; set d [read $fd]; close $fd
	P "  saved bytes=[string length $d] lines=[llength [split $d \n]]"
	rename tk_getOpenFile ::orig_gof; proc tk_getOpenFile {args} {return $::tsvf}
	CatalogPanelSeparateLoad
	rename tk_getOpenFile {}; rename ::orig_gof tk_getOpenFile
    } {status alldata}
    feature ai_merge_run {CatalogPanelAIMerge; P "  keys bound: [bind . <Key-n>]"} {status ai,groups ai,total ai,current ai,active extract_param,detect-thresh}
    feature ai_merge_next {CatalogPanelAINext; CatalogPanelAINext} {status ai,current ai,total}
    feature ai_merge_prev {CatalogPanelAIPrev; CatalogPanelAIPrev} {status ai,current}
    feature ai_merge_reject {CatalogPanelAIReject} {status ai,groups ai,total ai,current}
    feature ai_merge_accept {CatalogPanelAIAccept; P "  rows after accept: [llength [split [::ogf::cat::tsv] \n]]"} {status ai,groups ai,total ai,current ai,active merge,list merge,active}
    feature ai_merge_done {CatalogPanelAIDone; P "  keys unbound: [expr {[bind . <Key-n>] eq {}}]"} {status ai,groups ai,total ai,active}
    feature markers_mark_all {CatalogPanelMarkAll} {status markall,on}
    feature markers_clear {CatalogPanelClearMarkers} {status markall,on}
    feature markers_ctrl_select {CatalogPanelCtrlSelect 2; CatalogPanelCtrlSelect 3; CatalogPanelCtrlSelect 2} {status merge,list merge,active}
    feature markers_show_visible {set catpanel(visible_mode) 1; CatalogPanelShowVisible; set catpanel(visible_mode) 0; CatalogPanelShowVisible} {status visible_mode}
    # ---------------------------------------------------------------- catalog core: table fill, filter, sort, select, save/load, export, frame state
    proc tbl_rows {} {
	update; global catpanel; set db $catpanel(tbldb); global $db
	set out {}
	for {set r 1} {$r < [$catpanel(tbl) cget -rows]} {incr r} {
	    if {[info exists ${db}($r,1)]} {lappend out [set ${db}($r,1)]}
	}
	return [join $out ,]
    }
    proc fileio {save open body} {
	foreach c {tk_getSaveFile tk_getOpenFile} {rename $c ::cat_orig_$c}
	set ::cat_save $save; set ::cat_open $open
	proc tk_getSaveFile {args} {return $::cat_save}
	proc tk_getOpenFile {args} {return $::cat_open}
	set rc [catch {uplevel #0 $body} e]
	foreach c {tk_getSaveFile tk_getOpenFile} {rename $c {}; rename ::cat_orig_$c $c}
	if {$rc} {P "  ERROR: [norm $e]"}
    }
    proc dumpfile {f} {
	if {![file exists $f]} {P "  FILE [norm $f]: absent"; return}
	set fd [open $f r]; set d [read $fd]; close $fd
	P "  FILE [file tail $f]: [string length $d] bytes"; foreach l [lrange [split $d \n] 0 3] {P "    [norm [string range $l 0 150]]"}
    }
    set ::synth $synth
    CatalogPanelLoadTSV $synth synth
    feature catalog_load_tsv {P "  rows: [tbl_rows]"} {status alldata delim}
    feature catalog_filter {set catpanel(search_var) 3; CatalogPanelFilter; P "  rows: [tbl_rows]"; set catpanel(search_var) {}; CatalogPanelFilter; P "  rows: [tbl_rows]"} {status search_var}
    feature catalog_sort {CatalogPanelSort MAG_AUTO descending; P "  rows: [tbl_rows]"; CatalogPanelSort NUMBER ascending; P "  rows: [tbl_rows]"} {status sort,* alldata}
    feature catalog_header_click {
	set tb $catpanel(tbl); set bb [$tb bbox 0,2]
	CatalogPanelTableClick [expr {[lindex $bb 0]+2}] [expr {[lindex $bb 1]+2}]
	CatalogPanelTableClick [expr {[lindex $bb 0]+2}] [expr {[lindex $bb 1]+2}]
	P "  rows: [tbl_rows]"
    } {status sort,*}
    feature catalog_select_row {
	CatalogPanelSort NUMBER ascending
	CatalogPanelUpdateSelInfo 3
	CatalogPanelSelectCmd 1,1 3,1
	after 300 {set ::sel_done 1}; vwait ::sel_done
    } {status sel,text sel,nums}
    feature catalog_save_tsv {fileio [file join $::HOME c1.tsv] {} {CatalogPanelSaveCatalog}; dumpfile [file join $::HOME c1.tsv]} {status}
    feature catalog_save_csv {fileio [file join $::HOME c1.csv] {} {CatalogPanelSaveCatalog}; dumpfile [file join $::HOME c1.csv]} {status}
    feature catalog_load_csv {fileio {} [file join $::HOME c1.csv] {CatalogPanelLoadCatalog}; P "  rows: [tbl_rows]"} {status alldata}
    feature catalog_load_missing {fileio {} [file join $::HOME nonexist.tsv] {CatalogPanelLoadCatalog}} {status}
    feature catalog_export_regions {fileio [file join $::HOME c1.reg] {} {CatalogPanelExportRegions}; dumpfile [file join $::HOME c1.reg]} {status}
    feature catalog_export_fits {fileio [file join $::HOME c1.fits] {} {CatalogPanelExportFITS}} {status}
    feature catalog_temp_catalog {P "  tmp: [norm [CatalogPanelSaveTempCatalog unit]]"; dumpfile [file join $::HOME .ds9 unit_catalog.tsv]} {}
    feature catalog_add_columns {
	CatalogPanelAddColumnsFromTSV "NUMBER\tNEWCOL\tOTHER\n1\t1.5\tx\n3\t2.5\ty\n99\t0\tz" {NEWCOL OTHER}
	P "  header: [lindex [split [::ogf::cat::tsv] \n] 0]"
	P "  row1: [lindex [split [::ogf::cat::tsv] \n] 1]"
    } {status alldata}
    feature catalog_frame_state {
	set f0 [lindex $ds9(frames) 0]
	set catpanel(search_var) keepme; set catpanel(sort,col) NUMBER
	CatalogPanelSaveFrameState $f0
	P "  saved: [lsort [array names ::catpanel_fdata]]"
	set catpanel(alldata) {}; set catpanel(search_var) {}; set catpanel(sort,col) {}
	CatalogPanelRestoreFrameState $f0
	P "  restored: search=$catpanel(search_var) sort=$catpanel(sort,col) rows=[tbl_rows]"
	CatalogPanelDeleteFrameState $f0
	P "  after delete: [llength [array names ::catpanel_fdata]] keys"
    } {status search_var sort,col alldata}
    feature catalog_clear {CatalogPanelClear; P "  rows: [tbl_rows]"} {status alldata filename}
    feature catalog_clear_all {CatalogPanelLoadTSV $synth synth; set catpanel(merge,list) {1}; set catpanel(ai,active) 0; CatalogPanelClearAll; P "  rows: [tbl_rows]"} {status alldata filename sort,* visible_mode markall,on add_objects_mode trim,active merge,* ai,* psf,has_psf icl,has_* lsbg,has_* search_var}
    feature moving_status {OGFMovStatus "linking 3 of 7"} {status}
    P "SUMMARY-DONE steps=[llength $ogfsess(steps)] exec=[llength $::EXEC]"
    close $::fh
    exit
}
after 4000 {if {[catch run err]} {P "ERROR $err\n$::errorInfo"; close $::fh; exit 1}}
