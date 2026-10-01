# Scripted ds9 GUI session for scripts/verify_session_replay.py.
# Runs under  ds9 <fits...> -source verify_gui_session.tcl  with env:
#   OGF_VERIFY_KIND = hudf | m51 ; OGF_VERIFY_OUT = directory for the GUI outputs
# Every step calls the proc the corresponding menu entry calls (no separate code path).
global current catpanel ogfband ogfmask
set kind $::env(OGF_VERIFY_KIND)
set out  $::env(OGF_VERIFY_OUT)
file mkdir $out
set ::vlog [open [file join $out gui_steps.log] w]
proc vlog {m} {puts $::vlog $m; flush $::vlog; puts "VERIFY: $m"}
proc wait_until {cond {ms 600000}} {
    set t0 [clock milliseconds]
    while {![uplevel #0 $cond]} {
	update
	after 50
	if {[clock milliseconds] - $t0 > $ms} {error "timeout waiting for $cond"}
    }
}
proc nrows {} {global catpanel; set n 0; foreach l [lrange [split $catpanel(alldata) \n] 1 end] {if {[string trim $l] ne {}} {incr n}}; return $n}

proc gui_session {} {
    global current catpanel ogfband ogfmask ds9 kind out
    set t0 [clock milliseconds]
    set fits [lsort [glob -nocomplain $::env(OGF_VERIFY_FITS)]]
    # frames: one per FITS
    set files $::env(OGF_VERIFY_FILES)
    set bandfile {}
    # --- load images (first one is already loaded by the command line)
    foreach f [lrange $files 1 end] {
	CreateFrame
	LoadFitsFile $f {} {}
    }
    update
    # --- extraction on the detection image
    if {$kind eq "hudf"} {
	set detidx [lsearch -glob $files *f160w*]
	GotoFrame [lindex $ds9(frames) $detidx]
	update
	set catpanel(param,detect-thresh) 3.0
	set catpanel(param,mag-zeropoint) 25.9462
	set catpanel(param,pixel-scale) 0.06
    } else {
	set catpanel(param,detect-thresh) 2.0
    }
    set cur [CatalogPanelGetFITS]
    vlog "detection image $cur"
    set ts [clock milliseconds]
    CatalogPanelExtract
    vlog "extract: [nrows] rows, [expr {[clock milliseconds]-$ts}] ms"
    if {$kind eq "hudf"} {
	# --- register 3 bands (the procs behind Bands > Register Current Frame as Band...)
	foreach f $files {
	    set fr {}
	    foreach cand $ds9(frames) {
		if {[string match *[file tail $f] [string trim [$cand get fits file name full] "{}"]]} {set fr $cand}
	    }
	    set info [OGFImageInfo $f]
	    CatalogPanelBandsRegisterFrame $fr [dict get $info FILTER] [dict get $info ZP_AB] {} $f
	}
	CatalogPanelBandsSetDetect F160W
	vlog "bands: $ogfband(names) detect=$ogfband(detect)"
	set ts [clock milliseconds]
	CatalogPanelBandsMeasure 2.0
	vlog "forced photometry: [expr {[clock milliseconds]-$ts}] ms, cols=[llength [split [lindex [split $catpanel(alldata) \n] 0] \t]]"
	GotoFrame $ogfband(F160W,frame)
	update
    }
    # --- auto mask (current frame = detection band)
    set ogfmask(p,detect-thresh) 3.0
    if {$kind eq "m51"} {set ogfmask(p,bright-star-mag-limit) 5.0}
    set ts [clock milliseconds]
    CatalogPanelMaskAuto 0
    vlog "auto mask: $ogfmask(stats) [expr {[clock milliseconds]-$ts}] ms"
    # --- hand drawn regions: one circle (add), one box (erase)
    set fr $current(frame)
    if {$kind eq "hudf"} {
	set ::reg "image\ncircle(1800,1800,150)\n"
	set ::reg2 "image\nbox(1500,1500,200,120,30)\n"
    } else {
	set ::reg "image\ncircle(300,300,40)\n"
	set ::reg2 "image\nbox(150,150,60,40,20)\n"
    }
    $fr marker delete all
    $fr marker command ds9 var ::reg
    CatalogPanelMaskRegions 0
    vlog "mask add circle: $ogfmask(stats)"
    $fr marker delete all
    $fr marker command ds9 var ::reg2
    CatalogPanelMaskRegions 1
    vlog "mask erase box: $ogfmask(stats)"
    $fr marker delete all
    if {$kind eq "hudf"} {
	CatalogPanelMaskCopyToBands
	vlog "mask copied to bands: $catpanel(status)"
    }
    # --- catalog filter (Trim) + save
    if {$kind eq "m51"} {
	set catpanel(trim,MAG_AUTO,min) {}
	set catpanel(trim,MAG_AUTO,max) 14
    } else {
	set catpanel(trim,MAG_AUTO,min) 18
	set catpanel(trim,MAG_AUTO,max) 29
    }
    set catpanel(trim,active) 1
    CatalogPanelTrimApply
    vlog "trim: [nrows] rows"
    CatalogPanelSaveCatalogTo [file join $out catalog_gui.tsv]
    vlog "catalog saved"
    # extra python back ends on the small test
    if {$kind eq "m51"} {
	CatalogPanelSersicFit
	vlog "sersic: [nrows] rows cols=[llength [split [lindex [split $catpanel(alldata) \n] 0] \t]]"
	CatalogPanelMorphometry
	vlog "morphometry done"
	CatalogPanelSaveCatalogTo [file join $out catalog_gui_final.tsv]
	# --- ICL (background -> BCG centre -> profile -> measure); the centre is picked as the
	#     brightest catalogue source, exactly what the Set BCG Center dialog proposes
	set ts [clock milliseconds]
	CatalogPanelICLBackground polynomial
	vlog "ICL background: $catpanel(status) [expr {[clock milliseconds]-$ts}] ms"
	set ::icl_bcg_entry_id 1
	# like the Set BCG Center dialog: catalogue of the FIRST frame (ICL Background opened a new frame)
	global catpanel_fdata
	set ff [lindex $ds9(frames) 0]
	set cd $catpanel_fdata($ff,alldata)
	set hdr [split [lindex [split $cd \n] 0] \t]
	CatalogPanelICLSetCenterByID .none $cd [lsearch $hdr NUMBER] [lsearch $hdr X_IMAGE] [lsearch $hdr Y_IMAGE]
	vlog "ICL centre: $catpanel(icl,center_x) $catpanel(icl,center_y) ($catpanel(status))"
	set catpanel(icl,param,rmax) 150
	CatalogPanelICLProfile
	vlog "ICL profile: $catpanel(status)"
	CatalogPanelSaveCatalogTo [file join $out icl_profile_gui.tsv]
	CatalogPanelICLMeasure
	vlog "ICL measure: $catpanel(status)"
	CatalogPanelSaveCatalogTo [file join $out icl_measure_gui.tsv]
    }
    # --- copy GUI data products
    foreach g [glob -nocomplain ~/.ds9/mask_*.fits] {file copy -force $g [file join $out gui_[file tail $g]]}
    # session export
    set sp [file join $out ogfinder_session.py]
    CatalogPanelSessionSave $sp
    set fd [open [file join $out gui_log.txt] w]
    puts $fd [OGFSessLogText]
    close $fd
    vlog "session exported: $sp  total [expr {([clock milliseconds]-$t0)/1000.0}] s"
}
after 3000 {
    if {[catch gui_session err]} {
	vlog "ERROR: $err\n$::errorInfo"
	exit 3
    }
    vlog DONE
    exit 0
}
