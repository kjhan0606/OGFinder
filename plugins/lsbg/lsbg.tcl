# lsbg/lsbg.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/lsbg/plugin.json.

proc CatalogPanelLSBGUpdateFiles {fn} {
    global catpanel
    if {$fn eq {}} return
    set base [CatalogPanelFitsBaseName $fn]
    if {$base eq {}} return
    # Skip if already set for this base
    if {[info exists catpanel(lsbg,fits_base)] &&
	$catpanel(lsbg,fits_base) eq $base} return
    set catpanel(lsbg,fits_base) $base
    set ds9dir [file join [file normalize ~] .ds9]
    set catpanel(lsbg,mask_file)    [file join $ds9dir "mask_${base}_bool.fits"]
    set catpanel(lsbg,masked_file)  [file join $ds9dir "mask_${base}_masked.fits"]
    set catpanel(lsbg,bkg_file)     [file join $ds9dir "lsbg_background_${base}.fits"]
    set catpanel(lsbg,cleaned_file) [file join $ds9dir "lsbg_cleaned_${base}.fits"]
    set catpanel(lsbg,segmap_file)  [file join $ds9dir "lsbg_segmap_${base}.fits"]
    set catpanel(lsbg,catalog_file) [file join $ds9dir "lsbg_catalog_${base}.tsv"]
    # Detect if previous results exist for this FITS
    set catpanel(lsbg,has_mask)    [file exists $catpanel(lsbg,mask_file)]
    set catpanel(lsbg,has_clean)   [file exists $catpanel(lsbg,cleaned_file)]
    set catpanel(lsbg,has_detect)  [file exists $catpanel(lsbg,segmap_file)]
    set catpanel(lsbg,has_catalog) [file exists $catpanel(lsbg,catalog_file)]
}

proc CatalogPanelLSBGParamLoad {} {
    global catpanel

    set preffile [file join [file normalize ~] .ds9 lsbg.prf]
    if {![file exists $preffile]} return
    if {[catch {set fd [open $preffile r]} err]} return
    while {[gets $fd line] >= 0} {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set parts [split $line]
	if {[llength $parts] >= 2} {
	    set key [lindex $parts 0]
	    set val [lindex $parts 1]
	    if {[info exists catpanel(lsbg,param,$key)]} {
		set catpanel(lsbg,param,$key) $val
	    }
	}
    }
    close $fd
}

proc CatalogPanelLSBGParamSave {} {
    global catpanel

    set prefdir [file join [file normalize ~] .ds9]
    if {![file isdirectory $prefdir]} {
	file mkdir $prefdir
    }
    set preffile [file join $prefdir lsbg.prf]
    if {[catch {set fd [open $preffile w]} err]} return
    foreach pname {mask-detect-thresh mask-detect-minarea mask-expand-factor \
		   max-dilate-radius \
		   bright-star-mag-limit bright-star-radius-scale \
		   mask-mag-threshold interp-method \
		   lsb-protect lsb-mu-threshold \
		   bkg-method bkg-mesh-size bkg-poly-order \
		   bkg-sigma-clip bkg-n-iterations bkg-refine-thresh \
		   bkg-rms-quantile bkg-convergence-tol \
		   detect-thresh detect-minarea detect-filter-kernel \
		   deblend-nthresh deblend-mincont \
		   multiscale multiscale-factors \
		   sersic-fit sersic-n-min sersic-n-max sersic-re-min \
		   sersic-cutout-scale sersic-max-nfev \
		   phot-apertures mag-zeropoint pixel-scale \
		   mu-eff-min mu-eff-max r-eff-min r-eff-max \
		   ellipticity-max min-snr \
		   sersic-n-filter-min sersic-n-filter-max sersic-chi2-max \
		   svm-classify svm-threshold svm-checkpoint} {
	puts $fd "$pname $catpanel(lsbg,param,$pname)"
    }
    close $fd
}

proc CatalogPanelLSBGMask {} {
    OGFMaskPipelineMask lsbg
}

proc CatalogPanelLSBGViewMask {} {
    global ogfmask
    set ogfmask(overlay) 1
    CatalogPanelMaskToggleOverlay
}

proc CatalogPanelLSBGSaveMask {} {
    CatalogPanelMaskSaveAs
}

proc CatalogPanelLSBGImportMask {} {
    CatalogPanelMaskImport
}

proc CatalogPanelLSBGClean {method} {
    global catpanel

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    if {![OGFMaskEnsure lsbg]} {
	set catpanel(status) "LSBG: could not obtain a mask - run Mask > Auto Mask first"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Iterative background ($method)..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode clean \
	--mask $catpanel(lsbg,mask_file) \
	--bkg-method $method \
	--bkg-mesh-size $catpanel(lsbg,param,bkg-mesh-size) \
	--bkg-poly-order $catpanel(lsbg,param,bkg-poly-order) \
	--bkg-sigma-clip $catpanel(lsbg,param,bkg-sigma-clip) \
	--bkg-n-iterations $catpanel(lsbg,param,bkg-n-iterations) \
	--bkg-refine-thresh $catpanel(lsbg,param,bkg-refine-thresh) \
	--bkg-rms-quantile $catpanel(lsbg,param,bkg-rms-quantile) \
	--bkg-convergence-tol $catpanel(lsbg,param,bkg-convergence-tol) \
	--interp-method $catpanel(lsbg,param,interp-method) \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--bkg-output $catpanel(lsbg,bkg_file) \
	--cleaned-output $catpanel(lsbg,cleaned_file) \
	--n-workers $catpanel(param,n-workers)]

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG clean error: $err"
	return
    }

    set catpanel(lsbg,has_clean) 1

    # Auto-display cleaned image in new frame
    CreateFrame
    if {[catch {LoadFitsFile $catpanel(lsbg,cleaned_file) {} {}} err]} {
	set catpanel(status) "LSBG: Clean done but could not display: $err"
    } else {
	global scale
	set scale(mode) zscale
	ChangeScaleMode
    }

    set catpanel(status) "LSBG: Background cleaned ($method) (new frame)"
}

proc CatalogPanelLSBGViewClean {} {
    global catpanel

    if {![file exists $catpanel(lsbg,cleaned_file)]} {
	set catpanel(status) "LSBG: No cleaned image — run Background Model first"
	return
    }

    CreateFrame
    if {[catch {LoadFitsFile $catpanel(lsbg,cleaned_file) {} {}} err]} {
	set catpanel(status) "LSBG: Error loading cleaned image: $err"
	return
    }
    global scale
    set scale(mode) zscale
    ChangeScaleMode
    set catpanel(status) "LSBG: Cleaned image loaded in new frame"
}

proc CatalogPanelLSBGDetect {} {
    global catpanel

    if {![file exists $catpanel(lsbg,cleaned_file)]} {
	set catpanel(status) "LSBG: No cleaned image — run Background Model first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Detecting candidates..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode detect \
	--cleaned $catpanel(lsbg,cleaned_file) \
	--detect-thresh $catpanel(lsbg,param,detect-thresh) \
	--detect-minarea $catpanel(lsbg,param,detect-minarea) \
	--detect-filter-kernel $catpanel(lsbg,param,detect-filter-kernel) \
	--deblend-nthresh $catpanel(lsbg,param,deblend-nthresh) \
	--deblend-mincont $catpanel(lsbg,param,deblend-mincont) \
	--multiscale-factors $catpanel(lsbg,param,multiscale-factors) \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--segmap-output $catpanel(lsbg,segmap_file) \
	--n-workers $catpanel(param,n-workers)]
    if {$catpanel(lsbg,param,multiscale)} {
	lappend args --multiscale
    } else {
	lappend args --no-multiscale
    }

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG detect error: $err"
	return
    }

    # Parse detection results into table
    set catpanel(lsbg,detect_data) $result
    set catpanel(lsbg,has_detect) 1

    # Auto-display segmentation map in new frame
    if {[file exists $catpanel(lsbg,segmap_file)]} {
	CreateFrame
	if {[catch {LoadFitsFile $catpanel(lsbg,segmap_file) {} {}} err]} {
	    set catpanel(status) "LSBG: Detect done but could not display segmap: $err"
	} else {
	    global scale
	    set scale(mode) zscale
	    ChangeScaleMode
	}
    }

    # Load into panel table
    set catpanel(alldata) $result
    CatalogPanelLoadTSV $catpanel(alldata) "lsbg"

    # Count detections
    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]
    set catpanel(status) "LSBG: $nsrc candidates detected (segmap in new frame)"
}

proc CatalogPanelLSBGPhotometry {} {
    global catpanel

    if {![file exists $catpanel(lsbg,cleaned_file)]} {
	set catpanel(status) "LSBG: No cleaned image — run Background Model first"
	return
    }

    if {![file exists $catpanel(lsbg,segmap_file)]} {
	set catpanel(status) "LSBG: No detections — run Detect first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Measuring photometry..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode photometry \
	--cleaned $catpanel(lsbg,cleaned_file) \
	--segmap $catpanel(lsbg,segmap_file) \
	--detect-thresh $catpanel(lsbg,param,detect-thresh) \
	--detect-minarea $catpanel(lsbg,param,detect-minarea) \
	--detect-filter-kernel $catpanel(lsbg,param,detect-filter-kernel) \
	--deblend-nthresh $catpanel(lsbg,param,deblend-nthresh) \
	--deblend-mincont $catpanel(lsbg,param,deblend-mincont) \
	--phot-apertures $catpanel(lsbg,param,phot-apertures) \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--n-workers $catpanel(param,n-workers)]

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG photometry error: $err"
	return
    }

    set catpanel(alldata) $result
    set catpanel(lsbg,has_catalog) 1
    CatalogPanelLoadTSV $catpanel(alldata) "lsbg"

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    CatalogPanelMarkAll
    set catpanel(status) "LSBG: $nsrc sources — photometry done"
}

proc CatalogPanelLSBGSersic {} {
    global catpanel

    if {![file exists $catpanel(lsbg,cleaned_file)]} {
	set catpanel(status) "LSBG: No cleaned image — run Background Model first"
	return
    }

    if {![file exists $catpanel(lsbg,segmap_file)]} {
	set catpanel(status) "LSBG: No detections — run Detect first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Fitting Sérsic profiles..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode sersic \
	--cleaned $catpanel(lsbg,cleaned_file) \
	--segmap $catpanel(lsbg,segmap_file) \
	--detect-thresh $catpanel(lsbg,param,detect-thresh) \
	--detect-minarea $catpanel(lsbg,param,detect-minarea) \
	--detect-filter-kernel $catpanel(lsbg,param,detect-filter-kernel) \
	--deblend-nthresh $catpanel(lsbg,param,deblend-nthresh) \
	--deblend-mincont $catpanel(lsbg,param,deblend-mincont) \
	--phot-apertures $catpanel(lsbg,param,phot-apertures) \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--sersic-n-min $catpanel(lsbg,param,sersic-n-min) \
	--sersic-n-max $catpanel(lsbg,param,sersic-n-max) \
	--sersic-re-min $catpanel(lsbg,param,sersic-re-min) \
	--sersic-cutout-scale $catpanel(lsbg,param,sersic-cutout-scale) \
	--sersic-max-nfev $catpanel(lsbg,param,sersic-max-nfev) \
	--n-workers $catpanel(param,n-workers)]

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG Sérsic fit error: $err"
	return
    }

    set catpanel(alldata) $result
    set catpanel(lsbg,has_catalog) 1
    CatalogPanelLoadTSV $catpanel(alldata) "lsbg"

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    CatalogPanelMarkAll
    set catpanel(status) "LSBG: $nsrc sources — Sérsic fit done"
}

proc CatalogPanelLSBGFilter {} {
    global catpanel

    if {![file exists $catpanel(lsbg,cleaned_file)]} {
	set catpanel(status) "LSBG: No cleaned image — run Background Model first"
	return
    }

    if {![file exists $catpanel(lsbg,segmap_file)]} {
	set catpanel(status) "LSBG: No detections — run Detect first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Filtering + grading candidates..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode filter \
	--cleaned $catpanel(lsbg,cleaned_file) \
	--segmap $catpanel(lsbg,segmap_file) \
	--detect-thresh $catpanel(lsbg,param,detect-thresh) \
	--detect-minarea $catpanel(lsbg,param,detect-minarea) \
	--detect-filter-kernel $catpanel(lsbg,param,detect-filter-kernel) \
	--deblend-nthresh $catpanel(lsbg,param,deblend-nthresh) \
	--deblend-mincont $catpanel(lsbg,param,deblend-mincont) \
	--phot-apertures $catpanel(lsbg,param,phot-apertures) \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--mu-eff-min $catpanel(lsbg,param,mu-eff-min) \
	--mu-eff-max $catpanel(lsbg,param,mu-eff-max) \
	--r-eff-min $catpanel(lsbg,param,r-eff-min) \
	--r-eff-max $catpanel(lsbg,param,r-eff-max) \
	--ellipticity-max $catpanel(lsbg,param,ellipticity-max) \
	--min-snr $catpanel(lsbg,param,min-snr) \
	--sersic-n-min $catpanel(lsbg,param,sersic-n-min) \
	--sersic-n-max $catpanel(lsbg,param,sersic-n-max) \
	--sersic-re-min $catpanel(lsbg,param,sersic-re-min) \
	--sersic-cutout-scale $catpanel(lsbg,param,sersic-cutout-scale) \
	--sersic-max-nfev $catpanel(lsbg,param,sersic-max-nfev) \
	--sersic-n-filter-min $catpanel(lsbg,param,sersic-n-filter-min) \
	--sersic-n-filter-max $catpanel(lsbg,param,sersic-n-filter-max) \
	--sersic-chi2-max $catpanel(lsbg,param,sersic-chi2-max) \
	--n-workers $catpanel(param,n-workers)]
    if {$catpanel(lsbg,param,sersic-fit)} {
	lappend args --sersic-fit
    } else {
	lappend args --no-sersic-fit
    }

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG filter error: $err"
	return
    }

    set catpanel(alldata) $result
    set catpanel(lsbg,has_catalog) 1
    CatalogPanelLoadTSV $catpanel(alldata) "lsbg"

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    CatalogPanelMarkAll
    set catpanel(status) "LSBG: $nsrc candidates passed filtering"
}

proc CatalogPanelLSBGSVMClassify {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "LSBG SVM: No catalog loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG SVM: ds9_lsbg.py not found"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG SVM: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set catpanel(status) "LSBG SVM: Classifying candidates..."
    update idletasks

    # Save current catalog to temp file
    set tmpcat [file join [file normalize ~] .ds9 lsbg_svm_tmp.tsv]
    if {[catch {
	set fd [open $tmpcat w]
	puts $fd $catpanel(alldata)
	close $fd
    } err]} {
	set catpanel(status) "LSBG SVM: Cannot write temp catalog"
	return
    }

    set args [list [OGFPython] $script $fn --mode svm-classify \
	--catalog $tmpcat \
	--svm-threshold $catpanel(lsbg,param,svm-threshold)]
    if {$catpanel(lsbg,param,svm-checkpoint) ne {}} {
	lappend args --svm-checkpoint $catpanel(lsbg,param,svm-checkpoint)
    }

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG SVM error: $err"
	catch {file delete $tmpcat}
	return
    }
    catch {file delete $tmpcat}

    # Add SVM columns to catalog
    CatalogPanelAddColumnsFromTSV $result

    # Color markers by SVM result
    CatalogPanelLSBGSVMColorMarkers $result

    # Count results
    set lines [split $result \n]
    set n_lsbg 0
    set n_total 0
    foreach line [lrange $lines 1 end] {
	if {[string trim $line] eq {}} continue
	incr n_total
	set fields [split $line "\t"]
	if {[llength $fields] >= 2 && [lindex $fields 1] eq "1"} {
	    incr n_lsbg
	}
    }
    set catpanel(status) "LSBG SVM: $n_lsbg LSBG / $n_total total"
}

proc CatalogPanelLSBGSVMColorMarkers {result_tsv} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    set frame $current(frame)

    # Parse SVM results: NUMBER -> SVM_LSBG
    array set svm_class {}
    set lines [split $result_tsv \n]
    foreach line [lrange $lines 1 end] {
	set fields [split $line "\t"]
	if {[llength $fields] < 2} continue
	set src_num [lindex $fields 0]
	set svm_lsbg [lindex $fields 1]
	set svm_class($src_num) $svm_lsbg
    }

    # Delete existing markers and recreate with SVM colors
    OGFMarkDelete $frame

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set alllines [split $catpanel(alldata) \n]
    if {[llength $alllines] < 2} return

    set headers [split [lindex $alllines 0] "\t"]
    set ncols [llength $headers]

    # Find needed columns
    set col_x -1
    set col_y -1
    set col_a -1
    set col_b -1
    set col_theta -1
    set col_ir -1
    set col_num -1
    for {set c 0} {$c < $ncols} {incr c} {
	set hdr [string trim [lindex $headers $c]]
	switch -- $hdr {
	    NUMBER      { set col_num $c }
	    X_IMAGE     { set col_x $c }
	    Y_IMAGE     { set col_y $c }
	    A_IMAGE     { set col_a $c }
	    B_IMAGE     { set col_b $c }
	    THETA_IMAGE { set col_theta $c }
	    ISO_RADIUS  { set col_ir $c }
	}
    }
    if {$col_x < 0 || $col_y < 0} return

    # Build region strings with SVM colors
    set batch_size 500
    set reg "image\n"
    set count 0
    set batch_count 0
    global sextract_all_reg

    for {set i 1} {$i < [llength $alllines]} {incr i} {
	set line [lindex $alllines $i]
	if {[string trim $line] eq {}} continue
	set fields [split $line "\t"]

	set x [string trim [lindex $fields $col_x]]
	set y [string trim [lindex $fields $col_y]]
	if {![string is double -strict $x] || ![string is double -strict $y]} continue

	set src_num [expr {$i}]
	if {$col_num >= 0} {
	    set nv [string trim [lindex $fields $col_num]]
	    if {$nv ne {}} { set src_num $nv }
	}

	# Ellipse parameters
	set iso_radius 5.0
	set a_image 0
	set b_image 0
	set theta 0

	if {$col_ir >= 0} {
	    set val [string trim [lindex $fields $col_ir]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} {
		set iso_radius $v
	    }
	}
	if {$col_a >= 0} {
	    set val [string trim [lindex $fields $col_a]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set a_image $v }
	}
	if {$col_b >= 0} {
	    set val [string trim [lindex $fields $col_b]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set b_image $v }
	}
	if {$col_theta >= 0} {
	    set val [string trim [lindex $fields $col_theta]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0} { set theta $v }
	}

	set semi_a $iso_radius
	set semi_b $iso_radius
	if {$a_image > 0 && $b_image > 0} {
	    set semi_b [expr {$iso_radius * $b_image / $a_image}]
	}

	# Color: LSBG=green, contaminant=red, unclassified=yellow
	set color yellow
	if {[info exists svm_class($src_num)]} {
	    if {$svm_class($src_num) eq "1"} {
		set color green
	    } else {
		set color red
	    }
	}

	append reg "ellipse($x $y ${semi_a}i ${semi_b}i $theta) # color=$color width=1 tag={sextract_all} tag={sextract_src.$src_num} select=0 edit=0 move=0 rotate=0 delete=1 highlite=1 callback=highlite CatalogPanelMarkerCB {$src_num} callback=unhighlite CatalogPanelMarkerUnCB {$src_num}\n"

	incr count
	if {$count >= $batch_size} {
	    set sextract_all_reg $reg
	    OGFMarkSend $frame
	    set reg "image\n"
	    set count 0
	    incr batch_count
	}
    }

    # Flush remaining
    if {$count > 0} {
	set sextract_all_reg $reg
	OGFMarkSend $frame
    }
}

proc CatalogPanelLSBGRunAll {} {
    global catpanel

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    OGFMaskEnsure lsbg
    set catpanel(status) "LSBG: Running full pipeline..."
    update idletasks

    # Reset cmdlog for new session
    set catpanel(lsbg,cmdlog) {}

    set args [list [OGFPython] $script $fn --mode run \
	--mask-detect-thresh $catpanel(lsbg,param,mask-detect-thresh) \
	--mask-detect-minarea $catpanel(lsbg,param,mask-detect-minarea) \
	--mask-expand-factor $catpanel(lsbg,param,mask-expand-factor) \
	--max-dilate-radius $catpanel(lsbg,param,max-dilate-radius) \
	--bright-star-mag-limit $catpanel(lsbg,param,bright-star-mag-limit) \
	--bright-star-radius-scale $catpanel(lsbg,param,bright-star-radius-scale) \
	--mask-mag-threshold $catpanel(lsbg,param,mask-mag-threshold) \
	--interp-method $catpanel(lsbg,param,interp-method) \
	--lsb-mu-threshold $catpanel(lsbg,param,lsb-mu-threshold) \
	--bkg-method $catpanel(lsbg,param,bkg-method) \
	--bkg-mesh-size $catpanel(lsbg,param,bkg-mesh-size) \
	--bkg-poly-order $catpanel(lsbg,param,bkg-poly-order) \
	--bkg-sigma-clip $catpanel(lsbg,param,bkg-sigma-clip) \
	--bkg-n-iterations $catpanel(lsbg,param,bkg-n-iterations) \
	--bkg-refine-thresh $catpanel(lsbg,param,bkg-refine-thresh) \
	--bkg-rms-quantile $catpanel(lsbg,param,bkg-rms-quantile) \
	--bkg-convergence-tol $catpanel(lsbg,param,bkg-convergence-tol) \
	--detect-thresh $catpanel(lsbg,param,detect-thresh) \
	--detect-minarea $catpanel(lsbg,param,detect-minarea) \
	--detect-filter-kernel $catpanel(lsbg,param,detect-filter-kernel) \
	--deblend-nthresh $catpanel(lsbg,param,deblend-nthresh) \
	--deblend-mincont $catpanel(lsbg,param,deblend-mincont) \
	--multiscale-factors $catpanel(lsbg,param,multiscale-factors) \
	--sersic-n-min $catpanel(lsbg,param,sersic-n-min) \
	--sersic-n-max $catpanel(lsbg,param,sersic-n-max) \
	--sersic-re-min $catpanel(lsbg,param,sersic-re-min) \
	--sersic-cutout-scale $catpanel(lsbg,param,sersic-cutout-scale) \
	--sersic-max-nfev $catpanel(lsbg,param,sersic-max-nfev) \
	--phot-apertures $catpanel(lsbg,param,phot-apertures) \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--mu-eff-min $catpanel(lsbg,param,mu-eff-min) \
	--mu-eff-max $catpanel(lsbg,param,mu-eff-max) \
	--r-eff-min $catpanel(lsbg,param,r-eff-min) \
	--r-eff-max $catpanel(lsbg,param,r-eff-max) \
	--ellipticity-max $catpanel(lsbg,param,ellipticity-max) \
	--min-snr $catpanel(lsbg,param,min-snr) \
	--sersic-n-filter-min $catpanel(lsbg,param,sersic-n-filter-min) \
	--sersic-n-filter-max $catpanel(lsbg,param,sersic-n-filter-max) \
	--sersic-chi2-max $catpanel(lsbg,param,sersic-chi2-max) \
	--mask-input $catpanel(lsbg,mask_file) \
	--mask-output [OGFMaskRefinedPath lsbg] \
	--masked-output $catpanel(lsbg,masked_file) \
	--bkg-output $catpanel(lsbg,bkg_file) \
	--cleaned-output $catpanel(lsbg,cleaned_file) \
	--segmap-output $catpanel(lsbg,segmap_file) \
	--catalog-output $catpanel(lsbg,catalog_file) \
	--n-workers $catpanel(param,n-workers)]
    if {$catpanel(lsbg,param,lsb-protect)} {
	lappend args --lsb-protect
    } else {
	lappend args --no-lsb-protect
    }
    if {$catpanel(lsbg,param,multiscale)} {
	lappend args --multiscale
    } else {
	lappend args --no-multiscale
    }
    if {$catpanel(lsbg,param,sersic-fit)} {
	lappend args --sersic-fit
    } else {
	lappend args --no-sersic-fit
    }
    if {$catpanel(lsbg,param,svm-classify)} {
	lappend args --svm-classify
	lappend args --svm-threshold $catpanel(lsbg,param,svm-threshold)
	if {$catpanel(lsbg,param,svm-checkpoint) ne {}} {
	    lappend args --svm-checkpoint $catpanel(lsbg,param,svm-checkpoint)
	}
    }

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG pipeline error: $err"
	return
    }

    # Update state
    set catpanel(lsbg,has_mask) 1
    set catpanel(lsbg,has_clean) 1
    set catpanel(lsbg,has_detect) 1
    set catpanel(lsbg,has_catalog) 1

    # Load cleaned image in new frame
    CreateFrame
    if {[catch {LoadFitsFile $catpanel(lsbg,cleaned_file) {} {}} err]} {
	set catpanel(status) "LSBG: Warning — could not load cleaned image"
    } else {
	global scale
	set scale(mode) zscale
	ChangeScaleMode
    }

    # Load catalog
    set catpanel(alldata) $result
    CatalogPanelLoadTSV $catpanel(alldata) "lsbg"

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    # Mark all LSBG candidates
    CatalogPanelMarkAll
    set catpanel(status) "LSBG: Pipeline complete — $nsrc candidates"
}

proc CatalogPanelLSBGForcedPhot {} {
    global catpanel current

    # Check that we have a catalog
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "LSBG: No catalog — run pipeline first"
	return
    }

    # Select the other band FITS file
    set band_fits [tk_getOpenFile \
	-title "Select Other Band FITS Image" \
	-filetypes {
	    {{FITS Files} {.fits .fit .fts .fits.gz .fit.gz}}
	    {{All Files} {*}}
	}]
    if {$band_fits eq {}} return

    # Ask for band name
    set w .lsbg_bandname
    catch {destroy $w}
    toplevel $w
    wm title $w "Band Name"
    wm geometry $w 300x100
    wm transient $w .
    ttk::label $w.l -text "Enter band name (e.g. F606W):"
    ttk::entry $w.e -textvariable catpanel(lsbg,tmp_bandname)
    set catpanel(lsbg,tmp_bandname) ""
    ttk::frame $w.btns
    ttk::button $w.btns.ok -text "OK" -command [list set catpanel(lsbg,band_dialog_done) 1]
    ttk::button $w.btns.cancel -text "Cancel" -command [list set catpanel(lsbg,band_dialog_done) 0]
    pack $w.l -padx 10 -pady 5
    pack $w.e -padx 10 -fill x
    pack $w.btns -pady 5
    pack $w.btns.ok $w.btns.cancel -side left -padx 10
    focus $w.e
    bind $w.e <Return> [list set catpanel(lsbg,band_dialog_done) 1]
    tkwait variable catpanel(lsbg,band_dialog_done)
    set band_name $catpanel(lsbg,tmp_bandname)
    catch {destroy $w}
    if {!$catpanel(lsbg,band_dialog_done) || $band_name eq {}} return

    # Save current catalog to temp file
    set tmpcat [CatalogPanelSaveTempCatalog lsbg_forced]
    if {$tmpcat eq {}} {
	set catpanel(status) "LSBG: Failed to save temp catalog"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Forced photometry ($band_name)..."
    update idletasks

    set args [list [OGFPython] $script $band_fits --mode forced \
	--catalog $tmpcat \
	--band-name $band_name \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--n-workers $catpanel(param,n-workers)]

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG forced phot error: $err"
	return
    }

    # Merge result columns into current catalog
    set col_names [list FLUX_$band_name FLUXERR_$band_name MAG_$band_name MAGERR_$band_name]
    CatalogPanelAddColumnsFromTSV $result $col_names

    # Load the other band image in a new frame
    CreateFrame
    if {[catch {LoadFitsFile $band_fits {} {}} err]} {
	set catpanel(status) "LSBG: Forced phot done but could not display band image"
    } else {
	global scale
	set scale(mode) zscale
	ChangeScaleMode
    }

    set catpanel(status) "LSBG: Forced photometry ($band_name) complete — columns added"
}

