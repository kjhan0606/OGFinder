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

proc CatalogPanelLSBGSettings {} {
    global catpanel
    global ed

    set w .lsbgsettings
    if {[winfo exists $w]} {
	raise $w
	return
    }

    toplevel $w
    wm title $w "LSBG Detection Settings"
    wm geometry $w 500x620

    # Copy current values
    foreach pname {mask-detect-thresh mask-detect-minarea mask-expand-factor \
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
	set ed(lsbg,$pname) $catpanel(lsbg,param,$pname)
    }

    ttk::notebook $w.nb
    pack $w.nb -fill both -expand true -padx 8 -pady 8

    # --- Tab 1: Masking ---
    set t1 [ttk::frame $w.nb.mask]
    $w.nb add $t1 -text "Masking"

    set r 0
    ttk::label $t1.lmt -text "Mask mag threshold:"
    ttk::entry $t1.emt -textvariable ed(lsbg,mask-mag-threshold) -width 10
    grid $t1.lmt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.emt -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lef -text "Expand factor:"
    ttk::entry $t1.eef -textvariable ed(lsbg,mask-expand-factor) -width 10
    grid $t1.lef -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.eef -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lmdr -text "Max dilate radius (px):"
    ttk::entry $t1.emdr -textvariable ed(lsbg,max-dilate-radius) -width 10
    grid $t1.lmdr -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.emdr -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lml -text "Bright star mag limit:"
    ttk::entry $t1.eml -textvariable ed(lsbg,bright-star-mag-limit) -width 10
    grid $t1.lml -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.eml -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lrs -text "Bright star radius scale:"
    ttk::entry $t1.ers -textvariable ed(lsbg,bright-star-radius-scale) -width 10
    grid $t1.lrs -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.ers -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lim -text "Interpolation method:"
    ttk::combobox $t1.eim -textvariable ed(lsbg,interp-method) -width 10 \
	-values {linear cubic nearest}
    grid $t1.lim -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.eim -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.ldt -text "Mask detect threshold:"
    ttk::entry $t1.edt -textvariable ed(lsbg,mask-detect-thresh) -width 10
    grid $t1.ldt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.edt -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::checkbutton $t1.clp -text "LSB structure protection" \
	-variable ed(lsbg,lsb-protect)
    grid $t1.clp -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady 4
    incr r

    ttk::label $t1.llm -text "LSB mu threshold:"
    ttk::entry $t1.elm -textvariable ed(lsbg,lsb-mu-threshold) -width 10
    grid $t1.llm -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.elm -row $r -column 1 -sticky w -padx 4 -pady 4

    # --- Tab 2: Background ---
    set t2 [ttk::frame $w.nb.bkg]
    $w.nb add $t2 -text "Background"

    set r 0
    ttk::label $t2.lbm -text "Method:"
    ttk::combobox $t2.ebm -textvariable ed(lsbg,bkg-method) -width 12 \
	-values {sep_large polynomial chebyshev}
    grid $t2.lbm -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.ebm -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lms -text "Mesh size (px):"
    ttk::entry $t2.ems -textvariable ed(lsbg,bkg-mesh-size) -width 10
    grid $t2.lms -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.ems -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lbo -text "Polynomial order:"
    ttk::entry $t2.ebo -textvariable ed(lsbg,bkg-poly-order) -width 10
    grid $t2.lbo -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.ebo -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lsc -text "Sigma clip:"
    ttk::entry $t2.esc -textvariable ed(lsbg,bkg-sigma-clip) -width 10
    grid $t2.lsc -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.esc -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lni -text "Iterations:"
    ttk::entry $t2.eni -textvariable ed(lsbg,bkg-n-iterations) -width 10
    grid $t2.lni -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.eni -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lrt -text "Refine threshold (sigma):"
    ttk::entry $t2.ert -textvariable ed(lsbg,bkg-refine-thresh) -width 10
    grid $t2.lrt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.ert -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lrq -text "RMS quantile:"
    ttk::entry $t2.erq -textvariable ed(lsbg,bkg-rms-quantile) -width 10
    grid $t2.lrq -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.erq -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lct -text "Convergence tolerance:"
    ttk::entry $t2.ect -textvariable ed(lsbg,bkg-convergence-tol) -width 10
    grid $t2.lct -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.ect -row $r -column 1 -sticky w -padx 4 -pady 4

    # --- Tab 3: Detection ---
    set t3 [ttk::frame $w.nb.det]
    $w.nb add $t3 -text "Detection"

    set r 0
    ttk::label $t3.ldt -text "Detect threshold (sigma):"
    ttk::entry $t3.edt -textvariable ed(lsbg,detect-thresh) -width 10
    grid $t3.ldt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.edt -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lma -text "Min area (px):"
    ttk::entry $t3.ema -textvariable ed(lsbg,detect-minarea) -width 10
    grid $t3.lma -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.ema -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lfk -text "Filter kernel:"
    ttk::combobox $t3.efk -textvariable ed(lsbg,detect-filter-kernel) -width 12 \
	-values {none gauss3x3 gauss5x5 gauss7x7 gauss9x9 tophat5 tophat7 mexhat}
    grid $t3.lfk -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.efk -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lnt -text "Deblend N thresholds:"
    ttk::entry $t3.ent -textvariable ed(lsbg,deblend-nthresh) -width 10
    grid $t3.lnt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.ent -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lmc -text "Deblend min contrast:"
    ttk::entry $t3.emc -textvariable ed(lsbg,deblend-mincont) -width 10
    grid $t3.lmc -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.emc -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::checkbutton $t3.cms -text "Multi-scale detection" \
	-variable ed(lsbg,multiscale)
    grid $t3.cms -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady 4
    incr r

    ttk::label $t3.lmf -text "Scale factors:"
    ttk::entry $t3.emf -textvariable ed(lsbg,multiscale-factors) -width 10
    grid $t3.lmf -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.emf -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::checkbutton $t3.csf -text "Sérsic profile fitting" \
	-variable ed(lsbg,sersic-fit)
    grid $t3.csf -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady 4

    # --- Tab 4: Calibration & Filter ---
    set t4 [ttk::frame $w.nb.cal]
    $w.nb add $t4 -text "Calibration & Filter"

    set r 0
    ttk::label $t4.lzp -text "Mag zeropoint:"
    ttk::entry $t4.ezp -textvariable ed(lsbg,mag-zeropoint) -width 10
    grid $t4.lzp -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.ezp -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lps -text "Pixel scale (arcsec/px):"
    ttk::entry $t4.eps -textvariable ed(lsbg,pixel-scale) -width 10
    grid $t4.lps -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.eps -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lmn -text "mu_eff min (mag/arcsec2):"
    ttk::entry $t4.emn -textvariable ed(lsbg,mu-eff-min) -width 10
    grid $t4.lmn -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.emn -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lmx -text "mu_eff max (mag/arcsec2):"
    ttk::entry $t4.emx -textvariable ed(lsbg,mu-eff-max) -width 10
    grid $t4.lmx -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.emx -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lrn -text "R_eff min (arcsec):"
    ttk::entry $t4.ern -textvariable ed(lsbg,r-eff-min) -width 10
    grid $t4.lrn -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.ern -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lrx -text "R_eff max (arcsec):"
    ttk::entry $t4.erx -textvariable ed(lsbg,r-eff-max) -width 10
    grid $t4.lrx -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.erx -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lex -text "Ellipticity max:"
    ttk::entry $t4.eex -textvariable ed(lsbg,ellipticity-max) -width 10
    grid $t4.lex -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.eex -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lsn -text "Min SNR:"
    ttk::entry $t4.esn -textvariable ed(lsbg,min-snr) -width 10
    grid $t4.lsn -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.esn -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lap -text "Phot apertures (px):"
    ttk::entry $t4.eap -textvariable ed(lsbg,phot-apertures) -width 18
    grid $t4.lap -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.eap -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::separator $t4.sep -orient horizontal
    grid $t4.sep -row $r -column 0 -columnspan 2 -sticky ew -padx 8 -pady 6
    incr r

    ttk::label $t4.lsnm -text "Sérsic n min:"
    ttk::entry $t4.esnm -textvariable ed(lsbg,sersic-n-filter-min) -width 10
    grid $t4.lsnm -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.esnm -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lsnx -text "Sérsic n max:"
    ttk::entry $t4.esnx -textvariable ed(lsbg,sersic-n-filter-max) -width 10
    grid $t4.lsnx -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.esnx -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lsc2 -text "Sérsic chi2 max:"
    ttk::entry $t4.esc2 -textvariable ed(lsbg,sersic-chi2-max) -width 10
    grid $t4.lsc2 -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.esc2 -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::separator $t4.sep2 -orient horizontal
    grid $t4.sep2 -row $r -column 0 -columnspan 2 -sticky ew -padx 8 -pady 6
    incr r

    ttk::checkbutton $t4.csvm -text "SVM Classification" \
	-variable ed(lsbg,svm-classify)
    grid $t4.csvm -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady 4
    incr r

    ttk::label $t4.lsvmt -text "SVM threshold:"
    ttk::entry $t4.esvmt -textvariable ed(lsbg,svm-threshold) -width 10
    grid $t4.lsvmt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.esvmt -row $r -column 1 -sticky w -padx 4 -pady 4

    # Buttons
    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8

    ttk::button $bf.apply -text "Apply" -command [list CatalogPanelLSBGSettingsApply $w]
    ttk::button $bf.defaults -text "Defaults" -command CatalogPanelLSBGSettingsDefaults
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.defaults -side right -padx 4
    pack $bf.apply -side right -padx 4
}

proc CatalogPanelLSBGSettingsApply {w} {
    global catpanel
    global ed

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
	set catpanel(lsbg,param,$pname) $ed(lsbg,$pname)
    }
    CatalogPanelLSBGParamSave
    set catpanel(status) "LSBG settings applied and saved"
}

proc CatalogPanelLSBGSettingsDefaults {} {
    global ed

    set ed(lsbg,mask-detect-thresh) 1.5
    set ed(lsbg,mask-detect-minarea) 5
    set ed(lsbg,mask-expand-factor) 1.5
    set ed(lsbg,max-dilate-radius) 30
    set ed(lsbg,bright-star-mag-limit) 18.0
    set ed(lsbg,bright-star-radius-scale) 12.0
    set ed(lsbg,mask-mag-threshold) 22.0
    set ed(lsbg,interp-method) linear
    set ed(lsbg,lsb-protect) 1
    set ed(lsbg,lsb-mu-threshold) 24.0
    set ed(lsbg,bkg-method) sep_large
    set ed(lsbg,bkg-mesh-size) 256
    set ed(lsbg,bkg-poly-order) 3
    set ed(lsbg,bkg-sigma-clip) 3.0
    set ed(lsbg,bkg-n-iterations) 3
    set ed(lsbg,bkg-refine-thresh) 2.0
    set ed(lsbg,bkg-rms-quantile) 0.25
    set ed(lsbg,bkg-convergence-tol) 0.01
    set ed(lsbg,detect-thresh) 0.8
    set ed(lsbg,detect-minarea) 50
    set ed(lsbg,detect-filter-kernel) gauss5x5
    set ed(lsbg,deblend-nthresh) 32
    set ed(lsbg,deblend-mincont) 0.005
    set ed(lsbg,multiscale) 1
    set ed(lsbg,multiscale-factors) 1,2,4
    set ed(lsbg,sersic-fit) 1
    set ed(lsbg,sersic-n-min) 0.2
    set ed(lsbg,sersic-n-max) 10.0
    set ed(lsbg,sersic-re-min) 0.5
    set ed(lsbg,sersic-cutout-scale) 5.0
    set ed(lsbg,sersic-max-nfev) 500
    set ed(lsbg,phot-apertures) 5,10,20,40
    set ed(lsbg,mag-zeropoint) 25.0
    set ed(lsbg,pixel-scale) 0.06
    set ed(lsbg,mu-eff-min) 24.0
    set ed(lsbg,mu-eff-max) 30.0
    set ed(lsbg,r-eff-min) 2.5
    set ed(lsbg,r-eff-max) 60.0
    set ed(lsbg,ellipticity-max) 0.7
    set ed(lsbg,min-snr) 2.0
    set ed(lsbg,sersic-n-filter-min) 0.3
    set ed(lsbg,sersic-n-filter-max) 6.0
    set ed(lsbg,sersic-chi2-max) 10.0
    set ed(lsbg,svm-classify) 0
    set ed(lsbg,svm-threshold) 0.5
    set ed(lsbg,svm-checkpoint) {}
}

