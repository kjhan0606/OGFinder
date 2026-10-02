# lsbg/lsbg.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/lsbg/plugin.json.

proc CatalogPanelLSBGUpdateFiles {fn} {
    if {$fn eq {}} return
    set base [CatalogPanelFitsBaseName $fn]
    if {$base eq {}} return
    # Skip if already set for this base
    if {[::ogf::cat::exists lsbg,fits_base] &&
	[::ogf::cat::get lsbg,fits_base] eq $base} return
    ::ogf::cat::set lsbg,fits_base $base
    set ds9dir [file join [file normalize ~] .ds9]
    ::ogf::cat::set lsbg,mask_file [file join $ds9dir "mask_${base}_bool.fits"]
    ::ogf::cat::set lsbg,masked_file [file join $ds9dir "mask_${base}_masked.fits"]
    ::ogf::cat::set lsbg,bkg_file [file join $ds9dir "lsbg_background_${base}.fits"]
    ::ogf::cat::set lsbg,cleaned_file [file join $ds9dir "lsbg_cleaned_${base}.fits"]
    ::ogf::cat::set lsbg,segmap_file [file join $ds9dir "lsbg_segmap_${base}.fits"]
    ::ogf::cat::set lsbg,catalog_file [file join $ds9dir "lsbg_catalog_${base}.tsv"]
    # Detect if previous results exist for this FITS
    ::ogf::cat::set lsbg,has_mask [file exists [::ogf::cat::get lsbg,mask_file]]
    ::ogf::cat::set lsbg,has_clean [file exists [::ogf::cat::get lsbg,cleaned_file]]
    ::ogf::cat::set lsbg,has_detect [file exists [::ogf::cat::get lsbg,segmap_file]]
    ::ogf::cat::set lsbg,has_catalog [file exists [::ogf::cat::get lsbg,catalog_file]]
}

proc CatalogPanelLSBGParamLoad {} {

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
	    if {[::ogf::cat::exists lsbg,param,$key]} {
		::ogf::cat::set lsbg,param,$key $val
	    }
	}
    }
    close $fd
}

proc CatalogPanelLSBGParamSave {} {

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
	puts $fd "$pname [::ogf::cat::get lsbg,param,$pname]"
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

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	::ogf::cat::set status "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    if {![OGFMaskEnsure lsbg]} {
	::ogf::cat::set status "LSBG: could not obtain a mask - run Mask > Auto Mask first"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	::ogf::cat::set status "LSBG: ds9_lsbg.py not found"
	return
    }

    ::ogf::cat::set status "LSBG: Iterative background ($method)..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode clean \
	--mask [::ogf::cat::get lsbg,mask_file] \
	--bkg-method $method \
	--bkg-mesh-size [::ogf::cat::get lsbg,param,bkg-mesh-size] \
	--bkg-poly-order [::ogf::cat::get lsbg,param,bkg-poly-order] \
	--bkg-sigma-clip [::ogf::cat::get lsbg,param,bkg-sigma-clip] \
	--bkg-n-iterations [::ogf::cat::get lsbg,param,bkg-n-iterations] \
	--bkg-refine-thresh [::ogf::cat::get lsbg,param,bkg-refine-thresh] \
	--bkg-rms-quantile [::ogf::cat::get lsbg,param,bkg-rms-quantile] \
	--bkg-convergence-tol [::ogf::cat::get lsbg,param,bkg-convergence-tol] \
	--interp-method [::ogf::cat::get lsbg,param,interp-method] \
	--mag-zeropoint [::ogf::cat::get lsbg,param,mag-zeropoint] \
	--pixel-scale [::ogf::cat::get lsbg,param,pixel-scale] \
	--bkg-output [::ogf::cat::get lsbg,bkg_file] \
	--cleaned-output [::ogf::cat::get lsbg,cleaned_file] \
	--n-workers [::ogf::cat::get param,n-workers]]

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	::ogf::cat::set status "LSBG clean error: $err"
	return
    }

    ::ogf::cat::set lsbg,has_clean 1

    # Auto-display cleaned image in new frame
    CreateFrame
    if {[catch {LoadFitsFile [::ogf::cat::get lsbg,cleaned_file] {} {}} err]} {
	::ogf::cat::set status "LSBG: Clean done but could not display: $err"
    } else {
	global scale
	set scale(mode) zscale
	ChangeScaleMode
    }

    ::ogf::cat::set status "LSBG: Background cleaned ($method) (new frame)"
}

proc CatalogPanelLSBGViewClean {} {

    if {![file exists [::ogf::cat::get lsbg,cleaned_file]]} {
	::ogf::cat::set status "LSBG: No cleaned image — run Background Model first"
	return
    }

    CreateFrame
    if {[catch {LoadFitsFile [::ogf::cat::get lsbg,cleaned_file] {} {}} err]} {
	::ogf::cat::set status "LSBG: Error loading cleaned image: $err"
	return
    }
    global scale
    set scale(mode) zscale
    ChangeScaleMode
    ::ogf::cat::set status "LSBG: Cleaned image loaded in new frame"
}

proc CatalogPanelLSBGDetect {} {

    if {![file exists [::ogf::cat::get lsbg,cleaned_file]]} {
	::ogf::cat::set status "LSBG: No cleaned image — run Background Model first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	::ogf::cat::set status "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	::ogf::cat::set status "LSBG: ds9_lsbg.py not found"
	return
    }

    ::ogf::cat::set status "LSBG: Detecting candidates..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode detect \
	--cleaned [::ogf::cat::get lsbg,cleaned_file] \
	--detect-thresh [::ogf::cat::get lsbg,param,detect-thresh] \
	--detect-minarea [::ogf::cat::get lsbg,param,detect-minarea] \
	--detect-filter-kernel [::ogf::cat::get lsbg,param,detect-filter-kernel] \
	--deblend-nthresh [::ogf::cat::get lsbg,param,deblend-nthresh] \
	--deblend-mincont [::ogf::cat::get lsbg,param,deblend-mincont] \
	--multiscale-factors [::ogf::cat::get lsbg,param,multiscale-factors] \
	--mag-zeropoint [::ogf::cat::get lsbg,param,mag-zeropoint] \
	--pixel-scale [::ogf::cat::get lsbg,param,pixel-scale] \
	--segmap-output [::ogf::cat::get lsbg,segmap_file] \
	--n-workers [::ogf::cat::get param,n-workers]]
    if {[::ogf::cat::get lsbg,param,multiscale]} {
	lappend args --multiscale
    } else {
	lappend args --no-multiscale
    }

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	::ogf::cat::set status "LSBG detect error: $err"
	return
    }

    # Parse detection results into table
    ::ogf::cat::set lsbg,detect_data $result
    ::ogf::cat::set lsbg,has_detect 1

    # Auto-display segmentation map in new frame
    if {[file exists [::ogf::cat::get lsbg,segmap_file]]} {
	CreateFrame
	if {[catch {LoadFitsFile [::ogf::cat::get lsbg,segmap_file] {} {}} err]} {
	    ::ogf::cat::set status "LSBG: Detect done but could not display segmap: $err"
	} else {
	    global scale
	    set scale(mode) zscale
	    ChangeScaleMode
	}
    }

    # Load into panel table
    ::ogf::cat::set alldata $result
    CatalogPanelLoadTSV [::ogf::cat::tsv] "lsbg"

    # Count detections
    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]
    ::ogf::cat::set status "LSBG: $nsrc candidates detected (segmap in new frame)"
}

proc CatalogPanelLSBGPhotometry {} {

    if {![file exists [::ogf::cat::get lsbg,cleaned_file]]} {
	::ogf::cat::set status "LSBG: No cleaned image — run Background Model first"
	return
    }

    if {![file exists [::ogf::cat::get lsbg,segmap_file]]} {
	::ogf::cat::set status "LSBG: No detections — run Detect first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	::ogf::cat::set status "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	::ogf::cat::set status "LSBG: ds9_lsbg.py not found"
	return
    }

    ::ogf::cat::set status "LSBG: Measuring photometry..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode photometry \
	--cleaned [::ogf::cat::get lsbg,cleaned_file] \
	--segmap [::ogf::cat::get lsbg,segmap_file] \
	--detect-thresh [::ogf::cat::get lsbg,param,detect-thresh] \
	--detect-minarea [::ogf::cat::get lsbg,param,detect-minarea] \
	--detect-filter-kernel [::ogf::cat::get lsbg,param,detect-filter-kernel] \
	--deblend-nthresh [::ogf::cat::get lsbg,param,deblend-nthresh] \
	--deblend-mincont [::ogf::cat::get lsbg,param,deblend-mincont] \
	--phot-apertures [::ogf::cat::get lsbg,param,phot-apertures] \
	--mag-zeropoint [::ogf::cat::get lsbg,param,mag-zeropoint] \
	--pixel-scale [::ogf::cat::get lsbg,param,pixel-scale] \
	--n-workers [::ogf::cat::get param,n-workers]]

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	::ogf::cat::set status "LSBG photometry error: $err"
	return
    }

    ::ogf::cat::set alldata $result
    ::ogf::cat::set lsbg,has_catalog 1
    CatalogPanelLoadTSV [::ogf::cat::tsv] "lsbg"

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    CatalogPanelMarkAll
    ::ogf::cat::set status "LSBG: $nsrc sources — photometry done"
}

proc CatalogPanelLSBGSersic {} {

    if {![file exists [::ogf::cat::get lsbg,cleaned_file]]} {
	::ogf::cat::set status "LSBG: No cleaned image — run Background Model first"
	return
    }

    if {![file exists [::ogf::cat::get lsbg,segmap_file]]} {
	::ogf::cat::set status "LSBG: No detections — run Detect first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	::ogf::cat::set status "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	::ogf::cat::set status "LSBG: ds9_lsbg.py not found"
	return
    }

    ::ogf::cat::set status "LSBG: Fitting Sérsic profiles..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode sersic \
	--cleaned [::ogf::cat::get lsbg,cleaned_file] \
	--segmap [::ogf::cat::get lsbg,segmap_file] \
	--detect-thresh [::ogf::cat::get lsbg,param,detect-thresh] \
	--detect-minarea [::ogf::cat::get lsbg,param,detect-minarea] \
	--detect-filter-kernel [::ogf::cat::get lsbg,param,detect-filter-kernel] \
	--deblend-nthresh [::ogf::cat::get lsbg,param,deblend-nthresh] \
	--deblend-mincont [::ogf::cat::get lsbg,param,deblend-mincont] \
	--phot-apertures [::ogf::cat::get lsbg,param,phot-apertures] \
	--mag-zeropoint [::ogf::cat::get lsbg,param,mag-zeropoint] \
	--pixel-scale [::ogf::cat::get lsbg,param,pixel-scale] \
	--sersic-n-min [::ogf::cat::get lsbg,param,sersic-n-min] \
	--sersic-n-max [::ogf::cat::get lsbg,param,sersic-n-max] \
	--sersic-re-min [::ogf::cat::get lsbg,param,sersic-re-min] \
	--sersic-cutout-scale [::ogf::cat::get lsbg,param,sersic-cutout-scale] \
	--sersic-max-nfev [::ogf::cat::get lsbg,param,sersic-max-nfev] \
	--n-workers [::ogf::cat::get param,n-workers]]

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	::ogf::cat::set status "LSBG Sérsic fit error: $err"
	return
    }

    ::ogf::cat::set alldata $result
    ::ogf::cat::set lsbg,has_catalog 1
    CatalogPanelLoadTSV [::ogf::cat::tsv] "lsbg"

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    CatalogPanelMarkAll
    ::ogf::cat::set status "LSBG: $nsrc sources — Sérsic fit done"
}

proc CatalogPanelLSBGFilter {} {

    if {![file exists [::ogf::cat::get lsbg,cleaned_file]]} {
	::ogf::cat::set status "LSBG: No cleaned image — run Background Model first"
	return
    }

    if {![file exists [::ogf::cat::get lsbg,segmap_file]]} {
	::ogf::cat::set status "LSBG: No detections — run Detect first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	::ogf::cat::set status "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	::ogf::cat::set status "LSBG: ds9_lsbg.py not found"
	return
    }

    ::ogf::cat::set status "LSBG: Filtering + grading candidates..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode filter \
	--cleaned [::ogf::cat::get lsbg,cleaned_file] \
	--segmap [::ogf::cat::get lsbg,segmap_file] \
	--detect-thresh [::ogf::cat::get lsbg,param,detect-thresh] \
	--detect-minarea [::ogf::cat::get lsbg,param,detect-minarea] \
	--detect-filter-kernel [::ogf::cat::get lsbg,param,detect-filter-kernel] \
	--deblend-nthresh [::ogf::cat::get lsbg,param,deblend-nthresh] \
	--deblend-mincont [::ogf::cat::get lsbg,param,deblend-mincont] \
	--phot-apertures [::ogf::cat::get lsbg,param,phot-apertures] \
	--mag-zeropoint [::ogf::cat::get lsbg,param,mag-zeropoint] \
	--pixel-scale [::ogf::cat::get lsbg,param,pixel-scale] \
	--mu-eff-min [::ogf::cat::get lsbg,param,mu-eff-min] \
	--mu-eff-max [::ogf::cat::get lsbg,param,mu-eff-max] \
	--r-eff-min [::ogf::cat::get lsbg,param,r-eff-min] \
	--r-eff-max [::ogf::cat::get lsbg,param,r-eff-max] \
	--ellipticity-max [::ogf::cat::get lsbg,param,ellipticity-max] \
	--min-snr [::ogf::cat::get lsbg,param,min-snr] \
	--sersic-n-min [::ogf::cat::get lsbg,param,sersic-n-min] \
	--sersic-n-max [::ogf::cat::get lsbg,param,sersic-n-max] \
	--sersic-re-min [::ogf::cat::get lsbg,param,sersic-re-min] \
	--sersic-cutout-scale [::ogf::cat::get lsbg,param,sersic-cutout-scale] \
	--sersic-max-nfev [::ogf::cat::get lsbg,param,sersic-max-nfev] \
	--sersic-n-filter-min [::ogf::cat::get lsbg,param,sersic-n-filter-min] \
	--sersic-n-filter-max [::ogf::cat::get lsbg,param,sersic-n-filter-max] \
	--sersic-chi2-max [::ogf::cat::get lsbg,param,sersic-chi2-max] \
	--n-workers [::ogf::cat::get param,n-workers]]
    if {[::ogf::cat::get lsbg,param,sersic-fit]} {
	lappend args --sersic-fit
    } else {
	lappend args --no-sersic-fit
    }

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	::ogf::cat::set status "LSBG filter error: $err"
	return
    }

    ::ogf::cat::set alldata $result
    ::ogf::cat::set lsbg,has_catalog 1
    CatalogPanelLoadTSV [::ogf::cat::tsv] "lsbg"

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    CatalogPanelMarkAll
    ::ogf::cat::set status "LSBG: $nsrc candidates passed filtering"
}

proc CatalogPanelLSBGSVMClassify {} {

    if {![::ogf::cat::has]} {
	::ogf::cat::set status "LSBG SVM: No catalog loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	::ogf::cat::set status "LSBG SVM: ds9_lsbg.py not found"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	::ogf::cat::set status "LSBG SVM: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    ::ogf::cat::set status "LSBG SVM: Classifying candidates..."
    update idletasks

    # Save current catalog to temp file
    set tmpcat [file join [file normalize ~] .ds9 lsbg_svm_tmp.tsv]
    if {[catch {
	set fd [open $tmpcat w]
	puts $fd [::ogf::cat::tsv]
	close $fd
    } err]} {
	::ogf::cat::set status "LSBG SVM: Cannot write temp catalog"
	return
    }

    set args [list [OGFPython] $script $fn --mode svm-classify \
	--catalog $tmpcat \
	--svm-threshold [::ogf::cat::get lsbg,param,svm-threshold]]
    if {[::ogf::cat::get lsbg,param,svm-checkpoint] ne {}} {
	lappend args --svm-checkpoint [::ogf::cat::get lsbg,param,svm-checkpoint]
    }

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	::ogf::cat::set status "LSBG SVM error: $err"
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
    ::ogf::cat::set status "LSBG SVM: $n_lsbg LSBG / $n_total total"
}

proc CatalogPanelLSBGSVMColorMarkers {result_tsv} {
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

    if {![::ogf::cat::has]} return

    set alllines [split [::ogf::cat::tsv] \n]
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

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	::ogf::cat::set status "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	::ogf::cat::set status "LSBG: ds9_lsbg.py not found"
	return
    }

    OGFMaskEnsure lsbg
    ::ogf::cat::set status "LSBG: Running full pipeline..."
    update idletasks

    # Reset cmdlog for new session
    ::ogf::cat::set lsbg,cmdlog {}

    set args [list [OGFPython] $script $fn --mode run \
	--mask-detect-thresh [::ogf::cat::get lsbg,param,mask-detect-thresh] \
	--mask-detect-minarea [::ogf::cat::get lsbg,param,mask-detect-minarea] \
	--mask-expand-factor [::ogf::cat::get lsbg,param,mask-expand-factor] \
	--max-dilate-radius [::ogf::cat::get lsbg,param,max-dilate-radius] \
	--bright-star-mag-limit [::ogf::cat::get lsbg,param,bright-star-mag-limit] \
	--bright-star-radius-scale [::ogf::cat::get lsbg,param,bright-star-radius-scale] \
	--mask-mag-threshold [::ogf::cat::get lsbg,param,mask-mag-threshold] \
	--interp-method [::ogf::cat::get lsbg,param,interp-method] \
	--lsb-mu-threshold [::ogf::cat::get lsbg,param,lsb-mu-threshold] \
	--bkg-method [::ogf::cat::get lsbg,param,bkg-method] \
	--bkg-mesh-size [::ogf::cat::get lsbg,param,bkg-mesh-size] \
	--bkg-poly-order [::ogf::cat::get lsbg,param,bkg-poly-order] \
	--bkg-sigma-clip [::ogf::cat::get lsbg,param,bkg-sigma-clip] \
	--bkg-n-iterations [::ogf::cat::get lsbg,param,bkg-n-iterations] \
	--bkg-refine-thresh [::ogf::cat::get lsbg,param,bkg-refine-thresh] \
	--bkg-rms-quantile [::ogf::cat::get lsbg,param,bkg-rms-quantile] \
	--bkg-convergence-tol [::ogf::cat::get lsbg,param,bkg-convergence-tol] \
	--detect-thresh [::ogf::cat::get lsbg,param,detect-thresh] \
	--detect-minarea [::ogf::cat::get lsbg,param,detect-minarea] \
	--detect-filter-kernel [::ogf::cat::get lsbg,param,detect-filter-kernel] \
	--deblend-nthresh [::ogf::cat::get lsbg,param,deblend-nthresh] \
	--deblend-mincont [::ogf::cat::get lsbg,param,deblend-mincont] \
	--multiscale-factors [::ogf::cat::get lsbg,param,multiscale-factors] \
	--sersic-n-min [::ogf::cat::get lsbg,param,sersic-n-min] \
	--sersic-n-max [::ogf::cat::get lsbg,param,sersic-n-max] \
	--sersic-re-min [::ogf::cat::get lsbg,param,sersic-re-min] \
	--sersic-cutout-scale [::ogf::cat::get lsbg,param,sersic-cutout-scale] \
	--sersic-max-nfev [::ogf::cat::get lsbg,param,sersic-max-nfev] \
	--phot-apertures [::ogf::cat::get lsbg,param,phot-apertures] \
	--mag-zeropoint [::ogf::cat::get lsbg,param,mag-zeropoint] \
	--pixel-scale [::ogf::cat::get lsbg,param,pixel-scale] \
	--mu-eff-min [::ogf::cat::get lsbg,param,mu-eff-min] \
	--mu-eff-max [::ogf::cat::get lsbg,param,mu-eff-max] \
	--r-eff-min [::ogf::cat::get lsbg,param,r-eff-min] \
	--r-eff-max [::ogf::cat::get lsbg,param,r-eff-max] \
	--ellipticity-max [::ogf::cat::get lsbg,param,ellipticity-max] \
	--min-snr [::ogf::cat::get lsbg,param,min-snr] \
	--sersic-n-filter-min [::ogf::cat::get lsbg,param,sersic-n-filter-min] \
	--sersic-n-filter-max [::ogf::cat::get lsbg,param,sersic-n-filter-max] \
	--sersic-chi2-max [::ogf::cat::get lsbg,param,sersic-chi2-max] \
	--mask-input [::ogf::cat::get lsbg,mask_file] \
	--mask-output [OGFMaskRefinedPath lsbg] \
	--masked-output [::ogf::cat::get lsbg,masked_file] \
	--bkg-output [::ogf::cat::get lsbg,bkg_file] \
	--cleaned-output [::ogf::cat::get lsbg,cleaned_file] \
	--segmap-output [::ogf::cat::get lsbg,segmap_file] \
	--catalog-output [::ogf::cat::get lsbg,catalog_file] \
	--n-workers [::ogf::cat::get param,n-workers]]
    if {[::ogf::cat::get lsbg,param,lsb-protect]} {
	lappend args --lsb-protect
    } else {
	lappend args --no-lsb-protect
    }
    if {[::ogf::cat::get lsbg,param,multiscale]} {
	lappend args --multiscale
    } else {
	lappend args --no-multiscale
    }
    if {[::ogf::cat::get lsbg,param,sersic-fit]} {
	lappend args --sersic-fit
    } else {
	lappend args --no-sersic-fit
    }
    if {[::ogf::cat::get lsbg,param,svm-classify]} {
	lappend args --svm-classify
	lappend args --svm-threshold [::ogf::cat::get lsbg,param,svm-threshold]
	if {[::ogf::cat::get lsbg,param,svm-checkpoint] ne {}} {
	    lappend args --svm-checkpoint [::ogf::cat::get lsbg,param,svm-checkpoint]
	}
    }

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	::ogf::cat::set status "LSBG pipeline error: $err"
	return
    }

    # Update state
    ::ogf::cat::set lsbg,has_mask 1
    ::ogf::cat::set lsbg,has_clean 1
    ::ogf::cat::set lsbg,has_detect 1
    ::ogf::cat::set lsbg,has_catalog 1

    # Load cleaned image in new frame
    CreateFrame
    if {[catch {LoadFitsFile [::ogf::cat::get lsbg,cleaned_file] {} {}} err]} {
	::ogf::cat::set status "LSBG: Warning — could not load cleaned image"
    } else {
	global scale
	set scale(mode) zscale
	ChangeScaleMode
    }

    # Load catalog
    ::ogf::cat::set alldata $result
    CatalogPanelLoadTSV [::ogf::cat::tsv] "lsbg"

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    # Mark all LSBG candidates
    CatalogPanelMarkAll
    ::ogf::cat::set status "LSBG: Pipeline complete — $nsrc candidates"
}

proc CatalogPanelLSBGForcedPhot {} {
    global current ogflsbg

    # Check that we have a catalog
    if {![::ogf::cat::has]} {
	::ogf::cat::set status "LSBG: No catalog — run pipeline first"
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
    ttk::entry $w.e -textvariable ogflsbg(bandname)
    set ogflsbg(bandname) ""
    ttk::frame $w.btns
    ttk::button $w.btns.ok -text "OK" -command [list set ogflsbg(band_done) 1]
    ttk::button $w.btns.cancel -text "Cancel" -command [list set ogflsbg(band_done) 0]
    pack $w.l -padx 10 -pady 5
    pack $w.e -padx 10 -fill x
    pack $w.btns -pady 5
    pack $w.btns.ok $w.btns.cancel -side left -padx 10
    focus $w.e
    bind $w.e <Return> [list set ogflsbg(band_done) 1]
    tkwait variable ogflsbg(band_done)
    set band_name $ogflsbg(bandname)
    catch {destroy $w}
    if {!$ogflsbg(band_done) || $band_name eq {}} return

    # Save current catalog to temp file
    set tmpcat [CatalogPanelSaveTempCatalog lsbg_forced]
    if {$tmpcat eq {}} {
	::ogf::cat::set status "LSBG: Failed to save temp catalog"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	::ogf::cat::set status "LSBG: ds9_lsbg.py not found"
	return
    }

    ::ogf::cat::set status "LSBG: Forced photometry ($band_name)..."
    update idletasks

    set args [list [OGFPython] $script $band_fits --mode forced \
	--catalog $tmpcat \
	--band-name $band_name \
	--mag-zeropoint [::ogf::cat::get lsbg,param,mag-zeropoint] \
	--pixel-scale [::ogf::cat::get lsbg,param,pixel-scale] \
	--n-workers [::ogf::cat::get param,n-workers]]

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	::ogf::cat::set status "LSBG forced phot error: $err"
	return
    }

    # Merge result columns into current catalog
    set col_names [list FLUX_$band_name FLUXERR_$band_name MAG_$band_name MAGERR_$band_name]
    CatalogPanelAddColumnsFromTSV $result $col_names

    # Load the other band image in a new frame
    CreateFrame
    if {[catch {LoadFitsFile $band_fits {} {}} err]} {
	::ogf::cat::set status "LSBG: Forced phot done but could not display band image"
    } else {
	global scale
	set scale(mode) zscale
	ChangeScaleMode
    }

    ::ogf::cat::set status "LSBG: Forced photometry ($band_name) complete — columns added"
}


# hooks of the headless "full pipeline" step (::ogf::step::run lsbg run_all headless)
proc OGFLsbgHeadlessBefore {} {
    set fn [CatalogPanelGetFITS]
    CatalogPanelLSBGUpdateFiles $fn
    OGFMaskEnsure lsbg
    ::ogf::cat::set lsbg,cmdlog {}
}
proc OGFLsbgHeadlessAfter {} {
    foreach k {has_mask has_clean has_detect has_catalog} {::ogf::cat::set lsbg,$k 1}
    set result [::ogf::cat::tsv]       ;# the step runner has loaded the catalog; a new frame clears the table, so (like CatalogPanelLSBGRunAll) load it again afterwards
    CreateFrame
    set cleaned [file join [OGFSessWorkDir] lsbg_cleaned_[CatalogPanelFitsBaseName [CatalogPanelGetFITS]].fits]
    if {![file exists $cleaned] || [catch {LoadFitsFile $cleaned {} {}} err]} {
	::ogf::cat::set status "LSBG: Warning \u2014 could not load cleaned image"
    } else {
	global scale
	set scale(mode) zscale
	ChangeScaleMode
    }
    ::ogf::cat::set alldata $result
    CatalogPanelLoadTSV [::ogf::cat::tsv] "lsbg"
    CatalogPanelMarkAll
    ::ogf::cat::set status "LSBG: Pipeline complete \u2014 [::ogf::cat::nrows] candidates"
}
