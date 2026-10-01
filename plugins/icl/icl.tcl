# icl/icl.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/icl/plugin.json.

proc CatalogPanelICLUpdateFiles {fn} {
    global catpanel
    if {$fn eq {}} return
    set base [CatalogPanelFitsBaseName $fn]
    if {$base eq {}} return
    # Skip if already set for this base
    if {[info exists catpanel(icl,fits_base)] &&
	$catpanel(icl,fits_base) eq $base} return
    set catpanel(icl,fits_base) $base
    set ds9dir [file join [file normalize ~] .ds9]
    # shared mask (ds9_mask.py): boolean view + on-demand interpolated image
    set catpanel(icl,mask_file)    [file join $ds9dir "mask_${base}_bool.fits"]
    set catpanel(icl,masked_file)  [file join $ds9dir "mask_${base}_masked.fits"]
    set catpanel(icl,bkg_file)     [file join $ds9dir "icl_background_${base}.fits"]
    set catpanel(icl,bgsub_file)   [file join $ds9dir "icl_bgsub_${base}.fits"]
    set catpanel(icl,profile_file) [file join $ds9dir "icl_profile_${base}.tsv"]
    # Detect if previous results exist for this FITS
    set catpanel(icl,has_mask)    [file exists $catpanel(icl,mask_file)]
    set catpanel(icl,has_bkg)     [file exists $catpanel(icl,bkg_file)]
    set catpanel(icl,has_profile) [file exists $catpanel(icl,profile_file)]
}

proc CatalogPanelICLParamLoad {} {
    global catpanel

    set preffile [file join [file normalize ~] .ds9 icl.prf]
    if {![file exists $preffile]} return
    if {[catch {set fd [open $preffile r]} err]} return
    while {[gets $fd line] >= 0} {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set parts [split $line]
	if {[llength $parts] >= 2} {
	    set key [lindex $parts 0]
	    set val [lindex $parts 1]
	    if {[info exists catpanel(icl,param,$key)]} {
		set catpanel(icl,param,$key) $val
	    }
	}
    }
    close $fd
}

proc CatalogPanelICLParamSave {} {
    global catpanel

    set prefdir [file join [file normalize ~] .ds9]
    if {![file isdirectory $prefdir]} {
	file mkdir $prefdir
    }
    set preffile [file join $prefdir icl.prf]
    if {[catch {set fd [open $preffile w]} err]} return
    foreach pname {expand-factor max-dilate-radius \
		   bright-star-mag-limit bright-star-radius-scale \
		   interp-method detect-thresh \
		   bkg-method bkg-order bkg-sigma-clip bkg-sep-mesh \
		   bkg-iterative bkg-n-iterations bkg-convergence-tol bkg-refine-thresh \
		   rmin rmax nsteps spacing ellipticity pa \
		   mag-zeropoint pixel-scale \
		   mu-threshold mu-levels measure-radius} {
	puts $fd "$pname $catpanel(icl,param,$pname)"
    }
    close $fd
}

proc CatalogPanelICLMask {} {
    OGFMaskPipelineMask icl
}

proc CatalogPanelICLViewMask {} {
    global ogfmask
    set ogfmask(overlay) 1
    CatalogPanelMaskToggleOverlay
}

proc CatalogPanelICLSaveMask {} {
    CatalogPanelMaskSaveAs
}

proc CatalogPanelICLImportMask {} {
    CatalogPanelMaskImport
}

proc CatalogPanelICLBackground {method} {
    global catpanel

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "ICL: No FITS file loaded"
	return
    }
    CatalogPanelICLUpdateFiles $fn

    # Shared mask -> interpolated image produced now (on demand); raw if no mask
    set input $fn
    if {[OGFMaskExists]} {
	set input [OGFMaskMaskedFor icl $fn]
    }

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Fitting background ($method)..."
    update idletasks

    set args [list [OGFPython] $script $input --mode background \
	--bkg-method $method \
	--bkg-order $catpanel(icl,param,bkg-order) \
	--bkg-sigma-clip $catpanel(icl,param,bkg-sigma-clip) \
	--bkg-sep-mesh $catpanel(icl,param,bkg-sep-mesh) \
	--bkg-output $catpanel(icl,bkg_file) \
	--bgsub-output $catpanel(icl,bgsub_file)]

    if {[OGFMaskExists]} {
	lappend args --mask [OGFMaskBoolPath]
    }

    # Iterative background refinement
    if {$catpanel(icl,param,bkg-iterative)} {
	lappend args --iterative \
	    --interp-method $catpanel(icl,param,interp-method) \
	    --bkg-n-iterations $catpanel(icl,param,bkg-n-iterations) \
	    --bkg-convergence-tol $catpanel(icl,param,bkg-convergence-tol) \
	    --bkg-refine-thresh $catpanel(icl,param,bkg-refine-thresh) \
	    --mask-output [OGFMaskRefinedPath icl]
    }

    CatalogPanelCmdLog icl $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL background error: $err"
	return
    }

    set catpanel(icl,has_bkg) 1
    if {$catpanel(icl,param,bkg-iterative) && [file exists [OGFMaskRefinedPath icl]]} {
	# iterative refinement result feeds the profile step; shared mask untouched
	set catpanel(icl,mask_file) [OGFMaskRefinedPath icl]
    }

    # Auto-display bgsub image in new frame.  (CreateFrame resets the
    # per-frame panel state, so grab the path first.)
    set _bgsub $catpanel(icl,bgsub_file)
    if {[file exists $_bgsub]} {
	CreateFrame
	if {![catch {LoadFitsFile $_bgsub {} {}}]} {
	    global scale
	    set scale(mode) zscale
	    ChangeScaleMode
	}
    }

    set catpanel(status) "ICL: Background model ($method) complete"
}

proc CatalogPanelICLViewBkg {} {
    global catpanel

    if {![file exists $catpanel(icl,bgsub_file)]} {
	set catpanel(status) "ICL: No background model — run Background Model first"
	return
    }

    CreateFrame
    if {[catch {LoadFitsFile $catpanel(icl,bgsub_file) {} {}} err]} {
	set catpanel(status) "ICL: Error loading background-subtracted image: $err"
	return
    }
    global scale
    set scale(mode) zscale
    ChangeScaleMode
    set catpanel(icl,has_bkg) 1
    set catpanel(status) "ICL: Background-subtracted image loaded in new frame"
}

proc CatalogPanelICLSetCenter {} {
    global catpanel catpanel_fdata ds9 current

    # Use the first frame's catalog (original image SExtractor result)
    set catalog_data {}
    set first_frame [lindex $ds9(frames) 0]
    if {$first_frame ne {} &&
	[info exists catpanel_fdata($first_frame,alldata)] &&
	$catpanel_fdata($first_frame,alldata) ne {}} {
	set catalog_data $catpanel_fdata($first_frame,alldata)
    } elseif {[info exists catpanel(alldata)] && $catpanel(alldata) ne {}} {
	set catalog_data $catpanel(alldata)
    }

    if {$catalog_data eq {}} {
	set catpanel(status) "ICL: No SExtractor catalog found — run SExtract first"
	return
    }

    # Parse headers to find NUMBER, X_IMAGE, Y_IMAGE, MAG_AUTO
    set lines [split $catalog_data \n]
    set headers [split [lindex $lines 0] "\t"]
    set numcol -1
    set xcol -1
    set ycol -1
    set magcol -1
    for {set c 0} {$c < [llength $headers]} {incr c} {
	set h [string trim [lindex $headers $c]]
	if {$h eq "NUMBER"} { set numcol $c }
	if {$h eq "X_IMAGE"} { set xcol $c }
	if {$h eq "Y_IMAGE"} { set ycol $c }
	if {$h eq "MAG_AUTO"} { set magcol $c }
    }
    if {$numcol < 0 || $xcol < 0 || $ycol < 0} {
	set catpanel(status) "ICL: Catalog missing NUMBER/X_IMAGE/Y_IMAGE"
	return
    }

    # Find brightest source as default suggestion
    set best_id {}
    set best_mag 99.0
    for {set i 1} {$i < [llength $lines]} {incr i} {
	set row [split [lindex $lines $i] "\t"]
	if {[llength $row] <= $xcol} continue
	if {$magcol >= 0 && [llength $row] > $magcol} {
	    set m [string trim [lindex $row $magcol]]
	    if {[string is double $m] && $m < $best_mag} {
		set best_mag $m
		set best_id [string trim [lindex $row $numcol]]
	    }
	}
    }

    # Show dialog
    set w .icl_bcg_id
    if {[winfo exists $w]} { destroy $w }
    toplevel $w
    wm title $w "Set BCG Center"
    wm geometry $w 320x120
    wm resizable $w 0 0

    set default_id $best_id
    if {$default_id eq {}} { set default_id 1 }

    ttk::label $w.lbl -text "Enter source NUMBER from SExtractor catalog:"
    pack $w.lbl -padx 10 -pady {10 2} -anchor w

    ttk::frame $w.ef
    pack $w.ef -padx 10 -pady 2 -fill x
    ttk::label $w.ef.idlbl -text "Source ID:"
    ttk::entry $w.ef.entry -width 10 -textvariable icl_bcg_entry_id
    set ::icl_bcg_entry_id $default_id
    pack $w.ef.idlbl -side left -padx {0 5}
    pack $w.ef.entry -side left

    if {$best_id ne {} && $magcol >= 0} {
	ttk::label $w.ef.hint -text "(brightest: #${best_id}, mag=[format %.1f $best_mag])" \
	    -foreground gray
	pack $w.ef.hint -side left -padx 5
    }

    ttk::frame $w.bf
    pack $w.bf -padx 10 -pady 8 -fill x
    ttk::button $w.bf.ok -text "OK" \
	-command [list CatalogPanelICLSetCenterByID $w [list $catalog_data] $numcol $xcol $ycol]
    ttk::button $w.bf.click -text "Click on Image" \
	-command [list CatalogPanelICLSetCenterClickMode $w]
    ttk::button $w.bf.cancel -text "Cancel" -command [list destroy $w]
    pack $w.bf.cancel -side right -padx 2
    pack $w.bf.click -side right -padx 2
    pack $w.bf.ok -side right -padx 2

    bind $w.ef.entry <Return> [list $w.bf.ok invoke]
    focus $w.ef.entry
    $w.ef.entry selection range 0 end
}

proc CatalogPanelICLSetCenterByID {w catalog_data numcol xcol ycol} {
    global catpanel current

    set src_id [string trim $::icl_bcg_entry_id]
    if {$src_id eq {} || ![string is integer $src_id]} {
	set catpanel(status) "ICL: Enter a valid source NUMBER"
	return
    }

    # Search catalog for this NUMBER
    set lines [split $catalog_data \n]
    for {set i 1} {$i < [llength $lines]} {incr i} {
	set row [split [lindex $lines $i] "\t"]
	if {[llength $row] <= $xcol || [llength $row] <= $ycol} continue
	set n [string trim [lindex $row $numcol]]
	if {$n eq $src_id} {
	    set ix [string trim [lindex $row $xcol]]
	    set iy [string trim [lindex $row $ycol]]

	    # Store as 0-indexed
	    set catpanel(icl,center_x) [expr {$ix - 1.0}]
	    set catpanel(icl,center_y) [expr {$iy - 1.0}]

	    # Draw cyan cross
	    set frame $current(frame)
	    if {$frame ne {}} {
		catch {$frame marker catalog icl_bcg delete}
		global icl_bcg_reg
		set icl_bcg_reg "image\ncross point([format %.1f $ix] [format %.1f $iy]) # color=cyan width=2 point=cross 20 tag={icl_bcg} select=0 edit=0 move=0 rotate=0 delete=1\n"
		catch {$frame marker catalog command ds9 var icl_bcg_reg}
	    }

	    set catpanel(status) "ICL: BCG center set to source #$src_id ([format %.1f $ix], [format %.1f $iy])"
	    destroy $w
	    return
	}
    }

    set catpanel(status) "ICL: Source #$src_id not found in catalog"
}

proc CatalogPanelICLSetCenterClickMode {w} {
    global catpanel
    destroy $w
    set catpanel(icl,click_mode) 1
    set catpanel(status) "ICL: Click on image to set BCG center..."
}

proc CatalogPanelICLClickSetCenter {frame x y} {
    global catpanel

    # Convert canvas coords to image coords (1-based)
    set imgc [$frame get coordinates $x $y image]
    set ix [lindex $imgc 0]
    set iy [lindex $imgc 1]

    # Store as 0-indexed (Python convention)
    set catpanel(icl,center_x) [expr {$ix - 1.0}]
    set catpanel(icl,center_y) [expr {$iy - 1.0}]

    # Draw a cyan cross at BCG center
    catch {$frame marker catalog icl_bcg delete}

    global icl_bcg_reg
    set icl_bcg_reg "image\ncross point([format %.1f $ix] [format %.1f $iy]) # color=cyan width=2 point=cross 20 tag={icl_bcg} select=0 edit=0 move=0 rotate=0 delete=1\n"
    catch {$frame marker catalog command ds9 var icl_bcg_reg}

    set catpanel(icl,click_mode) 0
    set catpanel(status) "ICL: BCG center set to ([format %.1f $ix], [format %.1f $iy])"
}

proc CatalogPanelICLProfile {} {
    global catpanel

    if {$catpanel(icl,center_x) eq {} || $catpanel(icl,center_y) eq {}} {
	set catpanel(status) "ICL: Set BCG center first"
	return
    }

    # Use bgsub image if available, otherwise raw
    set fn [CatalogPanelGetFITS]
    if {[file exists $catpanel(icl,bgsub_file)]} {
	set fn $catpanel(icl,bgsub_file)
    }
    if {$fn eq {}} {
	set catpanel(status) "ICL: No image available"
	return
    }
    CatalogPanelICLUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Measuring SB profile..."
    update idletasks

    set center "$catpanel(icl,center_x),$catpanel(icl,center_y)"

    set args [list [OGFPython] $script $fn --mode profile \
	--center $center \
	--rmin $catpanel(icl,param,rmin) \
	--rmax $catpanel(icl,param,rmax) \
	--nsteps $catpanel(icl,param,nsteps) \
	--spacing $catpanel(icl,param,spacing) \
	--ellipticity $catpanel(icl,param,ellipticity) \
	--pa $catpanel(icl,param,pa) \
	--mag-zeropoint $catpanel(icl,param,mag-zeropoint) \
	--pixel-scale $catpanel(icl,param,pixel-scale) \
	--profile-output $catpanel(icl,profile_file)]

    CatalogPanelCmdLog icl $args
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL profile error: $err"
	return
    }

    set catpanel(icl,has_profile) 1

    # Display profile in catalog table
    CatalogPanelLoadTSV $data "icl_profile"

    # Draw annulus markers
    CatalogPanelICLDrawAnnuli

    set catpanel(status) "ICL: SB profile measured"
}

proc CatalogPanelICLDrawAnnuli {} {
    global catpanel current

    set frame $current(frame)
    if {$frame eq {}} return

    # Delete existing annulus markers
    catch {$frame marker catalog icl_annulus delete}

    # BCG center (0-indexed → 1-indexed for ds9 markers)
    set cx [expr {$catpanel(icl,center_x) + 1.0}]
    set cy [expr {$catpanel(icl,center_y) + 1.0}]

    # Draw BCG center cross
    global icl_ann_reg
    set icl_ann_reg "image\ncross point($cx $cy) # color=cyan width=2 point=cross 20 tag={icl_annulus} select=0 edit=0 move=0 rotate=0 delete=1\n"
    catch {$frame marker catalog command ds9 var icl_ann_reg}

    # Draw annulus rings at 25%, 50%, 75%, 100% of rmax
    set rmin $catpanel(icl,param,rmin)
    set rmax $catpanel(icl,param,rmax)
    foreach frac {0.25 0.50 0.75 1.00} {
	set r [expr {$rmin + ($rmax - $rmin) * $frac}]
	set ri [expr {int($r)}]
	set icl_ann_reg "image\ncircle($cx $cy ${r}i) # color=green dash=1 width=1 tag={icl_annulus} select=0 edit=0 move=0 rotate=0 delete=1 text={r=${ri}}\n"
	catch {$frame marker catalog command ds9 var icl_ann_reg}
    }
}

proc CatalogPanelICLSectorProfile {} {
    global catpanel ed

    if {$catpanel(icl,center_x) eq {} || $catpanel(icl,center_y) eq {}} {
	set catpanel(status) "ICL: Set BCG center first"
	return
    }

    set w .iclsector
    if {[winfo exists $w]} {
	raise $w
	return
    }

    toplevel $w
    wm title $w "ICL Sector Profile"
    wm geometry $w 300x180

    set ed(icl,sector-pa) 0.0
    set ed(icl,sector-width) 90.0

    set f [ttk::frame $w.content]
    pack $f -fill both -expand true -padx 8 -pady 8

    set r 0
    ttk::label $f.lpa -text "Sector PA (deg):"
    ttk::entry $f.epa -textvariable ed(icl,sector-pa) -width 10
    grid $f.lpa -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $f.epa -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $f.lw -text "Sector Width (deg):"
    ttk::entry $f.ew -textvariable ed(icl,sector-width) -width 10
    grid $f.lw -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $f.ew -row $r -column 1 -sticky w -padx 4 -pady 4

    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8

    ttk::button $bf.run -text "Measure" -command [list CatalogPanelICLSectorProfileRun $w]
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.run -side right -padx 4
}

proc CatalogPanelICLSectorProfileRun {w} {
    global catpanel ed

    set fn [CatalogPanelGetFITS]
    if {[file exists $catpanel(icl,bgsub_file)]} {
	set fn $catpanel(icl,bgsub_file)
    }
    if {$fn eq {}} {
	set catpanel(status) "ICL: No image available"
	return
    }
    CatalogPanelICLUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Measuring sector profile..."
    update idletasks

    set center "$catpanel(icl,center_x),$catpanel(icl,center_y)"

    set args [list [OGFPython] $script $fn --mode profile \
	--center $center \
	--rmin $catpanel(icl,param,rmin) \
	--rmax $catpanel(icl,param,rmax) \
	--nsteps $catpanel(icl,param,nsteps) \
	--spacing $catpanel(icl,param,spacing) \
	--mag-zeropoint $catpanel(icl,param,mag-zeropoint) \
	--pixel-scale $catpanel(icl,param,pixel-scale) \
	--sector-pa $ed(icl,sector-pa) \
	--sector-width $ed(icl,sector-width) \
	--profile-output $catpanel(icl,profile_file)]

    CatalogPanelCmdLog icl $args
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL sector profile error: $err"
	return
    }

    destroy $w
    set catpanel(icl,has_profile) 1
    CatalogPanelLoadTSV $data "icl_sector_profile"
    set catpanel(status) "ICL: Sector profile measured (PA=$ed(icl,sector-pa), width=$ed(icl,sector-width))"
}

proc CatalogPanelICLMeasure {} {
    global catpanel

    if {!$catpanel(icl,has_profile) || ![file exists $catpanel(icl,profile_file)]} {
	set catpanel(status) "ICL: No profile available — run Measure Profile first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} { set fn "dummy.fits" }
    CatalogPanelICLUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Computing ICL measurements..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode measure \
	--profile-file $catpanel(icl,profile_file) \
	--mu-threshold $catpanel(icl,param,mu-threshold) \
	--mu-levels $catpanel(icl,param,mu-levels) \
	--pixel-scale $catpanel(icl,param,pixel-scale)]

    CatalogPanelCmdLog icl $args
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL measure error: $err"
	return
    }

    CatalogPanelLoadTSV $data "icl_measurements"

    # Draw isophotal radius markers
    CatalogPanelICLDrawIsophotes $data

    set catpanel(status) "ICL: Measurements complete"
}

proc CatalogPanelICLDrawIsophotes {data} {
    global catpanel current

    set frame $current(frame)
    if {$frame eq {}} return
    if {$catpanel(icl,center_x) eq {} || $catpanel(icl,center_y) eq {}} return

    # Delete existing annulus markers
    catch {$frame marker catalog icl_annulus delete}

    set cx [expr {$catpanel(icl,center_x) + 1.0}]
    set cy [expr {$catpanel(icl,center_y) + 1.0}]

    # BCG center cross
    global icl_ann_reg
    set icl_ann_reg "image\ncross point($cx $cy) # color=cyan width=2 point=cross 20 tag={icl_annulus} select=0 edit=0 move=0 rotate=0 delete=1\n"
    catch {$frame marker catalog command ds9 var icl_ann_reg}

    # Parse TSV for R_MU* columns
    set lines [split $data \n]
    if {[llength $lines] < 2} return
    set headers [split [lindex $lines 0] "\t"]
    set vals [split [lindex $lines 1] "\t"]

    for {set c 0} {$c < [llength $headers]} {incr c} {
	set h [string trim [lindex $headers $c]]
	if {[string match "R_MU*" $h] && ![string match "*ARCSEC" $h]} {
	    set v [string trim [lindex $vals $c]]
	    if {$v ne "NaN" && [string is double $v]} {
		set r [expr {double($v)}]
		if {$r > 0 && $r < 1e6} {
		    # Extract mu level from column name (R_MU26 → 26)
		    set mu_label [string range $h 4 end]
		    set icl_ann_reg "image\ncircle($cx $cy ${r}i) # color=magenta dash=1 width=1 tag={icl_annulus} select=0 edit=0 move=0 rotate=0 delete=1 text={mu=$mu_label}\n"
		    catch {$frame marker catalog command ds9 var icl_ann_reg}
		}
	    }
	}
    }
}

proc CatalogPanelICLMeasureMulti {} {
    global catpanel

    if {!$catpanel(icl,has_profile) || ![file exists $catpanel(icl,profile_file)]} {
	set catpanel(status) "ICL: No profile available — run Measure Profile first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} { set fn "dummy.fits" }
    CatalogPanelICLUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Computing multi-threshold ICL..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode measure-multi \
	--profile-file $catpanel(icl,profile_file) \
	--pixel-scale $catpanel(icl,param,pixel-scale)]

    CatalogPanelCmdLog icl $args
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL multi-threshold error: $err"
	return
    }

    CatalogPanelLoadTSV $data "icl_multi_threshold"
    set catpanel(status) "ICL: Multi-threshold measurements complete"
}

proc CatalogPanelICLDecompose {} {
    global catpanel ds9 catpanel_fdata

    if {$catpanel(icl,center_x) eq {} || $catpanel(icl,center_y) eq {}} {
	set catpanel(status) "ICL: Set BCG center first"
	return
    }

    # Use the ORIGINAL image (first frame) for decompose —
    # the masked image has BCG removed, making Sérsic fitting impossible
    set fn {}
    set first_frame [lindex $ds9(frames) 0]
    if {$first_frame ne {}} {
	if {[info exists catpanel_fdata($first_frame,filename)] &&
	    $catpanel_fdata($first_frame,filename) ne {}} {
	    set fn $catpanel_fdata($first_frame,filename)
	}
    }
    if {$fn eq {}} {
	set fn [CatalogPanelGetFITS]
    }
    if {$fn eq {}} {
	set catpanel(status) "ICL: No FITS image available"
	return
    }
    CatalogPanelICLUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: BCG+ICL decomposition (original image)..."
    update idletasks

    # Measure profile on original image and decompose in one step
    set center "$catpanel(icl,center_x),$catpanel(icl,center_y)"

    # First: measure profile on original (unmasked) image
    set prof_args [list [OGFPython] $script $fn --mode profile \
	--center $center \
	--rmin $catpanel(icl,param,rmin) \
	--rmax $catpanel(icl,param,rmax) \
	--nsteps $catpanel(icl,param,nsteps) \
	--spacing $catpanel(icl,param,spacing) \
	--ellipticity $catpanel(icl,param,ellipticity) \
	--pa $catpanel(icl,param,pa) \
	--mag-zeropoint $catpanel(icl,param,mag-zeropoint) \
	--pixel-scale $catpanel(icl,param,pixel-scale) \
	--profile-output $catpanel(icl,profile_file)]

    CatalogPanelCmdLog icl $prof_args
    if {[catch {exec {*}$prof_args 2>@stderr} err]} {
	set catpanel(status) "ICL decompose: profile error: $err"
	return
    }

    # Then: decompose using that profile
    set args [list [OGFPython] $script $fn --mode decompose \
	--profile-file $catpanel(icl,profile_file) \
	--pixel-scale $catpanel(icl,param,pixel-scale) \
	--mag-zeropoint $catpanel(icl,param,mag-zeropoint)]

    CatalogPanelCmdLog icl $args
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL decompose error: $err"
	return
    }

    CatalogPanelLoadTSV $data "icl_decomposition"
    set catpanel(status) "ICL: BCG+ICL decomposition complete"
}

proc CatalogPanelICLColorProfile {} {
    global catpanel ed

    set w .iclcolor
    if {[winfo exists $w]} {
	raise $w
	return
    }

    toplevel $w
    wm title $w "ICL Color Profile"
    wm geometry $w 450x280

    set f [ttk::frame $w.content]
    pack $f -fill both -expand true -padx 8 -pady 8

    ttk::label $f.linfo -text "Specify 2-4 bands for color profile measurement:"
    grid $f.linfo -row 0 -column 0 -columnspan 3 -sticky w -padx 4 -pady 4

    for {set i 1} {$i <= 4} {incr i} {
	set ed(icl,color_name_$i) {}
	set ed(icl,color_file_$i) {}

	ttk::label $f.ln$i -text "Band $i name:"
	ttk::entry $f.en$i -textvariable ed(icl,color_name_$i) -width 8
	ttk::entry $f.ef$i -textvariable ed(icl,color_file_$i) -width 25
	ttk::button $f.bb$i -text "Browse" -command [list CatalogPanelICLColorBrowse $i]

	grid $f.ln$i -row $i -column 0 -sticky w -padx 4 -pady 2
	grid $f.en$i -row $i -column 1 -sticky w -padx 2 -pady 2
	grid $f.ef$i -row $i -column 2 -sticky ew -padx 2 -pady 2
	grid $f.bb$i -row $i -column 3 -sticky w -padx 2 -pady 2
    }

    grid columnconfigure $f 2 -weight 1

    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8
    ttk::button $bf.run -text "Run" -command [list CatalogPanelICLColorProfileRun $w]
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.run -side right -padx 4
}

proc CatalogPanelICLColorBrowse {idx} {
    global ed

    set fname [tk_getOpenFile -filetypes {{{FITS} {.fits .fit .fts}} {{All} *}}]
    if {$fname ne {}} {
	set ed(icl,color_file_$idx) $fname
    }
}

proc CatalogPanelICLColorProfileRun {w} {
    global catpanel ed

    if {$catpanel(icl,center_x) eq {} || $catpanel(icl,center_y) eq {}} {
	set catpanel(status) "ICL: Set BCG center first"
	return
    }

    # Build band spec
    set bands {}
    for {set i 1} {$i <= 4} {incr i} {
	set name [string trim $ed(icl,color_name_$i)]
	set file [string trim $ed(icl,color_file_$i)]
	if {$name ne {} && $file ne {} && [file exists $file]} {
	    if {$bands ne {}} { append bands , }
	    append bands "$name:$file"
	}
    }

    if {$bands eq {}} {
	set catpanel(status) "ICL: Need at least 2 bands with name and file"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} { set fn "dummy.fits" }
    CatalogPanelICLUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Measuring color profile..."
    update idletasks

    set center "$catpanel(icl,center_x),$catpanel(icl,center_y)"

    set args [list [OGFPython] $script $fn --mode color \
	--center $center \
	--bands $bands \
	--rmin $catpanel(icl,param,rmin) \
	--rmax $catpanel(icl,param,rmax) \
	--nsteps $catpanel(icl,param,nsteps) \
	--spacing $catpanel(icl,param,spacing) \
	--mag-zeropoint $catpanel(icl,param,mag-zeropoint) \
	--pixel-scale $catpanel(icl,param,pixel-scale)]

    if {$catpanel(icl,has_mask) && [file exists $catpanel(icl,mask_file)]} {
	lappend args --mask $catpanel(icl,mask_file)
    }

    catch {OGFSessFromCmdLog icl $args}
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL color profile error: $err"
	return
    }

    destroy $w
    CatalogPanelLoadTSV $data "icl_color_profile"
    set catpanel(status) "ICL: Color profile complete"
}

proc CatalogPanelICLSaveProfile {} {
    global catpanel

    if {!$catpanel(icl,has_profile) || ![file exists $catpanel(icl,profile_file)]} {
	set catpanel(status) "ICL: No profile available"
	return
    }

    set fname [tk_getSaveFile -defaultextension .tsv \
	-filetypes {{{TSV} {.tsv}} {{All} *}} \
	-initialfile [file tail $catpanel(icl,profile_file)]]
    if {$fname eq {}} return

    if {[catch {file copy -force $catpanel(icl,profile_file) $fname} err]} {
	set catpanel(status) "ICL: Save error: $err"
	return
    }
    set catpanel(status) "ICL: Profile saved to $fname"
}

proc CatalogPanelICLLoadProfile {} {
    global catpanel

    set fname [tk_getOpenFile -filetypes {{{TSV} {.tsv}} {{All} *}}]
    if {$fname eq {} || ![file exists $fname]} return

    if {[catch {set fd [open $fname r]} err]} {
	set catpanel(status) "ICL: Load error: $err"
	return
    }
    set data [read $fd]
    close $fd

    # Copy to standard location
    if {[catch {file copy -force $fname $catpanel(icl,profile_file)} err]} {
	# Non-fatal
    }

    set catpanel(icl,has_profile) 1
    CatalogPanelLoadTSV [string trim $data] "icl_profile"
    set catpanel(status) "ICL: Profile loaded from $fname"
}

