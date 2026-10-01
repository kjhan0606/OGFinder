# photometry/photometry.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/photometry/plugin.json.

proc CatalogPanelSegmentationMap {} {
    global catpanel current ds9

    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_segmap.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_segmap.py"
	return
    }

    set catpanel(status) "Generating segmentation map..."
    update idletasks

    set args [list [OGFPython] $script $fn]
    if {[info exists catpanel(param,detect-thresh)]} {
	lappend args --detect-thresh $catpanel(param,detect-thresh)
    }
    if {[info exists catpanel(param,detect-minarea)]} {
	lappend args --detect-minarea $catpanel(param,detect-minarea)
    }

    OGFSessLog analysis.segmap auto $args -title {Segmentation map} 
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "Segmentation map error: $err"
	return
    }

    # Parse output: "OK N_SOURCES OUTPUT_PATH"
    set parts [split $data]
    if {[llength $parts] >= 3 && [lindex $parts 0] eq "OK"} {
	set nsrc [lindex $parts 1]
	set outpath [lindex $parts 2]

	# Load in new frame
	if {[catch {
	    CreateFrame
	    LoadFitsFile $outpath {} {}
	} err2]} {
	    set catpanel(status) "Error loading segmap: $err2"
	    return
	}

	set catpanel(status) "Segmentation map: $nsrc sources"
    } else {
	set catpanel(status) "Segmentation map: unexpected output"
    }
}

proc CatalogPanelPSFPhotometry {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    # Check for PSF file
    if {![info exists catpanel(psf,file)] || ![file exists $catpanel(psf,file)]} {
	set catpanel(status) "No PSF file — build PSF first (Reconstruction > PSF Generation)"
	return
    }

    set script [CatalogPanelGetScript ds9_psf_phot.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_psf_phot.py"
	return
    }

    set tmpcat [CatalogPanelSaveTempCatalog "psfphot"]
    if {$tmpcat eq {}} {
	set catpanel(status) "Failed to save temp catalog"
	return
    }

    set catpanel(status) "PSF photometry..."
    update idletasks

    set args [list [OGFPython] $script $fn --catalog $tmpcat --psf $catpanel(psf,file)]
    if {[info exists catpanel(param,mag-zeropoint)]} {
	lappend args --mag-zeropoint $catpanel(param,mag-zeropoint)
    }
    lappend args --n-workers $catpanel(param,n-workers)

    OGFSessLog analysis.psf_phot auto $args -title {PSF photometry} -post [dict create kind add cols_list {FLUX_PSF FLUXERR_PSF MAG_PSF MAGERR_PSF CHI2_PSF X_PSF Y_PSF}] -requires [list psf]
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "PSF photometry error: $err"
	return
    }

    CatalogPanelAddColumnsFromTSV $data \
	{FLUX_PSF FLUXERR_PSF MAG_PSF MAGERR_PSF CHI2_PSF X_PSF Y_PSF}
    set catpanel(status) "PSF photometry complete"
}

proc CatalogPanelMultiBand {} {
    global catpanel ds9

    set w {.multibandphot}
    set ed(ok) 0

    DialogCreate $w {Multi-Band Photometry} ed(ok)

    set f [ttk::frame $w.param]

    # Detection image
    ttk::label $f.ldet -text "Detection Image:" -anchor w
    ttk::entry $f.edet -textvariable ed(mb,detect) -width 40
    ttk::button $f.bdet -text "Browse..." -command {
	set ff [tk_getOpenFile -title "Detection Image" \
	    -filetypes {{{FITS Files} {.fits .fit .fts}} {{All Files} {*}}}]
	if {$ff ne {}} { set ed(mb,detect) $ff }
    }
    grid $f.ldet $f.edet $f.bdet -padx 4 -pady 2 -sticky w

    # Current frame as default
    set ed(mb,detect) [CatalogPanelGetFITS]

    # Band list
    ttk::label $f.lbands -text "Bands (name:file, one per line):" -anchor w
    text $f.tbands -width 50 -height 6
    grid $f.lbands -padx 4 -pady 2 -sticky w -columnspan 3
    grid $f.tbands -padx 4 -pady 2 -sticky we -columnspan 3

    # Buttons
    set bf [ttk::frame $w.buttons]
    ttk::button $bf.ok -text {Run} -command {set ed(ok) 1} -default active
    ttk::button $bf.cancel -text {Cancel} -command {set ed(ok) 0}
    pack $bf.ok $bf.cancel -side left -expand true -padx 2 -pady 4

    bind $w <Return> {set ed(ok) 1}
    ttk::separator $w.sep -orient horizontal
    pack $w.buttons $w.sep -side bottom -fill x
    pack $w.param -side top -fill both -expand true

    DialogWait $w ed(ok) $f.edet
    set detect_img $ed(mb,detect)
    set bands_text [$f.tbands get 1.0 end]
    destroy $w

    if {!$ed(ok)} { unset ed; return }
    unset ed

    if {$detect_img eq {} || ![file exists $detect_img]} {
	set catpanel(status) "Detection image not found"
	return
    }

    # Parse band lines
    set bands {}
    foreach line [split $bands_text \n] {
	set line [string trim $line]
	if {$line eq {}} continue
	lappend bands $line
    }
    if {[llength $bands] == 0} {
	set catpanel(status) "No bands specified"
	return
    }

    set script [CatalogPanelGetScript ds9_multiband.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_multiband.py"
	return
    }

    set catpanel(status) "Multi-band photometry ([llength $bands] bands)..."
    update idletasks

    set args [list [OGFPython] $script --detect-image $detect_img \
	--bands [join $bands ","]]

    # Use existing catalog if available
    if {[info exists catpanel(alldata)] && $catpanel(alldata) ne {}} {
	set tmpcat [CatalogPanelSaveTempCatalog "multiband"]
	if {$tmpcat ne {}} {
	    lappend args --catalog $tmpcat
	}
    }
    if {[info exists catpanel(param,mag-zeropoint)]} {
	lappend args --mag-zeropoint $catpanel(param,mag-zeropoint)
    }
    lappend args --n-workers $catpanel(param,n-workers)

    OGFSessLog analysis.multiband auto $args -title {Multi-band photometry (dialog, ds9_multiband.py)} 
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "Multi-band error: $err"
	return
    }

    # Load as new catalog (replaces current)
    CatalogPanelLoadTSV $data "multi-band"
    set catpanel(status) "Multi-band photometry complete"
}

proc CatalogPanelCrowdedPhot {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    if {![info exists catpanel(psf,file)] || ![file exists $catpanel(psf,file)]} {
	set catpanel(status) "No PSF file — build PSF first"
	return
    }

    set script [CatalogPanelGetScript ds9_crowded_phot.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_crowded_phot.py"
	return
    }

    set tmpcat [CatalogPanelSaveTempCatalog "crowded"]
    if {$tmpcat eq {}} {
	set catpanel(status) "Failed to save temp catalog"
	return
    }

    set catpanel(status) "Crowded field photometry..."
    update idletasks

    set args [list [OGFPython] $script $fn --catalog $tmpcat \
	--psf $catpanel(psf,file)]
    if {[info exists catpanel(param,mag-zeropoint)]} {
	lappend args --mag-zeropoint $catpanel(param,mag-zeropoint)
    }
    lappend args --n-workers $catpanel(param,n-workers)

    OGFSessLog analysis.crowded_phot auto $args -title {Crowded field photometry} -post [dict create kind add cols_list {FLUX_CROWD FLUXERR_CROWD MAG_CROWD X_CROWD Y_CROWD N_NEIGHBORS}] -requires [list psf]
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "Crowded phot error: $err"
	return
    }

    CatalogPanelAddColumnsFromTSV $data \
	{FLUX_CROWD FLUXERR_CROWD MAG_CROWD X_CROWD Y_CROWD N_NEIGHBORS}
    set catpanel(status) "Crowded field photometry complete"
}

proc CatalogPanelCrossMatch {} {
    global catpanel ed

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first"
	return
    }

    set w {.crossmatch}
    set ed(ok) 0

    DialogCreate $w {Cross-Match (VizieR)} ed(ok)

    set f [ttk::frame $w.param]

    # Catalog selection
    ttk::label $f.lcat -text "VizieR Catalog:" -anchor w
    ttk::combobox $f.ecat -textvariable ed(xm,catalog) -width 25 \
	-values {GAIA_DR3 2MASS PanSTARRS_DR1 SDSS_DR17 ALLWISE} -state readonly
    set ed(xm,catalog) GAIA_DR3
    grid $f.lcat $f.ecat -padx 4 -pady 2 -sticky w

    # Match radius
    ttk::label $f.lrad -text "Match Radius (arcsec):" -anchor w
    ttk::entry $f.erad -textvariable ed(xm,radius) -width 12
    set ed(xm,radius) 2.0
    grid $f.lrad $f.erad -padx 4 -pady 2 -sticky w

    # Buttons
    set bf [ttk::frame $w.buttons]
    ttk::button $bf.ok -text {Run} -command {set ed(ok) 1} -default active
    ttk::button $bf.cancel -text {Cancel} -command {set ed(ok) 0}
    pack $bf.ok $bf.cancel -side left -expand true -padx 2 -pady 4

    bind $w <Return> {set ed(ok) 1}
    ttk::separator $w.sep -orient horizontal
    pack $w.buttons $w.sep -side bottom -fill x
    pack $w.param -side top -fill both -expand true

    DialogWait $w ed(ok) $f.ecat
    set vizcat $ed(xm,catalog)
    set matchrad $ed(xm,radius)
    destroy $w

    if {!$ed(ok)} { unset ed; return }
    unset ed

    set script [CatalogPanelGetScript ds9_crossmatch.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_crossmatch.py"
	return
    }

    set tmpcat [CatalogPanelSaveTempCatalog "crossmatch"]
    if {$tmpcat eq {}} {
	set catpanel(status) "Failed to save temp catalog"
	return
    }

    set catpanel(status) "Cross-matching with $vizcat (r=${matchrad}\")..."
    update idletasks

    OGFSessLog analysis.crossmatch auto [list [OGFPython] $script --catalog $tmpcat --vizier-cat $vizcat --radius $matchrad] -network 1 -title "Cross-match with $vizcat" -post [dict create kind add cols_list {MATCH_DIST MATCH_ID}]
    if {[catch {
	set data [exec [OGFPython] $script --catalog $tmpcat \
	    --vizier-cat $vizcat --radius $matchrad 2>@stderr]
    } err]} {
	set catpanel(status) "Cross-match error: $err"
	return
    }

    CatalogPanelAddColumnsFromTSV $data {MATCH_DIST MATCH_ID}
    set catpanel(status) "Cross-match complete ($vizcat)"
}

proc CatalogPanelCompleteness {} {
    global catpanel ed

    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set w {.completeness}
    set ed(ok) 0

    DialogCreate $w {Completeness Simulation} ed(ok)

    set f [ttk::frame $w.param]

    set ed(comp,ninject) 500
    set ed(comp,magmin) 20.0
    set ed(comp,magmax) 28.0
    set ed(comp,nbins) 16

    set row 0
    foreach {vname vlabel} {
	comp,ninject {N Inject (per bin)}
	comp,magmin {Mag Min}
	comp,magmax {Mag Max}
	comp,nbins {N Bins}
    } {
	ttk::label $f.l$row -text "$vlabel:" -anchor w
	ttk::entry $f.e$row -textvariable ed($vname) -width 12
	grid $f.l$row $f.e$row -padx 4 -pady 2 -sticky w
	incr row
    }

    set bf [ttk::frame $w.buttons]
    ttk::button $bf.ok -text {Run} -command {set ed(ok) 1} -default active
    ttk::button $bf.cancel -text {Cancel} -command {set ed(ok) 0}
    pack $bf.ok $bf.cancel -side left -expand true -padx 2 -pady 4

    bind $w <Return> {set ed(ok) 1}
    ttk::separator $w.sep -orient horizontal
    pack $w.buttons $w.sep -side bottom -fill x
    pack $w.param -side top -fill both -expand true

    DialogWait $w ed(ok) $f.e0
    set ninject $ed(comp,ninject)
    set magmin $ed(comp,magmin)
    set magmax $ed(comp,magmax)
    set nbins $ed(comp,nbins)
    destroy $w

    if {!$ed(ok)} { unset ed; return }
    unset ed

    set script [CatalogPanelGetScript ds9_completeness.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_completeness.py"
	return
    }

    set catpanel(status) "Completeness simulation ($nbins bins)..."
    update idletasks

    set args [list [OGFPython] $script $fn \
	--n-inject $ninject --mag-min $magmin --mag-max $magmax --n-bins $nbins]
    if {[info exists catpanel(param,detect-thresh)]} {
	lappend args --detect-thresh $catpanel(param,detect-thresh)
    }
    if {[info exists catpanel(param,mag-zeropoint)]} {
	lappend args --mag-zeropoint $catpanel(param,mag-zeropoint)
    }
    lappend args --n-workers $catpanel(param,n-workers)

    OGFSessLog analysis.completeness auto $args -title {Completeness simulation (seed 42)} 
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "Completeness error: $err"
	return
    }

    # Load completeness results as a new catalog view
    CatalogPanelLoadTSV $data "completeness"
    set catpanel(status) "Completeness simulation complete"
}

