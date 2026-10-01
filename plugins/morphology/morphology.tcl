# morphology/morphology.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/morphology/plugin.json.

proc CatalogPanelMorphometry {} {
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

    set script [CatalogPanelGetScript ds9_morphometry.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_morphometry.py"
	return
    }

    # Save temp catalog
    set tmpcat [CatalogPanelSaveTempCatalog "morphometry"]
    if {$tmpcat eq {}} {
	set catpanel(status) "Failed to save temp catalog"
	return
    }

    set catpanel(status) "Measuring morphometry (CAS/Gini/M20)..."
    update idletasks

    OGFSessLog analysis.morphometry auto [list [OGFPython] $script $fn --catalog $tmpcat --n-workers $catpanel(param,n-workers)] -title {Non-parametric morphology (CAS/Gini/M20)} -post [dict create kind add cols_list {CONC ASYM GINI M20 R_PETRO}]
    if {[catch {
	set data [exec [OGFPython] $script $fn --catalog $tmpcat \
	    --n-workers $catpanel(param,n-workers) 2>@stderr]
    } err]} {
	set catpanel(status) "Morphometry error: $err"
	return
    }

    if {[string trim $data] eq {}} {
	set catpanel(status) "Morphometry: no output"
	return
    }

    # Add columns to alldata
    CatalogPanelAddColumnsFromTSV $data {CONC ASYM GINI M20 R_PETRO}
    set catpanel(status) "Morphometry complete (CAS/Gini/M20/Petrosian)"
}

proc CatalogPanelSersicFit {} {
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

    set script [CatalogPanelGetScript ds9_sersic.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_sersic.py"
	return
    }

    set tmpcat [CatalogPanelSaveTempCatalog "sersic"]
    if {$tmpcat eq {}} {
	set catpanel(status) "Failed to save temp catalog"
	return
    }

    set catpanel(status) "Sérsic profile fitting..."
    update idletasks

    set args [list [OGFPython] $script $fn --catalog $tmpcat]
    if {[info exists catpanel(param,mag-zeropoint)]} {
	lappend args --mag-zeropoint $catpanel(param,mag-zeropoint)
    }
    lappend args --n-workers $catpanel(param,n-workers)

    OGFSessLog analysis.sersic auto $args -title {Sersic fitting} -post [dict create kind add cols_list {SERSIC_N SERSIC_RE SERSIC_IE SERSIC_ELLIP SERSIC_THETA SERSIC_CHI2}]
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "Sérsic fit error: $err"
	return
    }

    CatalogPanelAddColumnsFromTSV $data \
	{SERSIC_N SERSIC_RE SERSIC_IE SERSIC_ELLIP SERSIC_THETA SERSIC_CHI2}
    set catpanel(status) "Sérsic fitting complete"
}

