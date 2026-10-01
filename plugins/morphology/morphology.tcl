# morphology/morphology.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/morphology/plugin.json.

proc CatalogPanelMorphometry {} {

    if {![::ogf::cat::has]} {
	::ogf::cat::set status "Extract sources first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	::ogf::cat::set status "No FITS image loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_morphometry.py]
    if {![file exists $script]} {
	::ogf::cat::set status "Script not found: ds9_morphometry.py"
	return
    }

    # Save temp catalog
    set tmpcat [CatalogPanelSaveTempCatalog "morphometry"]
    if {$tmpcat eq {}} {
	::ogf::cat::set status "Failed to save temp catalog"
	return
    }

    ::ogf::cat::set status "Measuring morphometry (CAS/Gini/M20)..."
    update idletasks

    OGFSessLog analysis.morphometry auto [list [OGFPython] $script $fn --catalog $tmpcat --n-workers [::ogf::cat::get param,n-workers]] -title {Non-parametric morphology (CAS/Gini/M20)} -post [dict create kind add cols_list {CONC ASYM GINI M20 R_PETRO}]
    if {[catch {
	set data [exec [OGFPython] $script $fn --catalog $tmpcat \
	    --n-workers [::ogf::cat::get param,n-workers] 2>@stderr]
    } err]} {
	::ogf::cat::set status "Morphometry error: $err"
	return
    }

    if {[string trim $data] eq {}} {
	::ogf::cat::set status "Morphometry: no output"
	return
    }

    # Add columns to alldata
    CatalogPanelAddColumnsFromTSV $data {CONC ASYM GINI M20 R_PETRO}
    ::ogf::cat::set status "Morphometry complete (CAS/Gini/M20/Petrosian)"
}

proc CatalogPanelSersicFit {} {

    if {![::ogf::cat::has]} {
	::ogf::cat::set status "Extract sources first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	::ogf::cat::set status "No FITS image loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_sersic.py]
    if {![file exists $script]} {
	::ogf::cat::set status "Script not found: ds9_sersic.py"
	return
    }

    set tmpcat [CatalogPanelSaveTempCatalog "sersic"]
    if {$tmpcat eq {}} {
	::ogf::cat::set status "Failed to save temp catalog"
	return
    }

    ::ogf::cat::set status "Sérsic profile fitting..."
    update idletasks

    set args [list [OGFPython] $script $fn --catalog $tmpcat]
    if {[::ogf::cat::exists param,mag-zeropoint]} {
	lappend args --mag-zeropoint [::ogf::cat::get param,mag-zeropoint]
    }
    lappend args --n-workers [::ogf::cat::get param,n-workers]

    OGFSessLog analysis.sersic auto $args -title {Sersic fitting} -post [dict create kind add cols_list {SERSIC_N SERSIC_RE SERSIC_IE SERSIC_ELLIP SERSIC_THETA SERSIC_CHI2}]
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	::ogf::cat::set status "Sérsic fit error: $err"
	return
    }

    CatalogPanelAddColumnsFromTSV $data \
	{SERSIC_N SERSIC_RE SERSIC_IE SERSIC_ELLIP SERSIC_THETA SERSIC_CHI2}
    ::ogf::cat::set status "Sérsic fitting complete"
}

