# morphology/bulge_disk.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/morphology/plugin.json.

proc CatalogPanelBDParamLoad {} {
    global catpanel
    set prf [file join [file normalize ~] .ds9 bulge_disk.prf]
    if {![file exists $prf]} return
    if {[catch {
	set fd [open $prf r]
	set data [read $fd]
	close $fd
    }]} return
    foreach line [split $data "\n"] {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set eq [string first "=" $line]
	if {$eq < 0} continue
	set key [string trim [string range $line 0 [expr {$eq-1}]]]
	set val [string trim [string range $line [expr {$eq+1}] end]]
	set catpanel(bd,param,$key) $val
    }
}

proc CatalogPanelBDParamSave {} {
    global catpanel
    set prf [file join [file normalize ~] .ds9 bulge_disk.prf]
    catch {file mkdir [file dirname $prf]}
    if {[catch {set fd [open $prf w]}]} return
    foreach key {max-sources free-bulge-n mag-zeropoint pixel-scale} {
	if {[info exists catpanel(bd,param,$key)]} {
	    puts $fd "$key=$catpanel(bd,param,$key)"
	}
    }
    close $fd
}

proc CatalogPanelBulgeDisk {} {
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

    set script [CatalogPanelGetScript ds9_bulge_disk.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_bulge_disk.py"
	return
    }

    set tmpcat [CatalogPanelSaveTempCatalog "bulgedisk"]
    if {$tmpcat eq {}} {
	set catpanel(status) "Failed to save temp catalog"
	return
    }

    set catpanel(status) "Running Bulge+Disk decomposition..."
    update idletasks

    set args [list [OGFPython] $script $fn --catalog $tmpcat]
    lappend args --mag-zeropoint $catpanel(bd,param,mag-zeropoint)
    lappend args --pixel-scale $catpanel(bd,param,pixel-scale)
    lappend args --max-sources $catpanel(bd,param,max-sources)

    if {$catpanel(bd,param,free-bulge-n)} {
	lappend args --free-bulge-n
    }

    # Use PSF if available
    if {[info exists catpanel(psf,file)] && [file exists $catpanel(psf,file)]} {
	lappend args --psf $catpanel(psf,file)
    }

    lappend args --n-workers $catpanel(param,n-workers)

    OGFSessLog analysis.bulge_disk auto $args -title {Bulge+Disk decomposition} -post [dict create kind add cols_list {BT_RATIO BULGE_RE BULGE_MAG DISK_RS DISK_MAG BD_CHI2 BD_FLAG}]
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "Bulge+Disk error: $err"
	return
    }

    if {[string trim $data] eq {}} {
	set catpanel(status) "Bulge+Disk: no output"
	return
    }

    CatalogPanelAddColumnsFromTSV $data \
	{BT_RATIO BULGE_RE BULGE_MAG DISK_RS DISK_MAG BD_CHI2 BD_FLAG}
    set catpanel(status) "Bulge+Disk decomposition complete"
}

