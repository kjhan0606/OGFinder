# ds10core plugin (stand-alone shell of the shared Astrafex Core): Tcl side of "Open offline bundle".
# The bundle is restored by the command line (ds10.py open ...); this file only finds the restored workspace and shows its image and regions in ds9.
# Nothing of the compute core is used here.

# folder the bundle is restored into: the parameter "open-dest", else ~/astrafex-offline/<bundle file name without .zip> (the same rule as plugins/ds10core/ds10.py)
proc OGFDs10coreOpenDest {} {
    set d [string trim [::ogf::params::get ds10core open-dest]]
    if {$d ne {}} {return [file normalize $d]}
    set b [string trim [::ogf::params::get ds10core bundle]]
    return [file normalize [file join ~ astrafex-offline [file rootname [file tail $b]]]]
}

proc OGFDs10coreOpenAfter {} {
    global current
    set ws [OGFDs10coreOpenDest]
    set wj [file join $ws workspace.json]
    if {![file exists $wj]} {::ogf::status "Astrafex Core: the bundle could not be opened (see the text window)"; return}
    set fd [open $wj r]; set txt [read $fd]; close $fd
    ::ogf::cat::set ds10core,workspace $ws
    set shown {}
    # first file of the bundle that is an image (the list is in the order of the web session)
    if {[regexp {"local_path":\s*"([^"]+)"} $txt -> rel]} {
	set img [file join $ws $rel]
	if {[file exists $img] && ![catch {LoadFitsFile $img {} {}} err]} {
	    set shown $rel
	    set reg [file join $ws regions regions.reg]
	    if {[file exists $reg]} {catch {MarkerLoadFile $reg $current(frame) ds9 image fk5}}
	} elseif {[info exists err]} {
	    ::ogf::log ERROR "ds10core: cannot load $img: $err"
	}
    }
    ::ogf::status "Offline bundle restored to $ws[expr {$shown ne {} ? "; showing $shown" : {}}]. Run steps with the script \"$ws/session_script.json\" (Run astrafex-script)."
}

proc OGFDs10coreCatalogAfter {} {
    global current
    # the step writes markers.reg into its work directory ({work} = OGFSessWorkDir)
    set reg [file join [OGFSessWorkDir] markers.reg]
    if {[file exists $reg] && [info exists current(frame)]} {
	catch {MarkerLoadFile $reg $current(frame) ds9 image fk5}
	::ogf::status "Catalogue markers loaded from $reg"
    }
}

# ---- star aperture photometry (steps aper-phot / aper-series)
# before: the regions on the image (the stars you clicked: point or circle regions) -> {work}/aper_picks.reg, in image coordinates
proc OGFDs10coreAperBefore {} {
    global current
    set f [file join [OGFSessWorkDir] aper_picks.reg]
    set txt {}
    if {[info exists current(frame)] && $current(frame) ne {}} {catch {set txt [$current(frame) marker list ds9 image fk5 degrees 0]}}
    file mkdir [file dirname $f]
    set fd [open $f w]; puts -nonewline $fd $txt; close $fd
    set n [regexp -all -line {^\s*(point|circle)\(} $txt]
    if {[::ogf::params::get ds10core aper-stars] eq "regions"} {
	::ogf::status [expr {$n ? "Aperture photometry: $n picked star(s)" : "Aperture photometry: no stars picked (click on stars first); detecting stars instead"}]
    }
}

# after: apertures and sky annuli replace the picks on the image; the light curve (time series) is shown in its own window
proc OGFDs10coreAperAfter {} {
    global current
    set wd [OGFSessWorkDir]
    set reg [file join $wd aper_apertures.reg]
    if {[file exists $reg] && [info exists current(frame)]} {
	catch {$current(frame) marker delete all}
	catch {MarkerLoadFile $reg $current(frame) ds9 image fk5}
    }
    set sj [file join $wd aperphot_summary.json]
    set png [file join $wd aperseries_lightcurve.png]
    set series [expr {[file exists $png] && [file mtime $png] >= [clock seconds] - 600}]
    if {$series} {set sj [file join $wd aperseries_summary.json]}
    set msg "Aperture photometry done"
    if {[file exists $sj]} {
	set fd [open $sj r]; set s [read $fd]; close $fd
	set n ?; set fw ?; set zp ?
	regexp {"n_stars":\s*([0-9]+)} $s -> n
	regexp {"fwhm_px":\s*([0-9.]+)} $s -> fw
	regexp {"zeropoint":\s*\{[^\}]*"zp":\s*([-0-9.]+)} $s -> zp
	set msg "Aperture photometry: $n stars, FWHM [format %.2f $fw] px, ZP [format %.3f $zp]; table = catalogue, apertures/annuli on the image ($reg)"
    }
    ::ogf::status $msg
    if {$series} {OGFDs10coreAperLightCurve $png}
}

proc OGFDs10coreAperLightCurve {png} {
    set w .ogfaperlc
    catch {destroy $w}
    toplevel $w
    wm title $w "Astrafex Core: light curve"
    catch {wm geometry $w +[expr {[winfo rootx .] + [winfo width .] / 2}]+[expr {[winfo rooty .] + 260}]}
    if {[catch {image create photo ogfaperlcimg -file $png} err]} {
	label $w.l -text "light curve: $png ($err)"
    } else {
	label $w.l -image ogfaperlcimg
    }
    pack $w.l -fill both -expand 1
    button $w.b -text Close -command [list destroy $w]
    pack $w.b -side bottom
}
