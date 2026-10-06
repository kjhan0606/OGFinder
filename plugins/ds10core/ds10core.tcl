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
