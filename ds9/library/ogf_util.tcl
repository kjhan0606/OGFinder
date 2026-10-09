#  OGFinder helpers shared by the multi-band, link and mask modules.
#  (new file; see docs/multiband_mask_link.md)

package provide DS9 1.0

# Python interpreter for the OGFinder helpers that need astropy/sep.
# Override with the environment variable OGFINDER_PYTHON.
proc OGFPython {} {
    if {[info exists ::env(OGFINDER_PYTHON)] && $::env(OGFINDER_PYTHON) ne {}} {
	return $::env(OGFINDER_PYTHON)
    }
    return python3
}

proc OGFSextractBin {} {
    return [file join [OGFSessRoot] ogfmeas sextract.py]
}

# Extend LD_LIBRARY_PATH the same way CatalogPanelExtract does (conda libs)
proc OGFPrepareLibPath {} {
    if {$::tcl_platform(os) eq "Windows NT"} return
    set var LD_LIBRARY_PATH
    if {$::tcl_platform(os) eq "Darwin"} {set var DYLD_LIBRARY_PATH}
    set libpaths {}
    if {[info exists ::env(CONDA_PREFIX)]} {lappend libpaths "$::env(CONDA_PREFIX)/lib"}
    set home_conda [file join [file normalize ~] miniconda3/lib]
    if {[file isdirectory $home_conda]} {lappend libpaths $home_conda}
    if {[info exists ::env($var)]} {lappend libpaths $::env($var)}
    if {[llength $libpaths] > 0} {set ::env($var) [join $libpaths :]}
}

# Modal form: fields = list of {key label default ?check?}; returns dict or {} if
# cancelled.  Optional buttons = list of {label script}; the script runs in
# the caller's namespace with the array ::ogfform_v available (keyed by key).
proc OGFForm {title fields {buttons {}} {note {}}} {
    global ogfform_v ogfform_done
    set w .ogfform
    catch {destroy $w}
    toplevel $w
    wm title $w $title
    wm transient $w .
    set ogfform_done {}
    array unset ogfform_v
    set r 0
    if {$note ne {}} {
	ttk::label $w.note -text $note -justify left -wraplength 380
	grid $w.note -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady 4
	incr r
    }
    foreach f $fields {
	lassign $f key label def kind
	set ogfform_v($key) $def
	ttk::label $w.l$r -text $label
	if {$kind eq "check"} {
	    ttk::checkbutton $w.e$r -variable ogfform_v($key)
	} else {
	    ttk::entry $w.e$r -textvariable ogfform_v($key) -width 16
	}
	grid $w.l$r -row $r -column 0 -sticky w -padx 8 -pady 2
	grid $w.e$r -row $r -column 1 -sticky we -padx 8 -pady 2
	incr r
    }
    if {[llength $buttons]} {
	set bf [ttk::frame $w.pre]
	foreach b $buttons {
	    lassign $b bl bs
	    ttk::button $bf.b[incr bi] -text $bl -command $bs
	    pack $bf.b$bi -side left -padx 2
	}
	grid $bf -row $r -column 0 -columnspan 2 -pady 4
	incr r
    }
    set bf2 [ttk::frame $w.bb]
    ttk::button $bf2.ok -text OK -command {set ogfform_done ok}
    ttk::button $bf2.cancel -text Cancel -command {set ogfform_done cancel}
    pack $bf2.ok $bf2.cancel -side left -padx 4
    grid $bf2 -row $r -column 0 -columnspan 2 -pady 6
    bind $w <Return> {set ogfform_done ok}
    bind $w <Escape> {set ogfform_done cancel}
    wm protocol $w WM_DELETE_WINDOW {set ogfform_done cancel}
    catch {grab $w}
    vwait ogfform_done
    set res {}
    if {$ogfform_done eq "ok"} {set res [array get ogfform_v]}
    catch {grab release $w}
    destroy $w
    return $res
}

# Read-only text window
proc OGFTextWindow {title text} {
    set w .ogftext
    catch {destroy $w}
    toplevel $w
    wm title $w $title
    text $w.t -width 90 -height 14 -font TkFixedFont -wrap none
    $w.t insert end $text
    $w.t configure -state disabled
    ttk::button $w.close -text Close -command [list destroy $w]
    pack $w.t -fill both -expand true
    pack $w.close -pady 4
}
