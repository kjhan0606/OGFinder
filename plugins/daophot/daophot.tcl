# DAOPHOT plugin (plugins/daophot): Tcl side = apply the outputs of the CLI steps (frames, text/plot windows, PSF star list editing).
# State is kept in catalog keys daophot,* through ::ogf::cat (registered in ogf_core.tcl); no direct catpanel access.

proc OGFDaophotFile {name} {
    return [file join [OGFSessWorkDir] daophot daophot_$name]
}

proc OGFDaophotShowText {title name} {
    set f [OGFDaophotFile $name]
    if {![file exists $f]} {::ogf::status "DAOPHOT: $name does not exist yet - run the step first"; return}
    set fd [open $f r]; set txt [read $fd]; close $fd
    OGFTextWindow "$title ($f)" $txt
}

proc OGFDaophotLoadFrame {fn} {
    CreateFrame
    if {[catch {LoadFitsFile $fn {} {}} err]} {::ogf::log ERROR "daophot: cannot load $fn: $err"; return 0}
    global scale
    set scale(mode) zscale
    ChangeScaleMode
    return 1
}

# open the given daophot_*.fits files in new frames, the original frame stays current
proc OGFDaophotOpenFrames {names} {
    global current
    set orig $current(frame)
    set n 0
    foreach nm $names {
	set f [OGFDaophotFile $nm]
	if {[file exists $f] && [OGFDaophotLoadFrame $f]} {incr n}
    }
    catch {GotoFrame $orig}
    return $n
}

proc OGFDaophotAfter {} {
    ::ogf::cat::set daophot,stars_file [OGFDaophotFile stars.tsv]
    ::ogf::cat::set daophot,diag_file [OGFDaophotFile diag.png]
    ::ogf::cat::set daophot,psf_file [OGFDaophotFile psf.json]
    ::ogf::cat::set daophot,resid_file [OGFDaophotFile resid.fits]
    if {[::ogf::params::get daophot show-frames]} {catch {OGFDaophotOpenFrames {resid.fits}}}
}
proc OGFDaophotAfterPsf {} {::ogf::cat::set daophot,psf_file [OGFDaophotFile psf.json]}
proc OGFDaophotAfterSub {} {if {[::ogf::params::get daophot show-frames]} {catch {OGFDaophotOpenFrames {sub.fits}}}}
proc OGFDaophotAfterAdd {} {if {[::ogf::params::get daophot show-frames]} {catch {OGFDaophotOpenFrames {add.fits}}}}
proc OGFDaophotAfterMatch {} {::ogf::cat::set daophot,cmd_file [OGFDaophotFile cmd.png]}

proc OGFDaophotShowFrames {} {
    set n [OGFDaophotOpenFrames {resid.fits sub.fits add.fits}]
    ::ogf::status "DAOPHOT: $n frame(s) opened (residual / subtracted / artificial-star image where they exist)"
}

proc OGFDaophotTable {} {OGFDaophotShowText "DAOPHOT star table" stars.tsv}

proc OGFDaophotImageWindow {w title png} {
    if {![file exists $png]} {::ogf::status "DAOPHOT: no plot yet - run the step first"; return}
    if {[winfo exists $w]} {destroy $w}
    toplevel $w
    wm title $w $title
    catch {image delete ${w}_img}
    image create photo ${w}_img -file $png
    label $w.l -image ${w}_img
    pack $w.l -fill both -expand 1
    ttk::button $w.close -text Close -command [list destroy $w]
    pack $w.close -pady 3
}
proc OGFDaophotDiag {} {OGFDaophotImageWindow .ogfdaodiag "DAOPHOT diagnostics" [OGFDaophotFile diag.png]}
proc OGFDaophotCmd {} {OGFDaophotImageWindow .ogfdaocmd "DAOPHOT matched CMD" [OGFDaophotFile cmd.png]}

# --- PSF star list editing from the catalog selection: the parameters psf-add / psf-remove ("x,y;x,y", 1-based) are ordinary
# recorded parameters, so a replayed session reproduces the manual choices.
proc OGFDaophotSelectedXY {} {
    set sel [::ogf::cat::selection]
    if {$sel eq {}} {return {}}
    set out {}
    foreach r [::ogf::cat::rows] {
	if {[dict get $r NUMBER] in $sel && [dict exists $r X_IMAGE] && [dict exists $r Y_IMAGE]} {
	    lappend out "[dict get $r X_IMAGE],[dict get $r Y_IMAGE]"
	}
    }
    return $out
}
proc OGFDaophotEditList {name dropname} {
    set xy [OGFDaophotSelectedXY]
    if {$xy eq {}} {::ogf::status "DAOPHOT: select catalog rows first"; return}
    set cur [string trim [::ogf::params::get daophot $name]]
    set lst [expr {$cur eq {} ? {} : [split $cur \;]}]
    set other [string trim [::ogf::params::get daophot $dropname]]
    set olst [expr {$other eq {} ? {} : [split $other \;]}]
    foreach p $xy {
	if {$p ni $lst} {lappend lst $p}
	set i [lsearch -exact $olst $p]
	if {$i >= 0} {set olst [lreplace $olst $i $i]}
    }
    ::ogf::params::put daophot $name [join $lst \;]
    ::ogf::params::put daophot $dropname [join $olst \;]
    catch {::ogf::params::save daophot}
    ::ogf::status "DAOPHOT: $name = [llength $lst] position(s), $dropname = [llength $olst] (run PICKPSF / PSF again)"
}
proc OGFDaophotAddPsf {} {OGFDaophotEditList psf-add psf-remove}
proc OGFDaophotRemPsf {} {OGFDaophotEditList psf-remove psf-add}
