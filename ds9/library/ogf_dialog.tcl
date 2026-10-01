#  OGFinder core: ONE declarative parameter dialog.
#
#  Every settings dialog is a list of parameter specs (plugin.json "params"): name, label, type
#  (int|float|string|bool|choice|file), default, min, max, choices, help, group, expert, step.
#  OGFParamDialog renders it: groups, Basic/Expert toggle, per-field help line, range validation,
#  Reset to defaults, preset save/load, Apply/OK/Cancel.  The values live in the plugin's parameter store
#  (::ogfparam, or a legacy array declared in "store" so old code keeps reading catpanel(...)).
#  See docs/plugins.md

package provide DS9 1.0

# OGFParamDialog PLUGIN ?-step S? ?-title T? ?-group G? ?-parent W?
proc OGFParamDialog {plugin args} {
    global ogfdlg
    array set o {-step {} -title {} -group {} -parent .}
    array set o $args
    set specs [::ogf::params::specs $plugin $o(-step)]
    if {[llength $specs] == 0} {
	::ogf::status "no parameters declared for $plugin"
	return {}
    }
    set w .ogfparam_[string map {. _ - _} $plugin]
    if {[winfo exists $w]} {raise $w; return $w}
    toplevel $w
    set m [::ogf::reg::get $plugin]
    set title [expr {$o(-title) ne {} ? $o(-title) : "[::ogf::json::get $m name $plugin] - Settings"}]
    wm title $w $title
    catch {wm transient $w $o(-parent)}
    set ogfdlg($w,plugin) $plugin
    set ogfdlg($w,specs) $specs
    set ogfdlg($w,expert) 0
    set ogfdlg($w,help) {}
    set ogfdlg($w,err) {}
    set ogfdlg($w,preset) {}
    set ogfdlg($w,widgets) {}
    foreach s $specs {
	set ogfdlg($w,v,[dict get $s name]) [::ogf::params::get $plugin [dict get $s name]]
    }
    # top: presets + expert toggle
    ttk::frame $w.top
    ttk::label $w.top.pl -text "Preset:"
    ttk::combobox $w.top.pc -textvariable ogfdlg($w,preset) -width 16 -values [::ogf::params::preset_list $plugin]
    ttk::button $w.top.pload -text Load -width 5 -command [list OGFParamPresetLoad $w]
    ttk::button $w.top.psave -text Save -width 5 -command [list OGFParamPresetSave $w]
    ttk::checkbutton $w.top.ex -text "Expert" -variable ogfdlg($w,expert) -command [list OGFParamLayout $w]
    pack $w.top.pl $w.top.pc $w.top.pload $w.top.psave -side left -padx 2
    pack $w.top.ex -side right -padx 6
    pack $w.top -side top -fill x -pady 4
    # body: groups in order of first appearance
    ttk::frame $w.body
    pack $w.body -side top -fill both -expand 1 -padx 6
    set groups {}
    foreach s $specs {
	set g [::ogf::json::get $s group General]
	if {$g ni $groups} {lappend groups $g}
    }
    set ogfdlg($w,groups) $groups
    set gi 0
    foreach g $groups {
	ttk::labelframe $w.body.g$gi -text $g
	set ogfdlg($w,gf,$g) $w.body.g$gi
	set r 0
	foreach s $specs {
	    if {[::ogf::json::get $s group General] ne $g} continue
	    set n [dict get $s name]
	    set lab $w.body.g$gi.l$r
	    set e $w.body.g$gi.e$r
	    ttk::label $lab -text "[::ogf::json::get $s label $n][expr {[::ogf::json::get $s unit] ne {} ? " ([::ogf::json::get $s unit])" : {}}]:" -anchor w
	    set var ogfdlg($w,v,$n)
	    switch -- [::ogf::json::get $s type string] {
		bool {ttk::checkbutton $e -variable $var}
		choice {ttk::combobox $e -textvariable $var -values [::ogf::json::get $s choices] -state readonly -width 14}
		int - float {
		    set mn [::ogf::json::get $s min]; set mx [::ogf::json::get $s max]
		    if {$mn ne {} && $mx ne {}} {
			ttk::spinbox $e -textvariable $var -from $mn -to $mx -increment [expr {[::ogf::json::get $s type] eq "int" ? 1 : ($mx-$mn)/100.0}] -width 12
		    } else {
			ttk::entry $e -textvariable $var -width 14
		    }
		}
		file {
		    ttk::frame $e
		    ttk::entry $e.e -textvariable $var -width 22
		    ttk::button $e.b -text ... -width 3 -command [list OGFParamBrowse $w $n]
		    pack $e.e -side left -fill x -expand 1
		    pack $e.b -side left
		}
		default {ttk::entry $e -textvariable $var -width 18}
	    }
	    set hp [::ogf::json::get $s help]
	    foreach x [list $lab $e] {
		bind $x <Enter> [list set ogfdlg($w,help) $hp]
		bind $x <FocusIn> [list set ogfdlg($w,help) $hp]
	    }
	    if {[winfo exists $e.e]} {bind $e.e <FocusIn> [list set ogfdlg($w,help) $hp]}
	    lappend ogfdlg($w,widgets) [list $n $lab $e $g [::ogf::json::get $s expert 0]]
	    incr r
	}
	incr gi
    }
    ttk::label $w.help -textvariable ogfdlg($w,help) -wraplength 420 -justify left -anchor w -foreground gray30
    ttk::label $w.err -textvariable ogfdlg($w,err) -foreground #b00020 -wraplength 420 -anchor w
    ttk::frame $w.bb
    ttk::button $w.bb.reset -text "Reset to defaults" -command [list OGFParamReset $w]
    ttk::button $w.bb.apply -text Apply -command [list OGFParamApply $w 0]
    ttk::button $w.bb.ok -text OK -command [list OGFParamApply $w 1]
    ttk::button $w.bb.cancel -text Cancel -command [list OGFParamClose $w]
    pack $w.bb.cancel $w.bb.ok $w.bb.apply -side right -padx 3
    pack $w.bb.reset -side left -padx 3
    pack $w.bb -side bottom -fill x -pady 6 -padx 6
    pack $w.err -side bottom -fill x -padx 8
    pack $w.help -side bottom -fill x -padx 8 -pady 2
    bind $w <Escape> [list OGFParamClose $w]
    bind $w <Return> [list OGFParamApply $w 1]
    wm protocol $w WM_DELETE_WINDOW [list OGFParamClose $w]
    OGFParamLayout $w
    return $w
}

# (re)grid the groups / fields according to the Basic/Expert toggle; groups with no visible field are hidden
proc OGFParamLayout {w} {
    global ogfdlg
    foreach g $ogfdlg($w,groups) {
	set gf $ogfdlg($w,gf,$g)
	pack forget $gf
	foreach c [winfo children $gf] {grid forget $c}
    }
    set vis {}
    foreach g $ogfdlg($w,groups) {
	set gf $ogfdlg($w,gf,$g)
	set r 0
	foreach wd $ogfdlg($w,widgets) {
	    lassign $wd n lab e wg expert
	    if {$wg ne $g} continue
	    if {$expert && !$ogfdlg($w,expert)} continue
	    grid $lab -in $gf -row $r -column 0 -sticky w -padx 6 -pady 2
	    grid $e -in $gf -row $r -column 1 -sticky we -padx 6 -pady 2
	    incr r
	}
	grid columnconfigure $gf 1 -weight 1
	if {$r > 0} {pack $gf -in $w.body -side top -fill x -pady 3; lappend vis $g}
    }
    set ogfdlg($w,visible_groups) $vis
}

# number of fields currently shown (for tests)
proc OGFParamVisibleFields {w} {
    global ogfdlg
    set n 0
    foreach wd $ogfdlg($w,widgets) {if {[winfo manager [lindex $wd 1]] eq "grid"} {incr n}}
    return $n
}

proc OGFParamBrowse {w name} {
    global ogfdlg
    set f [tk_getOpenFile -parent $w]
    if {$f ne {}} {set ogfdlg($w,v,$name) $f}
}

proc OGFParamCollect {w} {
    global ogfdlg
    set d [dict create]
    foreach s $ogfdlg($w,specs) {dict set d [dict get $s name] $ogfdlg($w,v,[dict get $s name])}
    return $d
}

proc OGFParamReset {w} {
    global ogfdlg
    set ogfdlg($w,err) {}
    foreach s $ogfdlg($w,specs) {
	set ogfdlg($w,v,[dict get $s name]) [::ogf::json::get $s default]
    }
}

# validate everything, write to the store, save, call the plugin's on_apply hook
proc OGFParamApply {w close} {
    global ogfdlg
    set plugin $ogfdlg($w,plugin)
    set errs {}
    foreach s $ogfdlg($w,specs) {
	set e [::ogf::params::validate $s $ogfdlg($w,v,[dict get $s name])]
	if {$e ne {}} {lappend errs $e}
    }
    if {[llength $errs]} {set ogfdlg($w,err) [join $errs "\n"]; return 0}
    set ogfdlg($w,err) {}
    foreach s $ogfdlg($w,specs) {
	::ogf::params::put $plugin [dict get $s name] $ogfdlg($w,v,[dict get $s name])
    }
    ::ogf::params::save $plugin
    set h [::ogf::json::get [::ogf::reg::get $plugin] on_apply]
    if {$h ne {}} {catch {uplevel #0 $h}}
    ::ogf::status "[::ogf::json::get [::ogf::reg::get $plugin] name $plugin]: settings applied and saved"
    if {$close} {OGFParamClose $w}
    return 1
}

proc OGFParamClose {w} {
    global ogfdlg
    foreach k [array names ogfdlg $w,*] {unset ogfdlg($k)}
    destroy $w
}

proc OGFParamPresetSave {w} {
    global ogfdlg
    set n [string trim $ogfdlg($w,preset)]
    if {$n eq {} || ![regexp {^[A-Za-z0-9_.-]+$} $n]} {set ogfdlg($w,err) "preset name: letters, digits, _ . - only"; return}
    ::ogf::params::preset_save $ogfdlg($w,plugin) $n [OGFParamCollect $w]
    $w.top.pc configure -values [::ogf::params::preset_list $ogfdlg($w,plugin)]
    set ogfdlg($w,err) {}
    set ogfdlg($w,help) "preset '$n' saved"
}

proc OGFParamPresetLoad {w} {
    global ogfdlg
    set n [string trim $ogfdlg($w,preset)]
    if {[catch {set d [::ogf::params::preset_load $ogfdlg($w,plugin) $n]} err]} {
	set ogfdlg($w,err) "cannot load preset '$n'"
	return
    }
    dict for {k v} $d {if {[info exists ogfdlg($w,v,$k)]} {set ogfdlg($w,v,$k) $v}}
    set ogfdlg($w,err) {}
    set ogfdlg($w,help) "preset '$n' loaded (press Apply to use it)"
}
