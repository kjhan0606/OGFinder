#  OGFinder: workflow GUI built from the plugin registry.
#
#  Right-hand catalog panel (top to bottom):
#     row 0   [Workflow] [Tools] [Tile]  ......  progress strip  Detect > Classify > Measure   [Stop]
#     tab bar Detect | Classify | Measure | Low-SB | Time-domain | Results
#     row     one chip per plugin of the selected tab:  [Name v] [>] [gear]
#                 Name v : menu with all steps of the plugin (everything the old top-level menus offered)
#                 >      : run the plugin's primary step         gear : settings dialog
#     filter row, status line, selection summary, then the catalog table (geometry unchanged).
#  Everything is generated from plugins/*/plugin.json (docs/plugins.md); nothing here names a feature.

package provide DS9 1.0

namespace eval ::ogf::ui {
    variable built 0
}

proc OGFUIInit {} {
    global ogfui
    array set ogfui {tab Detect tile 0 tilemode grid lockwcs 1 sharescale 0 busy 0 pending {} nb {} content {} frames {}}
}

# ---------------------------------------------------------------- menus
proc OGFUIStepCommand {id s} {
    if {[dict exists $s proc]} {return [dict get $s proc]}
    return [list ::ogf::step::run $id [dict get $s id]]
}

# dynamic band cascades keep their post commands
proc OGFUIAddDynamic {m s} {
    set p [dict get $s proc]
    switch -- [lindex $p 1] {
	detect {
	    menu $m.det -tearoff 0 -postcommand [list OGFBandsPostDetect $m.det]
	    $m add cascade -label [dict get $s label] -menu $m.det
	}
	remove {
	    menu $m.rm -tearoff 0 -postcommand [list OGFBandsPostRemove $m.rm]
	    $m add cascade -label [dict get $s label] -menu $m.rm
	}
    }
}

# fill menu M with the steps of plugin ID (in manifest order); returns number of entries
proc OGFUIFillPluginMenu {m id} {
    set pm [::ogf::reg::get $id]
    set n 0
    set has_settings_step 0
    set stp [::ogf::json::get $pm settings]
    foreach s [::ogf::json::get $pm steps] {
	if {[::ogf::json::get $s settings_only 0]} {set has_settings_step 1}
	if {$stp ne {} && [::ogf::json::get $s proc] eq $stp} {set has_settings_step 1}
    }
    foreach s [::ogf::json::get $pm steps] {
	set lab [dict get $s label]
	if {[string match OGFUIBandMenu* [::ogf::json::get $s proc]]} {
	    OGFUIAddDynamic $m $s
	} elseif {[dict exists $s variants]} {
	    set sub $m.v[dict get $s id]
	    menu $sub -tearoff 0
	    foreach v [dict get $s variants] {
		$sub add command -label [dict get $v label] -command [dict get $v proc]
	    }
	    $m add cascade -label $lab -menu $sub
	} elseif {[dict exists $s toggle]} {
	    set var [dict get $s toggle]
	    if {[::ogf::json::get $s toggle_only 0]} {
		$m add checkbutton -label $lab -variable $var
	    } else {
		$m add checkbutton -label $lab -variable $var -command [dict get $s proc]
	    }
	} else {
	    $m add command -label $lab -command [OGFUIStepCommand $id $s]
	}
	incr n
    }
    if {[::ogf::json::get $pm settings] ne {} && !$has_settings_step} {
	if {$n > 0} {$m add separator}
	$m add command -label "Settings..." -command [list OGFUISettings $id]
	incr n
    }
    return $n
}

proc OGFUISettings {id} {
    set pm [::ogf::reg::get $id]
    set st [::ogf::json::get $pm settings]
    if {$st eq "params"} {
	OGFParamDialog $id
    } elseif {$st ne {}} {
	uplevel #0 $st
    }
}

# the menu widget of plugin ID's chip (for tests / scripts)
proc OGFUIPluginMenu {id} {
    global ogfui
    if {[info exists ogfui(menu,$id)]} {return $ogfui(menu,$id)}
    return {}
}

# ---------------------------------------------------------------- build
proc OGFUIBuild {f} {
    global ogfui catpanel ds9
    OGFUIInit
    set ogfui(root) $f
    # --- row 0: menubar
    set mb [ttk::frame $f.menubar]
    set catpanel(menubar) $mb
    ttk::style layout CatMenu.TMenubutton {
	Menubutton.focus -sticky nswe -children {
	    Menubutton.padding -sticky we -children {
		Menubutton.label -side left -sticky {}
	    }
	}
    }
    ttk::style configure CatMenu.TMenubutton -relief flat -padding {3 3}
    ttk::style map CatMenu.TMenubutton -relief {pressed flat active flat}
    ttk::menubutton $mb.workflow -text Workflow -menu $mb.workflow.m -style CatMenu.TMenubutton
    menu $mb.workflow.m -tearoff 0
    ttk::menubutton $mb.tools -text Tools -menu $mb.tools.m -style CatMenu.TMenubutton
    menu $mb.tools.m -tearoff 0
    ttk::style configure CatTile.Toolbutton -padding {4 1}
    ttk::style configure CatApply.TButton -padding {2 0}
    ttk::style layout CatChip.TMenubutton {
	Menubutton.focus -sticky nswe -children {
	    Menubutton.padding -sticky we -children {
		Menubutton.label -side left -sticky {}
	    }
	}
    }
    ttk::style configure CatChip.TMenubutton -relief raised -padding {4 1}
    ttk::style configure CatIcon.TButton -padding {1 0}
    ttk::label $mb.strip -text {} -anchor e
    ttk::button $mb.stop -text Stop -width 5 -style CatApply.TButton -command ::ogf::job::cancel
    pack $mb.workflow $mb.tools -side left -padx {0 2}
    pack $mb.strip -side right -padx 4
    set ogfui(strip) $mb.strip
    # The old 11-menu bar requested 620 px, which fixed the width of the catalog pane (and so the image
    # area).  Keep that requested width so the main layout does not change.
    set ogfui(panelwidth) 561
    $mb configure -width $ogfui(panelwidth) -height 25
    pack propagate $mb 0
    # --- tab bar + content inside the info area (created before the other info widgets are packed)
    set info $f.info
    ttk::style configure CatTabs.TNotebook -tabmargins {0 0 0 0} -padding 0
    ttk::style configure CatTabs.TNotebook.Tab -padding {7 1}
    ttk::notebook $info.tabs -style CatTabs.TNotebook -height 1
    set ogfui(nb) $info.tabs
    foreach t $::ogf::tabs {
	set pg [frame $info.tabs.pg_[OGFUITabId $t] -height 1 -width 1]
	$info.tabs add $pg -text $t
    }
    bind $info.tabs <<NotebookTabChanged>> OGFUITabChanged
    ttk::frame $info.content -height 24
    pack propagate $info.content 0
    set ogfui(content) $info.content
    foreach t $::ogf::tabs {
	set cf [ttk::frame $info.content.[OGFUITabId $t]]
	set ogfui(tabframe,$t) $cf
	set col 0
	foreach id [::ogf::reg::on_tab $t] {
	    if {[::ogf::json::get [::ogf::reg::get $id] menu] ne {}} continue
	    OGFUIChip $cf $id $col
	    incr col
	}
	if {$t eq "Time-domain"} {catch {OGFUITimeDomainControls $cf}}
    }
    # plugin menus living in Workflow / Tools
    OGFUIBuildMainMenus
    trace add variable ::ogfsess(steps) write OGFUIProgressSoon
    OGFUIShowTab Detect
    OGFUIProgressRefresh
}

proc OGFUITabId {t} {return [string tolower [string map {- _ { } _} $t]]}

proc OGFUIChip {parent id col} {
    global ogfui
    set pm [::ogf::reg::get $id]
    set c [ttk::frame $parent.c_$id]
    set mbtn [ttk::menubutton $c.m -text "[::ogf::json::get $pm short [dict get $pm name]] \u25be" \
	-style CatChip.TMenubutton -menu $c.m.m]
    menu $c.m.m -tearoff 0
    OGFUIFillPluginMenu $c.m.m $id
    set ogfui(menu,$id) $c.m.m
    pack $mbtn -side left
    set prim [::ogf::json::get $pm primary]
    if {$prim ne {}} {
	set s [::ogf::reg::step $id $prim]
	set cmd [OGFUIStepCommand $id $s]
	ttk::button $c.run -text "\u25b6" -width 2 -style CatIcon.TButton -command $cmd
	pack $c.run -side left
	OGFUITip $c.run "Run: [dict get $s label]"
	set ogfui(run,$id) $c.run
    }
    if {[::ogf::json::get $pm settings] ne {}} {
	ttk::button $c.set -text "\u2699" -width 2 -style CatIcon.TButton -command [list OGFUISettings $id]
	pack $c.set -side left
	OGFUITip $c.set "Settings..."
	set ogfui(set,$id) $c.set
    }
    OGFUITip $mbtn [::ogf::json::get $pm description [dict get $pm name]]
    pack $c -side left -padx {0 6}
}

proc OGFUITip {w text} {
    catch {tooltip::tooltip $w $text}
}

proc OGFUIBuildMainMenus {} {
    global ogfui
    set mb $ogfui(root).menubar
    set wm $mb.workflow.m
    foreach t $::ogf::tabs {
	$wm add radiobutton -label $t -variable ogfui(tab) -value $t -command [list OGFUIShowTab $t]
    }
    # plugins that live in the Workflow menu (session recorder)
    foreach id [::ogf::reg::ids] {
	if {[::ogf::json::get [::ogf::reg::get $id] menu] eq "Workflow"} {
	    $wm add separator
	    OGFUIFillPluginMenu $wm $id
	}
    }
    set tm $mb.tools.m
    $tm add checkbutton -label "Detach Catalog Panel" -variable catpanel(detached) -command CatalogPanelToggleDetach
    $tm add separator
    $tm add checkbutton -label "Tile all frames" -variable ogfui(tile) -command OGFUITileToggle
    menu $tm.tl -tearoff 0
    foreach {v l} {grid Grid column Column row Row} {
	$tm.tl add radiobutton -label $l -variable ogfui(tilemode) -value $v -command OGFUITileOptions
    }
    $tm add cascade -label "Tile layout" -menu $tm.tl
    $tm add checkbutton -label "Lock pan/zoom (WCS) in tile" -variable ogfui(lockwcs) -command OGFUITileOptions
    $tm add checkbutton -label "Share scale limits in tile" -variable ogfui(sharescale) -command OGFUITileOptions
    $tm add separator
    menu $tm.par -tearoff 0 -postcommand [list OGFUIPostParams $tm.par]
    $tm add cascade -label "Plugin settings" -menu $tm.par
    $tm add command -label "Pipeline overview..." -command OGFUIOverview
    $tm add command -label "Plugins..." -command OGFUIPluginManager
    $tm add command -label "Plugin log..." -command OGFUIShowLog
    foreach id [::ogf::reg::ids] {
	if {[::ogf::json::get [::ogf::reg::get $id] menu] eq "Tools"} {
	    $tm add separator
	    OGFUIFillPluginMenu $tm $id
	}
    }
}

proc OGFUIPostParams {m} {
    $m delete 0 end
    foreach id [::ogf::reg::ids] {
	set pm [::ogf::reg::get $id]
	if {[::ogf::json::get $pm params] ne {}} {
	    $m add command -label [dict get $pm name] -command [list OGFParamDialog $id]
	}
    }
}

# ---------------------------------------------------------------- tabs
proc OGFUIShowTab {t} {
    global ogfui
    set ogfui(tab) $t
    foreach x $::ogf::tabs {
	catch {pack forget $ogfui(tabframe,$x)}
    }
    pack $ogfui(tabframe,$t) -in $ogfui(content) -side left -fill y -padx 2
    set i [lsearch -exact $::ogf::tabs $t]
    if {[$ogfui(nb) index current] != $i} {$ogfui(nb) select $i}
}

proc OGFUITabChanged {} {
    global ogfui
    set i [$ogfui(nb) index current]
    set t [lindex $::ogf::tabs $i]
    if {$t ne {} && $t ne $ogfui(tab)} {OGFUIShowTab $t}
}

# widget-tree dump of a tab: list of {plugin steps...} actually mapped (for the verification scripts)
proc OGFUITabWidgets {t} {
    global ogfui
    set out {}
    foreach c [winfo children $ogfui(tabframe,$t)] {
	foreach w [winfo children $c] {
	    set cls [winfo class $w]
	    if {$cls eq "TMenubutton"} {
		lappend out [list menubutton $w [$w cget -text]]
	    } elseif {$cls eq "TButton"} {
		lappend out [list button $w [$w cget -text]]
	    } else {
		lappend out [list $cls $w]
	    }
	}
    }
    return $out
}

proc OGFUIMenuLabels {m} {
    set out {}
    set n [$m index end]
    if {$n eq "none"} {return $out}
    for {set i 0} {$i <= $n} {incr i} {
	if {[$m type $i] in {separator tearoff}} continue
	lappend out [$m entrycget $i -label]
    }
    return $out
}

# ---------------------------------------------------------------- progress strip
proc OGFUIProgressSoon {args} {
    global ogfui
    if {$ogfui(pending) ne {}} return
    set ogfui(pending) [after idle OGFUIProgressRefresh]
}

proc OGFUIProgressRefresh {} {
    global ogfui
    set ogfui(pending) {}
    if {![info exists ogfui(strip)] || ![winfo exists $ogfui(strip)]} return
    set p [::ogf::step::progress]
    set txt {}
    set i 0
    foreach st $::ogf::stages {
	lassign [dict get $p $st] state sec n
	set name [string totitle $st]
	switch -- $state {
	    done {append txt [format "%s %.1fs" $name $sec]}
	    failed {append txt "$name FAILED"}
	    default {append txt $name}
	}
	incr i
	if {$i < [llength $::ogf::stages]} {append txt "  >  "}
    }
    $ogfui(strip) configure -text $txt
    set ogfui(strip,state) $p
}

proc OGFUIJobState {busy} {
    global ogfui
    set mb $ogfui(root).menubar
    if {$busy} {pack $mb.stop -side right -padx 2 -before $mb.strip} else {pack forget $mb.stop}
}

# tile / single display: see ogf_tile.tcl (OGFUIDisplay, OGFUITileToggle, OGFUITileOptions)
proc OGFUIAfterLayout {} {
    catch {CatalogPanelSyncInfoHeight}
}

# ---------------------------------------------------------------- overview / plugin manager / log
proc OGFUIOverview {} {
    set t "OGFinder workflow overview (generated from the plugin registry)\n"
    append t "session class: AUTO = replayed in pipeline mode, CONFIG = replayed in replay mode only, MANUAL = needs a person, none = not recorded\n"
    foreach tab $::ogf::tabs {
	append t "\n== $tab ==\n"
	foreach id [::ogf::reg::on_tab $tab] {
	    set pm [::ogf::reg::get $id]
	    append t "  [dict get $pm name]  ($id)\n"
	    foreach s [::ogf::json::get $pm steps] {
		set cls [::ogf::json::get $s session NONE]
		set kind [expr {[dict exists $s cli] ? "cli" : "tcl"}]
		append t [format "     %-44s %-7s %-4s %s\n" [dict get $s label] $cls $kind [join [::ogf::json::get $s records] ,]]
	    }
	}
    }
    OGFTextWindow "OGFinder pipeline overview" $t
    return $t
}

proc OGFUIShowLog {} {
    OGFTextWindow "OGFinder plugin log" [join $::ogf::logbuf \n]
}

proc OGFUIPluginManager {} {
    global ogfui
    set w .ogfplugins
    catch {destroy $w}
    toplevel $w
    wm title $w "OGFinder plugins"
    wm transient $w .
    set r 0
    foreach id [::ogf::reg::ids 0] {
	set pm [::ogf::reg::get $id]
	set ogfui(pm,$id) [dict get $pm _enabled]
	ttk::checkbutton $w.c$r -text "[dict get $pm name] ($id)  [expr {[::ogf::json::get $pm example 0] ? {- example} : {}}]" \
	    -variable ogfui(pm,$id)
	ttk::label $w.d$r -text "tab: [dict get $pm tab]   [::ogf::json::get $pm description]" -foreground gray30 -wraplength 520 -justify left
	grid $w.c$r -row [expr {2*$r}] -column 0 -sticky w -padx 8 -pady {4 0}
	grid $w.d$r -row [expr {2*$r+1}] -column 0 -sticky w -padx 28
	incr r
    }
    ttk::label $w.note -text "Changes take effect after restarting OGFinder." 
    ttk::frame $w.bb
    ttk::button $w.bb.ok -text Save -command OGFUIPluginManagerSave
    ttk::button $w.bb.c -text Close -command [list destroy $w]
    pack $w.bb.ok $w.bb.c -side left -padx 4
    grid $w.note -row [expr {2*$r}] -column 0 -sticky w -padx 8 -pady 4
    grid $w.bb -row [expr {2*$r+1}] -column 0 -pady 6
}

proc OGFUIPluginManagerSave {} {
    global ogfui
    set en {}; set dis {}
    foreach id [::ogf::reg::ids 0] {
	if {$ogfui(pm,$id)} {lappend en $id} else {lappend dis $id}
    }
    ::ogf::reg::save_prefs $en $dis
    ::ogf::status "plugin preferences saved (restart to apply)"
    destroy .ogfplugins
}

# ---------------------------------------------------------------- startup glue
# called at the start of CreateCatalogPanel: discover manifests and source plugin Tcl files
proc OGFCoreStart {} {
    ::ogf::reg::load_all
}
# called at the end of CreateCatalogPanel: plugin init hooks, parameters
proc OGFCoreReady {} {
    ::ogf::reg::init_all
}
