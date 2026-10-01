# GUI smoke test for "Analysis > AI Services" (run by scripts/verify_ai_gui.py under Xvfb):
#   ds9 m51.fits -source verify_ai_gui.tcl      env: OGF_AI_OUT (dir), OGF_AI_SERVICES_FILE (optional)
# Calls the procs the menu entries call; no separate code path.
global catpanel ogfai ogfsess ds9
# count the "data leaves your computer" confirmations; answer yes
set ::nconfirm 0
rename tk_messageBox OGF_real_messageBox
proc tk_messageBox {args} {
    incr ::nconfirm
    set ::lastmsg $args
    return yes
}
set out $::env(OGF_AI_OUT)
file mkdir $out
set ::vlog [open [file join $out ai_gui_steps.log] w]
proc vlog {m} {puts $::vlog $m; flush $::vlog; puts "AIVERIFY: $m"}
proc ncols {} {global catpanel; return [llength [split [lindex [split $catpanel(alldata) \n] 0] \t]]}
proc nrows {} {global catpanel; set n 0; foreach l [lrange [split $catpanel(alldata) \n] 1 end] {if {[string trim $l] ne {}} {incr n}}; return $n}
proc hdr {} {global catpanel; return [split [lindex [split $catpanel(alldata) \n] 0] \t]}

proc ai_session {} {
    global catpanel ogfai ogfsess ds9 out
    set catpanel(param,detect-thresh) 5.0
    CatalogPanelExtract
    set n0 [nrows]; set c0 [ncols]
    vlog "extract: $n0 rows, $c0 cols"
    # ---- menu structure
    # AI Services plugin menu (workflow UI: Results tab, chip "AI Services")
    set sub [OGFUIPluginMenu ai_services]
    set labels {}
    for {set i 0} {$i <= [$sub index end]} {incr i} {
	if {[$sub type $i] in {separator tearoff}} continue
	lappend labels [$sub entrycget $i -label]
    }
    vlog "menu AI Services entries: $labels"
    # ---- registry dialog
    OGFAIRegistry
    update
    set w .ogfai_reg
    vlog "registry dialog exists: [winfo exists $w] rows=[llength [$w.tv children {}]] items=[$w.tv children {}]"
    vlog "registry message: [$w.msg cget -text]"
    $w.tv selection set mock
    OGFAIRegTest
    vlog "registry test(mock): [$w.msg cget -text]"
    $w.tv selection set mock
    update
    destroy $w
    # ---- dry run on selected rows
    CatalogPanelLinkSelect 1 replace 0
    CatalogPanelLinkSelect 2 add 0
    set ok [OGFAIRunTask photoz mock -rows selected -dry 1]
    vlog "dry run selected rows ok=$ok; columns unchanged: [expr {[ncols] == $c0}]; sessions steps: [llength $ogfsess(steps)]"
    destroy .ogfai_log
    # ---- real run (mock) all rows: photo-z -> columns
    set ok [OGFAIRunTask photoz mock -rename {PHOTOZ=PHOTO_Z PHOTOZ_ERR=PHOTO_Z_ERR}]
    vlog "run photoz/mock ok=$ok rows=[nrows] cols=[ncols] new=[lrange [hdr] $c0 end]"
    set c1 [ncols]
    # ---- morphology with cutouts (size in arcsec) on the mock
    set ok [OGFAIRunTask morphology mock -size 20 -unit arcsec -norm zscale -fmt png]
    vlog "run morphology/mock ok=$ok cols=[ncols] new=[lrange [hdr] $c1 end]"
    # ---- selected rows subset -> MANUAL step
    CatalogPanelLinkSelect 3 replace 0
    CatalogPanelLinkSelect 4 add 0
    set ok [OGFAIRunTask star_galaxy mock -rows selected]
    vlog "run star_galaxy/mock on selected ok=$ok cols=[ncols]"
    # sort/filter work on new columns
    CatalogPanelSort PHOTO_Z descending
    vlog "sorted by PHOTO_Z; first rows: [lrange [split $catpanel(alldata) \n] 1 2]"
    # ---- recorder
    OGFAIShowLog
    update
    vlog "last log window: [winfo exists .ogfai_log] first line: [lindex [split $ogfai(lastlog) \n] 0]"
    # ---- network service (local test HTTP server, env OGF_AI_REST=1): confirmation once per session
    if {[info exists ::env(OGF_AI_REST)]} {
	set n0 $::nconfirm
	set c2 [ncols]
	set ok [OGFAIRunTask photoz localrest -rename {PHOTOZ=REST_Z PHOTOZ_ERR=REST_Z_ERR}]
	vlog "run photoz/localrest ok=$ok cols $c2 -> [ncols] new=[lrange [hdr] $c2 end]"
	set ok2 [OGFAIRunTask photoz localrest -rename {PHOTOZ=REST_Z PHOTOZ_ERR=REST_Z_ERR}]
	vlog "second localrest run ok=$ok2 (columns already exist -> updated in place) cols=[ncols]"
	vlog "confirmation dialogs shown for network service: [expr {$::nconfirm - $n0}] (expected 1); text: [string range [dict get [dict create {*}$::lastmsg] -message] 0 160]"
	vlog "registry test of network service after confirmation: no extra dialog"
	OGFAIRegistry; update
	.ogfai_reg.tv selection set localrest
	OGFAIRegTest
	vlog "test localrest: [.ogfai_reg.msg cget -text] ; dialogs now [expr {$::nconfirm - $n0}]"
	destroy .ogfai_reg
    }
    # ---- run-dialog widgets build
    after 100 {set ogfai(run,done) cancel}
    OGFAIRunDialog
    vlog "run dialog opened and cancelled"
    # ---- backend hook default = local: with backend external -> external service columns
    set ogfai(backend,photoz) mock
    set before [ncols]
    CatalogPanelPhotoZ
    vlog "backend hook photoz=mock: status=$catpanel(status) cols $before -> [ncols]"
    set ogfai(backend,photoz) local
    set sp [file join $out ogfinder_session.py]
    set cf [file join $out catalog_gui.tsv]
    CatalogPanelSaveCatalogTo $cf
    CatalogPanelSessionSave $sp
    vlog "session exported: $sp  (steps: [llength $ogfsess(steps)])"
    vlog "session log:\n[OGFSessLogText]"
}
after 2500 {
    if {[catch ai_session err]} {
	vlog "ERROR: $err\n$::errorInfo"
	exit 3
    }
    vlog DONE
    exit 0
}
