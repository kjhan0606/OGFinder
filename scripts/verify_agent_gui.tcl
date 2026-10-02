# GUI test of the agent-CLI backends (run by scripts/verify_agent_gui.py under Xvfb with FAKE agent CLIs on PATH):
#   ds9 m51.fits -source verify_agent_gui.tcl      env: OGF_AGENT_OUT (dir), OGF_AI_SERVICES_FILE
global catpanel ogfai ogfsess
set ::answers yes
set ::nconfirm 0
set ::msgs {}
rename tk_messageBox OGF_real_messageBox
proc tk_messageBox {args} {
    incr ::nconfirm
    set d [dict create {*}$args]
    lappend ::msgs [dict get $d -message]
    return $::answers
}
set out $::env(OGF_AGENT_OUT)
set ::vlog [open [file join $out agent_gui_steps.log] w]
proc vlog {m} {puts $::vlog $m; flush $::vlog; puts "AGENTVERIFY: $m"}
proc ncols {} {global catpanel; return [llength [split [lindex [split $catpanel(alldata) \n] 0] \t]]}
proc hdr {} {global catpanel; return [split [lindex [split $catpanel(alldata) \n] 0] \t]}

# open the Run Task dialog, set backend/task/images/service, optionally edit fields, press Run (or Cancel)
proc run_dialog {be task {images 0} {press ok} {svc {}} {pre {}}} {
    global ogfai
    set ::dlg {}
    after 300 [list dlg_drive $be $task $images $press $svc $pre]
    OGFAIRunDialog
}
proc dlg_drive {be task images press svc pre} {
    global ogfai
    set w .ogfai_run
    set ogfai(run,task) $task
    set ogfai(run,bshow) [OGFAIBackendLabel $be]
    OGFAIRunBackendPicked
    if {$svc ne {}} {set ogfai(run,service) $svc; OGFAIRunServiceChanged}
    set ogfai(run,images) $images
    update
    set ::dlg [dict create backends [$w.cb cget -values] backend $ogfai(run,backend) services [$w.cs cget -values] \
        service $ogfai(run,service) agframe [winfo ismapped $w.ag] ast [$w.ag.st cget -text] exe $ogfai(run,exe) \
        xargs $ogfai(run,xargs) imgbox $ogfai(run,images) runstate [$w.bb.ok state]]
    if {$pre ne {}} {uplevel #0 $pre}
    set ogfai(run,done) $press
}

# number of recorded fake-CLI calls that carried a prompt (argv contains -p / exec / --print / --prompt-file)
proc prompt_calls {name} {
    set f [file join $::env(OGF_AGENT_BIN) $name.calls.jsonl]
    if {![file exists $f]} {return 0}
    set fd [open $f]; set d [read $fd]; close $fd
    set n 0
    foreach l [split $d \n] {if {$l ne {} && ![string match {*"argv": \["--version"\]*} $l]} {incr n}}
    return $n
}
proc save_agent {} {
    global ogfai
    set ogfai(run,exe) [file join $::env(OGF_AGENT_BIN) claude]
    set ogfai(run,xargs) {["--model","haiku"]}
    OGFAIRunSaveAgent
    set ::savedmsg [.ogfai_run.ag.st cget -text]
}

proc agent_session {} {
    global catpanel ogfai ogfsess out
    set catpanel(param,detect-thresh) 5.0
    CatalogPanelExtract
    set c0 [ncols]
    vlog "extract: cols=$c0"

    # ---- detection table (fake claude + codex + gemini on PATH; agy and grok missing)
    set ai [OGFAIAgentInfo]
    foreach n {agent_codex agent_claude agent_agy agent_grok} {
	vlog "detect $n: installed=[dict get $ai $n installed] path=[dict get $ai $n path]"
    }
    OGFAIRegistry; update
    set w .ogfai_reg
    vlog "registry rows: [$w.tv children {}]"
    vlog "registry note agent_claude: [lindex [$w.tv item agent_claude -values] 6] / agent_grok: [lindex [$w.tv item agent_grok -values] 6]"
    OGFAIRegDetect; update
    vlog "registry detect text: [string map {\n |} [$w.msg cget -text]]"
    $w.tv selection set agent_claude
    OGFAIRegTest
    vlog "registry test claude: [string map {\n |} [$w.msg cget -text]]"
    vlog "confirmation dialogs after detect/test: $::nconfirm"
    destroy $w

    # ---- Run dialog: backend dropdown lists the four + generic, with installed / not found
    run_dialog claude star_galaxy 0 cancel
    vlog "dialog backends: [dict get $::dlg backends]"
    vlog "dialog claude: services=[dict get $::dlg services] agframe=[dict get $::dlg agframe] st=[dict get $::dlg ast]"
    run_dialog generic star_galaxy 0 cancel
    vlog "dialog generic: services=[dict get $::dlg services] agframe=[dict get $::dlg agframe]"
    run_dialog grok star_galaxy 0 cancel
    vlog "dialog grok: services=[dict get $::dlg services] st=[dict get $::dlg ast]"

    vlog "prompt-carrying calls of the fakes before any confirmed run: claude=[prompt_calls claude] codex=[prompt_calls codex] gemini=[prompt_calls gemini]"
    # ---- not installed: no run, status says so
    set n0 $::nconfirm
    run_dialog grok star_galaxy 0 ok
    vlog "grok run (not installed): status=$catpanel(status) confirmations=[expr {$::nconfirm - $n0}] cols=[ncols]"

    # ---- dry run: nothing executed, exact prompt in the log
    run_dialog claude star_galaxy 0 ok {} {set ogfai(run,dry) 1}
    update
    vlog "dry run: status=$catpanel(status) cols=[ncols] confirmations=[expr {$::nconfirm - $n0}]"
    vlog "dry log has prompt: [string match *OBJECTS* $ogfai(lastlog)] argv: [string match {*--output-format*} $ogfai(lastlog)]"
    destroy .ogfai_log

    # ---- decline the confirmation: nothing runs
    set ::answers no
    set n0 $::nconfirm
    run_dialog claude star_galaxy 0 ok
    vlog "declined: status=$catpanel(status) confirmations=[expr {$::nconfirm - $n0}] cols=[ncols]"
    set ::answers yes

    # ---- real (fake-CLI) run with Claude: confirmation text, columns, provenance
    set n0 $::nconfirm
    run_dialog claude star_galaxy 0 ok
    vlog "claude run: status=$catpanel(status) cols $c0 -> [ncols] new=[lrange [hdr] $c0 end]"
    vlog "claude confirmation shown: [expr {$::nconfirm - $n0}]"
    vlog "CONFIRMTEXT: [string map {\n |} [lindex $::msgs end]]"
    # second run: no new confirmation (once per backend and session)
    set n0 $::nconfirm
    set c1 [ncols]
    run_dialog claude star_galaxy 0 ok
    vlog "claude 2nd run: confirmations=[expr {$::nconfirm - $n0}] cols=[ncols] (was $c1)"

    # ---- codex with images unticked, then ticked (asks again, images named in the text)
    set c2 [ncols]
    run_dialog codex morphology 0 ok
    vlog "codex calls after the no-image run: [prompt_calls codex]"
    vlog "codex morphology no images: status=$catpanel(status) new=[lrange [hdr] $c2 end]"
    vlog "CONFIRMTEXT2: [string map {\n |} [lindex $::msgs end]]"
    set n0 $::nconfirm
    run_dialog codex morphology 1 ok {} {set ogfai(run,size) 24}
    vlog "codex morphology WITH images: confirmations=[expr {$::nconfirm - $n0}] status=$catpanel(status)"
    vlog "CONFIRMTEXT3: [string map {\n |} [lindex $::msgs end]]"

    # ---- agy service falls back to the gemini executable
    run_dialog agy real_bogus 0 ok
    vlog "agy(gemini) real_bogus: status=$catpanel(status)"

    # ---- save executable / extra args to the profile from the dialog
    run_dialog claude photoz 0 ok {} save_agent
    vlog "save to profile: $::savedmsg"
    vlog "profile file now: [string map {\n { }} [exec cat $::env(OGF_AI_SERVICES_FILE)]]"

    # ---- recorder
    set sp [file join $out ogfinder_session.py]
    CatalogPanelSaveCatalogTo [file join $out catalog_gui.tsv]
    CatalogPanelSessionSave $sp
    vlog "session exported: [llength $ogfsess(steps)] steps"
    vlog "session log:\n[OGFSessLogText]"
}
after 2500 {
    if {[catch agent_session err]} {
	vlog "ERROR: $err\n$::errorInfo"
	exit 3
    }
    vlog DONE
    exit 0
}
