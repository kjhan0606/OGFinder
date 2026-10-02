#  OGFinder: "Analysis > AI Services" - connection points for external AI / astronomy services.
#
#  The GUI only builds and runs a ds9_ai_bridge.py command line (python package ai_bridge/); the
#  command line is what the session recorder stores (step "ai.run"), so the exported Python
#  pipeline script can replay it.  No models, no secrets: profiles name environment variables,
#  the GUI only ever shows "set"/"unset".  See docs/ai_services.md.
#
#  Menu (3 entries):   Service Registry...   Run Task on Catalog...   Show Last Run Log

package provide DS9 1.0

# keep in sync with ai_bridge/contracts.py (checked by scripts/verify_ai_gui.tcl)
set ::OGFAI_TASKS {photoz sed_fit morphology star_galaxy real_bogus transient_classification
    moving_object_classification anomaly_detection embedding_similarity captioning image_retrieval generic}

# existing local-model steps that can use an external service instead:  step -> {task renames}
set ::OGFAI_BACKEND_STEPS {
    photoz     {photoz     {PHOTOZ=PHOTO_Z PHOTOZ_ERR=PHOTO_Z_ERR}   "Photo-z (AI)"}
    sed_fit    {sed_fit    {}                                        "SED Fitting (AI)"}
    morphology {morphology {}                                        "Galaxy morphology"}
    star       {star_galaxy {}                                       "AI star classification"}
}

proc OGFAIInit {} {
    global ogfai
    if {[info exists ogfai(init)]} return
    set ogfai(init) 1
    set ogfai(sf) {}
    set ogfai(lastlog) "No AI Services run yet in this session."
    set ogfai(confirmed) {}
    foreach s {photoz sed_fit morphology star} {set ogfai(backend,$s) local}
    set ogfai(last,backend) generic
    set ogfai(last,size) {}
    set ogfai(last,unit) pix
    set ogfai(last,norm) {(profile)}
    set ogfai(last,fmt) {(profile)}
    set prf [OGFAIPrefFile]
    if {[file exists $prf] && ![catch {set fd [open $prf r]}]} {
	set data [read $fd]; close $fd
	foreach line [split $data \n] {
	    set line [string trim $line]
	    if {$line eq {} || [string index $line 0] eq "#"} continue
	    set eq [string first = $line]
	    if {$eq < 0} continue
	    set k [string trim [string range $line 0 $eq-1]]
	    set v [string trim [string range $line $eq+1 end]]
	    if {$k eq "services_file"} {set ogfai(sf) $v} elseif {$k eq "run_backend" && $v in {generic codex claude agy grok}} {set ogfai(last,backend) $v} elseif {[string match backend,* $k]} {set ogfai($k) $v}
	}
    }
}

proc OGFAIPrefFile {} {return [file join [OGFSessWorkDir] ai_services.prf]}

proc OGFAIPrefSave {} {
    global ogfai
    catch {file mkdir [OGFSessWorkDir]}
    if {[catch {set fd [open [OGFAIPrefFile] w]}]} return
    puts $fd "services_file=$ogfai(sf)"
    puts $fd "run_backend=$ogfai(last,backend)"
    foreach s {photoz sed_fit morphology star} {puts $fd "backend,$s=$ogfai(backend,$s)"}
    close $fd
}

proc OGFAIBuildMenu {m} {
    OGFAIInit
    menu $m -tearoff 0
    $m add command -label "Service Registry..." -command OGFAIRegistry
    $m add command -label "Run Task on Catalog..." -command OGFAIRunDialog
    $m add command -label "Show Last Run Log" -command OGFAIShowLog
    return $m
}

# ------------------------------------------------------------------ bridge plumbing
proc OGFAIScript {} {return [CatalogPanelGetScript ds9_ai_bridge.py]}

# argv prefix: python script --mode MODE [--services-file F]
proc OGFAIArgv {mode} {
    global ogfai
    OGFAIInit
    set a [list [OGFPython] [OGFAIScript] --mode $mode]
    if {$ogfai(sf) ne {}} {lappend a --services-file $ogfai(sf)}
    return $a
}

proc OGFAIAppendLog {text} {
    global ogfai
    set ogfai(lastlog) $text
    catch {
	file mkdir [OGFSessWorkDir]
	set fd [open [file join [OGFSessWorkDir] ai_last_run.log] w]
	puts $fd $text
	close $fd
    }
}

# run an argv, stderr to a file -> {rc stdout stderr}
proc OGFAIExec {argv} {
    set errf [file join [OGFSessWorkDir] ai_stderr.txt]
    catch {file mkdir [OGFSessWorkDir]}
    set rc [catch {exec {*}$argv 2>$errf} out]
    set err {}
    catch {set fd [open $errf r]; set err [read $fd]; close $fd}
    return [list $rc $out $err]
}

# list of dicts: name task transport enabled valid env envstatus target sha note ; plus file path / error in the 2nd value
proc OGFAIServices {} {
    if {![file exists [OGFAIScript]]} {return [list {} {} "ds9_ai_bridge.py not found"]}
    lassign [OGFAIExec [OGFAIArgv list-services]] rc out err
    set rows {}
    set path {}
    set ferr {}
    foreach l [split $out \n] {
	if {[string match "# profile file: *" $l]} {set path [string range $l 16 end]; continue}
	if {[string match "# ERROR: *" $l]} {set ferr [string range $l 9 end]; continue}
	if {$l eq {} || [string index $l 0] eq "#" || [string match "name\ttask\t*" $l]} continue
	lassign [split $l \t] name task transport enabled valid env st target sha note
	lappend rows [dict create name $name task $task transport $transport enabled $enabled valid $valid \
	    env $env envstatus $st target $target sha $sha note $note]
    }
    if {$rc && $ferr eq {}} {set ferr [string trim $err]}
    return [list $rows $path $ferr]
}

proc OGFAIServiceInfo {name} {
    lassign [OGFAIServices] rows path err
    foreach r $rows {if {[dict get $r name] eq $name} {return $r}}
    return {}
}

proc OGFAIIsNetwork {info} {
    return [expr {[dict exists $info transport] && [dict get $info transport] in {https_json https_multipart tap_query}}]
}

proc OGFAIIsAgent {info} {
    return [expr {[dict exists $info transport] && [dict get $info transport] eq "agent_cli"}]
}

# provider behind each agent CLI (shown in the confirmation; the CLI forwards the prompt to its cloud model)
set ::OGFAI_AGENT_PROVIDER {codex OpenAI claude Anthropic agy Google gemini Google grok xAI}

# agent_cli services: dict name -> {backend label installed path executable extra_args} from list-services --json
proc OGFAIAgentInfo {} {
    set res {}
    if {![file exists [OGFAIScript]]} {return $res}
    lassign [OGFAIExec [concat [OGFAIArgv list-services] --json]] rc out err
    if {$rc || [catch {set doc [::ogf::json::parse $out]}]} {return $res}
    foreach sv [::ogf::json::get $doc services] {
	set ag [::ogf::json::get $sv agent]
	if {$ag eq {}} continue
	dict set res [dict get $sv name] [dict create backend [::ogf::json::get $ag backend] \
	    label [::ogf::json::get $ag label] installed [::ogf::json::get $ag installed 0] path [::ogf::json::get $ag path] \
	    executable [::ogf::json::get $ag executable] extra_args [::ogf::json::get $ag extra_args] \
	    login [::ogf::json::get $ag login] enabled [expr {[::ogf::json::get $sv enabled 1] ? "yes" : "no"}] \
	    valid [expr {[::ogf::json::get $sv valid 1] ? "yes" : "no"}]]
    }
    return $res
}

# what exactly leaves the machine, built from the bridge's own dry run (--summary-only: no prompt text, no rows)
# -> dict (empty on failure)
proc OGFAIAgentSummary {sumargv} {
    lassign [OGFAIExec $sumargv] rc out err
    if {$rc || [catch {set d [::ogf::json::parse $out]}]} {return [dict create error [string trim $err]]}
    return $d
}

# confirmation for an agent CLI: once per backend and session (a run with image cutouts asks once more)
proc OGFAIConfirmAgent {info sumargv sendimg} {
    global ogfai OGFAI_AGENT_PROVIDER
    OGFAIInit
    set name [dict get $info name]
    set key "$name:[expr {$sendimg ? {images} : {rows}}]"
    if {[lsearch -exact $ogfai(confirmed) $key] >= 0} {return 1}
    set d [OGFAIAgentSummary $sumargv]
    if {[dict exists $d error]} {
	tk_messageBox -type ok -icon error -title "AI Services: $name" \
	    -message "Cannot prepare the request, nothing was sent:\n[dict get $d error]"
	return 0
    }
    set dl [::ogf::json::get $d data_leaving]
    set be [::ogf::json::get $dl backend]
    set prov [expr {[dict exists $OGFAI_AGENT_PROVIDER $be] ? [dict get $OGFAI_AGENT_PROVIDER $be] : "the CLI's provider"}]
    set nimg [::ogf::json::get $d images_total 0]
    set msg "'$name' hands data to the agent CLI [::ogf::json::get $dl backend]:\n  [::ogf::json::get $dl executable]\n\n"
    append msg "The CLI forwards the prompt to $prov's cloud model under the CLI's own login, so this data LEAVES YOUR COMPUTER:\n"
    append msg "  - [::ogf::json::get $dl objects] catalog object(s), per object: [join [::ogf::json::get $dl record_fields] {, }]\n"
    append msg "  - the full catalog row text: [expr {[::ogf::json::get $dl catalog_row_sent 0] ? {YES (this task needs it)} : {no}}]\n"
    if {$sendimg && $nimg > 0} {
	append msg "  - IMAGE DATA: $nimg cutout file(s) (you ticked 'Send image cutouts')\n"
    } else {
	append msg "  - image data: NONE (cutouts are not sent unless ticked)\n"
    }
    append msg "  - about [expr {[::ogf::json::get $d prompt_bytes_total 0] / 1024 + 1}] kB of prompt in [::ogf::json::get $d requests 1] request(s)\n"
    set envs [::ogf::json::get $dl env_names_passed]
    append msg "\nEnvironment variables handed to the CLI (names only): [expr {[llength $envs] ? [join $envs {, }] : {none}}]\n"
    append msg "No credentials are stored by OGFinder. The CLI runs with tools disabled in an empty temporary folder.\n"
    append msg "To see the exact prompt and command first, cancel and use 'Dry run'.\n"
    append msg "\nContinue? (asked once per backend per session)"
    if {[tk_messageBox -type yesno -icon warning -title "AI Services: send data to $name" -message $msg] ne "yes"} {return 0}
    lappend ogfai(confirmed) $key
    return 1
}

# one confirmation per service and session before anything leaves the machine
proc OGFAIConfirm {info what} {
    global ogfai
    OGFAIInit
    set name [dict get $info name]
    if {![OGFAIIsNetwork $info]} {return 1}
    if {[lsearch -exact $ogfai(confirmed) $name] >= 0} {return 1}
    set tgt [dict get $info target]
    set msg "Service '$name' is an external network service.\n\n$what\n\nTarget: $tgt\n"
    if {[dict get $info env] ne {-}} {append msg "Credentials are read from the environment variable [dict get $info env] ([dict get $info envstatus]).\n"}
    append msg "\nThis data leaves your computer. Continue? (asked once per service per session)"
    if {[tk_messageBox -type yesno -icon warning -title "AI Services: send data to $name" -message $msg] ne "yes"} {return 0}
    lappend ogfai(confirmed) $name
    return 1
}

proc OGFAIRowCount {} {
    if {![::ogf::cat::has]} {return 0}
    set n 0
    foreach l [lrange [split [::ogf::cat::tsv] \n] 1 end] {if {[string trim $l] ne {}} {incr n}}
    return $n
}

# images for cutouts: registered bands (detection band first = pixel grid of X_IMAGE) or the current frame
proc OGFAIImages {} {
    global ogfband
    set names {}; set files {}
    if {[info exists ogfband(names)] && [llength $ogfband(names)]} {
	set det $ogfband(detect)
	set order [list $det]
	foreach b [OGFBandsSorted] {if {$b ne $det} {lappend order $b}}
	foreach b $order {
	    if {$b eq {} || ![info exists ogfband($b,file)] || ![file exists $ogfband($b,file)]} continue
	    lappend names $b
	    lappend files [file normalize $ogfband($b,file)]
	}
    }
    if {[llength $files] == 0} {
	set fn [CatalogPanelGetFITS]
	if {$fn ne {} && [file exists $fn]} {lappend names main; lappend files [file normalize $fn]}
    }
    return [list $names $files]
}

# ------------------------------------------------------------------ run a task
# OGFAIRunTask TASK SERVICE ?-rows all|selected -dry 0|1 -size N -unit pix|arcsec -norm X -fmt X -rename {A=B ...}
#                              -title T -params {k=v ...} -numbers {..}?
# Returns 1 when columns were added (or a dry run completed), 0 otherwise.
proc OGFAIRunTask {task service args} {
    global ogfai ogfsess
    OGFAIInit
    set seq 0
    array set o {-rows all -dry 0 -size {} -unit pix -norm {} -fmt {} -rename {} -title {} -params {} -numbers {} -images 0}
    array set o $args
    if {![::ogf::cat::has]} {
	::ogf::cat::set status "AI Services: no catalog - extract sources first"
	return 0
    }
    if {$service eq {}} {::ogf::cat::set status "AI Services: no service chosen"; return 0}
    set script [OGFAIScript]
    if {![file exists $script]} {::ogf::cat::set status "AI Services: ds9_ai_bridge.py not found"; return 0}
    set info [OGFAIServiceInfo $service]
    if {$info eq {}} {::ogf::cat::set status "AI Services: service '$service' not found in the registry"; return 0}
    set numbers $o(-numbers)
    if {$o(-rows) eq "selected"} {
	set numbers [::ogf::cat::selection]
	if {[llength $numbers] == 0} {::ogf::cat::set status "AI Services: no rows selected in the table"; return 0}
    }
    set nrows [expr {[llength $numbers] ? [llength $numbers] : [OGFAIRowCount]}]
    set net [OGFAIIsNetwork $info]
    set agent [OGFAIIsAgent $info]
    set sendimg [expr {$agent && [string is true -strict $o(-images)]}]
    lassign [OGFAIImages] bands files
    set catfile [CatalogPanelSaveTempCatalog ai]
    if {$catfile eq {}} {::ogf::cat::set status "AI Services: cannot write the temporary catalog"; return 0}
    set dry [string is true -strict $o(-dry)]
    set argv [OGFAIArgv [expr {$dry ? "dry-run" : "run"}]]
    lappend argv --service $service --task $task --catalog $catfile
    if {[llength $files]} {lappend argv --image [join $files ,] --bands [join $bands ,]}
    if {[llength $numbers]} {lappend argv --numbers [join $numbers ,]}
    if {$o(-size) ne {}} {lappend argv [expr {$o(-unit) eq "arcsec" ? "--size-arcsec" : "--size-pix"}] $o(-size)}
    if {$o(-norm) ne {} && $o(-norm) ne "(profile)"} {lappend argv --normalize $o(-norm)}
    if {$o(-fmt) ne {} && $o(-fmt) ne "(profile)"} {lappend argv --cutout-format $o(-fmt)}
    foreach r $o(-rename) {lappend argv --rename $r}
    foreach p $o(-params) {lappend argv --param $p}
    if {$sendimg} {lappend argv --send-images}
    set sumargv [concat [lreplace $argv 3 3 dry-run] --summary-only]
    if {!$dry} {
	lappend argv --provenance-output [file join [OGFSessWorkDir] ai_last_run.provenance.json]
	# the user confirmed the transfer; the recorded argv keeps the flag so replay mode can run it.
	# Pipeline mode still refuses network steps unless the script itself gets --allow-network.
	if {$net} {lappend argv --allow-network}
	if {$agent} {lappend argv --allow-agent-cli}
    }
    set what "$nrows object(s) of the catalog (positions, magnitudes, catalog row text[expr {[llength $files] ? {, image cutouts} : {}}] as the profile's request template specifies) will be sent."
    if {$agent && !$dry} {
	if {![dict exists [OGFAIAgentInfo] $service] || ![dict get [OGFAIAgentInfo] $service installed]} {
	    ::ogf::cat::set status "AI Services: the agent CLI of '$service' was not found on PATH (Run Task dialog > Backend: enter its path and press 'Save to profile')"
	    return 0
	}
	set confirmed [OGFAIConfirmAgent $info $sumargv $sendimg]
    } else {
	set confirmed [expr {$dry || [OGFAIConfirm $info $what]}]
    }
    if {!$dry && !$confirmed} {
	::ogf::cat::set status "AI Services: cancelled"
	return 0
    }
    set title [expr {$o(-title) ne {} ? $o(-title) : "AI service $service / $task ($nrows objects)"}]
    ::ogf::cat::set status "AI Services: $service / $task on $nrows object(s) ..."
    update idletasks
    if {!$dry} {
	# session recorder: whole catalog = AUTO (replayed by the pipeline), a row subset = MANUAL (replay only)
	foreach im $files {if {[lsearch -exact $ogfsess(images) $im] < 0} {lappend ogfsess(images) $im}}
	set class [expr {[llength $numbers] ? "manual" : "auto"}]
	set seq 0
	catch {set seq [OGFSessLog ai.run $class $argv -title $title -tool python -requires catalog \
	    -network [expr {($net || $agent) ? 1 : 0}] \
	    -payload [expr {$agent ? [dict create service $service task $task profile_sha256 [dict get $info sha] \
		transport [dict get $info transport] objects $nrows backend [dict get [OGFAIAgentInfo] $service backend] \
		send_images [expr {$sendimg ? 1 : 0}]] : [dict create service $service task $task profile_sha256 [dict get $info sha] \
		transport [dict get $info transport] objects $nrows]}] \
	    -note [expr {[llength $numbers] ? "row subset chosen by hand (NUMBER list); replay only" : {}}]]}
    }
    set t0 [clock milliseconds]
    set xargv $argv
    lassign [OGFAIExec $xargv] rc data err
    set ms [expr {[clock milliseconds] - $t0}]
    set log "AI Services run  [clock format [clock seconds] -format {%Y-%m-%d %H:%M:%S %Z}]\nservice: $service   task: $task   mode: [expr {$dry ? {DRY RUN (nothing sent)} : {run}}]   objects: $nrows   elapsed: ${ms} ms\n\$ [join [lmap a $xargv {expr {[regexp {[\s]} $a] ? "'$a'" : $a}}] { }]\nexit status: [expr {$rc ? {non-zero} : 0}]\n"
    if {$dry} {append log "\n--- request preview (stdout) ---\n$data\n"}
    append log "\n--- bridge messages (stderr) ---\n$err\n"
    if {!$dry} {append log "\nprovenance: [file join [OGFSessWorkDir] ai_last_run.provenance.json]\n"}
    OGFAIAppendLog $log
    if {$dry} {
	::ogf::cat::set status "AI Services: dry run done (nothing was sent) - see Show Last Run Log"
	OGFAIShowLog
	return [expr {!$rc}]
    }
    if {$rc} {
	set last [lindex [split [string trim $err] \n] end]
	::ogf::cat::set status "AI Services: failed - $last"
	return 0
    }
    set hdr [split [lindex [split $data \n] 0] \t]
    set cols [lrange $hdr 1 end]
    if {[lindex $hdr 0] ne "NUMBER" || [llength $cols] == 0} {
	::ogf::cat::set status "AI Services: service returned no columns"
	return 0
    }
    # CatalogPanelAddColumnsFromTSV updates in place when ANY requested column already exists (and then silently
    # skips the new ones), so new and existing columns are merged in two calls (new first).  The recorded post
    # step "ai" repeats exactly this split in the exported script.
    set have [split [lindex [split [::ogf::cat::tsv] \n] 0] \t]
    set newc {}; set oldc {}
    foreach c $cols {if {$c in $have} {lappend oldc $c} else {lappend newc $c}}
    if {[llength $newc]} {CatalogPanelAddColumnsFromTSV $data $newc}
    if {[llength $oldc]} {CatalogPanelAddColumnsFromTSV $data $oldc}
    if {$seq ne {} && $seq ne 0} {OGFSessSet $seq post [dict create kind ai cols_list $cols]}
    set summary [lindex [split [string trim $err] \n] end]
    ::ogf::cat::set status "AI Services: [llength $cols] column(s) added (service $service) - $summary"
    return 1
}

# ---------------------------------------------------------------- hook for existing local steps
# Called first thing in the local Photo-z / SED / morphology / star-classification entry points.
# Returns 1 when an external service handled the step (the local code must then return).
proc OGFAIBackendHook {step} {
    global ogfai OGFAI_BACKEND_STEPS
    OGFAIInit
    if {![info exists ogfai(backend,$step)]} {return 0}
    set svc $ogfai(backend,$step)
    if {$svc eq {local} || $svc eq {}} {return 0}
    lassign [dict get $OGFAI_BACKEND_STEPS $step] task renames label
    OGFAIRunTask $task $svc -rename $renames -title "$label via external service $svc ($task)"
    return 1
}

# ------------------------------------------------------------------ Service Registry dialog
proc OGFAIRegistry {} {
    global ogfai
    OGFAIInit
    set w .ogfai_reg
    if {[winfo exists $w]} {raise $w; OGFAIRegRefresh; return}
    toplevel $w
    wm title $w "AI Services - Service Registry"
    ttk::frame $w.top
    ttk::label $w.top.l -text "Profile file:"
    ttk::entry $w.top.e -textvariable ogfai(sf_edit) -width 60
    ttk::button $w.top.b -text "Browse..." -command OGFAIRegBrowse
    pack $w.top.l -side left -padx 4
    pack $w.top.e -side left -fill x -expand true
    pack $w.top.b -side left -padx 4
    bind $w.top.e <Return> OGFAIRegApplyPath
    ttk::treeview $w.tv -columns {task transport enabled env envstatus target note} -show {tree headings} -height 10 -selectmode browse
    $w.tv heading #0 -text Name
    $w.tv column #0 -width 150 -stretch 0
    foreach {c t wd} {task Task 120 transport Transport 100 enabled Enabled 60 env {Auth env var} 130 envstatus {Env} 50 target Target 260 note Note 130} {
	$w.tv heading $c -text $t
	$w.tv column $c -width $wd -stretch [expr {$c eq "target"}]
    }
    ttk::label $w.msg -text {} -wraplength 760 -justify left
    ttk::frame $w.btn
    ttk::button $w.btn.en -text "Enable / Disable" -command OGFAIRegToggle
    ttk::button $w.btn.test -text "Test Connection" -command OGFAIRegTest
    ttk::button $w.btn.mk -text "Create from Example" -command OGFAIRegCreate
    ttk::button $w.btn.close -text Close -command [list destroy $w]
    pack $w.btn.en $w.btn.test $w.btn.mk -side left -padx 4
    pack $w.btn.close -side right -padx 4
    ttk::button $w.btn.ag -text "Detect Agent CLIs" -command OGFAIRegDetect
    pack $w.btn.ag -side left -padx 4
    ttk::button $w.btn.set -text "Profile / Backends..." -command {OGFParamDialog ai_services}
    pack $w.btn.set -side left -padx 4
    ttk::label $w.sec -text "Keys are never stored or shown: a profile names an environment variable; this window shows only set / unset.\nAgent CLIs (agent_*) use their own login; 'Test Connection' on them only runs <cli> --version." -foreground gray30 -justify left
    pack $w.top -fill x -pady 4
    pack $w.tv -fill both -expand true -padx 4
    pack $w.msg -fill x -padx 6 -pady 2
    pack $w.btn -fill x -pady 4
    pack $w.sec -padx 6 -pady 4 -anchor w
    OGFAIRegRefresh
}

proc OGFAIRegRefresh {} {
    global ogfai
    set w .ogfai_reg
    if {![winfo exists $w]} return
    lassign [OGFAIServices] rows path err
    set ogfai(sf_edit) $path
    $w.tv delete [$w.tv children {}]
    set ok {}
    foreach r $rows {
	set note [dict get $r note]
	if {[dict get $r valid] ne "yes"} {set note "INVALID $note"}
	$w.tv insert {} end -id [dict get $r name] -text [dict get $r name] -values [list \
	    [dict get $r task] [dict get $r transport] [dict get $r enabled] [dict get $r env] [dict get $r envstatus] \
	    [dict get $r target] $note]
    }
    set msg "$path"
    if {![file exists $path]} {
	append msg "  (file does not exist - only the built-in MOCK service is available)"
	pack $w.btn.mk -side left -padx 4
    } else {
	pack forget $w.btn.mk
    }
    if {$err ne {}} {append msg "\nERROR: $err"}
    $w.msg configure -text $msg
}

# choices of a "Backend of the local steps" setting: local + the enabled, valid services of the matching task
# (plugin.json ai_services "choices_proc")
proc OGFAIBackendChoices {step} {
    global OGFAI_BACKEND_STEPS
    OGFAIInit
    set task [lindex [dict get $OGFAI_BACKEND_STEPS $step] 0]
    lassign [OGFAIServices] rows path err
    set vals local
    foreach r $rows {
	if {[dict get $r enabled] eq "yes" && [dict get $r valid] eq "yes" && [dict get $r task] in [list $task any]} {lappend vals [dict get $r name]}
    }
    return $vals
}

# installed / not found for the four agent CLIs (no model is contacted; only `--version`)
proc OGFAIRegDetect {} {
    set w .ogfai_reg
    lassign [OGFAIExec [concat [OGFAIArgv detect-agents] --json]] rc out err
    if {$rc || [catch {set doc [::ogf::json::parse $out]}]} {
	$w.msg configure -text "detection failed: [string trim $err]"
	return
    }
    set lines {}
    foreach r [::ogf::json::get $doc agent_clis] {
	lappend lines [format "%-16s %s" [::ogf::json::get $r label] [expr {[::ogf::json::get $r installed 0] ? "installed: [::ogf::json::get $r path]  ([::ogf::json::get $r version])" : "NOT FOUND on PATH ([join [::ogf::json::get $r tried] {, }])"}]]
    }
    $w.msg configure -text [join $lines "\n"]
    OGFAIAppendLog "Detect agent CLIs\n[join $lines \n]\n"
}

proc OGFAIRegApplyPath {} {
    global ogfai
    set p [string trim $ogfai(sf_edit)]
    set ogfai(sf) $p
    OGFAIPrefSave
    OGFAIRegRefresh
}

proc OGFAIRegBrowse {} {
    global ogfai
    set f [tk_getOpenFile -title "AI services profile (JSON)" -filetypes {{{JSON} {.json}} {{All files} *}} \
	-initialdir [file dirname [expr {$ogfai(sf_edit) ne {} ? $ogfai(sf_edit) : [OGFSessWorkDir]}]]]
    if {$f eq {}} return
    set ogfai(sf_edit) $f
    OGFAIRegApplyPath
}

proc OGFAIRegCreate {} {
    global ogfai
    set ex [file join [OGFSessRoot] ai_services.example.json]
    set dest [string trim $ogfai(sf_edit)]
    if {$dest eq {}} {set dest [file join [OGFSessWorkDir] ai_services.json]}
    if {![file exists $ex]} {.ogfai_reg.msg configure -text "example not found: $ex"; return}
    if {[file exists $dest]} return
    file mkdir [file dirname $dest]
    file copy $ex $dest
    OGFAIRegRefresh
    .ogfai_reg.msg configure -text "Created $dest from the example (all template services are disabled and use placeholder URLs; edit it with a text editor)."
}

proc OGFAIRegSelected {} {
    set w .ogfai_reg
    set s [$w.tv selection]
    if {$s eq {}} {$w.msg configure -text "Select a service first."; return {}}
    return [lindex $s 0]
}

proc OGFAIRegToggle {} {
    set n [OGFAIRegSelected]
    if {$n eq {}} return
    set w .ogfai_reg
    set cur [lindex [$w.tv item $n -values] 2]
    lassign [OGFAIExec [concat [OGFAIArgv set-enabled] [list --service $n --enabled [expr {$cur eq "yes" ? "no" : "yes"}]]]] rc out err
    OGFAIRegRefresh
    if {$rc} {$w.msg configure -text "[string trim $err]"} else {$w.msg configure -text [string trim $out]}
    catch {$w.tv selection set $n}
}

proc OGFAIRegTest {} {
    set n [OGFAIRegSelected]
    if {$n eq {}} return
    set w .ogfai_reg
    set info [OGFAIServiceInfo $n]
    # agent CLI: only checks that the executable exists and prints its version; no prompt, no data, no cloud call
    if {[OGFAIIsNetwork $info] && ![OGFAIConfirm $info "A test request (HTTP GET of the service's health URL; no catalog data) will be sent."]} {
	$w.msg configure -text "Test cancelled."
	return
    }
    set argv [concat [OGFAIArgv check-service] [list --service $n]]
    if {[OGFAIIsNetwork $info]} {lappend argv --allow-network}
    $w.msg configure -text "Testing $n ..."
    update idletasks
    lassign [OGFAIExec $argv] rc out err
    set txt [string trim "$out\n$err"]
    $w.msg configure -text "[expr {$rc ? {NOT OK} : {OK}}]: $txt"
    OGFAIAppendLog "Test connection $n\n\$ [join $argv { }]\n$txt\n"
    catch {$w.tv selection set $n}
}

# ------------------------------------------------------------------ Run Task dialog
proc OGFAIRunDialog {} {
    global ogfai OGFAI_TASKS
    OGFAIInit
    if {![::ogf::cat::has]} {
	::ogf::cat::set status "AI Services: extract sources first"
	return
    }
    set w .ogfai_run
    catch {destroy $w}
    toplevel $w
    wm title $w "AI Services - Run Task on Catalog"
    wm transient $w .
    lassign [OGFAIServices] rows path err
    set ogfai(run,rows) $rows
    set ogfai(run,task) photoz
    set ogfai(run,service) {}
    set ogfai(run,agents) [OGFAIAgentInfo]
    set ogfai(run,backend) $ogfai(last,backend)
    set ogfai(run,bshow) {}
    set ogfai(run,images) 0
    set ogfai(run,exe) {}
    set ogfai(run,xargs) {}
    set ogfai(run,which) all
    set ogfai(run,dry) 0
    set ogfai(run,size) $ogfai(last,size)
    set ogfai(run,unit) $ogfai(last,unit)
    set ogfai(run,norm) $ogfai(last,norm)
    set ogfai(run,fmt) $ogfai(last,fmt)
    set ogfai(run,done) {}
    set nsel [llength [::ogf::cat::selection]]
    set nall [OGFAIRowCount]
    set r 0
    ttk::label $w.lt -text "Task"
    ttk::combobox $w.ct -textvariable ogfai(run,task) -state readonly -values $OGFAI_TASKS -width 34
    ttk::label $w.lb -text "Backend"
    ttk::combobox $w.cb -textvariable ogfai(run,bshow) -state readonly -width 44 -values [OGFAIBackendLabels]
    ttk::label $w.ls -text "Service"
    ttk::combobox $w.cs -textvariable ogfai(run,service) -state readonly -width 34
    ttk::labelframe $w.ag -text "Agent CLI (runs on your machine, uses its own login; data goes to the provider's cloud model)"
    ttk::label $w.ag.l1 -text "Executable"
    ttk::entry $w.ag.e1 -textvariable ogfai(run,exe) -width 44
    ttk::button $w.ag.b1 -text "Browse..." -command OGFAIRunBrowseExe
    ttk::label $w.ag.l2 -text "Extra args (JSON list)"
    ttk::entry $w.ag.e2 -textvariable ogfai(run,xargs) -width 44
    ttk::button $w.ag.b2 -text "Save to profile" -command OGFAIRunSaveAgent
    ttk::checkbutton $w.ag.img -text "Send image cutouts to the CLI (default OFF: only catalog numbers leave the machine)" -variable ogfai(run,images)
    ttk::label $w.ag.st -text {} -wraplength 480 -justify left
    grid $w.ag.l1 $w.ag.e1 $w.ag.b1 -padx 4 -pady 2 -sticky w
    grid $w.ag.l2 $w.ag.e2 $w.ag.b2 -padx 4 -pady 2 -sticky w
    grid $w.ag.img -columnspan 3 -padx 4 -pady 2 -sticky w
    grid $w.ag.st -columnspan 3 -padx 4 -pady 2 -sticky w
    ttk::label $w.lr -text "Rows"
    ttk::frame $w.fr
    ttk::radiobutton $w.fr.a -text "all ($nall)" -variable ogfai(run,which) -value all
    ttk::radiobutton $w.fr.s -text "selected ($nsel)" -variable ogfai(run,which) -value selected
    pack $w.fr.a $w.fr.s -side left -padx 4
    ttk::labelframe $w.cut -text "Cutouts (used when the service takes images)"
    ttk::label $w.cut.l1 -text "Size"
    ttk::entry $w.cut.e1 -textvariable ogfai(run,size) -width 7
    ttk::combobox $w.cut.u -textvariable ogfai(run,unit) -values {pix arcsec} -state readonly -width 7
    ttk::label $w.cut.l2 -text "Normalise"
    ttk::combobox $w.cut.n -textvariable ogfai(run,norm) -values {(profile) asinh zscale linear none} -state readonly -width 9
    ttk::label $w.cut.l3 -text "Format"
    ttk::combobox $w.cut.f -textvariable ogfai(run,fmt) -values {(profile) png fits npy} -state readonly -width 9
    grid $w.cut.l1 $w.cut.e1 $w.cut.u -padx 4 -pady 2 -sticky w
    grid $w.cut.l2 $w.cut.n -padx 4 -pady 2 -sticky w
    grid $w.cut.l3 $w.cut.f -padx 4 -pady 2 -sticky w
    ttk::checkbutton $w.dry -text "Dry run (print the request, send nothing)" -variable ogfai(run,dry)
    ttk::label $w.note -text {} -wraplength 380 -justify left -foreground gray30
    ttk::frame $w.bb
    ttk::button $w.bb.ok -text Run -command {set ogfai(run,done) ok}
    ttk::button $w.bb.cancel -text Cancel -command {set ogfai(run,done) cancel}
    pack $w.bb.ok $w.bb.cancel -side left -padx 4
    grid $w.lt $w.ct -padx 6 -pady 3 -sticky w
    grid $w.lb $w.cb -padx 6 -pady 3 -sticky w
    grid $w.ls $w.cs -padx 6 -pady 3 -sticky w
    grid $w.ag -columnspan 2 -padx 6 -pady 4 -sticky we
    grid $w.lr $w.fr -padx 6 -pady 3 -sticky w
    grid $w.cut -columnspan 2 -padx 6 -pady 4 -sticky we
    grid $w.dry -columnspan 2 -padx 6 -pady 2 -sticky w
    grid $w.note -columnspan 2 -padx 6 -pady 2 -sticky w
    grid $w.bb -columnspan 2 -pady 6
    bind $w.ct <<ComboboxSelected>> OGFAIRunTaskChanged
    bind $w.cb <<ComboboxSelected>> OGFAIRunBackendPicked
    bind $w.cs <<ComboboxSelected>> OGFAIRunServiceChanged
    set ogfai(run,bshow) [OGFAIBackendLabel $ogfai(run,backend)]
    OGFAIRunTaskChanged
    bind $w <Escape> {set ogfai(run,done) cancel}
    wm protocol $w WM_DELETE_WINDOW {set ogfai(run,done) cancel}
    catch {grab $w}
    vwait ogfai(run,done)
    catch {grab release $w}
    set go [expr {$ogfai(run,done) eq "ok"}]
    set task $ogfai(run,task); set svc $ogfai(run,service)
    set which $ogfai(run,which); set dry $ogfai(run,dry)
    foreach k {size unit norm fmt} {set ogfai(last,$k) $ogfai(run,$k)}
    set ogfai(last,backend) $ogfai(run,backend)
    set imgs $ogfai(run,images)
    OGFAIPrefSave
    destroy $w
    if {!$go} return
    if {$svc eq {}} {::ogf::cat::set status "AI Services: no enabled service for task $task (see Service Registry)"; return}
    OGFAIRunTask $task $svc -rows $which -dry $dry -size [string trim $ogfai(last,size)] -unit $ogfai(last,unit) \
	-norm $ogfai(last,norm) -fmt $ogfai(last,fmt) -images $imgs
}

# Backend dropdown: the four agent CLIs (with installed / not found) + the generic services of the registry
set ::OGFAI_BACKENDS {codex "Codex CLI (codex exec)" claude "Claude Code (claude -p)" agy "agy / Gemini CLI (agy, gemini -p)" grok "Grok (xAI CLI)"}

proc OGFAIBackendInstalled {be rows} {
    # rows: dict name -> info (OGFAIAgentInfo)
    set found 0
    dict for {n i} $rows {
	if {[dict get $i backend] eq $be && [dict get $i installed]} {set found 1}
    }
    return $found
}

proc OGFAIBackendLabel {be} {
    global ogfai OGFAI_BACKENDS
    if {$be eq "generic"} {return "Generic services (REST / local command / python / mock)"}
    set rows [expr {[info exists ogfai(run,agents)] ? $ogfai(run,agents) : {}}]
    return "[dict get $OGFAI_BACKENDS $be]  -  [expr {[OGFAIBackendInstalled $be $rows] ? {installed} : {NOT FOUND on PATH}}]"
}

proc OGFAIBackendLabels {} {
    global OGFAI_BACKENDS
    set l {}
    foreach be [dict keys $OGFAI_BACKENDS] {lappend l [OGFAIBackendLabel $be]}
    lappend l [OGFAIBackendLabel generic]
    return $l
}

proc OGFAIRunBackendPicked {} {
    global ogfai OGFAI_BACKENDS
    set sel $ogfai(run,bshow)
    set ogfai(run,backend) generic
    foreach be [dict keys $OGFAI_BACKENDS] {if {$sel eq [OGFAIBackendLabel $be]} {set ogfai(run,backend) $be}}
    set ogfai(run,service) {}
    OGFAIRunTaskChanged
}

proc OGFAIRunBrowseExe {} {
    global ogfai
    set f [tk_getOpenFile -title "Agent CLI executable" -parent .ogfai_run]
    if {$f ne {}} {set ogfai(run,exe) $f}
}

# write executable / extra args of the chosen agent service into the profile file (bridge --mode set-agent)
proc OGFAIRunSaveAgent {} {
    global ogfai
    set w .ogfai_run
    set svc $ogfai(run,service)
    if {$svc eq {} || ![dict exists $ogfai(run,agents) $svc]} return
    set xa [string trim $ogfai(run,xargs)]
    if {$xa eq {}} {set xa "\[\]"}
    lassign [OGFAIExec [concat [OGFAIArgv set-agent] [list --service $svc --executable [string trim $ogfai(run,exe)] --extra-args-json $xa]]] rc out err
    if {$rc} {
	$w.ag.st configure -text "NOT saved: [string trim $err]"
	return
    }
    set ogfai(run,agents) [OGFAIAgentInfo]
    $w.cb configure -values [OGFAIBackendLabels]
    set ogfai(run,bshow) [OGFAIBackendLabel $ogfai(run,backend)]
    OGFAIRunServiceChanged
    $w.ag.st configure -text "Saved: [string trim $out]"
}

# fill the agent frame for the chosen service
proc OGFAIRunServiceChanged {} {
    global ogfai
    set w .ogfai_run
    if {![winfo exists $w]} return
    set svc $ogfai(run,service)
    if {$ogfai(run,backend) eq "generic" || $svc eq {} || ![dict exists $ogfai(run,agents) $svc]} return
    set i [dict get $ogfai(run,agents) $svc]
    set ogfai(run,exe) [dict get $i executable]
    set ogfai(run,xargs) [expr {[llength [dict get $i extra_args]] ? "\[[join [lmap a [dict get $i extra_args] {format {"%s"} $a}] ,]\]" : {}}]
    if {[dict get $i installed]} {
	$w.ag.st configure -text "Found: [dict get $i path]    Login: [dict get $i login]"
    } else {
	$w.ag.st configure -text "NOT FOUND on PATH. Install the CLI, or enter its full path above and press 'Save to profile'.    Login: [dict get $i login]"
    }
}

proc OGFAIRunTaskChanged {} {
    global ogfai
    set w .ogfai_run
    set vals {}
    set be $ogfai(run,backend)
    foreach r $ogfai(run,rows) {
	if {[dict get $r enabled] ne "yes" || [dict get $r valid] ne "yes"} continue
	set isag [expr {[dict get $r transport] eq "agent_cli"}]
	if {$be eq "generic"} {
	    if {$isag} continue
	    if {[dict get $r task] in [list $ogfai(run,task) any]} {lappend vals [dict get $r name]}
	} else {
	    if {!$isag || ![dict exists $ogfai(run,agents) [dict get $r name]]} continue
	    if {[dict get [dict get $ogfai(run,agents) [dict get $r name]] backend] ne $be} continue
	    if {[dict get $r task] in [list $ogfai(run,task) any]} {lappend vals [dict get $r name]}
	}
    }
    $w.cs configure -values $vals
    if {$ogfai(run,service) ni $vals} {set ogfai(run,service) [lindex $vals 0]}
    if {$be eq "generic"} {
	grid remove $w.ag
    } else {
	grid $w.ag
	OGFAIRunServiceChanged
    }
    if {[llength $vals]} {
	if {$be eq "generic"} {
	    $w.note configure -text "Result columns are added to the catalog table; mock results are tagged SERVICE=mock."
	} else {
	    $w.note configure -text "Result columns are added to the catalog table. Before the first run you see exactly what leaves the machine; 'Dry run' shows the prompt and the command."
	}
	$w.bb.ok state !disabled
    } else {
	$w.note configure -text "No enabled, valid service for this backend / task. Open Service Registry... to enable one or create the profile file."
	$w.bb.ok state disabled
    }
}

# ------------------------------------------------------------------ last run log
proc OGFAIShowLog {} {
    global ogfai
    OGFAIInit
    set w .ogfai_log
    catch {destroy $w}
    toplevel $w
    wm title $w "AI Services - Last Run Log"
    text $w.t -width 110 -height 28 -font TkFixedFont -wrap none -yscrollcommand [list $w.sb set]
    ttk::scrollbar $w.sb -command [list $w.t yview]
    $w.t insert end $ogfai(lastlog)
    $w.t configure -state disabled
    ttk::button $w.close -text Close -command [list destroy $w]
    pack $w.close -side bottom -pady 4
    pack $w.sb -side right -fill y
    pack $w.t -fill both -expand true
}
