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
	    if {$k eq "services_file"} {set ogfai(sf) $v} elseif {[string match backend,* $k]} {set ogfai($k) $v}
	}
    }
}

proc OGFAIPrefFile {} {return [file join [OGFSessWorkDir] ai_services.prf]}

proc OGFAIPrefSave {} {
    global ogfai
    catch {file mkdir [OGFSessWorkDir]}
    if {[catch {set fd [open [OGFAIPrefFile] w]}]} return
    puts $fd "services_file=$ogfai(sf)"
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
    global catpanel
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {return 0}
    set n 0
    foreach l [lrange [split $catpanel(alldata) \n] 1 end] {if {[string trim $l] ne {}} {incr n}}
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
    global catpanel ogfai ogfsess
    OGFAIInit
    set seq 0
    array set o {-rows all -dry 0 -size {} -unit pix -norm {} -fmt {} -rename {} -title {} -params {} -numbers {}}
    array set o $args
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "AI Services: no catalog - extract sources first"
	return 0
    }
    if {$service eq {}} {set catpanel(status) "AI Services: no service chosen"; return 0}
    set script [OGFAIScript]
    if {![file exists $script]} {set catpanel(status) "AI Services: ds9_ai_bridge.py not found"; return 0}
    set info [OGFAIServiceInfo $service]
    if {$info eq {}} {set catpanel(status) "AI Services: service '$service' not found in the registry"; return 0}
    set numbers $o(-numbers)
    if {$o(-rows) eq "selected"} {
	set numbers $catpanel(sel,nums)
	if {[llength $numbers] == 0} {set catpanel(status) "AI Services: no rows selected in the table"; return 0}
    }
    set nrows [expr {[llength $numbers] ? [llength $numbers] : [OGFAIRowCount]}]
    set net [OGFAIIsNetwork $info]
    lassign [OGFAIImages] bands files
    set catfile [CatalogPanelSaveTempCatalog ai]
    if {$catfile eq {}} {set catpanel(status) "AI Services: cannot write the temporary catalog"; return 0}
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
    if {!$dry} {
	lappend argv --provenance-output [file join [OGFSessWorkDir] ai_last_run.provenance.json]
	# the user confirmed the transfer; the recorded argv keeps the flag so replay mode can run it.
	# Pipeline mode still refuses network steps unless the script itself gets --allow-network.
	if {$net} {lappend argv --allow-network}
    }
    set what "$nrows object(s) of the catalog (positions, magnitudes, catalog row text[expr {[llength $files] ? {, image cutouts} : {}}] as the profile's request template specifies) will be sent."
    if {!$dry && ![OGFAIConfirm $info $what]} {
	set catpanel(status) "AI Services: cancelled"
	return 0
    }
    set title [expr {$o(-title) ne {} ? $o(-title) : "AI service $service / $task ($nrows objects)"}]
    set catpanel(status) "AI Services: $service / $task on $nrows object(s) ..."
    update idletasks
    if {!$dry} {
	# session recorder: whole catalog = AUTO (replayed by the pipeline), a row subset = MANUAL (replay only)
	foreach im $files {if {[lsearch -exact $ogfsess(images) $im] < 0} {lappend ogfsess(images) $im}}
	set class [expr {[llength $numbers] ? "manual" : "auto"}]
	set seq 0
	catch {set seq [OGFSessLog ai.run $class $argv -title $title -tool python -requires catalog \
	    -network [expr {$net ? 1 : 0}] \
	    -payload [dict create service $service task $task profile_sha256 [dict get $info sha] \
		transport [dict get $info transport] objects $nrows] \
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
	set catpanel(status) "AI Services: dry run done (nothing was sent) - see Show Last Run Log"
	OGFAIShowLog
	return [expr {!$rc}]
    }
    if {$rc} {
	set last [lindex [split [string trim $err] \n] end]
	set catpanel(status) "AI Services: failed - $last"
	return 0
    }
    set hdr [split [lindex [split $data \n] 0] \t]
    set cols [lrange $hdr 1 end]
    if {[lindex $hdr 0] ne "NUMBER" || [llength $cols] == 0} {
	set catpanel(status) "AI Services: service returned no columns"
	return 0
    }
    # CatalogPanelAddColumnsFromTSV updates in place when ANY requested column already exists (and then silently
    # skips the new ones), so new and existing columns are merged in two calls (new first).  The recorded post
    # step "ai" repeats exactly this split in the exported script.
    set have [split [lindex [split $catpanel(alldata) \n] 0] \t]
    set newc {}; set oldc {}
    foreach c $cols {if {$c in $have} {lappend oldc $c} else {lappend newc $c}}
    if {[llength $newc]} {CatalogPanelAddColumnsFromTSV $data $newc}
    if {[llength $oldc]} {CatalogPanelAddColumnsFromTSV $data $oldc}
    if {$seq ne {} && $seq ne 0} {OGFSessSet $seq post [dict create kind ai cols_list $cols]}
    set summary [lindex [split [string trim $err] \n] end]
    set catpanel(status) "AI Services: [llength $cols] column(s) added (service $service) - $summary"
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
    ttk::treeview $w.tv -columns {task transport enabled env envstatus target note} -show {tree headings} -height 9 -selectmode browse
    $w.tv heading #0 -text Name
    $w.tv column #0 -width 150 -stretch 0
    foreach {c t wd} {task Task 120 transport Transport 100 enabled Enabled 60 env {Auth env var} 130 envstatus {Env} 50 target Target 260 note Note 90} {
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
    ttk::button $w.btn.set -text "Profile / Backends..." -command {OGFParamDialog ai_services}
    pack $w.btn.set -side left -padx 4
    ttk::label $w.sec -text "Keys are never stored or shown: a profile names an environment variable; this window shows only set / unset." -foreground gray30
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
    global ogfai catpanel OGFAI_TASKS
    OGFAIInit
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "AI Services: extract sources first"
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
    set ogfai(run,which) all
    set ogfai(run,dry) 0
    set ogfai(run,size) $ogfai(last,size)
    set ogfai(run,unit) $ogfai(last,unit)
    set ogfai(run,norm) $ogfai(last,norm)
    set ogfai(run,fmt) $ogfai(last,fmt)
    set ogfai(run,done) {}
    set nsel [llength $catpanel(sel,nums)]
    set nall [OGFAIRowCount]
    set r 0
    ttk::label $w.lt -text "Task"
    ttk::combobox $w.ct -textvariable ogfai(run,task) -state readonly -values $OGFAI_TASKS -width 34
    ttk::label $w.ls -text "Service"
    ttk::combobox $w.cs -textvariable ogfai(run,service) -state readonly -width 34
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
    grid $w.ls $w.cs -padx 6 -pady 3 -sticky w
    grid $w.lr $w.fr -padx 6 -pady 3 -sticky w
    grid $w.cut -columnspan 2 -padx 6 -pady 4 -sticky we
    grid $w.dry -columnspan 2 -padx 6 -pady 2 -sticky w
    grid $w.note -columnspan 2 -padx 6 -pady 2 -sticky w
    grid $w.bb -columnspan 2 -pady 6
    bind $w.ct <<ComboboxSelected>> OGFAIRunTaskChanged
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
    destroy $w
    if {!$go} return
    if {$svc eq {}} {set catpanel(status) "AI Services: no enabled service for task $task (see Service Registry)"; return}
    OGFAIRunTask $task $svc -rows $which -dry $dry -size [string trim $ogfai(last,size)] -unit $ogfai(last,unit) \
	-norm $ogfai(last,norm) -fmt $ogfai(last,fmt)
}

proc OGFAIRunTaskChanged {} {
    global ogfai
    set w .ogfai_run
    set vals {}
    foreach r $ogfai(run,rows) {
	if {[dict get $r enabled] eq "yes" && [dict get $r valid] eq "yes" && [dict get $r task] in [list $ogfai(run,task) any]} {
	    lappend vals [dict get $r name]
	}
    }
    $w.cs configure -values $vals
    if {$ogfai(run,service) ni $vals} {set ogfai(run,service) [lindex $vals 0]}
    if {[llength $vals]} {
	$w.note configure -text "Result columns are added to the catalog table; mock results are tagged SERVICE=mock."
	$w.bb.ok state !disabled
    } else {
	$w.note configure -text "No enabled, valid service for this task. Open Service Registry... to enable one or create the profile file."
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
