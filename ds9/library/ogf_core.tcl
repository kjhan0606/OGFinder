#  OGFinder core layer.
#
#  Services that plugins must use instead of touching each other's globals, the plugin registry
#  (OGFRegisterPlugin / plugin.json loader) and the step runner.  See docs/architecture.md and docs/plugins.md.
#
#    ::ogf::cat::*       catalog table access (columns, values, selection, add columns, temp file, load)
#    ::ogf::mask::*      mask service (paths, run ds9_mask.py through the recorder, overlay)
#    ::ogf::bands::*     band manager (names, files, frames, detection band)
#    ::ogf::session::*   session recorder registration (wraps OGFSessLog)
#    ::ogf::ai::*        AI service bridge (wraps OGFAIRunTask)
#    ::ogf::job::*       job runner: async argv with status, elapsed time, cancel, recorder integration
#    ::ogf::params::*    parameter store (plugin specs, values, JSON persistence, presets)
#    ::ogf::log          logging ring buffer (+ stderr when OGF_DEBUG is set)
#    ::ogf::reg::*       plugin registry

package provide DS9 1.0

namespace eval ::ogf {
    variable tabs {Detect Classify Measure Low-SB Time-domain Results}
    variable stages {detect classify measure}
    variable logbuf {}
    variable plugins [dict create]
    variable order {}
    variable steprun
}

# ---------------------------------------------------------------- logging
proc ::ogf::log {level msg} {
    variable logbuf
    set line "[clock format [clock seconds] -format %H:%M:%S] $level $msg"
    lappend logbuf $line
    if {[llength $logbuf] > 500} {set logbuf [lrange $logbuf end-499 end]}
    if {[info exists ::env(OGF_DEBUG)] || $level eq "ERROR"} {catch {puts stderr "OGF: $line"}}
}

proc ::ogf::status {msg} {
    catch {set ::catpanel(status) $msg}
    catch {update idletasks}
}

# ================================================================ catalog service
namespace eval ::ogf::cat {}

proc ::ogf::cat::tsv {} {
    if {![info exists ::catpanel(alldata)]} {return {}}
    return $::catpanel(alldata)
}
proc ::ogf::cat::has {} {return [expr {[tsv] ne {}}]}
proc ::ogf::cat::columns {} {
    set d [tsv]
    if {$d eq {}} {return {}}
    return [lmap h [split [lindex [split $d \n] 0] \t] {string trim $h}]
}
proc ::ogf::cat::nrows {} {
    set n 0
    foreach l [lrange [split [tsv] \n] 1 end] {if {[string trim $l] ne {}} {incr n}}
    return $n
}
# values of one column (all rows, catalog order)
proc ::ogf::cat::values {col} {
    set i [lsearch -exact [columns] $col]
    if {$i < 0} {return {}}
    set out {}
    foreach l [lrange [split [tsv] \n] 1 end] {
	if {[string trim $l] eq {}} continue
	lappend out [string trim [lindex [split $l \t] $i]]
    }
    return $out
}
# list of dicts, one per row
proc ::ogf::cat::rows {} {
    set cols [columns]
    set out {}
    foreach l [lrange [split [tsv] \n] 1 end] {
	if {[string trim $l] eq {}} continue
	set d {}
	foreach c $cols v [split $l \t] {dict set d $c [string trim $v]}
	lappend out $d
    }
    return $out
}
proc ::ogf::cat::selection {} {
    if {![info exists ::catpanel(sel,nums)]} {return {}}
    return $::catpanel(sel,nums)
}
# mode: replace | add | toggle
proc ::ogf::cat::select {nums {mode replace} {pan 1}} {
    set first 1
    foreach n $nums {
	CatalogPanelLinkSelect $n [expr {$first ? $mode : "add"}] [expr {$pan && $first}]
	set first 0
    }
}
proc ::ogf::cat::clear_selection {} {CatalogPanelClearSelection}
# result TSV (NUMBER + new columns) merged into the catalog; names = columns to take
proc ::ogf::cat::add_columns {result_tsv names} {CatalogPanelAddColumnsFromTSV $result_tsv $names}
proc ::ogf::cat::load_tsv {tsv name} {CatalogPanelLoadTSV $tsv $name}
proc ::ogf::cat::temp_file {suffix} {return [CatalogPanelSaveTempCatalog $suffix]}
proc ::ogf::cat::image_file {} {return [CatalogPanelGetFITS]}
proc ::ogf::cat::status {msg} {::ogf::status $msg}

# ================================================================ mask service
namespace eval ::ogf::mask {}
proc ::ogf::mask::paths {{fn {}}} {return [OGFMaskPaths $fn]}
proc ::ogf::mask::exists {} {return [OGFMaskExists]}
proc ::ogf::mask::bool_path {} {return [OGFMaskBoolPath]}
# run ds9_mask.py in MODE for image FN (default current frame); recorded by the session recorder
proc ::ogf::mask::run {mode {fn {}} {extra {}}} {return [OGFMaskRun $mode $fn $extra]}
proc ::ogf::mask::overlay {on} {set ::ogfmask(overlay) $on; CatalogPanelMaskToggleOverlay}

# ================================================================ band manager
namespace eval ::ogf::bands {}
proc ::ogf::bands::names {} {return [expr {[info exists ::ogfband(names)] ? $::ogfband(names) : {}}]}
proc ::ogf::bands::detect {} {return [expr {[info exists ::ogfband(detect)] ? $::ogfband(detect) : {}}]}
proc ::ogf::bands::path {name} {return [expr {[info exists ::ogfband($name,file)] ? $::ogfband($name,file) : {}}]}
proc ::ogf::bands::frames {} {return [OGFBandFrames 1]}
proc ::ogf::bands::sorted {} {return [OGFBandsSorted]}
proc ::ogf::bands::register {frame name zp {fwhm {}} {fn {}}} {CatalogPanelBandsRegisterFrame $frame $name $zp $fwhm $fn}
proc ::ogf::bands::set_detect {name} {CatalogPanelBandsSetDetect $name}
proc ::ogf::bands::pos {frame num x y} {return [OGFBandPos $frame $num $x $y]}

# ================================================================ session recorder
namespace eval ::ogf::session {}
# log STEP CLASS ARGV ?options of OGFSessLog? ; CLASS = auto|config|manual
proc ::ogf::session::log {step class argv args} {return [OGFSessLog $step $class $argv {*}$args]}
proc ::ogf::session::set_field {seq key val} {OGFSessSet $seq $key $val}
proc ::ogf::session::steps {} {return [expr {[info exists ::ogfsess(steps)] ? $::ogfsess(steps) : {}}]}
proc ::ogf::session::crc {s} {return [OGFSessCrc $s]}

# ================================================================ AI bridge
namespace eval ::ogf::ai {}
proc ::ogf::ai::run {task service args} {return [OGFAIRunTask $task $service {*}$args]}
proc ::ogf::ai::services {} {return [OGFAIServices]}

# ================================================================ job runner
namespace eval ::ogf::job {
    variable cur
    array set cur {busy 0 fd {} pid {} t0 0 title {} seq 0 step {} out {} errf {} done {} cancelled 0 plugin {}}
}
proc ::ogf::job::busy {} {variable cur; return $cur(busy)}
proc ::ogf::job::elapsed {} {
    variable cur
    return [expr {$cur(busy) ? ([clock milliseconds] - $cur(t0)) / 1000.0 : 0.0}]
}

# run ARGV asynchronously.  Options:
#   -step NAME -class auto|config|manual -title T -network 0|1 -outputs L -payload D -post D -requires L
#   -done SCRIPT   (called with  ok output ms)       -plugin ID
# The recorder entry (when -step is given) is written BEFORE the process starts, as the existing
# steps do; ms / failed / stdout_crc are filled in at the end.
proc ::ogf::job::run {argv args} {
    variable cur
    array set o {-step {} -class auto -title {} -network 0 -outputs {} -payload {} -post {} -requires {} -done {} -plugin {}}
    array set o $args
    if {$cur(busy)} {::ogf::status "another job is still running: $cur(title)"; return 0}
    set seq 0
    if {$o(-step) ne {}} {
	set opts [list -title $o(-title) -network $o(-network) -outputs $o(-outputs) -tool python]
	if {$o(-payload) ne {}} {lappend opts -payload $o(-payload)}
	if {$o(-post) ne {}} {lappend opts -post $o(-post)}
	if {$o(-requires) ne {}} {lappend opts -requires $o(-requires)}
	catch {set seq [OGFSessLog $o(-step) $o(-class) $argv {*}$opts]}
    }
    set errf [file join [OGFSessWorkDir] ogf_job_stderr.txt]
    catch {file mkdir [OGFSessWorkDir]}
    array set cur [list busy 1 t0 [clock milliseconds] title $o(-title) seq $seq step $o(-step) out {} errf $errf \
	done $o(-done) cancelled 0 plugin $o(-plugin) pid {}]
    ::ogf::status "$o(-title) ..."
    ::ogf::log INFO "job start: [join $argv { }]"
    catch {OGFPrepareLibPath}
    if {[catch {set fd [open |[concat $argv [list 2>$errf]] r]} err]} {
	set cur(busy) 0
	::ogf::status "cannot start: [string range $err 0 100]"
	::ogf::log ERROR "cannot start job: $err"
	return 0
    }
    set cur(fd) $fd
    set cur(pid) [pid $fd]
    fconfigure $fd -blocking 0 -buffering line -translation binary -encoding utf-8
    fileevent $fd readable ::ogf::job::_readable
    catch {OGFUIJobState 1}
    return 1
}

proc ::ogf::job::cancel {} {
    variable cur
    if {!$cur(busy)} return
    set cur(cancelled) 1
    catch {OGFsess_exec kill $cur(pid)}
    ::ogf::status "cancelling $cur(title) ..."
}

proc ::ogf::job::_readable {} {
    variable cur
    set fd $cur(fd)
    if {[catch {set chunk [read $fd]}]} {set chunk {}}
    append cur(out) $chunk
    if {![eof $fd]} return
    fconfigure $fd -blocking 1
    set rc [catch {close $fd} err]
    set cur(busy) 0
    set ms [expr {[clock milliseconds] - $cur(t0)}]
    set out $cur(out)
    # exec strips one trailing newline; the recorder and the exported script's checksums assume that
    if {[string index $out end] eq "\n"} {set out [string range $out 0 end-1]}
    set stderr_txt {}
    catch {set f [open $cur(errf) r]; set stderr_txt [read $f]; close $f}
    if {$cur(seq) > 0} {
	catch {OGFSessSet $cur(seq) ms $ms}
	if {$rc || $cur(cancelled)} {
	    catch {OGFSessSet $cur(seq) failed 1}
	    if {$cur(cancelled)} {catch {OGFSessSet $cur(seq) note "cancelled by the user"}}
	} else {
	    catch {OGFSessSet $cur(seq) stdout_crc [OGFSessCrc $out]}
	    catch {OGFSessSet $cur(seq) stdout_len [string length $out]}
	}
    }
    set ok [expr {!$rc && !$cur(cancelled)}]
    if {$cur(cancelled)} {
	::ogf::status "$cur(title) cancelled"
    } elseif {$rc} {
	set last [lindex [split [string trim $stderr_txt] \n] end]
	::ogf::status "$cur(title) failed: [string range [expr {$last ne {} ? $last : $err}] 0 140]"
	::ogf::log ERROR "job failed: $cur(title): $err"
    } else {
	::ogf::status "$cur(title) done ([format %.1f [expr {$ms/1000.0}]] s)"
    }
    set done $cur(done)
    set plugin $cur(plugin)
    catch {OGFUIJobState 0}
    if {$done ne {}} {catch {{*}$done $ok $out $ms} derr; if {$derr ne {}} {::ogf::log ERROR "done handler: $derr"}}
    catch {OGFUIProgressRefresh}
}

# ================================================================ parameter store
# Specs:  list of dicts  name label type default min max choices help group expert  (+ unit)
# Values live in  ::ogfparam(PLUGIN,NAME)  unless the plugin declares a "store" binding
#   "store": {"array": "catpanel", "key": "param,%s", "save": "CatalogPanelParamSave"}   (legacy variables)
namespace eval ::ogf::params {}

proc ::ogf::params::specs {plugin {step {}}} {
    set m [::ogf::reg::get $plugin]
    set l [::ogf::json::get $m params]
    if {$step eq {}} {return $l}
    return [lmap s $l {if {[::ogf::json::get $s step] ne {} && [::ogf::json::get $s step] ne $step} continue; set s}]
}
proc ::ogf::params::spec {plugin name} {
    foreach s [specs $plugin] {if {[dict get $s name] eq $name} {return $s}}
    return {}
}
proc ::ogf::params::var {plugin name} {
    set m [::ogf::reg::get $plugin]
    set st [::ogf::json::get $m store]
    if {$st ne {}} {
	set a [dict get $st array]
	if {![string match ::* $a]} {set a ::$a}
	return [list $a [format [dict get $st key] $name]]
    }
    return [list ::ogfparam $plugin,$name]
}
proc ::ogf::params::get {plugin name} {
    lassign [var $plugin $name] a k
    if {[info exists ${a}($k)]} {return [set ${a}($k)]}
    set s [spec $plugin $name]
    return [::ogf::json::get $s default]
}
proc ::ogf::params::put {plugin name val} {
    lassign [var $plugin $name] a k
    set ${a}($k) $val
}
proc ::ogf::params::defaults {plugin} {
    set d [dict create]
    foreach s [specs $plugin] {dict set d [dict get $s name] [::ogf::json::get $s default]}
    return $d
}
# make sure every declared parameter has a value (first load); then read the persisted file
proc ::ogf::params::init {plugin} {
    foreach s [specs $plugin] {
	lassign [var $plugin [dict get $s name]] a k
	if {![info exists ${a}($k)]} {set ${a}($k) [::ogf::json::get $s default]}
    }
    restore $plugin
}
proc ::ogf::params::prefpath {plugin} {
    return [file join [file normalize ~] .ds9 ogf_params $plugin.json]
}
proc ::ogf::params::legacy_store {plugin} {
    return [expr {[::ogf::json::get [::ogf::reg::get $plugin] store] ne {}}]
}
# persistence: plugins with a legacy store keep their own .prf (save proc), others use JSON
proc ::ogf::params::save {plugin} {
    set m [::ogf::reg::get $plugin]
    set st [::ogf::json::get $m store]
    if {$st ne {}} {
	set p [::ogf::json::get $st save]
	if {$p ne {}} {{*}$p}
	return
    }
    set d [dict create]
    foreach s [specs $plugin] {dict set d [dict get $s name] [get $plugin [dict get $s name]]}
    set f [prefpath $plugin]
    catch {::file mkdir [::file dirname $f]}
    if {![catch {set fd [open $f w]}]} {
	fconfigure $fd -encoding utf-8
	puts -nonewline $fd [::ogf::json::write_flat $d]
	close $fd
    }
}
proc ::ogf::params::restore {plugin} {
    if {[legacy_store $plugin]} return
    set f [prefpath $plugin]
    if {![::file exists $f]} return
    if {[catch {set fd [open $f r]; fconfigure $fd -encoding utf-8; set txt [read $fd]; close $fd; set d [::ogf::json::parse $txt]}]} return
    dict for {k v} $d {if {[spec $plugin $k] ne {}} {put $plugin $k $v}}
}
proc ::ogf::params::preset_dir {plugin} {return [::file join [file normalize ~] .ds9 ogf_params $plugin.presets]}
proc ::ogf::params::preset_list {plugin} {
    return [lsort [lmap f [glob -nocomplain -directory [preset_dir $plugin] *.json] {::file rootname [::file tail $f]}]]
}
proc ::ogf::params::preset_save {plugin name values} {
    ::file mkdir [preset_dir $plugin]
    set fd [open [::file join [preset_dir $plugin] $name.json] w]
    fconfigure $fd -encoding utf-8
    puts -nonewline $fd [::ogf::json::write_flat $values]
    close $fd
}
proc ::ogf::params::preset_load {plugin name} {
    set fd [open [::file join [preset_dir $plugin] $name.json] r]
    fconfigure $fd -encoding utf-8
    set txt [read $fd]; close $fd
    return [::ogf::json::parse $txt]
}
# check one value against its spec; returns {} when fine, else an error text
proc ::ogf::params::validate {s val} {
    set t [::ogf::json::get $s type string]
    set lab [::ogf::json::get $s label [dict get $s name]]
    switch -- $t {
	int {
	    if {![string is integer -strict $val]} {return "$lab: integer expected"}
	}
	float {
	    if {![string is double -strict $val]} {return "$lab: number expected"}
	}
	choice {
	    if {$val ni [::ogf::json::get $s choices]} {return "$lab: must be one of [join [::ogf::json::get $s choices] {, }]"}
	}
    }
    if {$t in {int float}} {
	set mn [::ogf::json::get $s min]; set mx [::ogf::json::get $s max]
	if {$mn ne {} && $val < $mn} {return "$lab: must be >= $mn"}
	if {$mx ne {} && $val > $mx} {return "$lab: must be <= $mx"}
    }
    return {}
}

# ================================================================ plugin registry
namespace eval ::ogf::reg {}

proc ::ogf::reg::get {id} {
    variable ::ogf::plugins
    if {![dict exists $::ogf::plugins $id]} {return {}}
    return [dict get $::ogf::plugins $id]
}
proc ::ogf::reg::ids {{only_enabled 1}} {
    set out {}
    foreach id $::ogf::order {
	if {$only_enabled && ![dict get [get $id] _enabled]} continue
	lappend out $id
    }
    return $out
}
proc ::ogf::reg::step {id sid} {
    foreach s [::ogf::json::get [get $id] steps] {if {[dict get $s id] eq $sid} {return $s}}
    return {}
}
# (plugin id, tab) -> ordered list of enabled plugin ids on a tab
proc ::ogf::reg::on_tab {tab} {
    set out {}
    foreach id [ids] {
	if {[::ogf::json::get [get $id] tab] eq $tab} {lappend out $id}
    }
    return $out
}

# user preference: ~/.ds9/ogf_plugins.json  {"enabled": [...], "disabled": [...]}  and env OGFINDER_PLUGINS="a,-b"
proc ::ogf::reg::prefs {} {
    set en {}; set dis {}
    set f [file join [file normalize ~] .ds9 ogf_plugins.json]
    if {[file exists $f]} {
	if {![catch {set fd [open $f r]; set txt [read $fd]; close $fd; set d [::ogf::json::parse $txt]}]} {
	    set en [::ogf::json::get $d enabled]; set dis [::ogf::json::get $d disabled]
	}
    }
    if {[info exists ::env(OGFINDER_PLUGINS)]} {
	foreach t [split $::env(OGFINDER_PLUGINS) ,] {
	    set t [string trim $t]
	    if {$t eq {}} continue
	    if {[string index $t 0] eq "-"} {lappend dis [string range $t 1 end]} else {lappend en $t}
	}
    }
    return [list $en $dis]
}
proc ::ogf::reg::save_prefs {en dis} {
    set f [file join [file normalize ~] .ds9 ogf_plugins.json]
    catch {file mkdir [file dirname $f]}
    set fd [open $f w]
    puts -nonewline $fd "\{\n  \"enabled\": [OGFJList $en],\n  \"disabled\": [OGFJList $dis]\n\}\n"
    close $fd
}

proc ::ogf::reg::enabled_default {m en dis} {
    set id [dict get $m id]
    if {[::ogf::json::get $m required 0]} {return 1}     ;# infrastructure plugins cannot be disabled
    if {$id in $dis} {return 0}
    if {$id in $en} {return 1}
    if {[::ogf::json::get $m example 0] || ![::ogf::json::get $m enabled 1]} {return 0}
    return 1
}

# validation: returns list of problems (empty = ok)
proc ::ogf::reg::validate {m} {
    set bad {}
    foreach k {id name tab} {if {![dict exists $m $k] || [dict get $m $k] eq {}} {lappend bad "missing '$k'"}}
    if {[dict exists $m tab] && [dict get $m tab] ni $::ogf::tabs} {
	lappend bad "unknown tab '[dict get $m tab]' (tabs: $::ogf::tabs)"
    }
    set seen {}
    foreach s [::ogf::json::get $m steps] {
	if {![dict exists $s id] || ![dict exists $s label]} {lappend bad "step without id/label"; continue}
	if {[dict get $s id] in $seen} {lappend bad "duplicate step id [dict get $s id]"}
	lappend seen [dict get $s id]
	set cls [string toupper [::ogf::json::get $s session AUTO]]
	if {$cls ni {AUTO CONFIG MANUAL NONE}} {lappend bad "step [dict get $s id]: session must be AUTO|CONFIG|MANUAL|NONE"}
	set n 0
	foreach k {cli proc} {if {[dict exists $s $k]} {incr n}}
	if {$n > 1} {lappend bad "step [dict get $s id]: give either 'cli' or 'proc'"}
    }
    set names {}
    foreach p [::ogf::json::get $m params] {
	if {![dict exists $p name]} {lappend bad "parameter without name"; continue}
	if {[dict get $p name] in $names} {lappend bad "duplicate parameter [dict get $p name]"}
	lappend names [dict get $p name]
	set t [::ogf::json::get $p type string]
	if {$t ni {int float string bool choice file}} {lappend bad "parameter [dict get $p name]: unknown type $t"}
    }
    return $bad
}

# THE registration entry point.  m = dict (the plugin.json content).  Options: -dir DIR (plugin directory).
# Returns 1 when registered.  Plugin Tcl files may call it directly with a dict they build.
proc ::ogf::reg::register {m args} {
    array set o {-dir {}}
    array set o $args
    set bad [validate $m]
    if {[llength $bad]} {
	::ogf::log ERROR "plugin [::ogf::json::get $m id ?] rejected: [join $bad {; }]"
	return 0
    }
    set id [dict get $m id]
    lassign [prefs] en dis
    dict set m _dir $o(-dir)
    dict set m _enabled [enabled_default $m $en $dis]
    dict set m _loaded 0
    if {![dict exists $::ogf::plugins $id]} {lappend ::ogf::order $id}
    dict set ::ogf::plugins $id $m
    ::ogf::log INFO "plugin registered: $id ([expr {[dict get $m _enabled] ? {enabled} : {disabled}}]) tab=[dict get $m tab]"
    return 1
}
proc OGFRegisterPlugin {m args} {return [::ogf::reg::register $m {*}$args]}

# directories scanned for plugin.json (later ones override earlier ones with the same id)
proc ::ogf::reg::search_dirs {} {
    set dirs {}
    catch {lappend dirs [file join $::ds9(root) plugins]}          ;# embedded in bin/ds9
    catch {lappend dirs [file join [OGFSessRoot] plugins]}         ;# next to bin/: <root>/plugins
    lappend dirs [file join [file normalize ~] .ds9 plugins]       ;# user plugins
    if {[info exists ::env(OGFINDER_PLUGIN_PATH)]} {
	foreach d [split $::env(OGFINDER_PLUGIN_PATH) :] {if {$d ne {}} {lappend dirs $d}}
    }
    return $dirs
}

proc ::ogf::reg::discover {} {
    set found [dict create]
    foreach d [search_dirs] {
	foreach f [lsort [glob -nocomplain -directory $d -types f */plugin.json]] {
	    if {[catch {
		set fd [open $f r]; fconfigure $fd -encoding utf-8
		set txt [read $fd]; close $fd
		set m [::ogf::json::parse $txt]
	    } err]} {
		::ogf::log ERROR "bad plugin.json $f: $err"
		continue
	    }
	    if {![dict exists $m id]} {::ogf::log ERROR "$f: no id"; continue}
	    dict set found [dict get $m id] [list $m [file dirname $f]]
	}
    }
    return $found
}

# read all manifests, register them (in the order: by "order" field, then id), source the Tcl files
proc ::ogf::reg::load_all {} {
    set found [discover]
    set items {}
    dict for {id v} $found {
	lassign $v m dir
	lappend items [list [::ogf::json::get $m order 100] $id $m $dir]
    }
    foreach it [lsort -integer -index 0 [lsort -index 1 $items]] {
	lassign $it ord id m dir
	register $m -dir $dir
    }
    foreach id [ids 0] {load_tcl $id}
}

proc ::ogf::reg::load_tcl {id} {
    set m [get $id]
    if {[dict get $m _loaded]} return
    # "tcl": "a.tcl" or ["a.tcl","b.tcl"].  "tcl_always": 1 = also source when the plugin is disabled: the file holds
    # procs that core code (layout.tcl, frame.tcl, the recorder, other plugins) calls by name; disabling the plugin then
    # only hides its workflow entries.
    set tfs [::ogf::json::get $m tcl]
    if {$tfs ne {} && ([dict get $m _enabled] || [::ogf::json::get $m tcl_always 0])} {
	foreach tf $tfs {
	    set path [file join [dict get $m _dir] $tf]
	    if {[catch {uplevel #0 [list source -encoding utf-8 $path]} err]} {
		::ogf::log ERROR "plugin $id: cannot source $tf: $err"
		catch {puts stderr "OGF: plugin $id: cannot source $tf: $err"}
		dict set ::ogf::plugins $id _enabled 0
		return
	    }
	}
    }
    dict set ::ogf::plugins $id _loaded 1
}

# run the plugin's "init" proc (once, when the panel exists) and load persisted parameters
proc ::ogf::reg::init_all {} {
    foreach id [ids] {
	set m [get $id]
	if {[::ogf::json::get $m params] ne {}} {::ogf::params::init $id}
	set p [::ogf::json::get $m init]
	if {$p ne {} && [info commands [lindex $p 0]] ne {}} {
	    if {[catch {uplevel #0 $p} err]} {::ogf::log ERROR "plugin $id init: $err"}
	}
    }
}

# requirements: {python modules, binaries}; returns list of missing items (python modules are checked lazily)
proc ::ogf::reg::missing {id {check_python 0}} {
    set m [get $id]
    set req [::ogf::json::get $m requires]
    set miss {}
    foreach b [::ogf::json::get $req binaries] {
	if {[string match */* $b]} {
	    if {![file executable $b]} {lappend miss "binary $b"}
	} elseif {[auto_execok $b] eq {} && ![file executable [file join [file dirname [info nameofexecutable]] $b]]} {
	    lappend miss "binary $b"
	}
    }
    if {$check_python} {
	foreach mod [::ogf::json::get $req python] {
	    if {[catch {OGFsess_exec [OGFPython] -c "import $mod"} e]} {lappend miss "python module $mod"}
	}
    }
    return $miss
}

# ================================================================ step execution
# argv template tokens:  {python} {plugin_dir} {work} {image} {catalog} {root}  and  {PARAM} (parameter store)
# a "cli" element is a string, or an object {"if": "param", "argv": [...]} (included when the parameter is
# true / non-empty / non-zero).
namespace eval ::ogf::step {}

# CLI drivers must be real files: prefer <root>/plugins/<id> on disk; a plugin that only exists inside the
# embedded zip is extracted once per run into ~/.ds9/plugin_cache/<id>.
proc ::ogf::step::plugin_dir {id} {
    set m [::ogf::reg::get $id]
    set d [dict get $m _dir]
    if {![string match {*zipfs*} $d] && ![string match {//*} $d]} {return $d}
    set disk [file join [OGFSessRoot] plugins $id]
    if {[file isdirectory $disk]} {return $disk}
    set cache [file join [OGFSessWorkDir] plugin_cache $id]
    file mkdir $cache
    foreach f [glob -nocomplain -directory $d *] {
	if {[file isfile $f]} {catch {file copy -force $f $cache}}
    }
    return $cache
}

proc ::ogf::step::context {id {catname {}}} {
    set m [::ogf::reg::get $id]
    set ctx [dict create python [OGFPython] plugin_dir [::ogf::step::plugin_dir $id] work [OGFSessWorkDir] root [OGFSessRoot]]
    dict set ctx image [CatalogPanelGetFITS]
    if {$catname ne {}} {
	dict set ctx catalog [::ogf::cat::temp_file $catname]
    }
    return $ctx
}

proc ::ogf::step::expand {id str ctx} {
    set out {}
    set i 0
    while {[regexp -start $i -indices {\{([A-Za-z0-9_.-]+)\}} $str m g]} {
	lassign $m a b; lassign $g c d
	append out [string range $str $i [expr {$a-1}]]
	set key [string range $str $c $d]
	if {[dict exists $ctx $key]} {
	    append out [dict get $ctx $key]
	} elseif {[::ogf::params::spec $id $key] ne {}} {
	    append out [::ogf::params::get $id $key]
	} else {
	    error "unknown token {$key} in argv template of plugin $id"
	}
	set i [expr {$b+1}]
    }
    append out [string range $str $i end]
    return $out
}

proc ::ogf::step::is_cond {e} {
    if {[catch {dict exists $e argv} r]} {return 0}
    return $r
}

proc ::ogf::step::build_argv {id step ctx} {
    set argv {}
    foreach e [::ogf::json::get $step cli] {
	if {[is_cond $e]} {
	    set cond [::ogf::json::get $e if]
	    set v [expr {[::ogf::params::spec $id $cond] ne {} ? [::ogf::params::get $id $cond] : 0}]
	    if {$v eq {} || $v eq "0" || [string is false -strict $v]} continue
	    foreach a [dict get $e argv] {lappend argv [expand $id $a $ctx]}
	} else {
	    lappend argv [expand $id $e $ctx]
	}
    }
    return $argv
}

# does the step need image / catalog?  returns an error text or {}
proc ::ogf::step::check_needs {step} {
    set needs [::ogf::json::get $step needs]
    foreach n $needs {
	switch -- $n {
	    image {if {[CatalogPanelGetFITS] eq {}} {return "load an image first"}}
	    catalog {if {![::ogf::cat::has]} {return "extract sources first (no catalog)"}}
	}
    }
    return {}
}

# run one step of plugin ID.  Tcl-proc steps run synchronously through the same procs the old menus called;
# cli steps go through the job runner and the session recorder.
proc ::ogf::step::run {id sid} {
    set m [::ogf::reg::get $id]
    set step [::ogf::reg::step $id $sid]
    if {$step eq {}} {::ogf::log ERROR "no step $sid in plugin $id"; return 0}
    set why [check_needs $step]
    if {$why ne {}} {::ogf::status "[dict get $step label]: $why"; return 0}
    set miss [::ogf::reg::missing $id]
    if {[llength $miss]} {::ogf::status "[dict get $step label]: missing [join $miss {, }]"; return 0}
    if {[dict exists $step proc]} {
	set cmd [dict get $step proc]
	::ogf::log INFO "step $id.$sid: $cmd"
	set n0 [llength [::ogf::session::steps]]
	set t0 [clock milliseconds]
	set rc [catch {uplevel #0 $cmd} res opts]
	catch {OGFUIProgressRefresh}
	if {$rc} {
	    ::ogf::log ERROR "step $id.$sid failed: $res"
	    ::ogf::status "[dict get $step label]: $res"
	    return 0
	}
	return 1
    }
    if {![dict exists $step cli]} {::ogf::log ERROR "step $id.$sid has neither cli nor proc"; return 0}
    set catname [::ogf::json::get $step catalog_tmp]
    if {$catname eq {} && [lsearch -glob [::ogf::json::get $step cli] *\{catalog\}*] >= 0} {set catname $sid}
    if {[catch {
	set ctx [context $id $catname]
	set argv [build_argv $id $step $ctx]
    } err]} {
	::ogf::status "[dict get $step label]: $err"
	::ogf::log ERROR "step $id.$sid: $err"
	return 0
    }
    set cls [string tolower [::ogf::json::get $step session AUTO]]
    set title [::ogf::json::get $step title [dict get $step label]]
    set out [::ogf::json::get $step output]
    set net [expr {[::ogf::json::get $step network [::ogf::json::get $m network 0]] ? 1 : 0}]
    set recname $id.$sid
    set opts [list -step $recname -class [expr {$cls eq "none" ? "auto" : $cls}] -title $title -network $net -plugin $id \
	-done [list ::ogf::step::done $id $sid]]
    if {$out ne {} && [::ogf::json::get $out mode] eq "add_columns"} {
	lappend opts -post [dict create kind add cols_list [::ogf::json::get $out columns]] -requires catalog
    } elseif {$out ne {} && [::ogf::json::get $out mode] eq "set"} {
	lappend opts -post [dict create kind set]
    }
    return [::ogf::job::run $argv {*}$opts]
}

# job finished: apply the declared output
proc ::ogf::step::done {id sid ok output ms} {
    if {!$ok} return
    set step [::ogf::reg::step $id $sid]
    set out [::ogf::json::get $step output]
    set mode [::ogf::json::get $out mode]
    switch -- $mode {
	add_columns {::ogf::cat::add_columns $output [::ogf::json::get $out columns]}
	set {::ogf::cat::load_tsv $output "[dict get $step label]"}
	text {OGFTextWindow "[dict get $step label]" $output}
	default {}
    }
    set p [::ogf::json::get $step after]
    if {$p ne {}} {catch {uplevel #0 $p}}
}

# stage progress for the strip: {stage {state elapsed_s}}; derived from the session recorder
proc ::ogf::step::stage_of_record {recstep} {
    # records globs come from the step manifests
    foreach id [::ogf::reg::ids] {
	foreach s [::ogf::json::get [::ogf::reg::get $id] steps] {
	    set st [::ogf::json::get $s stage]
	    if {$st eq {}} continue
	    foreach g [::ogf::json::get $s records] {
		if {[string match $g $recstep]} {return $st}
	    }
	    if {[::ogf::json::get $s proc] eq {} && $recstep eq "$id.[dict get $s id]"} {return $st}
	}
    }
    return {}
}
proc ::ogf::step::progress {} {
    set res [dict create]
    foreach st $::ogf::stages {dict set res $st [list pending 0 0]}
    foreach rec [::ogf::session::steps] {
	set st [stage_of_record [dict get $rec step]]
	if {$st eq {}} continue
	lassign [dict get $res $st] state sec n
	set ms [dict get $rec ms]
	if {[dict get $rec failed]} {
	    dict set res $st [list failed $sec $n]
	} else {
	    dict set res $st [list done [expr {$sec + ($ms eq {} ? 0 : $ms/1000.0)}] [expr {$n+1}]]
	}
    }
    return $res
}
