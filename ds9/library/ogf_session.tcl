#  OGFinder: session recorder ("Analysis > Save Session as Python Script...").
#
#  Every computational step that produces or changes a data product calls
#  OGFSessLog (directly, or through CatalogPanelCmdLog / OGFMaskRun /
#  CatalogPanelBands*).  A record is a dict:
#     seq step title class argv argv_t tool inputs outputs params payload
#     requires post band cat_in cat_after stdout_crc stdout_len ms failed ...
#  class: auto | manual | config | note
#  The records are exported into a self-contained Python 3 script
#  (template: ogf_session_template.py).  See docs/session_python_script.md

package provide DS9 1.0

proc OGFSessInit {} {
    global ogfsess
    set ogfsess(steps) {}
    set ogfsess(images) {}
    set ogfsess(maskfiles) {}
    set ogfsess(maskknown) {}
    set ogfsess(open) 0
    set ogfsess(t0) [clock seconds]
    set ogfsess(enabled) 1
    set ogfsess(last_note) {}
}

# ---------------------------------------------------------------- helpers
proc OGFSessRoot {} {
    return [file dirname [file dirname [info nameofexecutable]]]
}

proc OGFSessWorkDir {} {
    return [file join [file normalize ~] .ds9]
}

proc OGFSessCrc {s} {
    if {[catch {zlib crc32 [encoding convertto utf-8 $s]} c]} {return 0}
    return $c
}

proc OGFSessCatCrc {} {
    global catpanel
    if {![info exists catpanel(alldata)]} {return 0}
    return [OGFSessCrc $catpanel(alldata)]
}

# the review columns (plugins/report: REVIEW, REVIEW_NOTE, REVIEW_TIME) are hand annotations that the exported script never
# reproduces: leave them out of the catalog fingerprint so that "catalog identical to the GUI session" still holds in replay mode.
# Without such columns the text is returned untouched.
proc OGFSessPlainTSV {d} {
    set lines [split $d \n]
    set h [split [lindex $lines 0] \t]
    set drop {}
    set i 0
    foreach c $h {if {[string trim $c] in {REVIEW REVIEW_NOTE REVIEW_TIME}} {lappend drop $i}; incr i}
    if {![llength $drop]} {return $d}
    set out {}
    foreach l $lines {
	set f [split $l \t]
	foreach j [lreverse $drop] {if {$j < [llength $f]} {set f [lreplace $f $j $j]}}
	lappend out [join $f \t]
    }
    return [join $out \n]
}

proc OGFSessCatInfo {} {
    global catpanel
    set d {}
    if {[info exists catpanel(alldata)]} {set d [OGFSessPlainTSV $catpanel(alldata)]}
    if {$d eq {}} {return [list crc 0 rows 0 cols 0]}
    set lines [split $d \n]
    set rows 0
    foreach l [lrange $lines 1 end] {if {[string trim $l] ne {}} {incr rows}}
    return [list crc [OGFSessCrc $d] rows $rows cols [llength [split [lindex $lines 0] \t]]]
}

proc OGFSessIsImageFile {f} {
    if {$f eq {} || ![file isfile $f]} {return 0}
    if {[string first [OGFSessWorkDir] [file normalize $f]] == 0} {return 0}
    return [regexp -nocase {\.(fits|fit|fts)(\.gz)?$} $f]
}

# exposures of a multi-exposure tool: the values after --files (Moving Objects steps).  They are images of the session
# even when they live below the work directory (MAST cache ~/.ds9/mast_cache).
proc OGFSessFilesArg {argv} {
    set i [lsearch -exact $argv --files]
    if {$i < 0} {return {}}
    set out {}
    foreach a [lrange $argv [expr {$i+1}] end] {
	if {[string match -* $a]} break
	if {[file isfile $a] && [regexp -nocase {\.(fits|fit|fts)(\.gz)?$} $a]} {lappend out [file normalize $a]}
    }
    return $out
}

# image files used by an argv: first positional (python tools / sextract),
# --measure-image, --target-image
proc OGFSessArgvImages {argv} {
    set out {}
    set a1 [lindex $argv 1]
    if {[OGFSessIsImageFile $a1]} {lappend out [file normalize $a1]}
    set a2 [lindex $argv 2]
    if {[OGFSessIsImageFile $a2] && ![string match --* $a2]} {lappend out [file normalize $a2]}
    for {set i 0} {$i < [llength $argv]-1} {incr i} {
	if {[lindex $argv $i] in {--measure-image --target-image}} {
	    set v [lindex $argv [expr {$i+1}]]
	    if {[OGFSessIsImageFile $v]} {lappend out [file normalize $v]}
	}
    }
    foreach f [OGFSessFilesArg $argv] {lappend out $f}
    return [lsort -unique $out]
}

# argv -> params dict (flag -> value, bare flags -> 1)
proc OGFSessParams {argv} {
    set p {}
    set n [llength $argv]
    for {set i 1} {$i < $n} {incr i} {
	set a [lindex $argv $i]
	if {![string match --* $a]} continue
	set nx [lindex $argv [expr {$i+1}]]
	if {$i+1 < $n && ![regexp {^--[A-Za-z]} $nx]} {
	    dict set p $a $nx
	    incr i
	} else {
	    dict set p $a 1
	}
    }
    return $p
}

# log-time templating (image paths are mapped at export time)
proc OGFSessTemplate {argv {useroutputs {}}} {
    global catpanel
    set work [OGFSessWorkDir]
    set root [OGFSessRoot]
    set sbin [OGFSextractBin]
    set out {}
    set i 0
    set fl [OGFSessFilesArg $argv]
    foreach a $argv {
	set t $a
	if {[llength $fl] && [lsearch -exact $fl [file normalize $a]] >= 0 && [file isfile $a]} {
	    # exposure of a multi-exposure step: mapped to @{IMG:key} at export (the script substitutes the field images)
	    set t $a
	} elseif {[lsearch -exact $useroutputs $a] >= 0} {
	    set t @\{OUT\}/[file tail $a]
	} elseif {$i == 0 && $a eq [OGFPython]} {
	    set t @\{PY\}
	} elseif {$a eq $sbin || [file tail $a] eq "ds9_sextract"} {
	    set t @\{SEXTRACT\}
	} elseif {[string match *.py $a] && [file tail $a] ne $a && [string match ds9_* [file tail $a]]} {
	    set t @\{LIB\}/[file tail $a]
	} elseif {[string first $work/ $a] == 0} {
	    set rel [string range $a [string length $work/] end]
	    # temp catalogs are re-written by the script from the running catalog
	    if {[regexp {^(.*)_catalog\.tsv$} $rel -> nm] && [file exists $a] && [info exists catpanel(alldata)]} {
		set fd [open $a r]; fconfigure $fd -encoding utf-8
		set c [read $fd]; close $fd
		if {$c eq $catpanel(alldata)} {
		    set t @\{CAT:$nm\}
		} elseif {$c eq "$catpanel(alldata)\n"} {
		    set t @\{CATNL:$nm\}
		} else {
		    set t @\{WORK\}/$rel
		}
	    } else {
		set t @\{WORK\}/$rel
	    }
	} elseif {[string first $root/ $a] == 0} {
	    set t @\{ROOT\}/[string range $a [string length $root/] end]
	} elseif {[string first -- $a] != 0 && [string first $work $a] > 0} {
	    # "name:file" style arguments etc. -> left as is
	}
	lappend out $t
	incr i
    }
    return $out
}

# ------------------------------------------------------------ central log
# OGFSessLog STEP CLASS ARGV ?-title T -tool K -payload D -requires L -post D
#                          -outputs L -network 1 -note S -band K -noexec 1?
proc OGFSessLog {step class argv args} {
    global ogfsess catpanel
    if {![info exists ogfsess(enabled)] || !$ogfsess(enabled)} {return 0}
    array set o {-title {} -tool {} -payload {} -requires {} -post {} -outputs {} -network 0 -note {} -band {} -useroutputs {}}
    array set o $args
    OGFSessFinalize
    set tool $o(-tool)
    if {$tool eq {}} {
	set a0 [lindex $argv 0]
	if {[file tail $a0] eq "ds9_sextract"} {set tool sextract} elseif {$a0 ne {}} {set tool python} else {set tool internal}
    }
    set images [OGFSessArgvImages $argv]
    foreach im $images {
	if {[lsearch -exact $ogfsess(images) $im] < 0} {lappend ogfsess(images) $im}
    }
    set targv [expr {[llength $argv] ? [OGFSessTemplate $argv $o(-useroutputs)] : {}}]
    set inputs {}
    foreach a $argv {
	if {[string match -* $a]} continue
	if {[catch {file isfile $a} isf] || !$isf} continue
	lappend inputs $a
    }
    set winputs {}
    set nn [llength $argv]
    for {set i 1} {$i < $nn} {incr i} {
	set a [lindex $argv $i]
	if {[string match -* $a]} continue
	if {[string first [OGFSessWorkDir]/ $a] != 0 || ![file exists $a]} continue
	set pf [lindex $argv [expr {$i-1}]]
	if {[regexp {output|--catalog$|--regions$|--workdir$} $pf]} continue
	if {[lsearch -exact [OGFSessFilesArg $argv] [file normalize $a]] >= 0} continue
	if {[regexp {_catalog\.tsv$} $a]} continue
	if {[regexp {mask_regions\.reg$} $a]} continue
	if {[lindex $argv 1] ne {} && [file tail [lindex $argv 1]] eq "ds9_mask.py"} {
	    # mask file is an input only for modes that need an existing mask
	    set mi [lsearch -exact $argv --mode]
	    set mm [expr {$mi >= 0 ? [lindex $argv [expr {$mi+1}]] : {}}]
	    if {$mm ni {grow shrink masked export reproject undo redo}} continue
	}
	lappend winputs $i
    }
    set outs $o(-outputs)
    set n [llength $argv]
    for {set i 1} {$i < $n-1} {incr i} {
	if {[regexp {^--([a-z-]*output|output)$} [lindex $argv $i]]} {lappend outs [lindex $argv [expr {$i+1}]]}
    }
    set seq [expr {[llength $ogfsess(steps)] + 1}]
    set rec [dict create seq $seq step $step title $o(-title) class $class \
	argv $argv argv_t $targv tool $tool inputs $inputs outputs $outs \
	params [OGFSessParams $argv] payload $o(-payload) requires $o(-requires) \
	post $o(-post) band $o(-band) network $o(-network) note $o(-note) undone 0 winputs {} \
	winput_idx $winputs cat_in [OGFSessCatInfo] cat_after {} stdout_crc {} stdout_len {} ms {} failed 0 \
	wall [clock format [clock seconds] -format "%Y-%m-%dT%H:%M:%S%z"] \
	images $images]
    lappend ogfsess(steps) $rec
    set ogfsess(open) $seq
    return $seq
}

proc OGFSessSet {seq key val} {
    global ogfsess
    if {$seq < 1 || $seq > [llength $ogfsess(steps)]} return
    set i [expr {$seq-1}]
    set rec [lindex $ogfsess(steps) $i]
    dict set rec $key $val
    lset ogfsess(steps) $i $rec
}

proc OGFSessLast {} {
    global ogfsess
    if {![info exists ogfsess(steps)] || [llength $ogfsess(steps)] == 0} {return {}}
    return [lindex $ogfsess(steps) end]
}

# explicit post-processing declaration for the step that is currently open
proc OGFSessPost {post} {
    global ogfsess
    if {![info exists ogfsess(open)] || !$ogfsess(open)} return
    OGFSessSet $ogfsess(open) post $post
}

# --------------------------------------------- exec wrapper: timing, failure
rename exec OGFsess_exec
proc exec {args} {
    global ogfsess
    set t0 [clock milliseconds]
    set rc [catch {uplevel 1 [linsert $args 0 OGFsess_exec]} res opts]
    catch {OGFSessExecDone $args $rc $res $t0}
    return -options $opts $res
}

proc OGFSessExecDone {args rc res t0} {
    global ogfsess
    if {![info exists ogfsess(open)] || !$ogfsess(open)} return
    set seq $ogfsess(open)
    set rec [lindex $ogfsess(steps) [expr {$seq-1}]]
    set argv [dict get $rec argv]
    if {[llength $argv] == 0 || [dict get $rec ms] ne {}} return
    if {[lrange $args 0 [expr {[llength $argv]-1}]] ne $argv} return
    OGFSessSet $seq ms [expr {[clock milliseconds] - $t0}]
    if {$rc} {
	OGFSessSet $seq failed 1
	OGFSessSet $seq note "GUI run failed: [string range $res 0 200]"
    } else {
	OGFSessSet $seq stdout_crc [OGFSessCrc $res]
	OGFSessSet $seq stdout_len [string length $res]
    }
}

# The catalog fingerprint after a step is captured lazily, when the next step
# is logged (or at export): every catalog-changing action is logged explicitly,
# so nothing changes in between.
proc OGFSessFinalize {} {
    global ogfsess
    if {![info exists ogfsess(open)] || !$ogfsess(open)} return
    if {[llength $ogfsess(steps)] < $ogfsess(open)} return
    OGFSessSet $ogfsess(open) cat_after [OGFSessCatInfo]
    set ogfsess(open) 0
}

# CatalogPanelLoadTSV hook: link the step that is open to the catalog it produced.
# kind=set when the table now equals the step's stdout (so the script can rebuild it).
proc OGFSessOnCatalogLoad {source_name} {
    global ogfsess catpanel
    if {![info exists ogfsess(enabled)] || !$ogfsess(enabled)} return
    if {$source_name eq "visible"} {
	OGFSessLog catalog.unrecorded manual {} -tool internal \
	    -title "Show Visible Only replaced the catalog by the sources in view" \
	    -note "view-dependent subset (depends on zoom/pan); not replayable"
	return
    }
    if {![info exists ogfsess(open)] || !$ogfsess(open)} return
    set seq $ogfsess(open)
    set rec [lindex $ogfsess(steps) [expr {$seq-1}]]
    if {[dict get $rec post] ne {} || [dict get $rec stdout_crc] eq {}} return
    set crc [OGFSessCrc $catpanel(alldata)]
    if {$crc == [dict get $rec stdout_crc]} {
	OGFSessSet $seq post [dict create kind set]
    }
}

# CatalogPanelAddColumnsFromTSV hook (called with its raw arguments)
proc OGFSessOnAddColumns {result_data col_names} {
    global ogfsess
    if {![info exists ogfsess(enabled)] || !$ogfsess(enabled)} return
    if {![info exists ogfsess(open)] || !$ogfsess(open)} return
    set seq $ogfsess(open)
    set rec [lindex $ogfsess(steps) [expr {$seq-1}]]
    if {[dict get $rec post] ne {} || [dict get $rec stdout_crc] eq {}} return
    if {[OGFSessCrc $result_data] == [dict get $rec stdout_crc]} {
	OGFSessSet $seq post [dict create kind add cols_list $col_names]
    }
}

# ----------------------------------------------- ICL / LSBG command log
# Called from CatalogPanelCmdLog.  Mask commands are logged by OGFMaskRun.
proc OGFSessFromCmdLog {pipeline argv} {
    if {[llength $argv] < 2} return
    set script [file tail [lindex $argv 1]]
    if {$script eq "ds9_mask.py"} return
    if {![regexp {^ds9_(icl|lsbg)\.py$} $script]} return
    set mode {}
    set i [lsearch -exact $argv --mode]
    if {$i >= 0} {set mode [lindex $argv [expr {$i+1}]]}
    set payload {}
    set class auto
    set req {}
    set j [lsearch -exact $argv --center]
    if {$j >= 0} {
	dict set payload center [lindex $argv [expr {$j+1}]]
	lappend req extract_catalog
    }
    if {$pipeline eq "lsbg" && $mode in {forced}} {set class manual}
    set titles {icl,background "ICL background model" icl,profile "ICL surface-brightness profile"
	icl,measure "ICL measurements" icl,measure-multi "ICL multi-threshold" icl,decompose "ICL decomposition"
	icl,color "ICL colour profile" lsbg,clean "LSBG background/clean" lsbg,detect "LSBG detect"
	lsbg,photometry "LSBG photometry" lsbg,sersic "LSBG Sersic" lsbg,filter "LSBG filter"
	lsbg,svm-classify "LSBG SVM classify" lsbg,run "LSBG full pipeline" lsbg,forced "LSBG forced photometry"}
    set title $pipeline:$mode
    if {[dict exists $titles $pipeline,$mode]} {set title [dict get $titles $pipeline,$mode]}
    if {$pipeline eq "lsbg" && $mode eq "svm-classify"} {lappend req catalog}
    set post {}
    if {$pipeline eq "lsbg" && $mode in {detect photometry sersic filter run}} {
	set post [dict create kind set force 1]
    }
    OGFSessLog $pipeline.$mode $class $argv -title $title -payload $payload -requires $req -post $post
}

# ------------------------------------------------------------ JSON output
proc OGFJStr {s} {
    set s [string map [list \\ \\\\ \" \\\" \n \\n \r \\r \t \\t \b \\b \f \\f] $s]
    set out {}
    foreach ch [split $s {}] {
	scan $ch %c c
	if {$c < 32} {
	    append out [format {\u%04x} $c]
	} elseif {$c > 126} {
	    if {$c > 0xFFFF} {
		set c2 [expr {$c - 0x10000}]
		append out [format {\u%04x\u%04x} [expr {0xD800 + ($c2 >> 10)}] [expr {0xDC00 + ($c2 & 0x3FF)}]]
	    } else {
		append out [format {\u%04x} $c]
	    }
	} else {
	    append out $ch
	}
    }
    return "\"$out\""
}

proc OGFJList {l} {
    set parts {}
    foreach x $l {lappend parts [OGFJStr $x]}
    return "\[[join $parts {, }]\]"
}

# dict -> JSON object.  key suffix _list -> array of strings, _json -> raw JSON
proc OGFJDict {d {ind {}}} {
    if {[dict size $d] == 0} {return "\{\}"}
    set parts {}
    dict for {k v} $d {
	if {[string match *_list $k]} {
	    lappend parts "[OGFJStr [string range $k 0 end-5]]: [OGFJList $v]"
	} elseif {[string match *_json $k]} {
	    lappend parts "[OGFJStr [string range $k 0 end-5]]: $v"
	} else {
	    lappend parts "[OGFJStr $k]: [OGFJStr $v]"
	}
    }
    return "\{[join $parts {, }]\}"
}

proc OGFSessRecordJSON {rec} {
    set lines {}
    foreach k {seq step title class band tool requires network} {
	if {$k eq "requires"} {
	    lappend lines "  [OGFJStr $k]: [OGFJList [dict get $rec $k]]"
	} else {
	    lappend lines "  [OGFJStr $k]: [OGFJStr [dict get $rec $k]]"
	}
    }
    foreach k {argv argv_t inputs outputs winputs} {
	lappend lines "  [OGFJStr $k]: [OGFJList [dict get $rec $k]]"
    }
    foreach k {params payload post cat_in cat_after} {
	lappend lines "  [OGFJStr $k]: [OGFJDict [dict get $rec $k]]"
    }
    foreach k {stdout_crc stdout_len ms failed wall note undone} {
	lappend lines "  [OGFJStr $k]: [OGFJStr [dict get $rec $k]]"
    }
    return "\{\n[join $lines ",\n"]\n\}"
}

# --------------------------------------------------------------- export
# band-key assignment, image-path / basename tokens
proc OGFSessKeys {steps} {
    global ogfband ogfsess
    set keys [dict create]       ;# path -> key
    # registered band names win
    foreach rec $steps {
	if {[dict get $rec step] ne "bands.register"} continue
	set pl [dict get $rec payload]
	set f [file normalize [dict get $pl file]]
	if {![dict exists $keys $f]} {dict set keys $f [dict get $pl name]}
    }
    set n 0
    foreach im $ogfsess(images) {
	if {[dict exists $keys $im]} continue
	incr n
	dict set keys $im [expr {$n == 1 ? "main" : "img$n"}]
	if {$n == 1 && [lsearch -exact [dict values $keys] main] < 0} {}
    }
    return $keys
}

proc OGFSessMapArgs {argv_t keys} {
    # image paths (longest first), then basenames
    set paths [lsort -command {apply {{a b} {expr {[string length $b] - [string length $a]}}}} [dict keys $keys]]
    set m1 {}
    foreach p $paths {lappend m1 $p "@\{IMG:[dict get $keys $p]\}"}
    set m2 {}
    set bases {}
    foreach p $paths {lappend bases [list [CatalogPanelFitsBaseName $p] [dict get $keys $p]]}
    set bases [lsort -command {apply {{a b} {expr {[string length [lindex $b 0]] - [string length [lindex $a 0]]}}}} $bases]
    foreach b $bases {lappend m2 [lindex $b 0] "@\{BASE:[lindex $b 1]\}"}
    set out {}
    foreach a $argv_t {
	set r [string map $m1 $a]
	if {$r eq $a && ([string match @\{WORK\}/* $a] || [string match @\{OUT\}/* $a])} {
	    set r [string map $m2 $a]
	}
	lappend out $r
    }
    return $out
}

proc OGFSessExport {path} {
    global ogfsess ds9 catpanel
    set steps $ogfsess(steps)
    if {[llength $steps] == 0} {return 0}
    OGFSessFinalize
    set steps $ogfsess(steps)
    set keys [OGFSessKeys $steps]
    # detection key
    set detkey {}
    foreach rec $steps {
	if {[dict get $rec step] eq "extract" && [llength [dict get $rec images]]} {
	    set detkey [dict get $keys [lindex [dict get $rec images] 0]]
	    break
	}
    }
    if {$detkey eq {} && [llength [dict keys $keys]]} {set detkey [lindex [dict values $keys] 0]}
    set bdet {}
    global ogfband
    if {[info exists ogfband(detect)] && $ogfband(detect) ne {} && [dict exists $keys [file normalize $ogfband($ogfband(detect),file)]]} {
	set bdet [dict get $keys [file normalize $ogfband($ogfband(detect),file)]]
    }
    if {$bdet ne {}} {set detkey $bdet}
    # images table
    set imgs {}
    dict for {p k} $keys {
	set info [OGFImageInfo $p]
	set d [dict create key $k path $p base [CatalogPanelFitsBaseName $p]]
	foreach f {FILTER PIVOT ZP_AB PIXSCALE NAXIS1 NAXIS2 INSTRUME} {
	    if {[dict exists $info $f]} {dict set d $f [dict get $info $f]}
	}
	if {[file exists $p]} {dict set d size [file size $p]}
	lappend imgs [OGFJDict $d]
    }
    # steps
    set js {}
    set summary {}
    set k 0
    foreach rec $steps {
	incr k
	dict set rec argv_t [OGFSessMapArgs [dict get $rec argv_t] $keys]
	set wi {}
	foreach ix [dict get $rec winput_idx] {lappend wi [lindex [dict get $rec argv_t] $ix]}
	dict set rec winputs $wi
	set im [dict get $rec images]
	set bk {}
	if {[llength $im]} {set bk [dict get $keys [lindex $im 0]]}
	if {[dict get $rec band] eq {}} {dict set rec band $bk}
	# the literal argv keeps GUI paths; mask them in the exported file
	lappend js [OGFSessRecordJSON $rec]
	set cls [dict get $rec class]
	set tag [expr {[dict get $rec failed] ? "FAILED-IN-GUI" : [string toupper $cls]}]
	lappend summary [format "#   %3d  %-14s %-24s %s" $k $tag [dict get $rec step] [dict get $rec title]]
    }
    set root [OGFSessRoot]
    set git {}
    catch {set git [string trim [OGFsess_exec git -C $root rev-parse HEAD]]}
    set dirty {}
    catch {set dirty [expr {[string trim [OGFsess_exec git -C $root status --porcelain -uno -- ds9/library icl lsbg sersic_fit ai_merge photo_z psf_deconv star_finder sed_fit bulge_disk morphometry galaxy_morph parallel]] ne {}}]}
    set meta [dict create ogfinder_root $root ds9_version [expr {[info exists ds9(version)] ? $ds9(version) : {}}] \
	ds9_exe [info nameofexecutable] git_head $git git_dirty_sources $dirty \
	python [OGFPython] exported [clock format [clock seconds] -format "%Y-%m-%dT%H:%M:%S%z"] \
	detection_key $detkey workdir [OGFSessWorkDir] ephem_dir [expr {[file isdirectory [file join [OGFSessWorkDir] ephem]] ? [file join [OGFSessWorkDir] ephem] : {}}] n_workers_gui [expr {[info exists catpanel(param,n-workers)] ? $catpanel(param,n-workers) : {}}] \
	session_started [clock format $ogfsess(t0) -format "%Y-%m-%dT%H:%M:%S%z"]]
    set meta_json [OGFJDict $meta]
    set tplf [CatalogPanelGetScript ogf_session_template.py]
    if {![file exists $tplf]} {error "ogf_session_template.py not found"}
    set fd [open $tplf r]; fconfigure $fd -encoding utf-8
    set tpl [read $fd]; close $fd
    set header "# Recorded steps (AUTO = runs in --mode pipeline; MANUAL = only with --include-manual / --mode replay):\n[join $summary \n]\n"
    set data "SESSION_META = $meta_json\n\nIMAGES = \[\n  [join $imgs ",\n  "]\n\]\n\nSTEPS = \[\n[join $js ",\n"]\n\]\n"
    set tpl [string map [list @@HEADER@@ $header @@DATA@@ $data @@ROOT@@ $root] $tpl]
    set fd [open $path w]; fconfigure $fd -encoding utf-8
    puts -nonewline $fd $tpl
    close $fd
    catch {file attributes $path -permissions 0755}
    return [llength $steps]
}

# ----------------------------------------------------- menu procs (GUI)
proc CatalogPanelSessionSave {{path {}}} {
    global catpanel ogfsess
    if {![info exists ogfsess(steps)] || [llength $ogfsess(steps)] == 0} {
	set catpanel(status) "Session: nothing recorded yet"
	return
    }
    if {$path eq {}} {
	set path [tk_getSaveFile -title "Save Session as Python Script..." \
	    -initialfile ogfinder_session.py -defaultextension .py \
	    -filetypes {{{Python script} {.py}} {{All files} *}}]
    }
    if {$path eq {}} return
    if {[catch {set n [OGFSessExport $path]} err]} {
	set catpanel(status) "Session export failed: $err"
	return
    }
    set catpanel(status) "Session: $n steps written to [file tail $path]"
    return $path
}

proc CatalogPanelSessionReset {{ask 1}} {
    global catpanel ogfsess
    if {$ask && [tk_messageBox -type yesno -icon question -title "Reset Session Log" \
	    -message "Discard all [llength $ogfsess(steps)] recorded steps?"] ne "yes"} return
    OGFSessInit
    set catpanel(status) "Session log reset"
}

proc OGFSessLogText {} {
    global ogfsess
    set t "OGFinder session log - [llength $ogfsess(steps)] recorded step(s)\n"
    append t "(class AUTO = pipeline mode, MANUAL = replay only; argv is what the GUI executed)\n\n"
    foreach rec $ogfsess(steps) {
	set cls [string toupper [dict get $rec class]]
	if {[dict get $rec failed]} {set cls FAILED}
	append t [format "%3d  %-8s %-22s %s" [dict get $rec seq] $cls [dict get $rec step] [dict get $rec title]]
	if {[dict get $rec ms] ne {}} {append t [format "   (%.1f s)" [expr {[dict get $rec ms]/1000.0}]]}
	append t "\n"
	if {[llength [dict get $rec argv]]} {append t "       argv: [join [lrange [dict get $rec argv] 1 end] { }]\n"}
	if {[dict size [dict get $rec payload]]} {
	    set pl [dict get $rec payload]
	    foreach k [dict keys $pl] {append t "       payload.$k: [string range [dict get $pl $k] 0 200]\n"}
	}
	if {[dict get $rec note] ne {}} {append t "       note: [dict get $rec note]\n"}
    }
    return $t
}

proc CatalogPanelSessionShow {} {
    set w .ogfsesslog
    catch {destroy $w}
    toplevel $w
    wm title $w "Session Log"
    text $w.t -width 120 -height 30 -font TkFixedFont -wrap none -yscrollcommand [list $w.sb set]
    ttk::scrollbar $w.sb -command [list $w.t yview]
    $w.t insert end [OGFSessLogText]
    $w.t configure -state disabled
    ttk::button $w.close -text Close -command [list destroy $w]
    pack $w.close -side bottom -pady 4
    pack $w.sb -side right -fill y
    pack $w.t -fill both -expand true
}

# ------------------------------------------------- manual-step recorders
proc OGFSessNoteSeen {key} {
    global ogfsess
    if {[lsearch -exact $ogfsess(maskknown) $key] >= 0} {return 1}
    lappend ogfsess(maskknown) $key
    return 0
}

proc OGFSessInitMasks {} {
    global ogfsess
    if {![info exists ogfsess(maskfiles)]} {set ogfsess(maskfiles) {}}
}
