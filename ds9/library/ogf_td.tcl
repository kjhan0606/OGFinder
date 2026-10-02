#  OGFinder core service: the time-domain table.
#
#  ONE table (the catalog table of the panel) shows rows of several kinds:  galaxy | moving | transient
#  (+ "detection", a plugin-defined kind that is only shown on its own).  A kind filter (Galaxies / Moving /
#  Transients / All) selects what the table and the image markers show.
#
#    * galaxies   = the SExtractor catalog in catpanel(alldata).  It is never modified by this layer; in the
#                   Galaxies view the table is exactly what it always was (plus a trailing `kind` column once any
#                   other kind has data).
#    * other kinds are registered by plugins (::ogf::td::register_kind) and filled with
#                   ::ogf::td::set_rows KIND COLS ROWS.   Raw columns must contain  id ra dec  (degrees); `mag` is
#                   optional; a column named _pos holds "ra1 dec1 ra2 dec2 ..." extra marker positions (a track).
#    * the view always has the common columns  NUMBER kind X_IMAGE Y_IMAGE ALPHA_J2000 DELTA_J2000 MAG_AUTO
#      (NUMBER = key: galaxy "123", moving "M5", transient "T7"; X/Y are pixels of the detection grid).
#      With a single kind selected its own columns follow.  A plugin may add computed columns with a `decorate`
#      proc (host_id / host_sep for transients, overlap for moving objects use ::ogf::td::galaxy_near).
#    * sort, filter, save, hover, click, multi-select all go through the same procs as for galaxies (the link
#      code reads the displayed TSV through OGFViewTSV).  Sort / save of a non-galaxy view are recorded in the
#      session as informational ("note") steps: the galaxy catalog of the exported script is not touched.
#    * markers are drawn in sky coordinates (fk5) into EVERY frame, so they appear in all tiles and in every band;
#      a click on one (any tile) selects the table row (OGFTDKeyFromTags, used by CatalogPanelMarkerClick and
#      OGFTileClick).
#  See docs/plugins.md ("Time-domain table").

package provide DS9 1.0

namespace eval ::ogf::td {}

proc OGFTDInit {} {
    global ogftd
    array unset ogftd
    array set ogftd {
	kind galaxies radio galaxies selcbs {} selframes {} kinds {galaxies} stale 1 galstale 1 sortcol {} sortdir {} sorted {}
	view,tsv {} view,cols {} view,rows {} hostradius 5.0 pending {} hidden 0 gridkey {} gridinfo {}
    }
    set ogftd(def,galaxies) [dict create label Galaxies prefix {} color yellow columns {} decorate {} order 0 all 1 point {}]
    set ogftd(order) {galaxies}
    trace add variable ::catpanel(alldata) write OGFTDAlldataChanged
}

# ------------------------------------------------------------------ registry
# options: -label -prefix -color -columns {raw cols shown in the single-kind view} -decorate PROC -order N
#          -show_in_all 0|1 (default 1)
proc ::ogf::td::register_kind {kind args} {
    global ogftd
    array set o {-label {} -prefix {} -color white -columns {} -decorate {} -order 50 -show_in_all 1 -point {point=circle 12}}
    array set o $args
    if {$o(-label) eq {}} {set o(-label) [string totitle $kind]}
    if {$o(-prefix) eq {}} {set o(-prefix) [string toupper [string index $kind 0]]}
    set ogftd(def,$kind) [dict create label $o(-label) prefix $o(-prefix) color $o(-color) columns $o(-columns) \
	decorate $o(-decorate) order $o(-order) all $o(-show_in_all) point $o(-point)]
    if {$kind ni $ogftd(kinds)} {lappend ogftd(kinds) $kind}
    set ogftd(kinds) [lsort -command ::ogf::td::_cmp_order $ogftd(kinds)]
    OGFTDControlsRefresh
}
proc ::ogf::td::_cmp_order {a b} {
    global ogftd
    return [expr {[dict get $ogftd(def,$a) order] - [dict get $ogftd(def,$b) order]}]
}
proc ::ogf::td::kinds {} {global ogftd; return $ogftd(kinds)}
proc ::ogf::td::label {kind} {global ogftd; return [dict get $ogftd(def,$kind) label]}
proc ::ogf::td::color {kind} {global ogftd; return [dict get $ogftd(def,$kind) color]}
proc ::ogf::td::kind {} {global ogftd; return $ogftd(kind)}
proc ::ogf::td::active {} {global ogftd; return [expr {$ogftd(kind) ne "galaxies"}]}
proc ::ogf::td::has_data {kind} {global ogftd; return [info exists ogftd(raw,$kind)]}
proc ::ogf::td::any_data {} {
    global ogftd
    foreach k $ogftd(kinds) {if {$k ne "galaxies" && [info exists ogftd(raw,$k)]} {return 1}}
    return 0
}

# kinds visible in the current filter
proc ::ogf::td::visible_kinds {} {
    global ogftd
    switch -- $ogftd(kind) {
	galaxies {return {}}
	all {
	    set l {}
	    foreach k $ogftd(kinds) {
		if {$k eq "galaxies"} continue
		if {[dict get $ogftd(def,$k) all] && [info exists ogftd(raw,$k)]} {lappend l $k}
	    }
	    return $l
	}
	default {return [list $ogftd(kind)]}
    }
}
proc OGFTDHidesGalaxies {} {
    global ogftd
    return [expr {[info exists ogftd(kind)] && $ogftd(kind) ni {galaxies all}}]
}

# ------------------------------------------------------------------ data in
proc ::ogf::td::set_rows {kind cols rows {show 0}} {
    global ogftd
    if {![info exists ogftd(def,$kind)]} {error "time-domain kind $kind is not registered"}
    set ogftd(raw,$kind) [list $cols $rows]
    set ogftd(stale) 1
    catch {review_prune $kind}
    OGFTDControlsRefresh
    if {$show} {show $kind} elseif {[active]} {refresh}
}
proc ::ogf::td::clear {kind} {
    global ogftd
    unset -nocomplain ogftd(raw,$kind)
    if {[info exists ogftd(review)]} {
	set pre [dict get $ogftd(def,$kind) prefix]
	foreach k [dict keys $ogftd(review)] {if {[string match ${pre}* $k]} {dict unset ogftd(review) $k}}
    }
    set ogftd(stale) 1
    if {$ogftd(kind) eq $kind} {show galaxies} else {refresh}
    OGFTDControlsRefresh
}

# raw row of a key (dict col->value) or {}
proc ::ogf::td::row_of {key} {
    global ogftd
    foreach k $ogftd(kinds) {
	if {$k eq "galaxies" || ![info exists ogftd(raw,$k)]} continue
	set pre [dict get $ogftd(def,$k) prefix]
	if {![string match ${pre}* $key]} continue
	set id [string range $key [string length $pre] end]
	lassign $ogftd(raw,$k) cols rows
	set ic [lsearch -exact $cols id]
	foreach r $rows {
	    if {[lindex $r $ic] eq $id} {
		return [dict create kind $k id $id {*}[concat {*}[lmap c $cols v $r {list $c $v}]]]
	    }
	}
    }
    return {}
}
proc ::ogf::td::key_kind {key} {
    global ogftd
    foreach k $ogftd(kinds) {
	if {$k eq "galaxies"} continue
	set pre [dict get $ogftd(def,$k) prefix]
	if {[string match ${pre}* $key] && [info exists ogftd(raw,$k)]} {
	    if {![string is integer -strict [string range $key [string length $pre] end]]} continue
	    return $k
	}
    }
    return galaxies
}

# ------------------------------------------------------------------ galaxy catalog (sky index)
# nearest galaxies around (ra,dec): list of {number sep_arcsec extent_arcsec} within `rad` arcsec, sorted by sep
proc ::ogf::td::_gal_build {} {
    global ogftd catpanel
    if {!$ogftd(galstale)} return
    set ogftd(galstale) 0
    foreach k {num ra dec ext bins} {set ogftd(gal,$k) {}}
    set ogftd(gal,n) 0
    set ogftd(gal,maxext) 0.0
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return
    set lines [split $catpanel(alldata) \n]
    set h [split [lindex $lines 0] \t]
    set cn [lsearch -exact $h NUMBER]; set ca [lsearch -exact $h ALPHA_J2000]; set cd [lsearch -exact $h DELTA_J2000]
    set ci [lsearch -exact $h ISO_RADIUS]
    if {$ca < 0 || $cd < 0} return
    set ps [pixscale]
    set nums {}; set ras {}; set decs {}; set exts {}
    set i 0
    foreach l [lrange $lines 1 end] {
	incr i
	if {[string trim $l] eq {}} continue
	set f [split $l \t]
	set ra [string trim [lindex $f $ca]]; set de [string trim [lindex $f $cd]]
	if {![string is double -strict $ra] || ![string is double -strict $de]} continue
	set n $i
	if {$cn >= 0} {set nv [string trim [lindex $f $cn]]; if {$nv ne {}} {set n $nv}}
	set iso 5.0
	if {$ci >= 0} {set v [string trim [lindex $f $ci]]; if {[string is double -strict $v] && $v > 0} {set iso $v}}
	lappend nums $n; lappend ras $ra; lappend decs $de; lappend exts [expr {$iso * $ps}]
    }
    if {[llength $nums] == 0} return
    set d0 0.0
    foreach d $decs {set d0 [expr {$d0 + $d}]}
    set d0 [expr {$d0 / [llength $decs]}]
    set ogftd(gal,cosd) [expr {cos($d0 * 0.017453292519943295)}]
    set cell 30.0
    set ogftd(gal,cell) $cell
    set bins [dict create]
    set k 0; set mx 0.0
    foreach ra $ras de $decs e $exts {
	set x [expr {$ra * $ogftd(gal,cosd) * 3600.0}]; set y [expr {$de * 3600.0}]
	dict lappend bins "[expr {int(floor($x/$cell))}],[expr {int(floor($y/$cell))}]" $k
	if {$e > $mx} {set mx $e}
	incr k
    }
    set ogftd(gal,num) $nums; set ogftd(gal,ra) $ras; set ogftd(gal,dec) $decs; set ogftd(gal,ext) $exts
    set ogftd(gal,bins) $bins; set ogftd(gal,n) [llength $nums]; set ogftd(gal,maxext) $mx
}

proc ::ogf::td::have_galaxies {} {
    global ogftd
    _gal_build
    return [expr {$ogftd(gal,n) > 0}]
}
proc ::ogf::td::galaxy_max_extent {} {global ogftd; _gal_build; return $ogftd(gal,maxext)}

proc ::ogf::td::galaxy_near {ra dec rad} {
    global ogftd
    _gal_build
    if {$ogftd(gal,n) == 0} {return {}}
    set cell $ogftd(gal,cell); set cd $ogftd(gal,cosd)
    set x [expr {$ra * $cd * 3600.0}]; set y [expr {$dec * 3600.0}]
    set i0 [expr {int(floor(($x-$rad)/$cell))}]; set i1 [expr {int(floor(($x+$rad)/$cell))}]
    set j0 [expr {int(floor(($y-$rad)/$cell))}]; set j1 [expr {int(floor(($y+$rad)/$cell))}]
    set out {}
    set r2 [expr {$rad*$rad}]
    for {set i $i0} {$i <= $i1} {incr i} {
	for {set j $j0} {$j <= $j1} {incr j} {
	    if {![dict exists $ogftd(gal,bins) "$i,$j"]} continue
	    foreach k [dict get $ogftd(gal,bins) "$i,$j"] {
		set dx [expr {([lindex $ogftd(gal,ra) $k] - $ra) * $cd * 3600.0}]
		set dy [expr {([lindex $ogftd(gal,dec) $k] - $dec) * 3600.0}]
		set d2 [expr {$dx*$dx + $dy*$dy}]
		if {$d2 <= $r2} {lappend out [list [lindex $ogftd(gal,num) $k] [expr {sqrt($d2)}] [lindex $ogftd(gal,ext) $k]]}
	    }
	}
    }
    return [lsort -real -index 1 $out]
}

# ------------------------------------------------------------------ image grid (X/Y of non-galaxy rows)
# header info of the detection grid: the detection band if bands are registered, else the current frame's image
proc ::ogf::td::grid_info {} {
    global ogftd ogfband
    if {[info commands OGFBandsActive] ne {} && [OGFBandsActive] && $ogfband(detect) ne {} &&
	[info exists ogfband($ogfband(detect),info)]} {
	return $ogfband($ogfband(detect),info)
    }
    set fn {}
    catch {set fn [CatalogPanelGetFITS]}
    if {$fn eq {} || ![file exists $fn]} {return {}}
    if {$ogftd(gridkey) ne $fn} {
	set ogftd(gridkey) $fn
	set ogftd(gridinfo) [OGFImageInfo $fn]
    }
    return $ogftd(gridinfo)
}
proc ::ogf::td::pixscale {} {
    global catpanel
    set d [grid_info]
    if {[catch {
	set det [expr {abs([dict get $d CD1_1]*[dict get $d CD2_2] - [dict get $d CD1_2]*[dict get $d CD2_1])}]
	set s [expr {sqrt($det)*3600.0}]
    }] || $s <= 0} {
	set s 1.0
	catch {set s $catpanel(param,pixel-scale)}
    }
    return $s
}
proc ::ogf::td::world2pix {ra dec} {
    set d [grid_info]
    if {$d eq {}} {return {{} {}}}
    if {[catch {
	if {[dict get $d WCSVALID] != 1} {error nowcs}
	set r [expr {3.141592653589793/180.0}]
	lassign [OGFWorld2Pix $d [expr {$ra*$r}] [expr {$dec*$r}]] x y
    }]} {return {{} {}}}
    return [list [format %.2f $x] [format %.2f $y]]
}

# ------------------------------------------------------------------ the view
set ::ogf::td::common {NUMBER kind X_IMAGE Y_IMAGE ALPHA_J2000 DELTA_J2000 MAG_AUTO}

proc ::ogf::td::_fmt {v} {
    if {[string is double -strict $v] && ![string is integer -strict $v]} {return [format %.6g $v]}
    return $v
}

# rows of one kind as {cols rows} in view form: common columns, then the kind's own (when `own`)
proc ::ogf::td::_kind_view {kind own} {
    global ogftd
    lassign $ogftd(raw,$kind) cols rows
    set def $ogftd(def,$kind)
    set dec [dict get $def decorate]
    if {$dec ne {} && [info commands [lindex $dec 0]] ne {}} {
	if {![catch {uplevel #0 [list {*}$dec $cols $rows]} res]} {lassign $res cols rows}
    }
    set pre [dict get $def prefix]
    set ci [lsearch -exact $cols id]; set cr [lsearch -exact $cols ra]; set cd [lsearch -exact $cols dec]
    set cm [lsearch -exact $cols mag]
    set extra {}
    if {$own} {
	foreach c [dict get $def columns] {
	    set j [lsearch -exact $cols $c]
	    if {$j >= 0} {lappend extra $j}
	}
    }
    set outrows {}
    foreach r $rows {
	set ra [lindex $r $cr]; set de [lindex $r $cd]
	if {![string is double -strict $ra] || ![string is double -strict $de]} continue
	lassign [world2pix $ra $de] x y
	set m [expr {$cm >= 0 ? [lindex $r $cm] : {}}]
	if {[string is double -strict $m]} {set m [format %.3f $m]}
	set o [list $pre[lindex $r $ci] $kind $x $y [_fmt $ra] [_fmt $de] $m]
	foreach j $extra {lappend o [_fmt [lindex $r $j]]}
	lappend outrows $o
    }
    set ocols $::ogf::td::common
    foreach j $extra {lappend ocols [lindex $cols $j]}
    return [list $ocols $outrows]
}

proc ::ogf::td::build {} {
    global ogftd
    if {!$ogftd(stale)} return
    set ogftd(stale) 0
    set kind $ogftd(kind)
    set cols $::ogf::td::common
    set rows {}
    if {$kind eq "galaxies"} {
	# the view of the galaxy kind is alldata itself (displayed by the normal table code)
	set ogftd(view,cols) {}; set ogftd(view,rows) {}; set ogftd(view,tsv) {}
	return
    }
    set vis [visible_kinds]
    if {$kind eq "all"} {
	# galaxies first (common columns of the current galaxy catalog)
	_gal_view_rows rows
	foreach k $vis {
	    lassign [_kind_view $k 0] c r
	    lappend rows {*}$r
	}
    } else {
	lassign [_kind_view $kind 1] cols rows
    }
    # review columns (only once a decision exists or a REVIEW filter is on), then sort
    lassign [_review_view $cols $rows] cols rows
    if {$ogftd(sortcol) ne {}} {
	set ix [lsearch -exact $cols $ogftd(sortcol)]
	if {$ix >= 0} {
	    set num 1
	    foreach r $rows {
		set v [lindex $r $ix]
		if {$v ne {}} {if {![string is double -strict $v]} {set num 0}; break}
	    }
	    set cmp [expr {$num ? {-real} : {-dictionary}}]
	    set rows [lsort {*}$cmp -index $ix -[expr {$ogftd(sortdir) eq "descending" ? "decreasing" : "increasing"}] $rows]
	}
    }
    set ogftd(view,cols) $cols; set ogftd(view,rows) $rows
    set t [join $cols \t]
    foreach r $rows {append t \n [join $r \t]}
    set ogftd(view,tsv) $t
}

# galaxies in common-column form (for the All view)
proc ::ogf::td::_gal_view_rows {rv} {
    upvar 1 $rv rows
    global catpanel
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return
    set lines [split $catpanel(alldata) \n]
    set h [split [lindex $lines 0] \t]
    foreach {nm var} {NUMBER cn X_IMAGE cx Y_IMAGE cy ALPHA_J2000 ca DELTA_J2000 cd MAG_AUTO cm} {set $var [lsearch -exact $h $nm]}
    set i 0
    foreach l [lrange $lines 1 end] {
	incr i
	if {[string trim $l] eq {}} continue
	set f [split $l \t]
	set n $i
	if {$cn >= 0} {set nv [string trim [lindex $f $cn]]; if {$nv ne {}} {set n $nv}}
	set o [list $n galaxy]
	foreach c [list $cx $cy $ca $cd $cm] {lappend o [expr {$c >= 0 ? [string trim [lindex $f $c]] : {}}]}
	lappend rows $o
    }
}

# the TSV the table shows now (the link code, hover cache, save ... read this instead of catpanel(alldata))
proc OGFViewTSV {} {
    global catpanel ogftd
    if {[info exists ogftd(kind)] && $ogftd(kind) ne "galaxies"} {
	::ogf::td::build
	return $ogftd(view,tsv)
    }
    if {[info exists catpanel(alldata)]} {return $catpanel(alldata)}
    return {}
}

# ------------------------------------------------------------------ review state of non-galaxy rows
# Moving / transient / detection rows are generated views, so their accept / reject / uncertain decisions cannot live in catalog cells.
# They are kept here, keyed by the row key (M5, T7, D12) together with a signature of the row (its ra / dec strings): when a new run
# replaces the rows of a kind, only the decisions whose key AND signature still match are kept (review_prune); the rest are dropped,
# because a tracklet number of a new run is a different object.  The decisions appear as the columns REVIEW REVIEW_NOTE REVIEW_TIME at
# the right end of the view (and so in Save, the REVIEW filter and the report) once at least one decision exists.
#   ::ogf::td::review_set KEYS STATUS ?NOTE SETNOTE TIME?   -> number of rows changed (STATUS accept|reject|uncertain|clear)
#   ::ogf::td::review_get KEY -> {status note time} or {}      ::ogf::td::review_count -> number of stored decisions
#   ::ogf::td::review_export FILE / review_import FILE         TSV  key kind ra dec REVIEW REVIEW_NOTE REVIEW_TIME (matched by key + position)
set ::ogf::td::reviewcols {REVIEW REVIEW_NOTE REVIEW_TIME}
proc ::ogf::td::_sig {raw} {
    set ra [expr {[dict exists $raw ra] ? [dict get $raw ra] : {}}]
    set de [expr {[dict exists $raw dec] ? [dict get $raw dec] : {}}]
    return "$ra $de"
}
proc ::ogf::td::review_count {} {global ogftd; return [expr {[info exists ogftd(review)] ? [dict size $ogftd(review)] : 0}]}
proc ::ogf::td::review_get {key} {
    global ogftd
    if {[info exists ogftd(review)] && [dict exists $ogftd(review) $key]} {
	set d [dict get $ogftd(review) $key]
	return [list [dict get $d status] [dict get $d note] [dict get $d time]]
    }
    return {}
}
proc ::ogf::td::review_set {keys status {note {}} {setnote 0} {time {}}} {
    global ogftd
    if {![info exists ogftd(review)]} {set ogftd(review) [dict create]}
    if {$time eq {}} {set time [clock format [clock seconds] -format "%Y-%m-%d %H:%M:%S"]}
    set n 0
    foreach key $keys {
	set raw [row_of $key]
	if {$raw eq {}} continue
	incr n
	if {$status eq "clear"} {dict unset ogftd(review) $key; continue}
	set old [expr {[dict exists $ogftd(review) $key] ? [dict get $ogftd(review) $key] : {}}]
	set on [expr {[dict exists $old note] ? [dict get $old note] : {}}]
	dict set ogftd(review) $key [dict create status $status note [expr {$setnote ? $note : $on}] time $time sig [_sig $raw]]
    }
    set ogftd(stale) 1
    return $n
}
proc ::ogf::td::review_note {keys note} {
    global ogftd
    if {![info exists ogftd(review)]} {set ogftd(review) [dict create]}
    set n 0
    foreach key $keys {
	set raw [row_of $key]
	if {$raw eq {}} continue
	incr n
	if {[dict exists $ogftd(review) $key]} {
	    dict set ogftd(review) $key note $note
	} else {
	    dict set ogftd(review) $key [dict create status {} note $note time {} sig [_sig $raw]]
	}
    }
    set ogftd(stale) 1
    return $n
}
# after the rows of KIND were replaced: keep only decisions of that kind whose row is still there unchanged
proc ::ogf::td::review_prune {kind} {
    global ogftd
    if {![info exists ogftd(review)] || [dict size $ogftd(review)] == 0} {return 0}
    set pre [dict get $ogftd(def,$kind) prefix]
    set drop {}
    dict for {key d} $ogftd(review) {
	if {![string match ${pre}* $key] || ![string is integer -strict [string range $key [string length $pre] end]]} continue
	set raw [row_of $key]
	if {$raw eq {} || [_sig $raw] ne [dict get $d sig]} {lappend drop $key}
    }
    foreach k $drop {dict unset ogftd(review) $k}
    if {[llength $drop]} {::ogf::log INFO "review: dropped [llength $drop] decision(s) of kind $kind whose rows changed in the new run"}
    return [llength $drop]
}
# add the REVIEW columns to a view {cols rows}; galaxy rows of the All view take the decisions of the galaxy catalog
proc ::ogf::td::_review_view {cols rows} {
    global ogftd
    set galrev {}
    set ki [lsearch -exact $cols kind]
    if {$ki >= 0 && [::ogf::cat::has] && [lsearch -exact [::ogf::cat::columns] REVIEW] >= 0} {
	foreach d [::ogf::cat::rows] {
	    if {[dict exists $d NUMBER]} {dict set galrev [dict get $d NUMBER] $d}
	}
    }
    set need [expr {[review_count] > 0 || [dict size $galrev] > 0 && $ki >= 0 && [lsearch -exact [lmap r $rows {lindex $r $ki}] galaxy] >= 0}]
    if {!$need && [::ogf::cat::filter_active] && [dict exists [::ogf::cat::filters] REVIEW]} {set need 1}
    if {!$need} {return [list $cols $rows]}
    set out {}
    foreach r $rows {
	set key [lindex $r 0]
	set v [review_get $key]
	if {$v eq {} && [dict exists $galrev $key]} {
	    set g [dict get $galrev $key]
	    set v [list [expr {[dict exists $g REVIEW] ? [dict get $g REVIEW] : {}}] [expr {[dict exists $g REVIEW_NOTE] ? [dict get $g REVIEW_NOTE] : {}}] \
		[expr {[dict exists $g REVIEW_TIME] ? [dict get $g REVIEW_TIME] : {}}]]
	}
	if {$v eq {}} {set v {{} {} {}}}
	lappend out [concat $r $v]
    }
    return [list [concat $cols $::ogf::td::reviewcols] $out]
}
proc ::ogf::td::review_export {fn} {
    global ogftd
    set lines [list [join {key kind ra dec REVIEW REVIEW_NOTE REVIEW_TIME} \t]]
    if {[info exists ogftd(review)]} {
	dict for {key d} $ogftd(review) {
	    lassign [split [dict get $d sig] { }] ra de
	    set raw [row_of $key]
	    lappend lines [join [list $key [expr {$raw eq {} ? {} : [dict get $raw kind]}] $ra $de [dict get $d status] [dict get $d note] [dict get $d time]] \t]
	}
    }
    set fd [open $fn w]; fconfigure $fd -encoding utf-8; puts $fd [join $lines \n]; close $fd
    return [expr {[llength $lines] - 1}]
}
# rows of a review_export file are applied where key and position match the rows now loaded; returns {applied skipped}
proc ::ogf::td::review_import {fn} {
    global ogftd
    set fd [open $fn r]; fconfigure $fd -encoding utf-8; set txt [read $fd]; close $fd
    if {![info exists ogftd(review)]} {set ogftd(review) [dict create]}
    set ok 0; set skip 0
    foreach l [lrange [split $txt \n] 1 end] {
	if {[string trim $l] eq {}} continue
	lassign [split $l \t] key kind ra de st note time
	set raw [row_of $key]
	if {$raw eq {} || [_sig $raw] ne "$ra $de" || $st ni {accept reject uncertain {}}} {incr skip; continue}
	dict set ogftd(review) $key [dict create status $st note $note time $time sig [_sig $raw]]
	incr ok
    }
    set ogftd(stale) 1
    if {[active]} {refresh}
    return [list $ok $skip]
}

# ------------------------------------------------------------------ table output
proc ::ogf::td::fill {} {
    global catpanel ogftd
    build
    set cols $ogftd(view,cols)
    set nc [llength $cols]
    ::ogf::cat::table_begin
    ::ogf::cat::table_put 0 $cols
    set catpanel(cache,dirty) 1
    set pat [expr {[info exists catpanel(search_var)] ? $catpanel(search_var) : {}}]
    # column-value filters (::ogf::cat::filter_set, e.g. the review filter) apply to these rows too
    set fspecs {}
    if {[::ogf::cat::filter_active]} {set fspecs [::ogf::cat::filter_specs $cols]}
    set row 1
    foreach r $ogftd(view,rows) {
	if {$pat ne {} && ![string match -nocase "*${pat}*" [join $r \t]]} continue
	if {[llength $fspecs] && ![::ogf::cat::filter_row_ok $fspecs $r]} continue
	::ogf::cat::table_put $row $r
	incr row
    }
    ::ogf::cat::table_end $nc $row -state disabled
    catch {::ogf::cat::_table_filled}
    set n [expr {$row-1}]
    set tot [llength $ogftd(view,rows)]
    set why {}
    if {$pat ne {}} {lappend why "filter '$pat'"}
    if {[llength $fspecs]} {lappend why [::ogf::cat::filter_text]}
    set catpanel(status) "[expr {$ogftd(kind) eq {all} ? {All kinds} : [label $ogftd(kind)]}]: $n[expr {[llength $why] ? " of $tot ([join $why {; }])" : {}}] rows"
}

proc ::ogf::td::refresh {} {
    global ogftd
    set ogftd(stale) 1
    if {![active]} return
    fill
    draw_markers
}

proc ::ogf::td::_schedule {} {
    global ogftd
    if {$ogftd(pending) ne {}} return
    set ogftd(pending) [after idle ::ogf::td::_run_pending]
}
proc ::ogf::td::_run_pending {} {
    global ogftd
    set ogftd(pending) {}
    refresh
}

# switch the kind filter
proc ::ogf::td::show {kind} {
    global ogftd catpanel current
    if {$kind ne "galaxies" && $kind ne "all" && ![info exists ogftd(raw,$kind)]} {
	set catpanel(status) "No [string tolower [label $kind]] yet: run the corresponding Time-domain step first"
	set ogftd(kind) $ogftd(kind)
	set ogftd(radio) $ogftd(kind)
	return
    }
    set old $ogftd(kind)
    set ogftd(kind) $kind
    set ogftd(radio) $kind
    set ogftd(stale) 1
    set ogftd(sortcol) {}; set ogftd(sortdir) {}
    catch {CatalogPanelClearSelection}
    if {$kind eq "galaxies"} {
	delete_markers
	if {[info exists catpanel(alldata)] && $catpanel(alldata) ne {}} {
	    CatalogPanelFilterReload
	} else {
	    CatalogPanelClearTable
	}
	if {[info exists catpanel(markall,on)] && $catpanel(markall,on)} {catch {CatalogPanelCreateAllMarkers}}
    } else {
	fill
	if {$kind eq "all"} {
	    if {[have_galaxies]} {catch {CatalogPanelCreateAllMarkers}}
	} else {
	    # galaxy markers are hidden in the single-kind views (they come back with Galaxies / All)
	    foreach fr [_frames] {catch {$fr marker catalog sextract_all delete}}
	}
	draw_markers
	set n [llength [lsearch -all -not -exact [list {*}$ogftd(view,rows)] __none]]
	set catpanel(status) "[expr {$kind eq {all} ? {All kinds} : [label $kind]}]: $n rows"
    }
    OGFCacheDirtyAll
}

# a galaxy-table producer (CatalogPanelLoadTSV / Clear) took over the table: back to the galaxy view
proc OGFTDGalaxyLoaded {} {
    global ogftd
    if {![info exists ogftd(kind)]} return
    set ogftd(galstale) 1
    set ogftd(stale) 1
    if {$ogftd(kind) ne "galaxies"} {
	set ogftd(kind) galaxies
	set ogftd(radio) galaxies
	::ogf::td::delete_markers
    }
}
# kind column in the galaxy view (table only; catpanel(alldata) is untouched)
proc OGFTDAppendKindColumn {ncols nrows} {
    global catpanel ogftd
    if {![::ogf::td::any_data]} return
    set nc [expr {$ncols+1}]
    ::ogf::cat::cell_set 0 $nc kind
    for {set r 1} {$r < $nrows} {incr r} {::ogf::cat::cell_set $r $nc galaxy}
    $catpanel(tbl) configure -cols $nc
}

# alldata changed (new extraction, frame switch, merge ...): host / overlap columns must be recomputed
proc OGFTDAlldataChanged {args} {
    global ogftd
    set ogftd(galstale) 1
    set ogftd(stale) 1
    if {$ogftd(kind) ne "galaxies"} {::ogf::td::_schedule}
}

proc OGFCacheDirtyAll {} {
    set ::catpanel(cache,dirty) 1
}

# ------------------------------------------------------------------ sort / filter / save in a non-galaxy view
proc ::ogf::td::sort {col dir} {
    global ogftd
    set ogftd(sortcol) $col; set ogftd(sortdir) $dir
    set ogftd(stale) 1
    catch {OGFSessLog td.sort note {} -tool internal -title "Sort time-domain table ($ogftd(kind)) by $col $dir" \
	-note "view of the time-domain table only; the catalog of the exported script is not changed"}
    fill
    set ::catpanel(status) "Sorted by $col $dir ($ogftd(kind))"
}

proc ::ogf::td::save {fn} {
    global ogftd catpanel
    build
    set ext [string tolower [file extension $fn]]
    set cols $ogftd(view,cols)
    set outl [list $cols {*}$ogftd(view,rows)]
    set lines {}
    foreach r $outl {
	if {$ext eq ".csv"} {
	    lappend lines [join [lmap f $r {expr {[string match *,* $f] || [string match *\"* $f] ? "\"[string map {\" \"\"} $f]\"" : $f}}] ,]
	} else {
	    lappend lines [join $r \t]
	}
    }
    if {[catch {set fd [open $fn w]; puts -nonewline $fd [join $lines \n]; close $fd} err]} {
	set catpanel(status) "Save error: $err"
	return
    }
    catch {OGFSessLog td.save note {} -tool internal -title "Save time-domain table ($ogftd(kind)) as [file tail $fn]" \
	-note "time-domain view saved by hand ([file tail $fn]); the galaxy catalog of the exported script is unaffected"}
    set catpanel(status) "Saved [llength $ogftd(view,rows)] [string tolower [label $ogftd(kind)]] rows to [file tail $fn]"
}

# ------------------------------------------------------------------ selection info
proc ::ogf::td::selinfo {row} {
    global catpanel ogftd
    if {![active]} {return 0}
    set key [OGFNumberOfRow $row]
    if {$key eq {}} {return 0}
    set raw [row_of $key]
    if {$raw eq {}} {return 0}
    set k [dict get $raw kind]
    set g {}
    for {set c 1} {$c <= [::ogf::cat::table_ncols]} {incr c} {
	if {[::ogf::cat::cell_exists 0 $c] && [::ogf::cat::cell_exists $row $c]} {dict set g [::ogf::cat::cell 0 $c] [::ogf::cat::cell $row $c]}
    }
    proc _g {g k} {expr {[dict exists $g $k] && [dict get $g $k] ne {} ? [dict get $g $k] : "-"}}
    set l1 [format "%s %s   x,y = %s, %s" $key [string tolower [label $k]] [_g $g X_IMAGE] [_g $g Y_IMAGE]]
    set l2 [format "RA,Dec = %s, %s   mag = %s" [_g $g ALPHA_J2000] [_g $g DELTA_J2000] [_g $g MAG_AUTO]]
    set l3 [OGFTDKindSummary $k [dict merge $raw $g]]
    set catpanel(sel,text) "$l1\n$l2\n$l3"
    return 1
}
proc OGFTDKindSummary {k raw} {
    global ogftd
    set parts {}
    foreach c [dict get $ogftd(def,$k) columns] {
	if {[dict exists $raw $c] && [dict get $raw $c] ne {}} {lappend parts "$c=[::ogf::td::_fmt [dict get $raw $c]]"}
    }
    set s [join $parts "  "]
    if {[string length $s] > 84} {set s "[string range $s 0 80]..."}
    return $s
}

# ------------------------------------------------------------------ markers
proc ::ogf::td::_frames {} {
    global ds9
    set out {}
    foreach fr $ds9(frames) {
	if {![catch {$fr has fits} h] && $h} {lappend out $fr}
    }
    return $out
}

proc ::ogf::td::delete_markers {{frames {}}} {
    if {$frames eq {}} {set frames [_frames]}
    foreach fr $frames {catch {$fr marker catalog ogf_td delete}}
}

proc ::ogf::td::_regions {} {
    global ogftd
    set vis [visible_kinds]
    set batches {}
    set reg "fk5\n"
    set n 0
    foreach k $vis {
	set def $ogftd(def,$k)
	set pre [dict get $def prefix]; set col [dict get $def color]; set pt [dict get $def point]
	lassign $ogftd(raw,$k) cols rows
	set ci [lsearch -exact $cols id]; set cr [lsearch -exact $cols ra]; set cd [lsearch -exact $cols dec]
	set cp [lsearch -exact $cols _pos]
	foreach r $rows {
	    set id [lindex $r $ci]
	    set key $pre$id
	    set pos [list [lindex $r $cr] [lindex $r $cd]]
	    if {$cp >= 0 && [lindex $r $cp] ne {}} {set pos [lindex $r $cp]}
	    foreach {ra de} $pos {
		if {![string is double -strict $ra] || ![string is double -strict $de]} continue
		append reg "point($ra,$de) # $pt color=$col width=2 tag={ogf_td} tag={ogf_td.$key} select=0 edit=0 move=0 rotate=0 delete=1\n"
		if {[incr n] % 500 == 0} {lappend batches $reg; set reg "fk5\n"}
	    }
	}
    }
    if {$reg ne "fk5\n"} {lappend batches $reg}
    return $batches
}

proc ::ogf::td::draw_markers {} {
    global ogftd_reg
    delete_markers
    if {![active]} return
    set batches [_regions]
    foreach fr [_frames] {
	foreach b $batches {
	    set ogftd_reg $b
	    catch {$fr marker catalog command ds9 var ogftd_reg}
	}
    }
}
proc ::ogf::td::draw_frame {fr} {
    global ogftd_reg
    if {![active]} return
    catch {$fr marker catalog ogf_td delete}
    foreach b [_regions] {
	set ogftd_reg $b
	catch {$fr marker catalog command ds9 var ogftd_reg}
    }
}

# key of a clicked time-domain marker from its tag list ("" when it is not one)
proc OGFTDKeyFromTags {tags} {
    foreach t $tags {
	if {[string match "ogf_td.*" $t]} {return [string range $t 7 end]}
    }
    return {}
}

# frame switched: markers must exist in the new frame (frames loaded later), the host columns follow the catalog
proc OGFTDFrameChanged {old new} {
    global ogftd
    if {![info exists ogftd(kind)] || $ogftd(kind) eq "galaxies"} return
    ::ogf::td::_schedule
}

# the kind we were in before a frame switch is restored afterwards (RestoreFrameState reloads the galaxy table)
proc OGFTDKeep {} {global ogftd; return [expr {[info exists ogftd(kind)] ? $ogftd(kind) : {galaxies}}]}
proc OGFTDRestoreKind {k} {
    if {$k eq "galaxies"} return
    ::ogf::td::show $k
}

# ------------------------------------------------------------------ controls (Time-domain tab)
proc OGFUITimeDomainControls {parent} {
    global ogftd ogfui
    set f [ttk::frame $parent.td]
    ttk::label $f.l -text "Show:"
    pack $f.l -side left -padx {8 2}
    pack $f -side left -padx 2
    set ogftd(controls) $f
    OGFTDControlsRefresh
}
proc OGFTDControlsRefresh {} {
    global ogftd
    if {![info exists ogftd(controls)] || ![winfo exists $ogftd(controls)]} return
    set f $ogftd(controls)
    foreach w [winfo children $f] {if {$w ne "$f.l"} {destroy $w}}
    set choices {galaxies}
    foreach k $ogftd(kinds) {if {$k ne "galaxies"} {lappend choices $k}}
    lappend choices all
    # the detection kind is only offered once it has data; the others always
    foreach k $choices {
	if {$k ni {galaxies all} && ![dict get $ogftd(def,$k) all] && ![info exists ogftd(raw,$k)]} continue
	set lab [expr {$k eq "all" ? "All" : [::ogf::td::label $k]}]
	ttk::radiobutton $f.k_$k -text $lab -variable ogftd(radio) -value $k -style OGFTab.Toolbutton \
	    -command [list ::ogf::td::show $k]
	pack $f.k_$k -side left -padx 1
    }
}

# ------------------------------------------------------------------ selection of a non-galaxy row
proc ::ogf::td::on_select {proc} {global ogftd; if {$proc ni $ogftd(selcbs)} {lappend ogftd(selcbs) $proc}}

# called by CatalogPanelGotoSource for every row selection; returns 1 when the row was a time-domain row
proc ::ogf::td::goto {row {pan 1}} {
    global ogftd catpanel current
    if {![active]} {return 0}
    set key [OGFNumberOfRow $row]
    if {$key eq {}} {return 0}
    set raw [row_of $key]
    if {$raw eq {}} {return 0}
    set g [dict create]
    for {set c 1} {$c <= [::ogf::cat::table_ncols]} {incr c} {
	if {[::ogf::cat::cell_exists 0 $c] && [::ogf::cat::cell_exists $row $c]} {dict set g [::ogf::cat::cell 0 $c] [::ogf::cat::cell $row $c]}
    }
    set full [dict merge $raw $g]
    set ra [dict get $g ALPHA_J2000]; set de [dict get $g DELTA_J2000]
    clear_selection_markers
    if {[string is double -strict $ra] && [string is double -strict $de]} {
	global ogftd_sel_reg
	set ogftd_sel_reg "fk5\npoint($ra,$de) # point=cross 18 color=cyan width=2 tag={sextract_sel} select=0 edit=0 move=0 rotate=0 delete=1\ncircle($ra,$de,1.2\") # color=cyan width=1 tag={sextract_sel} select=0 edit=0 move=0 rotate=0 delete=1\n"
	foreach fr [_frames] {
	    catch {$fr marker catalog command ds9 var ogftd_sel_reg}
	    lappend ogftd(selframes) $fr
	}
	if {$pan && $current(frame) ne {}} {
	    catch {$current(frame) pan to wcs fk5 degrees $ra $de; UpdatePan $current(frame)}
	}
    }
    foreach cb $ogftd(selcbs) {catch {$cb $key $full}}
    return 1
}
proc ::ogf::td::clear_selection_markers {} {
    global ogftd
    foreach fr [lsort -unique $ogftd(selframes)] {catch {$fr marker catalog sextract_sel delete}}
    set ogftd(selframes) {}
}
proc OGFTDClearSel {} {::ogf::td::clear_selection_markers}

# helpers used when the view goes back to the galaxy catalog
proc CatalogPanelFilterReload {} {CatalogPanelFilter}
proc CatalogPanelClearTable {} {
    global catpanel
    global $catpanel(tbldb)
    $catpanel(tbl) configure -variable {}
    unset -nocomplain $catpanel(tbldb)
    $catpanel(tbl) configure -variable $catpanel(tbldb) -cols 19 -rows 20
}
