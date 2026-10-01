#  OGFinder core: minimal JSON reader (plugin.json, parameter presets).
#  Objects -> dict, arrays -> list, strings/numbers -> Tcl string, true/false -> 1/0, null -> {}.
#  (No tcllib in the embedded runtime.)  See docs/plugins.md

package provide DS9 1.0

namespace eval ::ogf::json {}

proc ::ogf::json::parse {text} {
    set i 0
    set v [_value $text i]
    _ws $text i
    if {$i < [string length $text]} {error "JSON: trailing characters at offset $i"}
    return $v
}

proc ::ogf::json::_ws {t iv} {
    upvar 1 $iv i
    set n [string length $t]
    while {$i < $n && [string is space [string index $t $i]]} {incr i}
}

proc ::ogf::json::_value {t iv} {
    upvar 1 $iv i
    _ws $t i
    set c [string index $t $i]
    switch -exact -- $c {
	\{ {
	    incr i
	    set d [dict create]
	    _ws $t i
	    if {[string index $t $i] eq "\}"} {incr i; return $d}
	    while 1 {
		_ws $t i
		if {[string index $t $i] ne "\""} {error "JSON: object key expected at offset $i"}
		set k [_string $t i]
		_ws $t i
		if {[string index $t $i] ne ":"} {error "JSON: ':' expected at offset $i"}
		incr i
		dict set d $k [_value $t i]
		_ws $t i
		set c [string index $t $i]
		incr i
		if {$c eq ","} continue
		if {$c eq "\}"} break
		error "JSON: ',' or '\}' expected at offset [expr {$i-1}]"
	    }
	    return $d
	}
	\[ {
	    incr i
	    set l {}
	    _ws $t i
	    if {[string index $t $i] eq "\]"} {incr i; return $l}
	    while 1 {
		lappend l [_value $t i]
		_ws $t i
		set c [string index $t $i]
		incr i
		if {$c eq ","} continue
		if {$c eq "\]"} break
		error "JSON: ',' or '\]' expected at offset [expr {$i-1}]"
	    }
	    return $l
	}
	\" {return [_string $t i]}
	default {
	    if {[regexp -start $i {\A(true|false|null)} $t m]} {
		incr i [string length $m]
		return [dict get {true 1 false 0 null {}} $m]
	    }
	    # Python's json module writes NaN / Infinity (not valid JSON, but our own tools produce them)
	    if {[regexp -start $i {\A(NaN|-?Infinity)} $t m]} {
		incr i [string length $m]
		return [dict get {NaN nan Infinity inf -Infinity -inf} $m]
	    }
	    if {[regexp -start $i {\A-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?} $t m]} {
		incr i [string length $m]
		return $m
	    }
	    error "JSON: unexpected character '$c' at offset $i"
	}
    }
}

proc ::ogf::json::_string {t iv} {
    upvar 1 $iv i
    incr i
    set out {}
    set n [string length $t]
    while {$i < $n} {
	set c [string index $t $i]
	if {$c eq "\""} {incr i; return $out}
	if {$c eq "\\"} {
	    incr i
	    set e [string index $t $i]
	    switch -exact -- $e {
		n {append out \n} t {append out \t} r {append out \r} b {append out \b}
		f {append out \f} / {append out /} \\ {append out \\} \" {append out \"}
		u {
		    scan [string range $t [expr {$i+1}] [expr {$i+4}]] %x cp
		    incr i 4
		    # surrogate pair
		    if {$cp >= 0xD800 && $cp < 0xDC00 && [string range $t [expr {$i+1}] [expr {$i+2}]] eq "\\u"} {
			scan [string range $t [expr {$i+3}] [expr {$i+6}]] %x lo
			set cp [expr {0x10000 + (($cp - 0xD800) << 10) + ($lo - 0xDC00)}]
			incr i 6
		    }
		    append out [format %c $cp]
		}
		default {error "JSON: bad escape \\$e"}
	    }
	    incr i
	    continue
	}
	append out $c
	incr i
    }
    error "JSON: unterminated string"
}

# dict/list/scalar helpers for schema-driven access
proc ::ogf::json::get {d key {default {}}} {
    if {[catch {dict exists $d $key} ex] || !$ex} {return $default}
    return [dict get $d $key]
}

# write a flat dict of strings as pretty JSON (presets)
proc ::ogf::json::write_flat {d} {
    set parts {}
    dict for {k v} $d {lappend parts "  [OGFJStr $k]: [OGFJStr $v]"}
    return "\{\n[join $parts ",\n"]\n\}\n"
}
