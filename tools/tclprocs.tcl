# usage: tclsh tclprocs.tcl file  -> "name start end" (1-based, start includes preceding # comment block)
set fd [open [lindex $argv 0]]; set L [split [read $fd] \n]; close $fd
set n [llength $L]; set i 0
while {$i < $n} {
    set l [lindex $L $i]
    if {[regexp {^proc\s+(\S+)\s} $l -> name]} {
        set buf $l; set j $i
        while {![info complete $buf]} {incr j; if {$j>=$n} {puts stderr "unterminated $name"; exit 1}; append buf \n [lindex $L $j]}
        set k [expr {$i-1}]
        while {$k >= 0 && [string index [lindex $L $k] 0] eq "#"} {incr k -1}
        puts "$name [expr {$k+2}] [expr {$j+1}]"
        set i [expr {$j+1}]
    } else {incr i}
}
