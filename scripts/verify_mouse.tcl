# REAL mouse / keyboard events (xdotool, XTEST) on the catalog panel widgets that the other GUI checks drive through procs:
# tab clicks, chip menus (post + choose an entry), table row click, Ctrl-click multi-select, header click sort, wheel scroll,
# the Tools menu "Detach Catalog Panel", and the review menu entries.  (verify_click_xevent.tcl covers clicks on the image markers.)
#   rm -rf ~/ds9.auto ~/ds9.auto.dir; DISPLAY=:77 OGF_MOUSE_OUT=/tmp/mouse.txt OGF_MOUSE_DIR=/tmp/mousedir bin/ds9 /workspace/fits/m51.fits -geometry 1300x950 -source scripts/verify_mouse.tcl
# Needs xdotool and an X server (no window manager is needed).  Every event travels X server -> Tk -> widget binding -> command.
set ::fh [open $::env(OGF_MOUSE_OUT) w]; set ::nf 0
file mkdir $::env(OGF_MOUSE_DIR)
proc R {tag ok {d {}}} {puts $::fh "[expr {$ok?{PASS}:{FAIL}}] $tag $d"; flush $::fh; if {!$ok} {incr ::nf}}
proc bgerror {m} {puts $::fh "BGERROR $m"; flush $::fh}
rename tk_messageBox ::orig_mb
proc tk_messageBox {args} {puts $::fh "MESSAGEBOX $args"; flush $::fh; return yes}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}
proc geom {} {
    update idletasks; update
    set tf [winfo parent $::catpanel(tbl)]
    return "[winfo rooty $tf] [winfo height $tf] [winfo height $::catpanel(infoarea)]"
}
proc rootxy {w x y} {return [list [expr {[winfo rootx $w] + int($x)}] [expr {[winfo rooty $w] + int($y)}]]}
proc mclick {rx ry {btn 1}} {
    exec xdotool mousemove --sync $rx $ry
    wait_idle 120
    exec xdotool click $btn
    wait_idle 300
}
proc wclick {w x y {btn 1}} {lassign [rootxy $w $x $y] rx ry; mclick $rx $ry $btn}
# centre of table cell r,c (widget coordinates)
proc cellxy {r c} {
    set bb [$::catpanel(tbl) bbox $r,$c]
    return [list [expr {[lindex $bb 0] + [lindex $bb 2]/2}] [expr {[lindex $bb 1] + [lindex $bb 3]/2}]]
}
proc tclick {r c {btn 1}} {lassign [cellxy $r $c] x y; wclick $::catpanel(tbl) $x $y $btn}
# post the chip menu / menubutton MB with a real click, then click entry INDEX of menu M with a real click
proc menu_choose {mb m index} {
    wclick $mb 6 6
    if {![winfo ismapped $m]} {return 0}
    update idletasks
    set y [$m yposition $index]
    set h [expr {[winfo reqheight $m]/([$m index end]+1)}]
    lassign [rootxy $m 20 [expr {$y + 6}]] rx ry
    mclick $rx $ry
    return 1
}
proc entry_index {m label} {
    set n [$m index end]
    for {set i 0} {$i <= $n} {incr i} {if {![catch {$m entrycget $i -label} l] && $l eq $label} {return $i}}
    return -1
}
proc run {} {
    global catpanel ds9 ogfui
    wait_idle 500
    if {[catch {exec xdotool version} v]} {R xdotool 0 $v; puts $::fh "SUMMARY failures=$::nf"; close $::fh; exit 1}
    R xdotool 1 $v
    set g0 [geom]
    R geometry_start [expr {$g0 eq "181 769 154"}] $g0
    # ---- real click on a notebook tab
    set nb $ogfui(nb)
    set want Classify
    set i [lsearch -exact $::ogf::tabs $want]
    set bb [$nb identify element 0 0]
    # the tab's bounding box: scan x for the tab with index i
    set tx {}
    for {set x 5} {$x < [winfo width $nb]} {incr x 6} {
        if {[catch {$nb index @$x,6} k]} continue
        if {$k == $i} {set tx $x; break}
    }
    R tab_found [expr {$tx ne {}}] "x=$tx"
    wclick $nb $tx 6
    R tab_click [expr {$ogfui(tab) eq $want}] "tab=$ogfui(tab)"
    wclick $nb 8 6   ;# first tab (Detect)
    R tab_click_back [expr {$ogfui(tab) eq {Detect}}] "tab=$ogfui(tab)"
    # ---- chip menu: Detect > SExtractor chip, real click on its menubutton and on a menu entry.  Use the Review plugin instead of a
    # long-running step: Results tab, "Review" chip, entry "Accept Selected" (needs a catalog -> extract first)
    CatalogPanelExtract
    set t0 [clock milliseconds]; while {![::ogf::cat::has] && [clock milliseconds]-$t0 < 90000} {update; after 100}
    wait_idle 500
    R extracted [expr {[::ogf::cat::nrows] > 50}] [::ogf::cat::nrows]
    set nums [::ogf::cat::values NUMBER]
    # ---- real click on a table row selects it (CatalogPanelSelectCmd via tktable's <<Select>>)
    tclick 3 1
    wait_idle 400
    R row_click_selects [expr {[::ogf::cat::selection] eq [list [lindex $nums 2]]}] "sel=[::ogf::cat::selection] want=[lindex $nums 2]"
    set txt [::ogf::cat::get sel,text {}]
    R row_click_summary [string match "Source #[lindex $nums 2]*" $txt] [string range $txt 0 60]
    # ---- Ctrl-click in the TABLE: measured behaviour is that it selects that one row (multi-selection is done by Ctrl-clicking markers on
    # the image, tested in verify_click_xevent.tcl; the table itself has no multi-select binding)
    exec xdotool keydown ctrl; wait_idle 150
    tclick 6 1
    exec xdotool keyup ctrl; wait_idle 300
    R ctrl_click_in_table_selects_row [expr {[::ogf::cat::selection] eq [list [lindex $nums 5]]}] "sel=[::ogf::cat::selection] want=[lindex $nums 5]"
    # ---- Results tab, Review chip, "Accept Selected" by real clicks on the menubutton and the entry
    OGFUIShowTab Results; wait_idle 300
    set m $ogfui(menu,report)
    set mb [winfo parent $m]
    R review_chip_mapped [winfo ismapped $mb] $mb
    set ix [entry_index $m "Accept Selected"]
    R review_entry_exists [expr {$ix >= 0}] $ix
    ::ogf::cat::select [lindex $nums 4]; wait_idle 300
    R menu_choose [menu_choose $mb $m $ix]
    wait_idle 500
    set rv [lindex [::ogf::cat::rows] 4]
    R review_accept_by_menu [expr {[dict exists $rv REVIEW] && [dict get $rv REVIEW] eq {accept}}] "row5=[expr {[dict exists $rv REVIEW] ? [dict get $rv REVIEW] : {-}}]"
    R menu_unposted [expr {![winfo ismapped $m]}]
    # ---- "Show by Review" submenu (cascade): real click on the cascade entry, then on "Accepted"
    set ix [entry_index $m "Show by Review"]
    wclick $mb 6 6
    lassign [rootxy $m 20 [expr {[$m yposition $ix] + 6}]] rx ry
    mclick $rx $ry
    set sub $m.vfilter
    R cascade_posted [winfo exists $sub] $sub
    if {[winfo exists $sub] && [winfo ismapped $sub]} {
        set j [entry_index $sub "Accepted"]
        lassign [rootxy $sub 20 [expr {[$sub yposition $j] + 6}]] rx ry
        mclick $rx $ry
        wait_idle 400
        R filter_by_menu [expr {[::ogf::cat::filter_active] && [llength [::ogf::cat::shown_numbers]] == 1}] "shown=[::ogf::cat::shown_numbers]"
    } else {
        R filter_by_menu 0 "submenu not mapped"
    }
    catch {::ogf::review::show all}
    catch {foreach w [list $m $m.vfilter] {if {[winfo exists $w]} {$w unpost}}}
    wait_idle 300
    # ---- header click sorts (real click on the MAG_AUTO header cell), twice -> descending
    set c [::ogf::cat::table_col MAG_AUTO]
    tclick 0 $c
    R header_click_sort [expr {[::ogf::cat::get sort,col] eq {MAG_AUTO} && [::ogf::cat::get sort,dir] eq {ascending}}] "sort=[::ogf::cat::get sort,col] [::ogf::cat::get sort,dir]"
    set mags {}
    foreach v [::ogf::cat::values MAG_AUTO] {if {[string is double -strict $v]} {lappend mags $v}}
    set asc 1; set p [lindex $mags 0]; foreach v $mags {if {$v < $p} {set asc 0; break}; set p $v}
    R sorted_ascending $asc
    tclick 0 [::ogf::cat::table_col MAG_AUTO]
    R header_click_toggle [expr {[::ogf::cat::get sort,dir] eq {descending}}] [::ogf::cat::get sort,dir]
    # ---- wheel: the binding is "natural" (Button-4 scrolls the table DOWN by 3 rows, Button-5 back up; layout.tcl)
    set tbl $catpanel(tbl)
    set top0 [lindex [$tbl yview] 0]
    lassign [cellxy 5 2] x y
    wclick $tbl $x $y 4; wclick $tbl $x $y 4
    set top1 [lindex [$tbl yview] 0]
    wclick $tbl $x $y 5; wclick $tbl $x $y 5
    set top2 [lindex [$tbl yview] 0]
    R wheel_scrolls [expr {$top1 > $top0 && $top2 < $top1}] "top: $top0 -> $top1 -> $top2"
    # ---- Tools menu: Detach Catalog Panel (check button) by real clicks, then back
    set tm $ogfui(root).menubar.tools.m
    set tmb $ogfui(root).menubar.tools
    set w0 [winfo width .]
    R tools_menu_click [menu_choose $tmb $tm [entry_index $tm "Detach Catalog Panel"]]
    wait_idle 800
    update idletasks
    R detach_by_menu [expr {$catpanel(detached) == 1 && [winfo width .] == 736}] "detached=$catpanel(detached) width=[winfo width .] (was $w0)"
    R geometry_detached [expr {[geom] eq $g0}] [geom]
    # the detached panel is a toplevel: its Tools menu is reached the same way
    menu_choose $tmb $tm [entry_index $tm "Detach Catalog Panel"]
    wait_idle 800
    R reattach_by_menu [expr {$catpanel(detached) == 0 && [winfo width .] == 1300}] "detached=$catpanel(detached) width=[winfo width .]"
    R geometry_end [expr {[geom] eq $g0}] "[geom] (start $g0)"
    puts $::fh "SUMMARY failures=$::nf"; close $::fh; exit
}
after 3000 {if {[catch run err]} {puts $::fh "ERROR $err $::errorInfo"; close $::fh; exit 1}}
