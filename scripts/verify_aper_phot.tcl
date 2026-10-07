# GUI test of the shared-core aperture photometry steps (plugins/ds10core: aper-phot, aper-series).  Run through scripts/verify_aper_phot.sh.
global current ds9
set ::fh [open $::env(OGF_DS_OUT) w]; set ::nf 0
set ::dir $::env(OGF_DS_DIR)
set ::shots $::env(OGF_DS_SHOTS)
proc R {tag ok {d {}}} {puts $::fh "[expr {$ok?{PASS}:{FAIL}}] $tag $d"; flush $::fh; if {!$ok} {incr ::nf}}
proc I {tag d} {puts $::fh "INFO $tag $d"; flush $::fh}
proc bgerror {m} {puts $::fh "BGERROR $m"; flush $::fh; incr ::nf}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}
proc wait_job {{ms 300000}} {set t [clock milliseconds]; while {[::ogf::job::busy] && [clock milliseconds]-$t < $ms} {update; after 50}; wait_idle 500}
proc steps_since {n} {lmap r [lrange [::ogf::session::steps] $n end] {list [dict get $r step] [dict get $r class]}}
proc run_step {plugin step} {
    set s0 [llength [::ogf::session::steps]]
    set ok [::ogf::step::run $plugin $step]
    wait_job
    return [list $ok [steps_since $s0]]
}
proc shot {name} {
    if {$::shots eq {}} return
    update idletasks; update; wait_idle 400
    file mkdir $::shots
    catch {exec $::env(OGFINDER_PYTHON) scripts/xshot.py $::env(OGF_DS_DISPLAY) [file join $::shots $name.png]} m
    puts $::fh "SHOT $name $m"; flush $::fh
}
proc nshape {shape} {
    global current
    set txt {}; catch {set txt [$current(frame) marker list ds9 image fk5 degrees 0]}
    return [regexp -all -line "^\\s*${shape}\\(" $txt]
}
proc load_points {pts} {
    global current
    set f [file join $::dir picks_in.reg]; set fd [open $f w]; puts $fd "image"
    foreach p $pts {lassign $p x y; puts $fd "point($x,$y) # point=x"}
    close $fd
    MarkerLoadFile $f $current(frame) ds9 image fk5
}
proc readtsv {f} {set fd [open $f r]; set l [split [string trim [read $fd]] \n]; close $fd; return $l}
# image pixel -> canvas point (ds9 frame coordinates are canvas-wide) -> root screen point for xdotool
proc img2canvas {fr ix iy} {
    global ds9
    set c [$ds9(canvas) bbox $fr]
    set x0 [expr {([lindex $c 0]+[lindex $c 2])/2 - 30}]; set y0 [expr {([lindex $c 1]+[lindex $c 3])/2 - 30}]
    lassign [$fr get coordinates $x0 $y0 image] ax ay
    lassign [$fr get coordinates [expr {$x0+60}] $y0 image] bx by
    lassign [$fr get coordinates $x0 [expr {$y0+60}] image] cx cy
    return [list [expr {$x0 + ($ix-$ax)*60.0/($bx-$ax)}] [expr {$y0 + ($iy-$ay)*60.0/($cy-$ay)}]]
}
proc xclick {cx cy} {
    global ds9
    set c $ds9(canvas)
    set rx [expr {int([winfo rootx $c] + $cx - [$c canvasx 0])}]
    set ry [expr {int([winfo rooty $c] + $cy - [$c canvasy 0])}]
    exec xdotool mousemove --sync $rx $ry; wait_idle 120
    exec xdotool click 1; wait_idle 300
}

proc run {} {
    global current ds9
    wait_idle 800
    set stars {}
    set fd [open [file join $::dir stars.txt] r]; foreach ln [split [string trim [read $fd]] \n] {lappend stars $ln}; close $fd
    R steps_registered [expr {[lsearch -index 0 [lmap s [dict get [::ogf::reg::get ds10core] steps] {list [dict get $s id]}] aper-phot] >= 0}]
    R param_defaults [expr {[::ogf::params::get ds10core aper-radii] eq "1,1.5,2,3" && [::ogf::params::get ds10core aper-stars] eq "regions"}]
    ::ogf::params::put ds10core aper-zeropoint auto
    # ---- the parameter dialog of the step (group "Aperture photometry")
    set dlg [OGFParamDialog ds10core -group {Aperture photometry} -title {Aperture photometry (shared core) - Settings}]
    R param_dialog [expr {$dlg ne {} && [winfo exists $dlg]}] $dlg
    if {$dlg ne {}} {catch {foreach t [$dlg.body.nb tabs] {if {[$dlg.body.nb tab $t -text] eq {Aperture photometry}} {$dlg.body.nb select $t}}}; catch {wm geometry $dlg 640x720+780+60}; wait_idle 800; shot aper_desktop_params_1440x900; catch {destroy $dlg}}
    # ---- pick 6 stars: real X clicks in region edit mode (xdotool), else point regions created in the frame
    catch {$current(frame) marker delete all}
    set fr $current(frame)
    set npick 6; set how none
    if {![catch {exec xdotool version}]} {
	catch {set current(mode) region; UpdateEditMenu}
	catch {set ::marker(shape) circle}
	foreach s [lrange $stars 0 [expr {$npick-1}]] {lassign $s x y; lassign [img2canvas $fr [expr {$x+1.2}] [expr {$y-0.8}]] cx cy; catch {xclick $cx $cy}}
	set got [expr {[nshape circle] + [nshape point]}]
	I xdotool_picks "regions created by real clicks: $got"
	if {$got == $npick} {set how xclick} else {catch {$current(frame) marker delete all}}
    }
    if {$how eq "none"} {
	load_points [lmap s [lrange $stars 0 [expr {$npick-1}]] {lassign $s x y; list [expr {$x+1.2}] [expr {$y-0.8}]}]
	set how marker_create
    }
    R picks_on_image [expr {[nshape circle] + [nshape point] == $npick}] "how=$how"
    # ---- run the step: picks -> aperture photometry (catalogue table + overlay)
    lassign [run_step ds10core aper-phot] ok recs
    set cat [file join [OGFSessWorkDir] aperphot_catalog.tsv]
    set rows {}; if {[file exists $cat]} {set rows [readtsv $cat]}
    set hdr [split [lindex $rows 0] \t]
    R aper_phot_ok [expr {$ok && [llength $rows] == $npick + 1 && [lsearch $hdr MAG] >= 0 && [lsearch $hdr FLAGS] >= 0}] "$recs rows=[expr {[llength $rows]-1}]"
    R aper_phot_recorded [expr {[lindex $recs 0 0] eq "ds10core.aper_phot"}] $recs
    # recentred on the true positions; flux within 3 % of the injected flux (aperture-corrected)
    set ix [lsearch $hdr X_IMAGE]; set iy [lsearch $hdr Y_IMAGE]; set ifl [lsearch $hdr FLUX]; set worst 0.0; set dmax 0.0
    foreach r [lrange $rows 1 end] s [lrange $stars 0 [expr {$npick-1}]] {
	set c [split $r \t]; lassign $s x y f
	set dmax [expr {max($dmax, hypot([lindex $c $ix]-$x, [lindex $c $iy]-$y))}]
	set worst [expr {max($worst, abs([lindex $c $ifl]/$f - 1))}]
    }
    R aper_phot_centroids [expr {$dmax < 0.2}] "max offset [format %.3f $dmax] px"
    R aper_phot_flux [expr {$worst < 0.03}] "worst |flux/true-1| [format %.4f $worst]"
    wait_idle 1000
    R aper_phot_table [expr {[::ogf::cat::nrows] == $npick && [lsearch [::ogf::cat::columns] MAGERR] >= 0}] "rows=[::ogf::cat::nrows]"
    R aper_overlay [expr {[nshape annulus] == $npick && [nshape circle] == 4 * $npick && [nshape point] == 0}] "annuli=[nshape annulus] circles=[nshape circle]"
    set st {}; catch {set st $::catpanel(status)}
    R aper_status [string match "*Aperture photometry*" $st] $st
    catch {destroy .ogftext}; wait_idle 800
    shot aper_desktop_table_1440x900
    lassign [lindex $stars 1] zx zy
    catch {set current(zoom) {4 4}; $current(frame) zoom to 4 4; UpdateZoomMenu}
    catch {$current(frame) pan to image $zx $zy; UpdatePan $current(frame)}
    wait_idle 800
    shot aper_desktop_overlay_1440x900
    catch {ZoomToFit}
    # ---- light curve over the 6-frame sequence: target = star 1 (8 % sinusoid), comparisons 2,3, check 4
    set fl {}; foreach k {1 2 3 4 5} {lappend fl [file join $::dir seq_$k.fits]}
    foreach {k v} [list aper-frames [join $fl ,] aper-target 1 aper-comps 2,3 aper-check 4 aper-track wcs] {::ogf::params::put ds10core $k $v}
    catch {$current(frame) marker delete all}
    load_points [lrange $stars 0 3]
    lassign [run_step ds10core aper-series] ok recs
    set lc [file join [OGFSessWorkDir] aperseries_lightcurve.tsv]
    set lrows {}; if {[file exists $lc]} {set lrows [readtsv $lc]}
    R aper_series_ok [expr {$ok && [llength $lrows] == 7}] "$recs frames=[expr {[llength $lrows]-1}]"
    set h [split [lindex $lrows 0] \t]; set ir [lsearch $h REL_FLUX]; set ic [lsearch $h DMAG_CHECK]
    set rel {}; set chk {}
    foreach r [lrange $lrows 1 end] {set c [split $r \t]; lappend rel [lindex $c $ir]; lappend chk [lindex $c $ic]}
    set dev 0.0
    if {[llength $rel] == 6} {
	set r0 [lindex $rel 0]
	foreach k {0 1 2 3 4 5} v $rel {set want [expr {1 + 0.08*sin(2*3.141592653589793*$k/6)}]; set dev [expr {max($dev, abs($v/$r0 - $want))}]}
    }
    R aper_series_variable [expr {[llength $rel] == 6 && $dev < 0.02}] "max |rel/rel0 - true| [format %.4f $dev]"
    set cm [expr {[llength $chk] ? [tcl::mathop::+ {*}$chk]/[llength $chk] : 0}]; set cr 0.0
    foreach v $chk {set cr [expr {max($cr, abs($v-$cm))}]}
    R aper_series_check_flat [expr {[llength $chk] == 6 && $cr < 0.02}] "check max dev [format %.4f $cr] mag"
    set aav [file join [OGFSessWorkDir] aperseries_aavso.txt]
    set a {}; if {[file exists $aav]} {set fd [open $aav r]; set a [read $fd]; close $fd}
    R aper_series_aavso [string match "#TYPE=EXTENDED*" $a]
    R aper_series_lc_window [winfo exists .ogfaperlc]
    catch {ZoomToFit}
    wait_idle 800
    shot aper_desktop_lightcurve_1440x900
    catch {destroy .ogfaperlc}
    puts $::fh "SUMMARY failures=$::nf"; close $::fh
    exit
}
after 3000 {if {[catch run err]} {puts $::fh "BGERROR $err"; puts $::fh "SUMMARY failures=99"; close $::fh; exit}}
