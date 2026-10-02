# GUI test of spectra "Fit Kinematics" + "Kinematics Viewer" (called by scripts/verify_spectra_kin.sh; synthetic slit object 1, cube object 2, no spectrum for object 3).
global catpanel current ds9
set ::fh [open $::env(OGF_SK_OUT) w]; set ::nf 0
proc R {tag ok {d {}}} {puts $::fh "[expr {$ok?{PASS}:{FAIL}}] $tag $d"; flush $::fh; if {!$ok} {incr ::nf}}
proc bgerror {m} {puts $::fh "BGERROR $m"; flush $::fh; incr ::nf}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}
proc geom {} {
    update idletasks; update
    set tf [winfo parent $::catpanel(tbl)]
    return "[winfo rooty $tf] [winfo height $tf] [winfo height $::catpanel(infoarea)] [winfo width .]x[winfo height .]"
}
proc wait_job {{ms 300000}} {set t [clock milliseconds]; while {[::ogf::job::busy] && [clock milliseconds]-$t < $ms} {update; after 50}; wait_idle 300}
proc readjson {f} {set fd [open $f r]; set t [read $fd]; close $fd; return [::ogf::json::parse $t]}
proc run {} {
    global catpanel current ds9
    R geometry_before [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    set fd [open $::env(OGF_SK_CAT) r]; set data [read $fd]; close $fd
    ::ogf::cat::load_tsv $data kin
    wait_idle 800
    R catalog_loaded [expr {[::ogf::cat::nrows] == 3}] "rows=[::ogf::cat::nrows]"
    ::ogf::params::put spectra spec-dir $::env(OGF_SK_SPEC)
    ::ogf::params::put spectra z-column Z
    ::ogf::params::put spectra cube-xy 1
    ::ogf::params::put spectra kin-inst-fwhm 2.0
    ::ogf::params::put spectra kin-window-kms 600
    ::ogf::params::put spectra kin-inc 0
    set s0 [llength [::ogf::session::steps]]
    set ok [::ogf::step::run spectra kin]
    wait_job
    R kin_ran $ok
    set recs [lrange [::ogf::session::steps] $s0 end]
    R kin_recorded [expr {[llength $recs] == 1 && [dict get [lindex $recs 0] step] eq "analysis.spectra_kin"}] [lmap r $recs {dict get $r step}]
    set cols [::ogf::cat::columns]
    foreach c {SP_KIN_VSINI SP_KIN_VC SP_KIN_PA SP_KIN_INC SP_KIN_SIGMA SP_KIN_RT} {R col_$c [expr {$c in $cols}]}
    set pa [lindex [::ogf::cat::values SP_KIN_PA] 1]; set inc [lindex [::ogf::cat::values SP_KIN_INC] 1]; set vc [lindex [::ogf::cat::values SP_KIN_VC] 1]
    R cube_pa [expr {$pa ne {} && abs($pa-130.0) < 1.5}] "PA=$pa (truth 130)"
    R cube_inc [expr {$inc ne {} && abs($inc-55.0) < 4}] "inc=$inc (truth 55)"
    R cube_vc [expr {$vc ne {} && abs($vc-200.0) < 8}] "vc=$vc (truth 200)"
    set vs [lindex [::ogf::cat::values SP_KIN_VSINI] 0]
    R slit_vsini [expr {$vs ne {} && abs($vs-220*sin(60*acos(-1)/180)) < 5}] "vsini=$vs (truth [format %.1f [expr {220*sin(60*acos(-1)/180)}]])"
    R none_for_obj3 [expr {[lindex [::ogf::cat::values SP_KIN_VSINI] 2] eq {}}]
    set w [file join [OGFSessWorkDir] spectra]
    foreach f {kin_1.png kin_2.png kin_2_maps.fits kin_2_vel.fits spectra_kin.json} {R file_$f [expr {[file exists [file join $w $f]] && [file size [file join $w $f]] > 100}]}
    # viewer
    set ok [::ogf::step::run spectra kinview]
    wait_idle 800
    R viewer_open [expr {[winfo exists .ogfkinview] && [winfo ismapped .ogfkinview]}]
    R viewer_image [expr {[image width ogfkinimg] > 200 && [image height ogfkinimg] > 100}] "[image width ogfkinimg]x[image height ogfkinimg]"
    OGFSpectraKinShow 2
    wait_idle 300
    set txt [.ogfkinview.t get 1.0 end]
    R viewer_text [expr {[string match "*object 2*cube*PA*" $txt]}] [string map {"\n" " | "} [string trim $txt]]
    set n0 [llength $ds9(frames)]
    set n [OGFSpectraKinFrames 2]
    wait_idle 1000
    R viewer_frames [expr {$n == 3 && [llength $ds9(frames)] == $n0 + 3}] "opened $n, frames [llength $ds9(frames)]"
    destroy .ogfkinview
    GotoFrame [lindex $ds9(frames) 0]
    wait_idle 500
    R geometry_after [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    puts $::fh "SUMMARY failures=$::nf"; close $::fh
    exit
}
after 1500 {if {[catch run err]} {puts $::fh "FAIL exception $err $::errorInfo"; puts $::fh "SUMMARY failures=99"; close $::fh; exit}}
