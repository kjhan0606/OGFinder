# GUI test of the ZOGY options of Moving Objects > Difference (item 2): default argv unchanged, parameters map to CLI flags,
# the recorded session step carries exactly the argv that is run.  OGFMovRun is wrapped so no Python is started.
#   rm -rf ~/ds9.auto ~/ds9.auto.dir; DISPLAY=:77 OGF_MOVOPT_OUT=/tmp/movopt.txt bin/ds9 /workspace/fits/m51.fits -geometry 1300x950 -source scripts/verify_moving_options.tcl
set ::fh [open $::env(OGF_MOVOPT_OUT) w]; set ::nf 0
proc R {tag ok {d {}}} {puts $::fh "[expr {$ok?{PASS}:{FAIL}}] $tag $d"; flush $::fh; if {!$ok} {incr ::nf}}
proc bgerror {m} {puts $::fh "BGERROR $m"; flush $::fh}
proc geom {} {
    update idletasks; update
    set tf [winfo parent $::catpanel(tbl)]
    return "[winfo rooty $tf] [winfo height $tf] [winfo height $::catpanel(infoarea)] [winfo width .]x[winfo height .]"
}
proc run {} {
    global ogfmov
    R geometry_before [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    R plugin_enabled [expr {"moving" in [::ogf::reg::ids]}]
    foreach p {source-noise astrom psf-tile template-psf} {R param_$p [expr {![catch {::ogf::params::get moving $p}]}] [catch {::ogf::params::get moving $p}]}
    R defaults [expr {[::ogf::params::get moving source-noise] == 0 && [::ogf::params::get moving astrom] eq "off" && [::ogf::params::get moving psf-tile] == 0 && [::ogf::params::get moving template-psf] eq "target"}]
    R default_args_empty [expr {[OGFMovZogyArgs] eq {}}] [OGFMovZogyArgs]
    # the session recorder gets the argv of OGFMovRun
    rename OGFMovRun OGFMovRun_real
    proc OGFMovRun {step title extra done {network 0} {outputs {}}} {set ::last [list $step $extra]}
    set ogfmov(files) {/tmp/a.fits /tmp/b.fits}; set ogfmov(ra) 150.1; set ogfmov(dec) 2.3
    OGFMovDifference
    lassign $::last st ex
    set base $ex
    R default_argv_unchanged [expr {$st eq "moving.difference" && $ex eq [list --mode difference --workdir [OGFMovWork] --half-pix $ogfmov(halfpix) --snr 5 --files /tmp/a.fits /tmp/b.fits --ra 150.1 --dec 2.3]}] $ex
    ::ogf::params::put moving source-noise 1; ::ogf::params::put moving astrom measure
    ::ogf::params::put moving psf-tile 256; ::ogf::params::put moving template-psf measure
    OGFMovDifference
    lassign $::last st ex
    R option_argv [expr {[lrange $ex end-6 end] eq {--source-noise --astrom measure --psf-tile 256 --template-psf measure}}] $ex
    R prefix_same [expr {[lrange $ex 0 [expr {[llength $base]-1}]] eq $base}]
    # the recorded step (real OGFMovRun_real path) would log this argv; check the recorder accepts and exports it
    set n0 [llength [::ogf::session::steps]]
    set seq [OGFSessLog moving.difference auto [list python ds9_moving.py {*}$ex] -title t -network 0 -outputs {} -tool python]
    set rec [lindex [::ogf::session::steps] end]
    R recorded [expr {[llength [::ogf::session::steps]] == $n0+1 && [dict get $rec step] eq "moving.difference"}]
    R recorded_argv [expr {"--astrom" in [dict get $rec argv] && "--psf-tile" in [dict get $rec argv]}] [dict get $rec argv]
    ::ogf::params::put moving astrom 0.3
    R astrom_value_passed_through [expr {[OGFMovZogyArgs] ne {}}] [OGFMovZogyArgs]
    ::ogf::params::put moving source-noise 0; ::ogf::params::put moving astrom off; ::ogf::params::put moving psf-tile 0; ::ogf::params::put moving template-psf target
    R back_to_default [expr {[OGFMovZogyArgs] eq {}}]
    R geometry_after [expr {[geom] eq "181 769 154 1300x950"}] [geom]
    puts $::fh "SUMMARY failures=$::nf"; close $::fh; exit
}
after 3000 run
