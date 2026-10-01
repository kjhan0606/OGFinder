# Unit test of the ::ogf::cat key accessor (get/set/exists/unset/append/lappend/keys/trace/registry) inside a running ds9.
#   OGF_CAT_API_OUT=/tmp/api.txt DISPLAY=:77 bin/ds9 /workspace/fits/m51.fits -source scripts/verify_cat_api.tcl
set ::fh [open $::env(OGF_CAT_API_OUT) w]; set ::nf 0
proc R {tag ok {d {}}} {puts $::fh "[expr {$ok?{PASS}:{FAIL}}] $tag $d"; if {!$ok} {incr ::nf}}
proc run {} {
    global catpanel
    R registry_nonempty [expr {[llength [::ogf::cat::registry]] > 20}] [llength [::ogf::cat::registry]]
    # round trip through the legacy array
    ::ogf::cat::set param,zz-test 7
    R set_visible_in_array [expr {$catpanel(param,zz-test) == 7}]
    set catpanel(param,zz-test) 8
    R get_sees_legacy_write [expr {[::ogf::cat::get param,zz-test] == 8}]
    R get_default [expr {[::ogf::cat::get param,nope 42] == 42}]
    R get_missing_errors [catch {::ogf::cat::get param,nope}]
    R exists [expr {[::ogf::cat::exists param,zz-test] && ![::ogf::cat::exists param,nope]}]
    ::ogf::cat::unset param,zz-test param,nope
    R unset [expr {![info exists catpanel(param,zz-test)]}]
    ::ogf::cat::unset_glob morph,*
    ::ogf::cat::set morph,1 a; ::ogf::cat::set morph,2 b
    R keys_glob [expr {[::ogf::cat::keys morph,*] eq {morph,1 morph,2}}] [::ogf::cat::keys morph,*]
    ::ogf::cat::unset_glob morph,*
    R unset_glob [expr {[::ogf::cat::keys morph,*] eq {}}]
    ::ogf::cat::unset morph,map
    ::ogf::cat::lappend morph,map 3; ::ogf::cat::lappend morph,map 4
    R lappend [expr {$catpanel(morph,map) eq {3 4}}]
    ::ogf::cat::set status {}; ::ogf::cat::append status ab cd
    R append [expr {$catpanel(status) eq {abcd}}]
    ::ogf::cat::unset morph,map
    # trace fires for accessor writes and for legacy writes, with key and value
    set ::seen {}
    proc tr {k v} {lappend ::seen [list $k $v]}
    ::ogf::cat::trace add param,zz-tr tr
    ::ogf::cat::set param,zz-tr 1
    set catpanel(param,zz-tr) 2
    R trace_fires [expr {$::seen eq {{param,zz-tr 1} {param,zz-tr 2}}}] $::seen
    R trace_info [expr {[::ogf::cat::trace info param,zz-tr] eq {tr}}]
    ::ogf::cat::trace remove param,zz-tr tr
    set catpanel(param,zz-tr) 3
    R trace_removed [expr {[llength $::seen] == 2}]
    ::ogf::cat::unset param,zz-tr
    # registry
    R describe_known [expr {[lindex [::ogf::cat::describe psf,param,rl-iterations] 2] eq {star_psf}}]
    R describe_unknown [expr {[::ogf::cat::describe totally,unknown] eq {}}]
    ::ogf::cat::set totally,unknown 1; ::ogf::cat::unset totally,unknown
    R unregistered_warns [expr {[lsearch -glob $::ogf::logbuf "*unregistered key \"totally,unknown\"*"] >= 0}]
    # builtin service procs still work and see the same store
    R existing_has [expr {[::ogf::cat::has] == ([info exists catpanel(alldata)] && $catpanel(alldata) ne {})}]
    # strict mode
    set ::env(OGF_CAT_STRICT) 1
    R strict_errors [catch {::ogf::cat::get totally,unknown2 1}]
    unset ::env(OGF_CAT_STRICT)
    puts $::fh "SUMMARY failures=$::nf"; close $::fh; exit
}
after 3000 {if {[catch run err]} {puts $::fh "ERROR $err $::errorInfo"; close $::fh; exit 1}}
