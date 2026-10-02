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
    # ---- table cell / per-frame accessors (item 4); these need a loaded catalog
    ::ogf::cat::load_tsv "NUMBER\tX_IMAGE\tY_IMAGE\tMAG_AUTO\n1\t10.5\t20.5\t18.2\n2\t30\t40\t19\n3\t5\t6\t20.1\n" apitest
    update
    set db $catpanel(tbldb); global $db
    R table_ncols [expr {[::ogf::cat::table_ncols] == 4}] [::ogf::cat::table_ncols]
    R table_nrows [expr {[::ogf::cat::table_nrows] == 3}] [::ogf::cat::table_nrows]
    R header [expr {[::ogf::cat::header 2] eq {X_IMAGE} && [::ogf::cat::header 99] eq {}}]
    R table_col [expr {[::ogf::cat::table_col MAG_AUTO] == 4 && [::ogf::cat::table_col NOPE] == -1}]
    R cell_matches_array [expr {[::ogf::cat::cell 2 3] eq [set ${db}(2,3)] && [::ogf::cat::cell 2 3] eq {40}}]
    R cell_default [expr {[::ogf::cat::cell 9 9] eq {} && [::ogf::cat::cell 9 9 zz] eq {zz} && ![::ogf::cat::cell_exists 9 9] && [::ogf::cat::cell_exists 1 1]}]
    R row_of [expr {[::ogf::cat::row_of 3] == 3 && [::ogf::cat::row_of 77] == -1}]
    R number_of [expr {[::ogf::cat::number_of 2] eq {2} && [::ogf::cat::number_of 50] eq {}}]
    # filtering keeps the accessors in step with the table
    set catpanel(search_var) 30; CatalogPanelFilter
    R filter_rows [expr {[::ogf::cat::table_nrows] == 1 && [::ogf::cat::row_of 2] == 1 && [::ogf::cat::row_of 1] == -1}] [::ogf::cat::table_nrows]
    set catpanel(search_var) {}; CatalogPanelFilter
    R filter_back [expr {[::ogf::cat::table_nrows] == 3}]
    # begin/put/end round trip: short rows are padded, long rows cut, fields trimmed
    ::ogf::cat::table_begin
    ::ogf::cat::table_put 0 {A B C}
    ::ogf::cat::table_put 1 { 1 2} 3
    ::ogf::cat::table_put 2 {1 2 3 4} 3
    ::ogf::cat::table_end 3 3
    R table_put [expr {[::ogf::cat::cell 1 1] eq {1} && [::ogf::cat::cell 1 3] eq {} && [::ogf::cat::cell 2 3] eq {3} && ![::ogf::cat::cell_exists 2 4]}]
    R table_nrows2 [expr {[::ogf::cat::table_nrows] == 2}]
    ::ogf::cat::table_reset
    R table_reset [expr {[::ogf::cat::table_ncols] == 19 && [::ogf::cat::table_nrows] == 19 && ![::ogf::cat::cell_exists 0 1]}]
    ::ogf::cat::load_tsv [::ogf::cat::tsv] apitest
    R reload_after_reset [expr {[::ogf::cat::table_nrows] == 3}]
    # bind_var / bump
    R bind_var [expr {[::ogf::cat::bind_var plot,logx] eq {::catpanel(plot,logx)}}]
    set b0 [::ogf::cat::get plot,counter 0]
    R bump [expr {[::ogf::cat::bump plot,counter] == $b0 + 1 && [::ogf::cat::bump plot,counter 2] == $b0 + 3}]
    set catpanel(plot,counter) $b0
    # per-frame snapshots
    ::ogf::cat::frame_set zz1 alldata ABC
    R frame_set_array [expr {$::catpanel_fdata(zz1,alldata) eq {ABC}}]
    R frame_get [expr {[::ogf::cat::frame_get zz1 alldata] eq {ABC} && [::ogf::cat::frame_get zz1 nope d] eq {d}}]
    R frame_get_missing_errors [catch {::ogf::cat::frame_get zz1 nope}]
    R frame_exists [expr {[::ogf::cat::frame_exists zz1 alldata] && ![::ogf::cat::frame_exists zz1 nope]}]
    ::ogf::cat::frame_set zz1 morph,map {1 2}; ::ogf::cat::frame_set zz2 alldata X
    CatalogPanelDeleteFrameState zz1
    R frame_delete [expr {![::ogf::cat::frame_exists zz1 alldata] && ![::ogf::cat::frame_exists zz1 morph,map] && [::ogf::cat::frame_exists zz2 alldata]}]
    ::ogf::cat::frame_delete zz2
    CatalogPanelClear
    puts $::fh "SUMMARY failures=$::nf"; close $::fh; exit
}
after 3000 {if {[catch run err]} {puts $::fh "ERROR $err $::errorInfo"; close $::fh; exit 1}}
