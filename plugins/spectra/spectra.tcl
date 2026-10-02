# Spectroscopy (plugins/spectra): link table / results windows and the spectrum viewer.

proc OGFSpectraFile {name} {return [file join [OGFSessWorkDir] spectra $name]}

proc OGFSpectraAfter {} {
    foreach {k f} {links_file spectra_links.tsv results_file spectra_results.json} {
	set p [OGFSpectraFile $f]
	if {[file exists $p]} {::ogf::cat::set spectra,$k $p}
    }
}

# command line of "spectra.py --task plot" from the current parameters
proc OGFSpectraPlotArgs {num out} {
    set P [list spec-dir --spec-dir spec-pattern --spec-pattern link-file --link-file kind --kind wave-unit --wave-unit row-column --row-column extract-mode --extract-mode \
	extract-halfwidth --extract-halfwidth sky-inner --sky-inner sky-outer --sky-outer aperture-radius --aperture-radius snr-min --snr-min z-column --z-column]
    set a [list --task plot --number $num --out $out --work [file join [OGFSessWorkDir] spectra]]
    foreach {p opt} $P {lappend a $opt [::ogf::params::get spectra $p]}
    if {[::ogf::params::get spectra cube-xy]} {lappend a --cube-xy 1}
    return $a
}

proc OGFSpectraRender {num} {
    set cat [::ogf::cat::temp_file .tsv]
    set out [OGFSpectraFile view_$num.png]
    file mkdir [file dirname $out]
    set script [file join [::ogf::step::plugin_dir spectra] spectra.py]
    if {[catch {exec [OGFPython] $script --catalog $cat {*}[OGFSpectraPlotArgs $num $out]} err]} {
	::ogf::status "Spectra: $num: [lindex [split $err \n] end]"
	return {}
    }
    return $out
}

proc OGFSpectraLinkedNumbers {} {
    set f [OGFSpectraFile spectra_links.tsv]
    set nums {}
    if {![file exists $f]} {return $nums}
    set fd [open $f r]; set lines [split [read $fd] \n]; close $fd
    foreach l [lrange $lines 1 end] {
	set t [split $l \t]
	if {[llength $t] >= 7 && [lindex $t 6] == 1} {lappend nums [lindex $t 0]}
    }
    return $nums
}

proc OGFSpectraViewer {{num {}}} {
    set nums [OGFSpectraLinkedNumbers]
    if {$nums eq {}} {::ogf::status "Spectra: no linked spectra - run Link Spectra to Catalog first"; return {}}
    if {$num eq {}} {
	set sel [::ogf::cat::selection]
	set num [expr {[llength $sel] ? [lindex $sel 0] : [lindex $nums 0]}]
    }
    set w .ogfspecview
    if {![winfo exists $w]} {
	toplevel $w
	wm title $w "Spectrum viewer"
	label $w.l
	pack $w.l -fill both -expand 1
	frame $w.b
	ttk::button $w.b.prev -text "< Prev" -command [list OGFSpectraStep -1]
	ttk::label $w.b.nl -text "NUMBER"
	ttk::combobox $w.b.num -width 8 -values $nums -textvariable ::ogf_spec_num
	ttk::button $w.b.go -text Show -command {OGFSpectraShow $::ogf_spec_num}
	ttk::button $w.b.next -text "Next >" -command [list OGFSpectraStep 1]
	ttk::button $w.b.close -text Close -command [list destroy $w]
	pack $w.b.prev $w.b.nl $w.b.num $w.b.go $w.b.next $w.b.close -side left -padx 3
	pack $w.b -pady 3
	bind $w.b.num <<ComboboxSelected>> {OGFSpectraShow $::ogf_spec_num}
	bind $w.b.num <Return> {OGFSpectraShow $::ogf_spec_num}
    } else {
	$w.b.num configure -values $nums
    }
    OGFSpectraShow $num
    return $w
}

proc OGFSpectraShow {num} {
    set ::ogf_spec_num $num
    set png [OGFSpectraRender $num]
    if {$png eq {} || ![file exists $png]} return
    catch {image delete ogfspecimg}
    image create photo ogfspecimg -file $png
    .ogfspecview.l configure -image ogfspecimg
    wm title .ogfspecview "Spectrum viewer - object $num"
}

proc OGFSpectraStep {d} {
    set nums [OGFSpectraLinkedNumbers]
    set i [lsearch -exact $nums $::ogf_spec_num]
    set i [expr {$i < 0 ? 0 : ($i + $d) % [llength $nums]}]
    OGFSpectraShow [lindex $nums $i]
}

proc OGFSpectraTable {} {
    set f [OGFSpectraFile spectra_links.tsv]
    if {![file exists $f]} {::ogf::status "Spectra: nothing yet - run Link Spectra to Catalog first"; return}
    set fd [open $f r]; set txt [read $fd]; close $fd
    OGFTextWindow "Spectra: link table" $txt
}

# ------------------------------------------------------------------ kinematics (2D slit rotation curves, IFU velocity / dispersion maps)
proc OGFSpectraKinAfter {} {
    set p [OGFSpectraFile spectra_kin.json]
    if {[file exists $p]} {::ogf::cat::set spectra,kin_file $p}
    if {[winfo exists .ogfkinview]} {catch {OGFSpectraKinShow $::ogf_kin_num}}
}

proc OGFSpectraKinNumbers {} {
    set f [OGFSpectraFile spectra_kin.json]
    if {![file exists $f]} {return {}}
    set fd [open $f r]; set t [read $fd]; close $fd
    if {[catch {::ogf::json::parse $t} d]} {return {}}
    set nums {}
    dict for {k v} $d {if {[file exists [OGFSpectraFile kin_$k.png]]} {lappend nums $k}}
    return [lsort -integer $nums]
}

proc OGFSpectraKinViewer {{num {}}} {
    set nums [OGFSpectraKinNumbers]
    if {$nums eq {}} {::ogf::status "Spectra: no kinematics yet - run Fit Kinematics on linked 2D spectra / cubes first"; return {}}
    if {$num eq {}} {
	set sel [::ogf::cat::selection]
	set num [expr {[llength $sel] && [lindex $sel 0] in $nums ? [lindex $sel 0] : [lindex $nums 0]}]
    }
    set w .ogfkinview
    if {![winfo exists $w]} {
	toplevel $w
	wm title $w "Kinematics viewer"
	label $w.l
	pack $w.l -fill both -expand 1
	text $w.t -height 8 -width 90 -font TkFixedFont -wrap none -state disabled
	pack $w.t -fill x
	frame $w.b
	ttk::label $w.b.nl -text "NUMBER"
	ttk::combobox $w.b.num -width 8 -values $nums -textvariable ::ogf_kin_num
	ttk::button $w.b.go -text Show -command {OGFSpectraKinShow $::ogf_kin_num}
	ttk::button $w.b.fr -text "Open maps in frames" -command {OGFSpectraKinFrames $::ogf_kin_num}
	ttk::button $w.b.close -text Close -command [list destroy $w]
	pack $w.b.nl $w.b.num $w.b.go $w.b.fr $w.b.close -side left -padx 3
	pack $w.b -pady 3
	bind $w.b.num <<ComboboxSelected>> {OGFSpectraKinShow $::ogf_kin_num}
    } else {
	$w.b.num configure -values $nums
    }
    OGFSpectraKinShow $num
    return $w
}

proc OGFSpectraKinShow {num} {
    set ::ogf_kin_num $num
    set png [OGFSpectraFile kin_$num.png]
    if {![file exists $png]} {::ogf::status "Spectra: no kinematics for object $num"; return}
    catch {image delete ogfkinimg}
    image create photo ogfkinimg -file $png
    .ogfkinview.l configure -image ogfkinimg
    wm title .ogfkinview "Kinematics viewer - object $num"
    set fd [open [OGFSpectraFile spectra_kin.json] r]; set t [read $fd]; close $fd
    set d [::ogf::json::parse $t]
    .ogfkinview.t configure -state normal
    .ogfkinview.t delete 1.0 end
    if {[dict exists $d $num]} {
	set o [dict get $d $num]
	.ogfkinview.t insert end "object $num: [dict get $o kind], line [dict get $o line], z = [format %.5f [dict get $o z]] ([dict get $o z_source])\n"
	if {[dict exists $o velocity_field]} {
	    set v [dict get $o velocity_field]
	    .ogfkinview.t insert end [format "  PA %.1f deg  inc %.1f deg  vc %.1f km/s  rt %.2f px  vsys %.1f km/s  chi2r %.2f  (%d pixels)\n" [dict get $v pa] [dict get $v inc] [dict get $v vc] [dict get $v rt] [dict get $v vsys] [dict get $v chi2r] [dict get $v n]]
	} elseif {[dict exists $o rotation_curve]} {
	    set v [dict get $o rotation_curve]
	    if {[dict get $v ok]} {.ogfkinview.t insert end [format "  v_obs %.1f km/s  rt %.2f px  vsys %.1f km/s  chi2r %.2f  (%d bins)\n" [dict get $v vflat_obs] [dict get $v rt] [dict get $v vsys] [dict get $v chi2r] [dict get $o n_bins]]}
	}
	if {[dict exists $o dispersion]} {.ogfkinview.t insert end "  intrinsic dispersion (flux-weighted mean) [format %.1f [dict get [dict get $o dispersion] mean]] km/s\n"}
    }
    .ogfkinview.t configure -state disabled
}

proc OGFSpectraKinFrames {num} {
    global current scale
    set n 0
    foreach {nm label} [list kin_${num}_vel.fits velocity kin_${num}_sigma.fits dispersion kin_${num}_flux.fits flux] {
	set f [OGFSpectraFile $nm]
	if {![file exists $f]} continue
	CreateFrame
	if {[catch {LoadFitsFile $f {} {}} err]} {::ogf::log ERROR "spectra kin: cannot load $f: $err"; continue}
	set scale(mode) minmax
	catch {ChangeScaleMode}
	incr n
    }
    ::ogf::status "Spectra: $n kinematic map(s) opened in new frames"
    return $n
}
