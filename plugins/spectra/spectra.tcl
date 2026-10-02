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
