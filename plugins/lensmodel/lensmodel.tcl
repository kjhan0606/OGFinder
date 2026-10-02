# Strong-lens model (plugins/lensmodel): results dialog (summary / images / counter-images / delays / curves), overlay on the current frame,
# map frames, selection helpers.  The numerical work is done by lensmodel.py (CLI steps of the manifest); everything here only reads its files.

proc OGFLensFile {name} {return [file join [OGFSessWorkDir] lensmodel $name]}

proc OGFLensRead {} {
    set f [OGFLensFile lens_model.json]
    if {![file exists $f]} {return {}}
    set fd [open $f r]; set t [read $fd]; close $fd
    if {[catch {::ogf::json::parse $t} d]} {::ogf::log ERROR "lensmodel: bad lens_model.json: $d"; return {}}
    return $d
}

# delete the regions of an earlier overlay (they form the group 'lensmodel'; other regions are untouched), then load the new file
proc OGFLensOverlay {} {
    global current
    set f [OGFLensFile lens_overlay.reg]
    if {![file exists $f]} {::ogf::status "Lens model: no overlay yet - run Fit Lens Model first"; return 0}
    if {$current(frame) eq {}} {return 0}
    catch {$current(frame) marker unselect all; $current(frame) marker {lensmodel} select; $current(frame) marker delete select}
    if {[catch {MarkerLoadFile $f $current(frame) ds9 image {}} err]} {::ogf::log ERROR "lensmodel overlay: $err"; return 0}
    ::ogf::status "Lens model: critical curves (red), caustics (cyan, source plane), predicted images drawn"
    return 1
}

proc OGFLensAfter {} {
    if {[::ogf::params::get lensmodel show-overlay]} {catch {OGFLensOverlay}}
    if {[winfo exists .ogflens]} {catch {OGFLensFill}}
}

proc OGFLensOpenFits {name label} {
    global current scale
    set f [OGFLensFile $name]
    if {![file exists $f]} {::ogf::status "Lens model: $name does not exist yet"; return 0}
    set orig $current(frame)
    CreateFrame
    if {[catch {LoadFitsFile $f {} {}} err]} {::ogf::log ERROR "lensmodel: cannot load $f: $err"; catch {GotoFrame $orig}; return 0}
    set scale(mode) zscale
    catch {ChangeScaleMode}
    ::ogf::status "Lens model: $label opened in a new frame"
    return 1
}
proc OGFLensMagFrame {} {OGFLensOpenFits lens_magnification.fits "magnification map"}
proc OGFLensSourceFrame {} {OGFLensOpenFits lens_source.fits "source-plane reconstruction"}

proc OGFLensPickImages {} {
    set sel [::ogf::cat::selection]
    if {[llength $sel] < 2} {::ogf::status "Lens model: select at least 2 table rows (the multiple images)"; return 0}
    ::ogf::params::put lensmodel image-numbers [join $sel ,]
    catch {::ogf::params::save lensmodel}
    ::ogf::status "Lens model: images = [join $sel ,]"
    return 1
}
proc OGFLensPickLens {} {
    set sel [::ogf::cat::selection]
    if {[llength $sel] != 1} {::ogf::status "Lens model: select exactly one table row (the lens galaxy)"; return 0}
    ::ogf::params::put lensmodel lens-number [lindex $sel 0]
    catch {::ogf::params::save lensmodel}
    ::ogf::status "Lens model: lens galaxy = NUMBER [lindex $sel 0]"
    return 1
}

# ------------------------------------------------------------------ results dialog
proc OGFLensDialog {} {
    set w .ogflens
    if {[winfo exists $w]} {raise $w; OGFLensFill; return $w}
    toplevel $w
    wm title $w "Lens model"
    ttk::notebook $w.nb
    foreach {tab label} {sum Summary img Images dly Curves} {
	ttk::frame $w.nb.$tab
	$w.nb add $w.nb.$tab -text $label
    }
    text $w.nb.sum.t -width 100 -height 24 -font TkFixedFont -wrap none -state disabled
    pack $w.nb.sum.t -fill both -expand 1
    ttk::treeview $w.nb.img.t -columns {kind x y mu parity delay res} -show headings -height 10
    foreach {c t wd} {kind Image 120 x X 80 y Y 80 mu mu 80 parity Parity 70 delay {Delay (d)} 90 res {Residual (")} 90} {
	$w.nb.img.t heading $c -text $t
	$w.nb.img.t column $c -width $wd -anchor center
    }
    pack $w.nb.img.t -fill both -expand 1
    ttk::treeview $w.nb.dly.t -columns {kind theta area} -show headings -height 6
    foreach {c t wd} {kind Curve 140 theta {theta_eff (")} 120 area {Area (arcsec^2)} 140} {
	$w.nb.dly.t heading $c -text $t
	$w.nb.dly.t column $c -width $wd -anchor center
    }
    pack $w.nb.dly.t -fill both -expand 1
    pack $w.nb -fill both -expand 1
    ttk::frame $w.b
    ttk::button $w.b.fit -text "Fit" -command {::ogf::step::run lensmodel fit}
    ttk::button $w.b.ov -text "Overlay" -command OGFLensOverlay
    ttk::button $w.b.mag -text "Magnification map" -command {::ogf::step::run lensmodel magmap}
    ttk::button $w.b.src -text "Source plane" -command {::ogf::step::run lensmodel source}
    ttk::button $w.b.pr -text "Predict" -command {::ogf::step::run lensmodel predict}
    ttk::button $w.b.cl -text Close -command [list destroy $w]
    pack $w.b.fit $w.b.ov $w.b.mag $w.b.src $w.b.pr $w.b.cl -side left -padx 3 -pady 4
    pack $w.b -fill x
    OGFLensFill
    return $w
}

proc OGFLensFill {} {
    set w .ogflens
    if {![winfo exists $w]} return
    set d [OGFLensRead]
    set t $w.nb.sum.t
    $t configure -state normal
    $t delete 1.0 end
    if {$d eq {}} {
	$t insert end "No lens model yet. Select the multiple images in the catalog table, press 'Use Selected Rows as Images' (and 'Use Selected Row as Lens'), then Fit Lens Model."
    } else {
	set p [dict get $d params]
	set e [expr {[dict exists $d errors] ? [dict get $d errors] : {}}]
	$t insert end [format "chi2 = %.3f for %d dof,  image-plane rms = %.4f arcsec (sigma_pos %.4f)\n" [dict get $d chi2] [dict get $d dof] [dict get $d rms_arcsec] [dict get $d sigma_pos]]
	foreach k {theta_E q phi gamma phi_g} {
	    set err [expr {[dict exists $e $k] ? [format " +- %.4f" [dict get $e $k]] : " (fixed)"}]
	    $t insert end [format "  %-8s = %9.4f%s\n" $k [dict get $p $k] $err]
	}
	lassign [dict get $d source_arcsec] sx sy
	lassign [dict get $d source_pixel] spx spy
	$t insert end [format "  source   = (%.4f, %.4f) arcsec  pixel (%.2f, %.2f)\n" $sx $sy $spx $spy]
	$t insert end [format "  theta_E (tangential critical curve) = %.4f arcsec\n" [dict get $d theta_E_critical_curve]]
	if {[dict exists $d pa_sky_deg] && [dict get $d pa_sky_deg] ne {}} {$t insert end [format "  sky PA of the mass major axis (E of N) = %.1f deg\n" [dict get $d pa_sky_deg]]}
	dict for {k v} [dict get $d physical] {$t insert end [format "  %s = %.4g\n" $k $v]}
	$t insert end [format "  predicted images: %d, counter-images not given: %d\n" [dict get $d n_images_predicted] [dict get $d n_counter_images]]
	foreach wmsg [dict get $d warnings] {$t insert end "  WARNING: $wmsg\n"}
    }
    $t configure -state disabled
    $w.nb.img.t delete [$w.nb.img.t children {}]
    $w.nb.dly.t delete [$w.nb.dly.t children {}]
    if {$d eq {}} return
    set nums [dict get $d image_numbers]
    set k 0
    foreach im [dict get $d images] {
	set m [dict get $im matched]
	set kind [expr {$m >= 0 ? "image [lindex $nums $m]" : "COUNTER-IMAGE"}]
	set res [expr {$m >= 0 ? [format %.4f [lindex [dict get $d residuals_arcsec] $m]] : ""}]
	set dl [expr {[dict exists $im delay_days] ? [format %.2f [dict get $im delay_days]] : ""}]
	$w.nb.img.t insert {} end -values [list $kind [format %.2f [dict get $im px]] [format %.2f [dict get $im py]] [format %+.2f [dict get $im mu]] [dict get $im parity] $dl $res]
    }
    set f [OGFLensFile lens_curves.json]
    if {[file exists $f]} {
	set fd [open $f r]; set tx [read $fd]; close $fd
	if {![catch {::ogf::json::parse $tx} cd]} {
	    foreach c [dict get $cd curves] {
		$w.nb.dly.t insert {} end -values [list [dict get $c kind] [format %.4f [dict get $c theta_eff]] [format %.3f [dict get $c area]]]
	    }
	}
    }
}
