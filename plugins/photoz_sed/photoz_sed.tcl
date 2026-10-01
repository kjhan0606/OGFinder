# photoz_sed/photoz_sed.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/photoz_sed/plugin.json.

proc CatalogPanelPhotoZParamLoad {} {
    global catpanel
    set prf [file join [file normalize ~] .ds9 photo_z.prf]
    if {![file exists $prf]} return
    if {[catch {
	set fd [open $prf r]
	set data [read $fd]
	close $fd
    }]} return
    foreach line [split $data "\n"] {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set eq [string first "=" $line]
	if {$eq < 0} continue
	set key [string trim [string range $line 0 [expr {$eq-1}]]]
	set val [string trim [string range $line [expr {$eq+1}] end]]
	set catpanel(photoz,param,$key) $val
    }
}

proc CatalogPanelPhotoZParamSave {} {
    global catpanel
    set prf [file join [file normalize ~] .ds9 photo_z.prf]
    catch {file mkdir [file dirname $prf]}
    if {[catch {set fd [open $prf w]}]} return
    foreach key {bands mag-columns checkpoint} {
	if {[info exists catpanel(photoz,param,$key)]} {
	    puts $fd "$key=$catpanel(photoz,param,$key)"
	}
    }
    close $fd
}

proc CatalogPanelPhotoZ {} {
    global catpanel
    if {[info commands OGFAIBackendHook] ne {} && [OGFAIBackendHook photoz]} return  ;# ogf_ai.tcl: backend local|external

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first"
	return
    }

    # Get column names for band mapping
    set lines [split $catpanel(alldata) "\n"]
    set header [split [lindex $lines 0] "\t"]

    set w .catphotoz
    if {[winfo exists $w]} { raise $w; return }
    toplevel $w
    wm title $w "Photo-z (AI)"
    wm geometry $w 400x280

    ttk::label $w.lbands -text "Bands (comma-sep):"
    ttk::entry $w.bands -width 30
    $w.bands insert 0 $catpanel(photoz,param,bands)

    ttk::label $w.lmags -text "Mag columns (comma-sep):"
    ttk::entry $w.mags -width 30
    # Auto-detect magnitude columns
    set magcols {}
    foreach c $header {
	if {[string match "MAG_*" $c] || [string match "mag_*" $c]} {
	    lappend magcols $c
	}
    }
    set default_mags [join $magcols ","]
    if {$catpanel(photoz,param,mag-columns) ne {}} {
	set default_mags $catpanel(photoz,param,mag-columns)
    }
    $w.mags insert 0 $default_mags

    ttk::label $w.lckpt -text "Checkpoint:"
    ttk::entry $w.ckpt -width 30
    $w.ckpt insert 0 $catpanel(photoz,param,checkpoint)

    ttk::frame $w.btns
    ttk::button $w.btns.run -text "Run Photo-z" \
	-command [list CatalogPanelPhotoZRun $w]
    ttk::button $w.btns.close -text "Close" -command [list destroy $w]
    pack $w.btns.run $w.btns.close -side left -padx 5

    grid $w.lbands $w.bands -padx 5 -pady 4 -sticky w
    grid $w.lmags  $w.mags  -padx 5 -pady 4 -sticky w
    grid $w.lckpt  $w.ckpt  -padx 5 -pady 4 -sticky w
    grid $w.btns   -        -padx 5 -pady 10
}

proc CatalogPanelPhotoZRun {dlg} {
    global catpanel

    set bands [$dlg.bands get]
    set magcols [$dlg.mags get]
    set ckpt [$dlg.ckpt get]

    # Save params
    set catpanel(photoz,param,bands) $bands
    set catpanel(photoz,param,mag-columns) $magcols
    set catpanel(photoz,param,checkpoint) $ckpt
    CatalogPanelPhotoZParamSave

    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_photo_z.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_photo_z.py"
	return
    }

    set tmpcat [CatalogPanelSaveTempCatalog "photoz"]
    if {$tmpcat eq {}} {
	set catpanel(status) "Failed to save temp catalog"
	return
    }

    set catpanel(status) "Running Photo-z estimation..."
    update idletasks

    set args [list [OGFPython] $script $fn --catalog $tmpcat]
    if {$bands ne {}} { lappend args --bands $bands }
    if {$magcols ne {}} { lappend args --mag-columns $magcols }
    if {$ckpt ne {} && [file exists $ckpt]} {
	lappend args --checkpoint $ckpt
    }
    lappend args --n-workers $catpanel(param,n-workers)

    OGFSessLog analysis.photo_z auto $args -title {Photo-z} -post [dict create kind add cols_list {PHOTO_Z PHOTO_Z_ERR PHOTO_Z_Q68 PHOTO_Z_OUTLIER}]
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "Photo-z error: $err"
	return
    }

    if {[string trim $data] eq {}} {
	set catpanel(status) "Photo-z: no output"
	return
    }

    CatalogPanelAddColumnsFromTSV $data \
	{PHOTO_Z PHOTO_Z_ERR PHOTO_Z_Q68 PHOTO_Z_OUTLIER}
    set catpanel(status) "Photo-z estimation complete"
    catch {destroy $dlg}
}

proc CatalogPanelSEDParamLoad {} {
    global catpanel
    set prf [file join [file normalize ~] .ds9 sed_fit.prf]
    if {![file exists $prf]} return
    if {[catch {
	set fd [open $prf r]
	set data [read $fd]
	close $fd
    }]} return
    foreach line [split $data "\n"] {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set eq [string first "=" $line]
	if {$eq < 0} continue
	set key [string trim [string range $line 0 [expr {$eq-1}]]]
	set val [string trim [string range $line [expr {$eq+1}] end]]
	set catpanel(sed,param,$key) $val
    }
}

proc CatalogPanelSEDParamSave {} {
    global catpanel
    set prf [file join [file normalize ~] .ds9 sed_fit.prf]
    catch {file mkdir [file dirname $prf]}
    if {[catch {set fd [open $prf w]}]} return
    foreach key {bands mag-columns photoz-column checkpoint-emulator checkpoint-inverse backend} {
	if {[info exists catpanel(sed,param,$key)]} {
	    puts $fd "$key=$catpanel(sed,param,$key)"
	}
    }
    close $fd
}

proc CatalogPanelSEDFit {} {
    global catpanel
    if {[info commands OGFAIBackendHook] ne {} && [OGFAIBackendHook sed_fit]} return  ;# ogf_ai.tcl: backend local|external

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first"
	return
    }

    # Check for PHOTO_Z column
    set lines [split $catpanel(alldata) "\n"]
    set header [split [lindex $lines 0] "\t"]
    set has_pz 0
    foreach c $header {
	if {$c eq "PHOTO_Z"} { set has_pz 1; break }
    }

    set w .catsedfit
    if {[winfo exists $w]} { raise $w; return }
    toplevel $w
    wm title $w "SED Fitting (AI)"
    wm geometry $w 420x380

    if {!$has_pz} {
	ttk::label $w.warn -text "Warning: PHOTO_Z column not found.\nRun Photo-z first for best results." \
	    -foreground red
	grid $w.warn - -padx 5 -pady 5
    }

    ttk::label $w.lbackend -text "SPS Backend:"
    ttk::combobox $w.backend -width 15 -state readonly \
	-values {auto analytic fsps bagpipes prospector cigale dense_basis}
    $w.backend set $catpanel(sed,param,backend)

    ttk::label $w.lbands -text "Bands (comma-sep):"
    ttk::entry $w.bands -width 30
    $w.bands insert 0 $catpanel(sed,param,bands)

    ttk::label $w.lmags -text "Mag columns (comma-sep):"
    ttk::entry $w.mags -width 30
    set magcols {}
    foreach c $header {
	if {[string match "MAG_*" $c] || [string match "mag_*" $c]} {
	    lappend magcols $c
	}
    }
    set default_mags [join $magcols ","]
    if {$catpanel(sed,param,mag-columns) ne {}} {
	set default_mags $catpanel(sed,param,mag-columns)
    }
    $w.mags insert 0 $default_mags

    ttk::label $w.lpzcol -text "Photo-z column:"
    ttk::entry $w.pzcol -width 20
    $w.pzcol insert 0 $catpanel(sed,param,photoz-column)

    ttk::label $w.lcke -text "Emulator checkpoint:"
    ttk::entry $w.cke -width 30
    $w.cke insert 0 $catpanel(sed,param,checkpoint-emulator)

    ttk::label $w.lcki -text "Inverse checkpoint:"
    ttk::entry $w.cki -width 30
    $w.cki insert 0 $catpanel(sed,param,checkpoint-inverse)

    ttk::frame $w.btns
    ttk::button $w.btns.run -text "Run SED Fit" \
	-command [list CatalogPanelSEDFitRun $w]
    ttk::button $w.btns.close -text "Close" -command [list destroy $w]
    pack $w.btns.run $w.btns.close -side left -padx 5

    grid $w.lbackend $w.backend -padx 5 -pady 4 -sticky w
    grid $w.lbands $w.bands -padx 5 -pady 4 -sticky w
    grid $w.lmags  $w.mags  -padx 5 -pady 4 -sticky w
    grid $w.lpzcol $w.pzcol -padx 5 -pady 4 -sticky w
    grid $w.lcke   $w.cke   -padx 5 -pady 4 -sticky w
    grid $w.lcki   $w.cki   -padx 5 -pady 4 -sticky w
    grid $w.btns   -        -padx 5 -pady 10
}

proc CatalogPanelSEDFitRun {dlg} {
    global catpanel

    set backend [$dlg.backend get]
    set bands [$dlg.bands get]
    set magcols [$dlg.mags get]
    set pzcol [$dlg.pzcol get]
    set cke [$dlg.cke get]
    set cki [$dlg.cki get]

    set catpanel(sed,param,backend) $backend
    set catpanel(sed,param,bands) $bands
    set catpanel(sed,param,mag-columns) $magcols
    set catpanel(sed,param,photoz-column) $pzcol
    set catpanel(sed,param,checkpoint-emulator) $cke
    set catpanel(sed,param,checkpoint-inverse) $cki
    CatalogPanelSEDParamSave

    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_sed_fit.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_sed_fit.py"
	return
    }

    set tmpcat [CatalogPanelSaveTempCatalog "sedfit"]
    if {$tmpcat eq {}} {
	set catpanel(status) "Failed to save temp catalog"
	return
    }

    set catpanel(status) "Running SED fitting ($backend)..."
    update idletasks

    set args [list [OGFPython] $script $fn --catalog $tmpcat]
    if {$backend ne {}} { lappend args --backend $backend }
    if {$bands ne {}} { lappend args --bands $bands }
    if {$magcols ne {}} { lappend args --mag-columns $magcols }
    if {$pzcol ne {}} { lappend args --photo-z-column $pzcol }
    if {$cke ne {} && [file exists $cke]} {
	lappend args --checkpoint-emulator $cke
    }
    if {$cki ne {} && [file exists $cki]} {
	lappend args --checkpoint-inverse $cki
    }
    lappend args --n-workers $catpanel(param,n-workers)

    OGFSessLog analysis.sed_fit auto $args -title {SED fitting} -post [dict create kind add cols_list {LOG_MASS LOG_MASS_ERR LOG_AGE LOG_AGE_ERR LOG_Z AV SFR}]
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "SED fit error: $err"
	return
    }

    if {[string trim $data] eq {}} {
	set catpanel(status) "SED fit: no output"
	return
    }

    CatalogPanelAddColumnsFromTSV $data \
	{LOG_MASS LOG_MASS_ERR LOG_AGE LOG_AGE_ERR LOG_Z AV SFR}
    set catpanel(status) "SED fitting complete ($backend)"
    catch {destroy $dlg}
}

