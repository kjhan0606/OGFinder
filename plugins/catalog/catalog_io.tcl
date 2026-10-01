# catalog/catalog_io.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/catalog/plugin.json.

# Load tab-separated catalog data into the panel
proc CatalogPanelLoadTSV {data source_name} {
    global catpanel

    global $catpanel(tbldb)

    # Unbind table from variable while modifying
    $catpanel(tbl) configure -variable {}

    unset -nocomplain $catpanel(tbldb)

    set lines [split $data \n]
    set nlines [llength $lines]

    if {$nlines < 2} {
	set catpanel(status) "No sources detected"
	$catpanel(tbl) configure -variable $catpanel(tbldb)
	return
    }

    # Store for filtering
    set catpanel(alldata) $data
    set catpanel(delim) "\t"
    catch {OGFSessOnCatalogLoad $source_name}

    # Parse header
    set headers [split [lindex $lines 0] "\t"]
    set ncols [llength $headers]

    # Fill header row
    for {set c 0} {$c < $ncols} {incr c} {
	set ${catpanel(tbldb)}(0,[expr {$c+1}]) \
	    [string trim [lindex $headers $c]]
    }

    # Fill data rows
    set row 1
    for {set i 1} {$i < $nlines} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	set fields [split $line "\t"]
	for {set c 0} {$c < $ncols} {incr c} {
	    set ${catpanel(tbldb)}($row,[expr {$c+1}]) \
		[string trim [lindex $fields $c]]
	}
	incr row
    }

    # Rebind table and configure dimensions to trigger full refresh
    $catpanel(tbl) configure -variable $catpanel(tbldb) \
	-cols $ncols -rows $row -state disabled

    set nobj [expr {$row - 1}]
    set catpanel(status) "$source_name: $nobj sources extracted"
    catch {OGFTDGalaxyLoaded}
    catch {OGFTDAppendKindColumn $ncols $row}
}

proc CatalogPanelClear {} {
    global catpanel
    global current

    # Delete all sextract markers
    if {$current(frame) != {}} {
	catch {$current(frame) marker catalog sextract_sel delete}
	catch {$current(frame) marker catalog sextract_all delete}
	catch {$current(frame) marker catalog sextract_merge delete}
    }

    catch {OGFSessLog catalog.clear manual {} -tool native -title "Clear catalog"}
    global $catpanel(tbldb)
    $catpanel(tbl) configure -variable {}
    unset -nocomplain $catpanel(tbldb)
    $catpanel(tbl) configure -variable $catpanel(tbldb) \
	-cols 19 -rows 20

    set catpanel(status) {Ready}
    catch {CatalogPanelClearSelection}
    set catpanel(sel,text) {No source selected}
    set catpanel(filename) {}
    set catpanel(alldata) {}
    catch {OGFTDGalaxyLoaded}

    # Reset merge state
    set catpanel(merge,list) {}
    set catpanel(merge,active) 0

    # Reset mark all state
    set catpanel(markall,on) 0

    # Reset visible mode
    set catpanel(visible_mode) 0

    # Reset add objects mode
    set catpanel(add_objects_mode) 0

    # Reset trim state
    set catpanel(trim,active) 0

    # Reset AI merge state
    if {$current(frame) != {}} {
	catch {$current(frame) marker catalog ai_merge delete}
    }
    set catpanel(ai,groups) {}
    set catpanel(ai,active) 0
    set catpanel(ai,total) 0
    set catpanel(ai,current) 0
    CatalogPanelAIUnbindKeys
}

proc CatalogPanelSaveCatalog {} {
    global catpanel

    if {[catch {::ogf::td::active} _tda]} {set _tda 0}
    if {!$_tda && (![info exists catpanel(alldata)] || $catpanel(alldata) eq {})} {
	set catpanel(status) "No catalog to save"
	return
    }

    set fn [tk_getSaveFile \
		-title "Save Catalog" \
		-defaultextension ".tsv" \
		-filetypes {
		    {{Tab-Separated Values} {.tsv}}
		    {{CSV Files} {.csv}}
		    {{All Files} {*}}
		}]
    if {$fn eq {}} return
    CatalogPanelSaveCatalogTo $fn
}

# Write the catalog to FN (.tsv or .csv).  Also what the session recorder replays.
proc CatalogPanelSaveCatalogTo {fn} {
    global catpanel
    if {[catch {::ogf::td::active} _tda]} {set _tda 0}
    if {$_tda} {::ogf::td::save $fn; return}
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "No catalog to save"
	return
    }
    OGFSessLog catalog.save auto {} -tool native -title "Save catalog as [file tail $fn]" \
	-payload [dict create name [file tail $fn]] -requires catalog
    set ext [string tolower [file extension $fn]]

    if {$ext eq ".csv"} {
	# Convert TSV to CSV
	set lines [split $catpanel(alldata) \n]
	set csvdata {}
	foreach line $lines {
	    if {[string trim $line] eq {}} continue
	    set fields [split $line \t]
	    set csvfields {}
	    foreach fld $fields {
		set fld [string trim $fld]
		if {[string match *,* $fld] || [string match *\"* $fld]} {
		    regsub -all {"} $fld {""} fld
		    set fld "\"$fld\""
		}
		lappend csvfields $fld
	    }
	    lappend csvdata [join $csvfields ,]
	}
	set outdata [join $csvdata \n]
    } else {
	set outdata $catpanel(alldata)
    }

    if {[catch {
	set fd [open $fn w]
	puts -nonewline $fd $outdata
	close $fd
    } err]} {
	set catpanel(status) "Save error: $err"
	return
    }

    set nlines [llength [split $catpanel(alldata) \n]]
    set nobj [expr {$nlines - 1}]
    set catpanel(status) "Saved $nobj sources to [file tail $fn]"
}

proc CatalogPanelLoadCatalog {} {
    global catpanel

    set fn [tk_getOpenFile \
		-title "Load Catalog" \
		-filetypes {
		    {{Tab-Separated Values} {.tsv}}
		    {{CSV Files} {.csv}}
		    {{All Files} {*}}
		}]
    if {$fn eq {}} return

    if {[catch {
	set fd [open $fn r]
	set rawdata [read $fd]
	close $fd
    } err]} {
	set catpanel(status) "Load error: $err"
	return
    }

    set rawdata [string trimright $rawdata \n]
    if {$rawdata eq {}} {
	set catpanel(status) "Empty file: [file tail $fn]"
	return
    }

    set ext [string tolower [file extension $fn]]

    if {$ext eq ".csv"} {
	# Convert CSV to TSV
	set lines [split $rawdata \n]
	set tsvlines {}
	foreach line $lines {
	    set line [string trimright $line \r]
	    if {$line eq {}} continue
	    # Simple CSV parse: split on comma, handle quoted fields
	    set fields {}
	    set cur {}
	    set inquote 0
	    for {set i 0} {$i < [string length $line]} {incr i} {
		set ch [string index $line $i]
		if {$inquote} {
		    if {$ch eq "\""} {
			if {$i+1 < [string length $line] && [string index $line [expr {$i+1}]] eq "\""} {
			    append cur "\""
			    incr i
			} else {
			    set inquote 0
			}
		    } else {
			append cur $ch
		    }
		} else {
		    if {$ch eq "\""} {
			set inquote 1
		    } elseif {$ch eq ","} {
			lappend fields $cur
			set cur {}
		    } else {
			append cur $ch
		    }
		}
	    }
	    lappend fields $cur
	    lappend tsvlines [join $fields \t]
	}
	set data [join $tsvlines \n]
    } else {
	set data $rawdata
    }

    OGFSessLog catalog.load manual {} -tool native -title "Load catalog [file tail $fn]" \
	-payload [dict create file [file normalize $fn]]
    CatalogPanelLoadTSV $data [file tail $fn]
}

proc CatalogPanelGetScript {scriptname} {
    set bindir [file dirname [info nameofexecutable]]
    set script [file join $bindir $scriptname]
    if {![file exists $script]} {
	set libdir [file join [file dirname $bindir] ds9 library]
	set script [file join $libdir $scriptname]
    }
    return $script
}

proc CatalogPanelGetFITS {} {
    global current
    set fn {}
    if {$current(frame) != {}} {
	catch {set fn [$current(frame) get fits file name full]}
    }
    set fn [string trim $fn "{}"]
    regsub {\[.*\]$} $fn {} fn
    return $fn
}

proc CatalogPanelFitsBaseName {fn} {
    set base [file tail $fn]
    # Strip .gz first if present
    if {[string match "*.gz" $base]} {
	set base [file rootname $base]
    }
    # Strip .fits/.fit/.fts
    set base [file rootname $base]
    return $base
}

proc CatalogPanelSaveTempCatalog {suffix} {
    global catpanel
    set tmpdir [file join [file normalize ~] .ds9]
    if {![file isdirectory $tmpdir]} { file mkdir $tmpdir }
    set tmpfile [file join $tmpdir "${suffix}_catalog.tsv"]
    if {[catch {
	set fd [open $tmpfile w]
	puts -nonewline $fd $catpanel(alldata)
	close $fd
    } err]} {
	return {}
    }
    return $tmpfile
}

proc CatalogPanelAddColumnsFromTSV {result_data col_names} {
    global catpanel
    catch {OGFSessOnAddColumns $result_data $col_names}

    # Parse result
    set rlines [split $result_data \n]
    set rheaders [split [lindex $rlines 0] "\t"]

    # Find NUMBER column in results
    set r_num_col -1
    for {set c 0} {$c < [llength $rheaders]} {incr c} {
	if {[string trim [lindex $rheaders $c]] eq "NUMBER"} {
	    set r_num_col $c
	    break
	}
    }
    if {$r_num_col < 0} return

    # Build lookup: number -> values
    array set rdata {}
    for {set i 1} {$i < [llength $rlines]} {incr i} {
	set line [lindex $rlines $i]
	if {[string trim $line] eq {}} continue
	set fields [split $line "\t"]
	set num [string trim [lindex $fields $r_num_col]]
	set vals {}
	foreach cn $col_names {
	    set cidx -1
	    for {set c 0} {$c < [llength $rheaders]} {incr c} {
		if {[string trim [lindex $rheaders $c]] eq $cn} {
		    set cidx $c
		    break
		}
	    }
	    if {$cidx >= 0 && $cidx < [llength $fields]} {
		lappend vals [string trim [lindex $fields $cidx]]
	    } else {
		lappend vals {}
	    }
	}
	set rdata($num) $vals
    }

    # Parse alldata
    set lines [split $catpanel(alldata) \n]
    set headers [split [lindex $lines 0] "\t"]
    set ncols [llength $headers]

    # Find NUMBER col in alldata
    set num_col -1
    for {set c 0} {$c < $ncols} {incr c} {
	if {[string trim [lindex $headers $c]] eq "NUMBER"} {
	    set num_col $c
	    break
	}
    }
    if {$num_col < 0} return

    # Check if columns already exist
    set existing 0
    foreach cn $col_names {
	for {set c 0} {$c < $ncols} {incr c} {
	    if {[string trim [lindex $headers $c]] eq $cn} {
		set existing 1
		break
	    }
	}
	if {$existing} break
    }

    # Build new data
    set newdata {}
    if {$existing} {
	# Update existing columns
	set col_indices {}
	foreach cn $col_names {
	    set cidx -1
	    for {set c 0} {$c < $ncols} {incr c} {
		if {[string trim [lindex $headers $c]] eq $cn} {
		    set cidx $c
		    break
		}
	    }
	    lappend col_indices $cidx
	}

	append newdata [lindex $lines 0]
	for {set i 1} {$i < [llength $lines]} {incr i} {
	    set line [lindex $lines $i]
	    if {[string trim $line] eq {}} continue
	    set fields [split $line "\t"]
	    set sn [string trim [lindex $fields $num_col]]
	    if {[info exists rdata($sn)]} {
		set vals $rdata($sn)
		for {set v 0} {$v < [llength $col_names]} {incr v} {
		    set ci [lindex $col_indices $v]
		    if {$ci >= 0} {
			lset fields $ci [lindex $vals $v]
		    }
		}
	    }
	    append newdata "\n" [join $fields "\t"]
	}
    } else {
	# Append new columns
	append newdata [lindex $lines 0]
	foreach cn $col_names {
	    append newdata "\t" $cn
	}
	for {set i 1} {$i < [llength $lines]} {incr i} {
	    set line [lindex $lines $i]
	    if {[string trim $line] eq {}} continue
	    set fields [split $line "\t"]
	    set sn [string trim [lindex $fields $num_col]]
	    append newdata "\n" $line
	    if {[info exists rdata($sn)]} {
		foreach v $rdata($sn) {
		    append newdata "\t" $v
		}
	    } else {
		foreach cn $col_names {
		    append newdata "\t"
		}
	    }
	}
    }

    set catpanel(alldata) $newdata
    CatalogPanelLoadTSV $catpanel(alldata) "analysis"
}

proc CatalogPanelExportRegions {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "No catalog to export"
	return
    }

    set fn [tk_getSaveFile \
	-title "Export DS9 Regions" \
	-defaultextension ".reg" \
	-filetypes {
	    {{DS9 Region Files} {.reg}}
	    {{All Files} {*}}
	}]
    if {$fn eq {}} return

    # Parse alldata
    set lines [split $catpanel(alldata) \n]
    set headers [split [lindex $lines 0] "\t"]
    set col_map {}
    for {set c 0} {$c < [llength $headers]} {incr c} {
	dict set col_map [string trim [lindex $headers $c]] $c
    }

    foreach needed {X_IMAGE Y_IMAGE A_IMAGE B_IMAGE THETA_IMAGE NUMBER} {
	if {![dict exists $col_map $needed]} {
	    set catpanel(status) "Missing column: $needed"
	    return
	}
    }

    set col_x [dict get $col_map X_IMAGE]
    set col_y [dict get $col_map Y_IMAGE]
    set col_a [dict get $col_map A_IMAGE]
    set col_b [dict get $col_map B_IMAGE]
    set col_t [dict get $col_map THETA_IMAGE]
    set col_n [dict get $col_map NUMBER]

    if {[catch {
	set fd [open $fn w]
	puts $fd "# Region file format: DS9 version 4.1"
	puts $fd {global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1}
	puts $fd "image"

	set nreg 0
	for {set i 1} {$i < [llength $lines]} {incr i} {
	    set line [lindex $lines $i]
	    if {[string trim $line] eq {}} continue
	    set fields [split $line "\t"]
	    set x [string trim [lindex $fields $col_x]]
	    set y [string trim [lindex $fields $col_y]]
	    set a [string trim [lindex $fields $col_a]]
	    set b [string trim [lindex $fields $col_b]]
	    set t [string trim [lindex $fields $col_t]]
	    set n [string trim [lindex $fields $col_n]]
	    puts $fd "ellipse($x,$y,$a,$b,$t) # text=\{$n\}"
	    incr nreg
	}
	close $fd
    } err]} {
	set catpanel(status) "Export error: $err"
	return
    }

    set catpanel(status) "Exported $nreg regions to [file tail $fn]"
}

proc CatalogPanelExportFITS {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "No catalog to export"
	return
    }

    set fn [tk_getSaveFile \
	-title "Export FITS Table" \
	-defaultextension ".fits" \
	-filetypes {
	    {{FITS Files} {.fits}}
	    {{All Files} {*}}
	}]
    if {$fn eq {}} return

    # Save temp TSV
    set tmpfile [CatalogPanelSaveTempCatalog "fits_export"]
    if {$tmpfile eq {}} {
	set catpanel(status) "Failed to save temp catalog"
	return
    }

    set script [CatalogPanelGetScript ds9_fits_export.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_fits_export.py"
	return
    }

    set catpanel(status) "Exporting FITS table..."
    update idletasks

    OGFSessLog catalog.export_fits auto [list [OGFPython] $script --input $tmpfile --output $fn] -title "Export FITS table [file tail $fn]" -useroutputs [list $fn]
    if {[catch {
	set data [exec [OGFPython] $script --input $tmpfile --output $fn 2>@stderr]
    } err]} {
	set catpanel(status) "FITS export error: $err"
	return
    }

    set catpanel(status) "Exported FITS table to [file tail $fn]"
}

proc CatalogPanelSaveFrameState {frame} {
    global catpanel catpanel_fdata

    if {$frame eq {}} return
    if {![info exists catpanel(alldata)]} return

    # Core catalog + display + merge + AI merge
    foreach key {
	alldata filename sort,col sort,dir
	visible_mode markall,on add_objects_mode trim,active
	merge,list merge,active
	ai,groups ai,current ai,total ai,active
	psf,stars psf,star_indices psf,file psf,has_psf
	icl,has_mask icl,has_bkg icl,has_profile icl,center_x icl,center_y
	icl,mask_file icl,masked_file icl,bkg_file icl,bgsub_file icl,profile_file
	lsbg,has_mask lsbg,has_clean lsbg,has_detect lsbg,has_catalog lsbg,detect_data
	lsbg,mask_file lsbg,masked_file lsbg,bkg_file lsbg,cleaned_file
	lsbg,segmap_file lsbg,catalog_file
	status search_var
    } {
	if {[info exists catpanel($key)]} {
	    set catpanel_fdata($frame,$key) $catpanel($key)
	}
    }

    # Morph data: save map + per-source entries
    if {[info exists catpanel(morph,map)]} {
	set catpanel_fdata($frame,morph,map) $catpanel(morph,map)
	foreach src_num $catpanel(morph,map) {
	    if {[info exists catpanel(morph,$src_num)]} {
		set catpanel_fdata($frame,morph,$src_num) $catpanel(morph,$src_num)
	    }
	}
    } else {
	set catpanel_fdata($frame,morph,map) {}
    }
}

proc CatalogPanelRestoreFrameState {frame} {
    global catpanel catpanel_fdata

    if {$frame eq {}} return

    # Unbind AI keys if active
    if {[info exists catpanel(ai,active)] && $catpanel(ai,active)} {
	CatalogPanelAIUnbindKeys
    }

    # Check if we have saved data for this frame
    if {![info exists catpanel_fdata($frame,alldata)]} {
	# No saved state — initialize to empty (without deleting markers)
	global $catpanel(tbldb)
	$catpanel(tbl) configure -variable {}
	unset -nocomplain $catpanel(tbldb)
	$catpanel(tbl) configure -variable $catpanel(tbldb) \
	    -cols 19 -rows 20

	set catpanel(alldata) {}
	set catpanel(filename) {}
	set catpanel(sort,col) {}
	set catpanel(sort,dir) {}
	set catpanel(visible_mode) 0
	set catpanel(markall,on) 0
	set catpanel(add_objects_mode) 0
	set catpanel(trim,active) 0
	set catpanel(merge,list) {}
	set catpanel(merge,active) 0
	set catpanel(ai,groups) {}
	set catpanel(ai,current) 0
	set catpanel(ai,total) 0
	set catpanel(ai,active) 0
	set catpanel(psf,stars) {}
	set catpanel(psf,star_indices) {}
	set catpanel(psf,file) [file join [file normalize ~] .ds9 psf_current.fits]
	set catpanel(psf,has_psf) 0
	set catpanel(icl,fits_base) {}
	set catpanel(icl,has_mask) 0
	set catpanel(icl,has_bkg) 0
	set catpanel(icl,has_profile) 0
	set catpanel(icl,center_x) {}
	set catpanel(icl,center_y) {}
	set catpanel(icl,click_mode) 0
	set catpanel(icl,mask_file) [file join [file normalize ~] .ds9 icl_mask.fits]
	set catpanel(icl,masked_file) [file join [file normalize ~] .ds9 icl_masked.fits]
	set catpanel(icl,bkg_file) [file join [file normalize ~] .ds9 icl_background.fits]
	set catpanel(icl,bgsub_file) [file join [file normalize ~] .ds9 icl_bgsub.fits]
	set catpanel(icl,profile_file) [file join [file normalize ~] .ds9 icl_profile.tsv]
	set catpanel(lsbg,fits_base) {}
	set catpanel(lsbg,has_mask) 0
	set catpanel(lsbg,has_clean) 0
	set catpanel(lsbg,has_detect) 0
	set catpanel(lsbg,has_catalog) 0
	set catpanel(lsbg,detect_data) {}
	set catpanel(lsbg,mask_file) [file join [file normalize ~] .ds9 lsbg_mask.fits]
	set catpanel(lsbg,masked_file) [file join [file normalize ~] .ds9 lsbg_masked.fits]
	set catpanel(lsbg,bkg_file) [file join [file normalize ~] .ds9 lsbg_background.fits]
	set catpanel(lsbg,cleaned_file) [file join [file normalize ~] .ds9 lsbg_cleaned.fits]
	set catpanel(lsbg,segmap_file) [file join [file normalize ~] .ds9 lsbg_segmap.fits]
	set catpanel(lsbg,catalog_file) [file join [file normalize ~] .ds9 lsbg_catalog.tsv]
	set catpanel(status) {Ready}
	set catpanel(search_var) {}

	# Clear morph state
	if {[info exists catpanel(morph,map)]} {
	    foreach src_num $catpanel(morph,map) {
		unset -nocomplain catpanel(morph,$src_num)
	    }
	}
	set catpanel(morph,map) {}

	return
    }

    # Restore saved state
    foreach key {
	alldata filename sort,col sort,dir
	visible_mode markall,on add_objects_mode trim,active
	merge,list merge,active
	ai,groups ai,current ai,total ai,active
	psf,stars psf,star_indices psf,file psf,has_psf
	icl,has_mask icl,has_bkg icl,has_profile icl,center_x icl,center_y
	icl,mask_file icl,masked_file icl,bkg_file icl,bgsub_file icl,profile_file
	lsbg,has_mask lsbg,has_clean lsbg,has_detect lsbg,has_catalog lsbg,detect_data
	lsbg,mask_file lsbg,masked_file lsbg,bkg_file lsbg,cleaned_file
	lsbg,segmap_file lsbg,catalog_file
	status search_var
    } {
	if {[info exists catpanel_fdata($frame,$key)]} {
	    set catpanel($key) $catpanel_fdata($frame,$key)
	}
    }

    # Restore morph data
    # First clear old morph entries
    if {[info exists catpanel(morph,map)]} {
	foreach src_num $catpanel(morph,map) {
	    unset -nocomplain catpanel(morph,$src_num)
	}
    }
    set catpanel(morph,map) {}
    if {[info exists catpanel_fdata($frame,morph,map)]} {
	set catpanel(morph,map) $catpanel_fdata($frame,morph,map)
	foreach src_num $catpanel(morph,map) {
	    if {[info exists catpanel_fdata($frame,morph,$src_num)]} {
		set catpanel(morph,$src_num) $catpanel_fdata($frame,morph,$src_num)
	    }
	}
    }

    # Reload the table from alldata
    if {$catpanel(alldata) ne {}} {
	CatalogPanelLoadTSV $catpanel(alldata) [file tail $catpanel(filename)]
    } else {
	global $catpanel(tbldb)
	$catpanel(tbl) configure -variable {}
	unset -nocomplain $catpanel(tbldb)
	$catpanel(tbl) configure -variable $catpanel(tbldb) \
	    -cols 19 -rows 20
    }
}

proc CatalogPanelFrameChanged {old_frame new_frame} {
    # registered band frames share one catalog: nothing to swap
    if {$old_frame ne {} && [OGFBandsShareCatalog $old_frame $new_frame]} {
	return
    }
    set _tdk [OGFTDKeep]
    if {$old_frame ne {} && $old_frame ne $new_frame} {
	CatalogPanelSaveFrameState $old_frame
    }
    CatalogPanelRestoreFrameState $new_frame
    catch {OGFTDRestoreKind $_tdk}
}

proc CatalogPanelDeleteFrameState {frame} {
    global catpanel_fdata

    # Remove all saved state for the deleted frame
    foreach key [array names catpanel_fdata "$frame,*"] {
	unset catpanel_fdata($key)
    }
}

proc CatalogPanelClearAll {} {
    global catpanel

    # Clear table
    if {[info exists catpanel(tbldb)] && [info exists catpanel(tbl)]} {
	global $catpanel(tbldb)
	$catpanel(tbl) configure -variable {}
	unset -nocomplain $catpanel(tbldb)
	$catpanel(tbl) configure -variable $catpanel(tbldb) \
	    -cols 19 -rows 20
    }

    set catpanel(alldata) {}
    set catpanel(filename) {}
    set catpanel(sort,col) {}
    set catpanel(sort,dir) {}
    set catpanel(visible_mode) 0
    set catpanel(markall,on) 0
    set catpanel(add_objects_mode) 0
    set catpanel(trim,active) 0
    set catpanel(merge,list) {}
    set catpanel(merge,active) 0
    set catpanel(ai,groups) {}
    set catpanel(ai,current) 0
    set catpanel(ai,total) 0
    set catpanel(ai,active) 0
    set catpanel(psf,stars) {}
    set catpanel(psf,star_indices) {}
    set catpanel(psf,has_psf) 0
    set catpanel(icl,has_mask) 0
    set catpanel(icl,has_bkg) 0
    set catpanel(icl,has_profile) 0
    set catpanel(icl,center_x) {}
    set catpanel(icl,center_y) {}
    set catpanel(icl,click_mode) 0
    set catpanel(lsbg,has_mask) 0
    set catpanel(lsbg,has_clean) 0
    set catpanel(lsbg,has_detect) 0
    set catpanel(lsbg,has_catalog) 0
    set catpanel(lsbg,detect_data) {}
    set catpanel(status) {Ready}
    set catpanel(search_var) {}

    # Clear morph state
    if {[info exists catpanel(morph,map)]} {
	foreach src_num $catpanel(morph,map) {
	    unset -nocomplain catpanel(morph,$src_num)
	}
    }
    set catpanel(morph,map) {}
}

