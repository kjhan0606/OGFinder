# catalog/cli_script.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/catalog/plugin.json.

proc CatalogPanelCmdLog {pipeline args_list} {
    global catpanel
    lappend catpanel($pipeline,cmdlog) $args_list
    # session recorder (python script export); the ICL/LSBG .sh export above is unchanged
    catch {OGFSessFromCmdLog $pipeline $args_list}
}

proc ShellQuote {s} {
    if {[string index $s 0] eq "\$"} {
	return $s
    }
    if {[regexp {[[:space:]'\"\\;&|<>()]} $s]} {
	return "'[string map {' '\\''} $s]'"
    }
    return $s
}

proc CatalogPanelExportCLIScript {pipeline} {
    global catpanel

    set types {{"Shell Script" {.sh}} {"All Files" *}}
    set outfile [tk_getSaveFile -title "Export $pipeline CLI Script" \
	-filetypes $types \
	-defaultextension .sh \
	-initialfile "${pipeline}_pipeline.sh"]
    if {$outfile eq {}} return

    set fits_file [CatalogPanelGetFITS]

    if {[llength $catpanel($pipeline,cmdlog)] > 0} {
	set commands $catpanel($pipeline,cmdlog)
    } else {
	set commands [CatalogPanel[string totitle $pipeline]GenerateScript]
    }

    if {[llength $commands] == 0} {
	set catpanel(status) "$pipeline: No commands to export"
	return
    }

    set script_dir [file dirname [CatalogPanelGetScript ds9_${pipeline}.py]]
    set ds9dir [file join [file normalize ~] .ds9]
    set PIPELINE [string toupper $pipeline]

    # --- Collect all --param value pairs across commands as shell variables ---
    set param_vars {}      ;# list of {VARNAME value flag_name}
    set param_var_map {}   ;# dict: value -> VARNAME (for dedup)
    set output_files {}    ;# list of {VARNAME basename flag_name}

    foreach cmd $commands {
	for {set i 0} {$i < [llength $cmd]} {incr i} {
	    set arg [lindex $cmd $i]
	    if {[string match "--*" $arg] && $i + 1 < [llength $cmd]} {
		set val [lindex $cmd [expr {$i + 1}]]
		# Skip positional args, python3, script path
		if {$val eq "python3" || [string match "--*" $val]} continue

		# Convert flag name to shell var: --bkg-method → BKG_METHOD
		set varname [string range $arg 2 end]
		set varname [string map {- _} $varname]
		set varname [string toupper $varname]

		if {[string match "--*-output" $arg]} {
		    # Output file → use $OUTDIR/basename
		    lappend output_files [list $varname [file tail $val] $arg]
		} elseif {[string match "${ds9dir}/*" $val] ||
			  [regexp {\.(fits|tsv)(\.gz)?$} $val]} {
		    # Intermediate file → use $OUTDIR/basename
		    lappend output_files [list $varname [file tail $val] $arg]
		} elseif {![dict exists $param_var_map $arg]} {
		    # Parameter value → shell variable
		    lappend param_vars [list $varname $val $arg]
		    dict set param_var_map $arg $varname
		}
	    }
	}
    }

    # --- Write script ---
    set fd [open $outfile w]
    puts $fd "#!/bin/bash"
    puts $fd "# [string toupper $pipeline] pipeline — exported from DS9 GUI"
    puts $fd "# Generated: [clock format [clock seconds] -format {%Y-%m-%d %H:%M:%S}]"
    puts $fd "#"
    puts $fd "# Usage: bash [file tail $outfile] <input.fits> \[output_dir\]"
    puts $fd ""
    puts $fd "set -euo pipefail"
    puts $fd ""

    # Paths
    puts $fd "SCRIPT_DIR=\"\$(cd \"\$(dirname \"\$0\")\" && pwd)\""
    puts $fd "${PIPELINE}_PY=\"$script_dir/ds9_${pipeline}.py\""
    puts $fd ""
    puts $fd "# --- Input / Output ---"
    puts $fd "FITS=\"\${1:?Usage: bash [file tail $outfile] <input.fits> \[output_dir\]}\""
    puts $fd "OUTDIR=\"\${2:-./output}\""
    puts $fd "mkdir -p \"\$OUTDIR\""
    puts $fd ""

    # Parameter variables (grouped)
    if {[llength $param_vars] > 0} {
	puts $fd "# --- Parameters ---"
	foreach pv $param_vars {
	    lassign $pv varname val flag
	    puts $fd "$varname=$val"
	}
	puts $fd ""
    }

    # Validate
    puts $fd "# --- Validate ---"
    puts $fd "if \[ ! -f \"\$FITS\" \]; then"
    puts $fd "    echo \"ERROR: Input FITS not found: \$FITS\" >&2"
    puts $fd "    exit 1"
    puts $fd "fi"
    puts $fd ""

    # Commands
    set step 0
    foreach cmd $commands {
	incr step

	# Extract --mode value for echo
	set mode "step$step"
	for {set i 0} {$i < [llength $cmd]} {incr i} {
	    if {[lindex $cmd $i] eq "--mode"} {
		set mode [lindex $cmd [expr {$i + 1}]]
		break
	    }
	}
	puts $fd "# --- Step $step: $mode ---"
	puts $fd "echo \">>> Step $step: $mode ...\" >&2"
	puts $fd ""

	# Build command with shell variables and line continuations
	set first_line {}
	set cont_lines {}
	set found_script 0
	set found_input 0

	for {set i 0} {$i < [llength $cmd]} {incr i} {
	    set arg [lindex $cmd $i]

	    # python3
	    if {$arg eq "python3"} {
		set first_line "python3"
		continue
	    }

	    # Script path
	    if {!$found_script && [string match "*ds9_${pipeline}.py" $arg]} {
		append first_line " \"\$${PIPELINE}_PY\""
		set found_script 1
		continue
	    }

	    # Input FITS (first positional arg)
	    if {$found_script && !$found_input && [string index $arg 0] ne "-"} {
		if {$fits_file ne {} && $arg eq $fits_file} {
		    append first_line " \"\$FITS\""
		} elseif {[string match "${ds9dir}/*" $arg]} {
		    append first_line " \"\$OUTDIR/[file tail $arg]\""
		} else {
		    append first_line " \"\$FITS\""
		}
		set found_input 1
		continue
	    }

	    # Flag + value pair
	    if {[string match "--*" $arg] && $i + 1 < [llength $cmd]} {
		set val [lindex $cmd [expr {$i + 1}]]

		if {[string match "--*" $val]} {
		    # Boolean flag (no value)
		    lappend cont_lines "    $arg"
		    continue
		}

		incr i ;# consume value

		# Output files → $OUTDIR/basename
		if {[string match "--*-output" $arg] ||
		    ([string match "${ds9dir}/*" $val] &&
		     [regexp {\.(fits|tsv)(\.gz)?$} $val])} {
		    lappend cont_lines "    $arg \"\$OUTDIR/[file tail $val]\""
		    continue
		}

		# Known intermediate file flags
		if {$arg in {--mask --cleaned --segmap --catalog
			     --profile-file --import-mask-file} &&
		    [regexp {\.(fits|tsv)(\.gz)?$} $val]} {
		    lappend cont_lines "    $arg \"\$OUTDIR/[file tail $val]\""
		    continue
		}

		# Parameter with shell variable
		set varname [string range $arg 2 end]
		set varname [string map {- _} $varname]
		set varname [string toupper $varname]

		if {[dict exists $param_var_map $arg]} {
		    lappend cont_lines "    $arg \$$varname"
		} else {
		    lappend cont_lines "    $arg [ShellQuote $val]"
		}
		continue
	    }

	    # Boolean flag (standalone)
	    if {[string match "--*" $arg]} {
		lappend cont_lines "    $arg"
		continue
	    }
	}

	# Write with line continuations
	if {[llength $cont_lines] > 0} {
	    puts $fd "$first_line \\"
	    for {set j 0} {$j < [llength $cont_lines]} {incr j} {
		set line [lindex $cont_lines $j]
		if {$j < [llength $cont_lines] - 1} {
		    puts $fd "$line \\"
		} else {
		    puts $fd $line
		}
	    }
	} else {
	    puts $fd $first_line
	}
	puts $fd ""
    }

    # Summary
    puts $fd "echo \"\" >&2"
    puts $fd "echo \"=== $PIPELINE Pipeline Complete ($step steps) ===\" >&2"
    puts $fd "echo \"Output: \$OUTDIR/\" >&2"
    puts $fd "ls -lh \"\$OUTDIR/\"*.fits \"\$OUTDIR/\"*.tsv 2>/dev/null | while read line; do"
    puts $fd "    echo \"  \$line\" >&2"
    puts $fd "done"

    close $fd
    file attributes $outfile -permissions 0755

    set catpanel(status) "$pipeline: Script exported → [file tail $outfile] ($step steps)"
}

proc CatalogPanelIclGenerateScript {} {
    global catpanel
    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} return
    CatalogPanelICLUpdateFiles $fn
    set script [CatalogPanelGetScript ds9_icl.py]
    set commands {}

    # Step 1: shared mask (ds9_mask.py) + interpolated image for the fit
    set mscript [CatalogPanelGetScript ds9_mask.py]
    set mp [OGFMaskPaths $fn]
    set args [list [OGFPython] $mscript $fn --mode auto \
	--mask [dict get $mp mask] \
	--expand-factor $catpanel(icl,param,expand-factor) \
	--max-dilate-radius $catpanel(icl,param,max-dilate-radius) \
	--bright-star-mag-limit $catpanel(icl,param,bright-star-mag-limit) \
	--bright-star-radius-scale $catpanel(icl,param,bright-star-radius-scale) \
	--detect-thresh $catpanel(icl,param,detect-thresh) \
	--mag-zeropoint $catpanel(param,mag-zeropoint)]
    lappend commands $args
    set args [list [OGFPython] $mscript $fn --mode masked \
	--mask [dict get $mp mask] \
	--interp-method $catpanel(icl,param,interp-method) \
	--masked-output $catpanel(icl,masked_file)]
    lappend commands $args

    # Step 2: Background
    set args [list [OGFPython] $script $catpanel(icl,masked_file) --mode background \
	--bkg-method $catpanel(icl,param,bkg-method) \
	--bkg-order $catpanel(icl,param,bkg-order) \
	--bkg-sigma-clip $catpanel(icl,param,bkg-sigma-clip) \
	--bkg-sep-mesh $catpanel(icl,param,bkg-sep-mesh) \
	--bkg-output $catpanel(icl,bkg_file) \
	--bgsub-output $catpanel(icl,bgsub_file) \
	--mask $catpanel(icl,mask_file)]
    lappend commands $args

    # Step 3: Profile (if center is set)
    if {$catpanel(icl,center_x) ne {} && $catpanel(icl,center_y) ne {}} {
	set args [list [OGFPython] $script $catpanel(icl,bgsub_file) --mode profile \
	    --center "$catpanel(icl,center_x),$catpanel(icl,center_y)" \
	    --rmin $catpanel(icl,param,rmin) \
	    --rmax $catpanel(icl,param,rmax) \
	    --nsteps $catpanel(icl,param,nsteps) \
	    --spacing $catpanel(icl,param,spacing) \
	    --ellipticity $catpanel(icl,param,ellipticity) \
	    --pa $catpanel(icl,param,pa) \
	    --mag-zeropoint $catpanel(icl,param,mag-zeropoint) \
	    --pixel-scale $catpanel(icl,param,pixel-scale) \
	    --mask $catpanel(icl,mask_file) \
	    --profile-output $catpanel(icl,profile_file)]
	lappend commands $args

	# Step 4: Measure
	set args [list [OGFPython] $script $fn --mode measure \
	    --profile-file $catpanel(icl,profile_file) \
	    --mu-threshold $catpanel(icl,param,mu-threshold) \
	    --mu-levels $catpanel(icl,param,mu-levels) \
	    --pixel-scale $catpanel(icl,param,pixel-scale)]
	lappend commands $args
    }

    return $commands
}

proc CatalogPanelLsbgGenerateScript {} {
    global catpanel
    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} return
    CatalogPanelLSBGUpdateFiles $fn
    set script [CatalogPanelGetScript ds9_lsbg.py]

    set args [list [OGFPython] $script $fn --mode run \
	--mask-detect-thresh $catpanel(lsbg,param,mask-detect-thresh) \
	--mask-detect-minarea $catpanel(lsbg,param,mask-detect-minarea) \
	--mask-expand-factor $catpanel(lsbg,param,mask-expand-factor) \
	--max-dilate-radius $catpanel(lsbg,param,max-dilate-radius) \
	--bright-star-mag-limit $catpanel(lsbg,param,bright-star-mag-limit) \
	--bright-star-radius-scale $catpanel(lsbg,param,bright-star-radius-scale) \
	--mask-mag-threshold $catpanel(lsbg,param,mask-mag-threshold) \
	--interp-method $catpanel(lsbg,param,interp-method) \
	--lsb-mu-threshold $catpanel(lsbg,param,lsb-mu-threshold) \
	--bkg-method $catpanel(lsbg,param,bkg-method) \
	--bkg-mesh-size $catpanel(lsbg,param,bkg-mesh-size) \
	--bkg-poly-order $catpanel(lsbg,param,bkg-poly-order) \
	--bkg-sigma-clip $catpanel(lsbg,param,bkg-sigma-clip) \
	--bkg-n-iterations $catpanel(lsbg,param,bkg-n-iterations) \
	--bkg-refine-thresh $catpanel(lsbg,param,bkg-refine-thresh) \
	--bkg-rms-quantile $catpanel(lsbg,param,bkg-rms-quantile) \
	--bkg-convergence-tol $catpanel(lsbg,param,bkg-convergence-tol) \
	--detect-thresh $catpanel(lsbg,param,detect-thresh) \
	--detect-minarea $catpanel(lsbg,param,detect-minarea) \
	--detect-filter-kernel $catpanel(lsbg,param,detect-filter-kernel) \
	--deblend-nthresh $catpanel(lsbg,param,deblend-nthresh) \
	--deblend-mincont $catpanel(lsbg,param,deblend-mincont) \
	--multiscale-factors $catpanel(lsbg,param,multiscale-factors) \
	--sersic-n-min $catpanel(lsbg,param,sersic-n-min) \
	--sersic-n-max $catpanel(lsbg,param,sersic-n-max) \
	--sersic-re-min $catpanel(lsbg,param,sersic-re-min) \
	--sersic-cutout-scale $catpanel(lsbg,param,sersic-cutout-scale) \
	--sersic-max-nfev $catpanel(lsbg,param,sersic-max-nfev) \
	--phot-apertures $catpanel(lsbg,param,phot-apertures) \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--mu-eff-min $catpanel(lsbg,param,mu-eff-min) \
	--mu-eff-max $catpanel(lsbg,param,mu-eff-max) \
	--r-eff-min $catpanel(lsbg,param,r-eff-min) \
	--r-eff-max $catpanel(lsbg,param,r-eff-max) \
	--ellipticity-max $catpanel(lsbg,param,ellipticity-max) \
	--min-snr $catpanel(lsbg,param,min-snr) \
	--sersic-n-filter-min $catpanel(lsbg,param,sersic-n-filter-min) \
	--sersic-n-filter-max $catpanel(lsbg,param,sersic-n-filter-max) \
	--sersic-chi2-max $catpanel(lsbg,param,sersic-chi2-max) \
	--mask-input $catpanel(lsbg,mask_file) \
	--mask-output [OGFMaskRefinedPath lsbg] \
	--masked-output $catpanel(lsbg,masked_file) \
	--bkg-output $catpanel(lsbg,bkg_file) \
	--cleaned-output $catpanel(lsbg,cleaned_file) \
	--segmap-output $catpanel(lsbg,segmap_file) \
	--catalog-output $catpanel(lsbg,catalog_file) \
	--n-workers $catpanel(param,n-workers)]
    if {$catpanel(lsbg,param,lsb-protect)} {
	lappend args --lsb-protect
    } else {
	lappend args --no-lsb-protect
    }
    if {$catpanel(lsbg,param,multiscale)} {
	lappend args --multiscale
    } else {
	lappend args --no-multiscale
    }
    if {$catpanel(lsbg,param,sersic-fit)} {
	lappend args --sersic-fit
    } else {
	lappend args --no-sersic-fit
    }
    if {$catpanel(lsbg,param,svm-classify)} {
	lappend args --svm-classify
	lappend args --svm-threshold $catpanel(lsbg,param,svm-threshold)
	if {$catpanel(lsbg,param,svm-checkpoint) ne {}} {
	    lappend args --svm-checkpoint $catpanel(lsbg,param,svm-checkpoint)
	}
    }

    return [list $args]
}

proc CatalogPanelICLExportScript {} {
    CatalogPanelExportCLIScript icl
}

proc CatalogPanelLSBGExportScript {} {
    CatalogPanelExportCLIScript lsbg
}

proc CatalogPanelImportCLIScript {pipeline} {
    global catpanel current

    # --- Step 0: Select script file first (before requiring FITS) ---
    set types {{"Shell Script" {.sh}} {"All Files" *}}
    set infile [tk_getOpenFile -title "Import $pipeline CLI Script" \
	-filetypes $types]
    if {$infile eq {}} return

    if {[catch {set fd [open $infile r]; set content [read $fd]; close $fd} err]} {
	set catpanel(status) "$pipeline: Cannot read script: $err"
	return
    }

    # --- Get current FITS, or auto-detect from script ---
    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	# Try to find the input FITS path from the script
	# Look for FITS=... variable assignment
	set script_fits {}
	foreach rawline [split $content \n] {
	    set line [string trim $rawline]
	    if {[regexp {^FITS=["']?([^"'\$]+\.fits[^"']*)["']?$} $line -> fpath]} {
		set script_fits $fpath
		break
	    }
	}
	# Also try first positional arg of first python3 command
	if {$script_fits eq {}} {
	    foreach rawline [split $content \n] {
		set line [string trim $rawline]
		if {[string match "python3 *" $line]} {
		    # python3 script.py INPUT.fits ...
		    set parts [split $line " "]
		    if {[llength $parts] >= 3} {
			set candidate [lindex $parts 2]
			if {[string match "*.fits*" $candidate] &&
			    ![string match {$*} $candidate] &&
			    ![string match {--*} $candidate]} {
			    set script_fits $candidate
			    break
			}
		    }
		}
	    }
	}
	# Resolve relative path from script directory
	if {$script_fits ne {} && [file pathtype $script_fits] ne "absolute"} {
	    set script_fits [file join [file dirname $infile] $script_fits]
	}
	if {$script_fits ne {} && [file exists $script_fits]} {
	    # Auto-load the FITS file
	    if {$current(frame) eq {}} { CreateFrame }
	    LoadFitsFile $script_fits {} {}
	    global scale
	    set scale(mode) zscale
	    ChangeScaleMode
	    set fn $script_fits
	} else {
	    # Ask user to select FITS file
	    set ftypes {{"FITS Files" {.fits .fit .fts .fits.gz}} {"All Files" *}}
	    set fn [tk_getOpenFile -title "Select input FITS file" \
		-filetypes $ftypes]
	    if {$fn eq {}} return
	    if {$current(frame) eq {}} { CreateFrame }
	    LoadFitsFile $fn {} {}
	    global scale
	    set scale(mode) zscale
	    ChangeScaleMode
	}
    }

    # Update intermediate file paths for this FITS
    if {$pipeline eq "icl"} {
	CatalogPanelICLUpdateFiles $fn
    } elseif {$pipeline eq "lsbg"} {
	CatalogPanelLSBGUpdateFiles $fn
    }

    # --- Phase 1: Parse shell variable assignments ---
    set shell_vars [dict create]
    foreach rawline [split $content \n] {
	set line [string trim $rawline]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	if {[regexp {^([A-Za-z_][A-Za-z_0-9]*)=(.+)$} $line -> vname vval]} {
	    # Skip command substitutions $(...) and backticks
	    if {[string match {*$(*} $vval] || [string match {*`*} $vval]} continue
	    # Strip surrounding quotes
	    if {([string index $vval 0] eq "\"" && [string index $vval end] eq "\"") ||
		([string index $vval 0] eq "'" && [string index $vval end] eq "'")} {
		set vval [string range $vval 1 end-1]
	    } else {
		# Remove inline comments: "0.263  # arcsec/pixel" → "0.263"
		# Only strip if not inside quotes
		set cidx [string first " #" $vval]
		if {$cidx >= 0} {
		    set vval [string trimright [string range $vval 0 $cidx-1]]
		}
	    }
	    # Expand already-known variables in value
	    dict for {k v} $shell_vars {
		set vval [string map [list "\$$k" $v "\${$k}" $v] $vval]
	    }
	    dict set shell_vars $vname $vval
	}
    }

    # --- Phase 2: Join continuation lines, extract python3 commands ---
    set commands {}
    set buf {}
    set cont 0
    foreach rawline [split $content \n] {
	set trimmed [string trim $rawline]
	if {!$cont} {
	    if {$trimmed eq {} || [string index $trimmed 0] eq "#"} continue
	    set buf $trimmed
	} else {
	    append buf " " $trimmed
	}
	if {[string index $buf end] eq "\\"} {
	    set buf [string range $buf 0 end-1]
	    set cont 1
	    continue
	}
	set cont 0
	if {[string match "python3 *" $buf] &&
	    ![string match "python3 -c *" $buf] &&
	    ![string match "python3 -c*" $buf]} {
	    # Expand shell variables
	    dict for {k v} $shell_vars {
		set buf [string map [list \
		    "\"\$$k\"" $v "\$$k" $v \
		    "\"\${$k}\"" $v "\${$k}" $v] $buf]
	    }
	    lappend commands $buf
	}
	set buf {}
    }

    if {[llength $commands] == 0} {
	set catpanel(status) "$pipeline: No python3 commands found in script"
	return
    }

    # --- Phase 2.5: Check if output files already exist ---
    # Scan commands for output files and check if they exist
    set existing_outputs {}
    foreach cmdline $commands {
	# Look for --output FILE, --*-output FILE, > FILE patterns
	set parts [split $cmdline " "]
	set nparts [llength $parts]
	for {set i 0} {$i < $nparts} {incr i} {
	    set p [lindex $parts $i]
	    set next_i [expr {$i + 1}]
	    if {$next_i < $nparts} {
		set nextarg [lindex $parts $next_i]
		# Strip quotes from nextarg
		set nextarg [string trim $nextarg "\"'"]
		if {([string match "--*-output" $p] || $p eq "--output" || $p eq ">") &&
		    [regexp {\.(tsv|fits)(\.gz)?$} $nextarg]} {
		    # Resolve path: substitute shell vars, check ds9dir
		    set resolved $nextarg
		    dict for {k v} $shell_vars {
			set resolved [string map [list "\$$k" $v "\${$k}" $v] $resolved]
		    }
		    # Also check ~/.ds9/ version (Import redirects outputs there)
		    set ds9ver [file join [file join [file normalize ~] .ds9] [file tail $resolved]]
		    if {[file exists $resolved] && [string match "*.tsv" $resolved]} {
			lappend existing_outputs $resolved
		    } elseif {[file exists $ds9ver] && [string match "*.tsv" $ds9ver]} {
			lappend existing_outputs $ds9ver
		    }
		}
	    }
	}
    }

    # If output TSV files exist, offer to load them directly
    if {[llength $existing_outputs] > 0} {
	set flist [join $existing_outputs "\n  "]
	set msg "Previously generated output files found:\n  $flist\n\nLoad existing results without re-running?"
	set answer [tk_messageBox -type yesnocancel -icon question \
	    -title "$pipeline: Existing Results Found" \
	    -message $msg \
	    -detail "Yes = Load existing files (fast)\nNo = Re-run all steps\nCancel = Abort"]
	if {$answer eq "cancel"} return
	if {$answer eq "yes"} {
	    # Load the best existing TSV (prefer one with X_IMAGE)
	    set best_tsv {}
	    foreach tsvfile $existing_outputs {
		if {[catch {set fd [open $tsvfile r]; set data [read $fd]; close $fd}]} continue
		if {$data eq {}} continue
		set hdr [lindex [split $data \n] 0]
		if {[string match "*X_IMAGE*" $hdr]} {
		    set best_tsv $data
		    break
		}
		if {$best_tsv eq {}} { set best_tsv $data }
	    }
	    if {$best_tsv ne {}} {
		set catpanel(alldata) $best_tsv
		CatalogPanelLoadTSV $best_tsv $pipeline
		CatalogPanelMarkAll

		# Update pipeline state from existing files
		if {$pipeline eq "icl"} {
		    if {[file exists $catpanel(icl,mask_file)]}    { set catpanel(icl,has_mask) 1 }
		    if {[file exists $catpanel(icl,bkg_file)]}     { set catpanel(icl,has_bkg) 1 }
		    if {[file exists $catpanel(icl,profile_file)]} { set catpanel(icl,has_profile) 1 }
		} elseif {$pipeline eq "lsbg"} {
		    if {[file exists $catpanel(lsbg,mask_file)]}    { set catpanel(lsbg,has_mask) 1 }
		    if {[file exists $catpanel(lsbg,cleaned_file)]} { set catpanel(lsbg,has_clean) 1 }
		    if {[file exists $catpanel(lsbg,segmap_file)]}  { set catpanel(lsbg,has_detect) 1 }
		    set catpanel(lsbg,has_catalog) 1
		}

		set n [expr {[llength [split $best_tsv \n]] - 1}]
		set catpanel(status) "$pipeline: Loaded existing results ($n sources)"
		return
	    }
	    set catpanel(status) "$pipeline: Could not read existing files, re-running..."
	}
    }

    # --- Phase 3: Execute commands ---
    set script_path [CatalogPanelGetScript ds9_${pipeline}.py]
    set script_dir [file dirname $script_path]
    set ds9dir [file join [file normalize ~] .ds9]

    # Reset cmdlog for new session
    set catpanel($pipeline,cmdlog) {}

    set total [llength $commands]
    set step 0
    set loaded_tsv 0
    set last_tsv_data {}
    set last_input_fits $fn
    foreach cmdline $commands {
	incr step
	set catpanel(status) "$pipeline: Running step $step/$total..."
	update idletasks

	# Standard variable substitution (Export-generated scripts)
	set cmdline [string map [list \
	    "\$INPUT" $fn "\${INPUT}" $fn \
	    "\$OUTDIR" $ds9dir "\${OUTDIR}" $ds9dir \
	    "\$SCRIPT_DIR" $script_dir "\${SCRIPT_DIR}" $script_dir] $cmdline]

	# Shell-split into args (handle single and double quotes)
	set args {}
	set in_sq 0
	set in_dq 0
	set cur {}
	foreach ch [split $cmdline {}] {
	    if {$ch eq "'" && !$in_dq} {
		set in_sq [expr {!$in_sq}]
	    } elseif {$ch eq "\"" && !$in_sq} {
		set in_dq [expr {!$in_dq}]
	    } elseif {$ch eq " " && !$in_sq && !$in_dq} {
		if {$cur ne {}} { lappend args $cur; set cur {} }
	    } else {
		append cur $ch
	    }
	}
	if {$cur ne {}} { lappend args $cur }

	if {[llength $args] == 0} continue

	# --- Smart path substitution ---
	set new_args {}
	set found_script 0
	set found_input 0
	set prev_flag {}
	set cmd_output_file {}

	foreach arg $args {
	    # 1. Replace script path (any form of ds9_<pipeline>.py)
	    if {!$found_script && [string match "*ds9_${pipeline}.py" $arg]} {
		lappend new_args $script_path
		set found_script 1
		set prev_flag {}
		continue
	    }
	    # Also recognize other known pipeline scripts
	    if {!$found_script && [string match "*.py" $arg] &&
		[string match "*ds9_*" $arg] || [string match "*sep_detect*" $arg]} {
		# Non-ds9 script: resolve relative to script_dir or keep as-is
		if {[file exists $arg]} {
		    lappend new_args $arg
		} else {
		    set resolved [file join [file dirname $script_path] [file tail $arg]]
		    if {[file exists $resolved]} {
			lappend new_args $resolved
		    } else {
			lappend new_args $arg
		    }
		}
		set found_script 1
		set prev_flag {}
		continue
	    }

	    # 2. Redirect --*-output and --output values to ~/.ds9/
	    if {[string match "--*-output" $prev_flag] ||
		$prev_flag eq "--output"} {
		set arg [file join $ds9dir [file tail $arg]]
		lappend new_args $arg
		if {[string match "*.tsv" $arg]} {
		    set cmd_output_file $arg
		}
		set prev_flag {}
		continue
	    }

	    # 3. Redirect shell > output to ~/.ds9/
	    if {$prev_flag eq ">"} {
		set arg [file join $ds9dir [file tail $arg]]
		lappend new_args $arg
		if {[string match "*.tsv" $arg]} {
		    set cmd_output_file $arg
		}
		set prev_flag {}
		continue
	    }

	    # 4. Redirect intermediate file inputs to ~/.ds9/
	    if {$prev_flag in {--mask --cleaned --segmap --catalog
			       --profile-file --import-mask-file
			       --detection --forced} &&
		[regexp {\.(fits|tsv|fit|fts)(\.gz)?$} $arg]} {
		set arg [file join $ds9dir [file tail $arg]]
		lappend new_args $arg
		set prev_flag {}
		continue
	    }

	    # 5. Input file (first positional arg after script)
	    if {$found_script && !$found_input &&
		[string index $arg 0] ne "-" &&
		$arg ne ">"} {
		if {[string match "${ds9dir}/*" $arg]} {
		    # Intermediate file in ~/.ds9/
		    lappend new_args $arg
		} elseif {[file exists $arg]} {
		    # Original file exists, keep it
		    lappend new_args $arg
		    set last_input_fits $arg
		} else {
		    # File doesn't exist, use current FITS
		    lappend new_args $fn
		}
		set found_input 1
		set prev_flag {}
		continue
	    }

	    # 6. Track flag for next arg
	    set prev_flag $arg
	    lappend new_args $arg
	}

	# Log and execute
	CatalogPanelCmdLog $pipeline $new_args

	if {[catch {set result [exec {*}$new_args 2>@stderr]} err]} {
	    set catpanel(status) "$pipeline import: Step $step failed: $err"
	    return
	}

	# Collect TSV: from stdout or from output file
	set tsv {}
	if {[string match "*\t*" $result] && [llength [split $result \n]] > 1} {
	    set tsv $result
	} elseif {$cmd_output_file ne {} && [file exists $cmd_output_file]} {
	    catch {set fd [open $cmd_output_file r]
		set tsv [read $fd]; close $fd}
	}
	# Keep the TSV with X_IMAGE (has coordinates); remember last TSV as fallback
	if {$tsv ne {} && [string match "*\t*" $tsv]} {
	    set hdr [lindex [split $tsv \n] 0]
	    if {[string match "*X_IMAGE*" $hdr]} {
		set last_tsv_data $tsv
	    } elseif {$last_tsv_data eq {}} {
		set last_tsv_data $tsv
	    }
	    set loaded_tsv 1
	}
    }

    # Load the best TSV collected during execution
    if {$last_tsv_data ne {}} {
	set catpanel(alldata) $last_tsv_data
	CatalogPanelLoadTSV $last_tsv_data $pipeline
	CatalogPanelMarkAll
    }

    # Update pipeline state
    if {$pipeline eq "icl"} {
	if {[file exists $catpanel(icl,mask_file)]}    { set catpanel(icl,has_mask) 1 }
	if {[file exists $catpanel(icl,bkg_file)]}     { set catpanel(icl,has_bkg) 1 }
	if {[file exists $catpanel(icl,profile_file)]} { set catpanel(icl,has_profile) 1 }
    } elseif {$pipeline eq "lsbg"} {
	if {[file exists $catpanel(lsbg,mask_file)]}    { set catpanel(lsbg,has_mask) 1 }
	if {[file exists $catpanel(lsbg,cleaned_file)]} { set catpanel(lsbg,has_clean) 1 }
	if {[file exists $catpanel(lsbg,segmap_file)]}  { set catpanel(lsbg,has_detect) 1 }
	set catpanel(lsbg,has_catalog) 1
    }

    # Display appropriate image
    if {$pipeline eq "icl"} {
	# ICL: show bgsub or masked image (better for ICL visualization)
	set display_file {}
	if {[file exists $catpanel(icl,bgsub_file)]} {
	    set display_file $catpanel(icl,bgsub_file)
	} elseif {[file exists $catpanel(icl,masked_file)]} {
	    set display_file $catpanel(icl,masked_file)
	} else {
	    set display_file $fn
	}
	CreateFrame
	if {![catch {LoadFitsFile $display_file {} {}}]} {
	    global scale
	    set scale(mode) zscale
	    ChangeScaleMode
	}
    } else {
	# LSBG: show first input FITS (detection image, coordinates match catalog)
	if {![catch {LoadFitsFile $fn {} {}}]} {
	    global scale
	    set scale(mode) zscale
	    ChangeScaleMode
	}
    }

    set catpanel(status) "$pipeline: Script imported — $total steps executed"
}

proc CatalogPanelICLImportScript {} {
    CatalogPanelImportCLIScript icl
}

proc CatalogPanelLSBGImportScript {} {
    CatalogPanelImportCLIScript lsbg
}

