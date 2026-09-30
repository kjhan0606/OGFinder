#  OGFinder: unified mask manager ("Mask" menu).
#  One bit-flag mask per image (~/.ds9/mask_<base>.fits) edited by
#  ds9_mask.py; displayed through DS9's native mask layer; consumed by the
#  ICL and LSBG steps as ~/.ds9/mask_<base>_bool.fits.
#  See docs/multiband_mask_link.md

package provide DS9 1.0

proc OGFMaskInit {} {
    global catpanel ogfmask
    set ogfmask(overlay) 0
    set ogfmask(color) red
    set ogfmask(transparency) 50
    set ogfmask(stats) {}
    set ogfmask(frame) {}
    # shared Auto Mask parameters.  Defaults = the ICL values; LSBG-specific
    # values are offered as a preset in the Auto Mask dialog.
    set ogfmask(p,detect-thresh) $catpanel(icl,param,detect-thresh)
    set ogfmask(p,minarea) 5
    set ogfmask(p,expand-factor) $catpanel(icl,param,expand-factor)
    set ogfmask(p,max-dilate-radius) $catpanel(icl,param,max-dilate-radius)
    set ogfmask(p,bright-star-mag-limit) $catpanel(icl,param,bright-star-mag-limit)
    set ogfmask(p,bright-star-radius-scale) $catpanel(icl,param,bright-star-radius-scale)
    set ogfmask(p,mag-threshold) 99
    set ogfmask(p,lsb-protect) 0
}

proc OGFMaskBase {fn} {
    return [CatalogPanelFitsBaseName $fn]
}

# Paths for the image in the current frame.  Returns {} when none.
proc OGFMaskPaths {{fn {}}} {
    global catpanel
    if {$fn eq {}} {set fn [CatalogPanelGetFITS]}
    if {$fn eq {}} {return {}}
    set base [OGFMaskBase $fn]
    set d [file join [file normalize ~] .ds9]
    if {![file isdirectory $d]} {file mkdir $d}
    return [dict create fits $fn base $base \
	mask [file join $d "mask_${base}.fits"] \
	bool [file join $d "mask_${base}_bool.fits"] \
	masked [file join $d "mask_${base}_masked.fits"]]
}

proc OGFMaskRun {mode fn extra} {
    global catpanel ogfmask
    set p [OGFMaskPaths $fn]
    if {$p eq {}} {
	set catpanel(status) "Mask: no FITS image loaded"
	return {}
    }
    set script [CatalogPanelGetScript ds9_mask.py]
    if {![file exists $script]} {
	set catpanel(status) "Mask: ds9_mask.py not found"
	return {}
    }
    set args [list [OGFPython] $script [dict get $p fits] --mode $mode \
	--mask [dict get $p mask] {*}$extra]
    if {[info exists catpanel(icl,cmdlog)]} {catch {CatalogPanelCmdLog icl $args}}
    if {[catch {set out [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "Mask $mode error: [string range $err 0 160]"
	return {}
    }
    foreach line [split $out \n] {
	if {[string match "#MASK_STATS*" $line]} {
	    OGFMaskParseStats $line
	    return $line
	}
    }
    return $out
}

proc OGFMaskParseStats {line} {
    global ogfmask catpanel
    foreach kv [lrange $line 1 end] {
	if {[regexp {^([A-Z_]+)=(.*)$} $kv -> k v]} {set ogfmask(s,$k) $v}
    }
    if {![info exists ogfmask(s,N_MASKED)]} return
    set frac [expr {100.0 * $ogfmask(s,FRACTION)}]
    set txt [format "Mask: %d px masked (%.2f%% of image)" $ogfmask(s,N_MASKED) $frac]
    if {[info exists ogfmask(s,FRAC_OF_VALID)]} {
	append txt [format ", %.2f%% of valid" [expr {100.0*$ogfmask(s,FRAC_OF_VALID)}]]
    }
    if {[info exists ogfmask(s,UNDO)]} {append txt "  (undo $ogfmask(s,UNDO) / redo $ogfmask(s,REDO))"}
    set ogfmask(stats) $txt
    set catpanel(status) $txt
}

# Put the shared mask on the native DS9 mask layer of the current frame.
proc OGFMaskRefreshOverlay {} {
    global ogfmask current catpanel
    if {!$ogfmask(overlay)} return
    set fr $current(frame)
    if {$fr eq {} || ![$fr has fits]} return
    set p [OGFMaskPaths]
    if {$p eq {} || ![file exists [dict get $p bool]]} {
	catch {$fr mask clear}
	return
    }
    catch {$fr mask clear}
    if {[catch {LoadFitsFile [dict get $p bool] mask {}} err]} {
	set catpanel(status) "Mask overlay: $err"
	return
    }
    OGFMaskApplyStyle $fr
    set ogfmask(frame) $fr
}

proc OGFMaskApplyStyle {fr} {
    global ogfmask
    catch {
	$fr mask color $ogfmask(color)
	$fr mask mark nonzero
	$fr mask transparency $ogfmask(transparency)
	$fr mask blend screen
    }
}

proc CatalogPanelMaskToggleOverlay {} {
    global ogfmask current catpanel
    set fr $current(frame)
    if {$ogfmask(overlay)} {
	set p [OGFMaskPaths]
	if {$p eq {} || ![file exists [dict get $p mask]]} {
	    set ogfmask(overlay) 0
	    set catpanel(status) "Mask: none yet - use Mask > Auto Mask..."
	    return
	}
	OGFMaskRefreshOverlay
    } else {
	catch {$fr mask clear}
    }
}

proc CatalogPanelMaskOverlaySettings {} {
    global ogfmask current
    set r [OGFForm "Mask Overlay" [list \
	{color "Colour (red, green, cyan, #rrggbb)" $ogfmask(color)} \
	{trans "Transparency 0-100 (%)" $ogfmask(transparency)}]]
    if {$r eq {}} return
    set ogfmask(color) [dict get $r color]
    set t [dict get $r trans]
    if {[string is integer -strict $t]} {set ogfmask(transparency) [expr {max(0,min(100,$t))}]}
    if {$ogfmask(overlay)} {OGFMaskApplyStyle $current(frame)}
}

# ------------------------------------------------------------- auto mask
proc OGFMaskAutoArgs {} {
    global ogfmask catpanel
    set a [list --detect-thresh $ogfmask(p,detect-thresh) \
	--minarea $ogfmask(p,minarea) \
	--expand-factor $ogfmask(p,expand-factor) \
	--max-dilate-radius $ogfmask(p,max-dilate-radius) \
	--bright-star-mag-limit $ogfmask(p,bright-star-mag-limit) \
	--bright-star-radius-scale $ogfmask(p,bright-star-radius-scale) \
	--mag-threshold $ogfmask(p,mag-threshold) \
	--mag-zeropoint $catpanel(param,mag-zeropoint) \
	--n-workers $catpanel(param,n-workers)]
    if {[info exists catpanel(lsbg,param,pixel-scale)]} {
	lappend a --pixel-scale $catpanel(lsbg,param,pixel-scale) \
	    --lsb-mu-threshold $catpanel(lsbg,param,lsb-mu-threshold)
    }
    if {$ogfmask(p,lsb-protect)} {lappend a --lsb-protect}
    return $a
}

proc OGFMaskPresetICL {} {
    global ogfmask catpanel
    set ogfmask(p,detect-thresh) $catpanel(icl,param,detect-thresh)
    set ogfmask(p,expand-factor) $catpanel(icl,param,expand-factor)
    set ogfmask(p,max-dilate-radius) $catpanel(icl,param,max-dilate-radius)
    set ogfmask(p,bright-star-mag-limit) $catpanel(icl,param,bright-star-mag-limit)
    set ogfmask(p,bright-star-radius-scale) $catpanel(icl,param,bright-star-radius-scale)
    set ogfmask(p,mag-threshold) 99
    set ogfmask(p,lsb-protect) 0
    OGFMaskFormRefresh
}
proc OGFMaskPresetLSBG {} {
    global ogfmask catpanel
    set ogfmask(p,detect-thresh) $catpanel(lsbg,param,mask-detect-thresh)
    set ogfmask(p,minarea) $catpanel(lsbg,param,mask-detect-minarea)
    set ogfmask(p,expand-factor) $catpanel(lsbg,param,mask-expand-factor)
    set ogfmask(p,max-dilate-radius) $catpanel(lsbg,param,max-dilate-radius)
    set ogfmask(p,bright-star-mag-limit) $catpanel(lsbg,param,bright-star-mag-limit)
    set ogfmask(p,bright-star-radius-scale) $catpanel(lsbg,param,bright-star-radius-scale)
    set ogfmask(p,mag-threshold) $catpanel(lsbg,param,mask-mag-threshold)
    set ogfmask(p,lsb-protect) $catpanel(lsbg,param,lsb-protect)
    OGFMaskFormRefresh
}
proc OGFMaskFormRefresh {} {
    global ogfmask ogfform_v
    foreach k {detect-thresh minarea expand-factor max-dilate-radius bright-star-mag-limit
	       bright-star-radius-scale mag-threshold lsb-protect} {
	set ogfform_v($k) $ogfmask(p,$k)
    }
}

proc CatalogPanelMaskAutoDialog {} {
    global ogfmask
    set fields [list \
	[list detect-thresh "Detection threshold (sigma)" $ogfmask(p,detect-thresh)] \
	[list minarea "Min area (px)" $ogfmask(p,minarea)] \
	[list expand-factor "Expansion factor" $ogfmask(p,expand-factor)] \
	[list max-dilate-radius "Max dilate radius (px)" $ogfmask(p,max-dilate-radius)] \
	[list bright-star-mag-limit "Bright-star mag limit" $ogfmask(p,bright-star-mag-limit)] \
	[list bright-star-radius-scale "Bright-star radius scale" $ogfmask(p,bright-star-radius-scale)] \
	[list mag-threshold "Mask only brighter than mag (99=all)" $ogfmask(p,mag-threshold)] \
	[list lsb-protect "Protect LSB structures" $ogfmask(p,lsb-protect) check]]
    set r [OGFForm "Auto Mask" $fields [list \
	{"ICL defaults" OGFMaskPresetICL} {"LSBG defaults" OGFMaskPresetLSBG}] \
	"Existing manual add/erase edits are kept (only source and star bits are regenerated)."]
    if {$r eq {}} {return 0}
    foreach k [dict keys $r] {set ogfmask(p,$k) [dict get $r $k]}
    return 1
}

proc CatalogPanelMaskAuto {{ask 1}} {
    global catpanel ogfmask
    if {$ask && ![CatalogPanelMaskAutoDialog]} return
    set catpanel(status) "Mask: running Auto Mask..."
    update idletasks
    set extra [OGFMaskAutoArgs]
    if {[info exists catpanel(alldata)] && $catpanel(alldata) ne {}} {
	set cf [CatalogPanelSaveTempCatalog mask]
	if {$cf ne {}} {lappend extra --catalog $cf}
    }
    set t0 [clock milliseconds]
    set r [OGFMaskRun auto {} $extra]
    if {$r eq {}} return
    set ogfmask(last_ms) [expr {[clock milliseconds] - $t0}]
    OGFMaskAfterEdit
}

# make the ICL/LSBG "has mask" flags and overlay consistent after any edit
proc OGFMaskAfterEdit {} {
    global catpanel ogfmask
    set p [OGFMaskPaths]
    if {$p ne {}} {
	OGFMaskSyncPipelines $p
    }
    OGFMaskRefreshOverlay
}

# Point the ICL / LSBG steps at the shared mask (boolean file).
proc OGFMaskSyncPipelines {p} {
    global catpanel
    foreach pl {icl lsbg} {
	set catpanel($pl,mask_file) [dict get $p bool]
	set catpanel($pl,has_mask) [file exists [dict get $p bool]]
	set catpanel($pl,fits_base_mask) [dict get $p base]
    }
}

# Guarantee a shared mask exists; generate with the current shared
# parameters (seeded from the calling pipeline) when missing.  Returns 1/0.
proc OGFMaskEnsure {pipeline} {
    global catpanel ogfmask
    set p [OGFMaskPaths]
    if {$p eq {}} {return 0}
    if {[file exists [dict get $p mask]] && [file exists [dict get $p bool]]} {
	OGFMaskSyncPipelines $p
	return 1
    }
    if {$pipeline eq "lsbg"} {OGFMaskPresetLSBG} else {OGFMaskPresetICL}
    CatalogPanelMaskAuto 0
    set p [OGFMaskPaths]
    return [file exists [dict get $p bool]]
}

# ----------------------------------------------- region add / erase, edits
proc OGFMaskRegionsFile {} {
    global current catpanel
    set fr $current(frame)
    if {$fr eq {} || ![$fr has fits]} {return {}}
    set f [file join [file normalize ~] .ds9 mask_regions.reg]
    set txt [$fr marker list ds9 image fk5 degrees 0]
    # Region shapes only (drop the DS9 catalog markers and OGFinder tags)
    set out "image\n"
    set n 0
    foreach line [split $txt \n] {
	if {[string match "#*" $line] || [string match "global*" $line]} continue
	if {[string match "*tag=\{sextract*" $line] || [string match "*tag=\{psf_star*" $line]} continue
	if {[regexp {^-?(circle|ellipse|box|polygon|annulus)\(} $line]} {
	    append out $line "\n"
	    incr n
	}
    }
    if {$n == 0} {return {}}
    set fd [open $f w]; puts -nonewline $fd $out; close $fd
    return [list $f $n]
}

proc CatalogPanelMaskRegions {erase} {
    global catpanel
    set rf [OGFMaskRegionsFile]
    if {$rf eq {}} {
	set catpanel(status) "Mask: draw circle/ellipse/box/polygon regions first (Region menu)"
	return
    }
    lassign $rf f n
    set mode [expr {$erase ? "erase" : "add"}]
    set r [OGFMaskRun $mode {} [list --regions $f]]
    if {$r eq {}} return
    append ::catpanel(status) "  ($n regions $mode)"
    OGFMaskAfterEdit
}

proc CatalogPanelMaskGrow {shrink} {
    global catpanel
    set r [OGFForm [expr {$shrink ? "Shrink Mask" : "Grow Mask"}] {{px "Pixels" 3}}]
    if {$r eq {}} return
    set px [dict get $r px]
    if {![string is double -strict $px] || $px <= 0} return
    if {[OGFMaskRun [expr {$shrink ? "shrink" : "grow"}] {} [list --pixels $px]] eq {}} return
    OGFMaskAfterEdit
}

proc CatalogPanelMaskSimple {mode} {
    if {[OGFMaskRun $mode {} {}] eq {}} return
    OGFMaskAfterEdit
}

proc CatalogPanelMaskStats {} {
    if {[OGFMaskRun stats {} {}] eq {}} return
}

proc CatalogPanelMaskSaveAs {} {
    global catpanel
    set p [OGFMaskPaths]
    if {$p eq {} || ![file exists [dict get $p mask]]} {
	set catpanel(status) "Mask: nothing to save"
	return
    }
    set fname [tk_getSaveFile -title "Save Mask As..." \
	-initialfile "[dict get $p base]_mask.fits" \
	-filetypes {{{FITS files} {.fits .fit}} {{All files} *}}]
    if {$fname eq {}} return
    set isbool [expr {[tk_messageBox -type yesno -icon question -title "Save Mask" \
	-message "Save as plain 0/1 mask?" \
	-detail "Yes = boolean mask (readable by any tool).\nNo = bit-flag mask (1 src, 2 star, 4 add, 8 erase, 16 import) - keeps causes."] eq "yes"}]
    set extra [list --file $fname]
    if {$isbool} {lappend extra --boolean}
    if {[OGFMaskRun export {} $extra] eq {}} return
    append ::catpanel(status) "  saved: $fname"
}

proc CatalogPanelMaskImport {} {
    global catpanel
    set fname [tk_getOpenFile -title "Import Mask FITS..." \
	-filetypes {{{FITS files} {.fits .fit}} {{All files} *}}]
    if {$fname eq {}} return
    if {[OGFMaskRun import {} [list --file $fname]] eq {}} return
    OGFMaskAfterEdit
}

# On-demand interpolated image (the only place the masked image is written)
proc OGFMaskMakeMasked {interp} {
    set p [OGFMaskPaths]
    if {$p eq {}} {return {}}
    if {[OGFMaskRun masked {} [list --masked-output [dict get $p masked] \
	--interp-method $interp]] eq {}} {return {}}
    return [dict get $p masked]
}

proc CatalogPanelMaskShowMasked {} {
    global catpanel
    set p [OGFMaskPaths]
    if {$p eq {} || ![file exists [dict get $p mask]]} {
	set catpanel(status) "Mask: none yet - use Mask > Auto Mask..."
	return
    }
    set catpanel(status) "Mask: interpolating masked image..."
    update idletasks
    set f [OGFMaskMakeMasked $catpanel(icl,param,interp-method)]
    if {$f eq {}} return
    CreateFrame
    if {![catch {LoadFitsFile $f {} {}}]} {
	global scale
	set scale(mode) zscale
	ChangeScaleMode
    }
    set catpanel(status) "Mask: masked image shown in a new frame (written on demand to $f)"
}

# Reuse the detection-band mask on every other registered band
proc CatalogPanelMaskCopyToBands {} {
    global catpanel ogfband
    if {![info exists ogfband] || $ogfband(detect) eq {}} {
	set catpanel(status) "Mask: register bands first (Bands menu)"
	return
    }
    set det $ogfband(detect)
    set dp [OGFMaskPaths $ogfband($det,file)]
    if {![file exists [dict get $dp mask]]} {
	set catpanel(status) "Mask: make a mask on the detection band ($det) first"
	return
    }
    set script [CatalogPanelGetScript ds9_mask.py]
    set n 0
    foreach b $ogfband(names) {
	if {$b eq $det} continue
	set bp [OGFMaskPaths $ogfband($b,file)]
	set args [list [OGFPython] $script $ogfband($det,file) --mode reproject \
	    --mask [dict get $dp mask] --target-image $ogfband($b,file) \
	    --output [dict get $bp mask]]
	if {[catch {exec {*}$args 2>@stderr} out]} {
	    set catpanel(status) "Mask copy to $b failed: [string range $out 0 120]"
	    return
	}
	incr n
    }
    set catpanel(status) "Mask from $det reprojected to $n other band(s) (~/.ds9/mask_<band>.fits)"
}

# ------------------------------------------- ICL / LSBG integration
proc OGFMaskExists {} {
    set p [OGFMaskPaths]
    if {$p eq {}} {return 0}
    return [expr {[file exists [dict get $p mask]] && [file exists [dict get $p bool]]}]
}

proc OGFMaskBoolPath {} {
    set p [OGFMaskPaths]
    if {$p eq {}} {return {}}
    return [dict get $p bool]
}

# File written by ICL/LSBG iterative refinement (never overwrites the shared mask)
proc OGFMaskRefinedPath {pipeline} {
    set p [OGFMaskPaths]
    if {$p eq {}} {return [file join [file normalize ~] .ds9 ${pipeline}_mask_refined.fits]}
    return [file join [file dirname [dict get $p mask]] "mask_[dict get $p base]_${pipeline}refined.fits"]
}

# Interpolated image for a downstream step; produced now from the shared mask.
proc OGFMaskMaskedFor {pipeline fn} {
    global catpanel
    set p [OGFMaskPaths $fn]
    if {$p eq {}} {return $fn}
    set interp $catpanel(icl,param,interp-method)
    if {$pipeline eq "lsbg"} {set interp $catpanel(lsbg,param,interp-method)}
    set out [dict get $p masked]
    set catpanel(status) "Mask: producing masked image for the $pipeline step..."
    update idletasks
    set f [OGFMaskMakeMasked $interp]
    if {$f eq {}} {return $fn}
    return $f
}

# ICL "1. Source Masking" / LSBG "1. Mask Bright Sources"
# Reuse the shared mask when it exists (keeps manual edits); otherwise build
# it with the calling method's parameters as defaults for the shared dialog.
proc OGFMaskPipelineMask {pipeline} {
    global catpanel ogfmask
    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "[string toupper $pipeline]: No FITS file loaded"
	return
    }
    CatalogPanelICLUpdateFiles $fn
    CatalogPanelLSBGUpdateFiles $fn
    if {[OGFMaskExists]} {
	OGFMaskSyncPipelines [OGFMaskPaths]
	set ogfmask(overlay) 1
	OGFMaskRefreshOverlay
	OGFMaskRun stats {} {}
	append catpanel(status) "  (existing shared mask reused; Mask > Auto Mask... to regenerate)"
	return
    }
    if {$pipeline eq "lsbg"} {OGFMaskPresetLSBG} else {OGFMaskPresetICL}
    if {![CatalogPanelMaskAutoDialog]} return
    set ogfmask(overlay) 1
    CatalogPanelMaskAuto 0
}
