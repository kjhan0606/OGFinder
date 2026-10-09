#  OGFinder: multi-band manager ("Bands" menu).
#  Band registry, detection band, forced photometry in all bands
#  (ds9_sextract --forced-catalog/--measure-image), colour columns,
#  tiled view and selection that follows the source across band frames.
#  See docs/multiband_mask_link.md

package provide DS9 1.0

proc OGFBandsInit {} {
    global ogfband
    set ogfband(names) {}
    set ogfband(detect) {}
    set ogfband(snrmin) 2.0
}

# ---------------------------------------------------------------- info
proc OGFImageInfo {fn} {
    set bin [OGFSextractBin]
    set d {}
    if {![file exists $bin]} {return $d}
    if {[catch {set out [exec [OGFPython] $bin --info $fn 2>@stderr]}]} {return $d}
    foreach line [split $out \n] {
	if {[regexp {^([A-Z0-9_]+)=(.*)$} $line -> k v]} {dict set d $k [string trim $v]}
    }
    return $d
}

proc OGFSameGrid {a b} {
    foreach k {NAXIS1 NAXIS2} {
	if {![dict exists $a $k] || ![dict exists $b $k] ||
	    [dict get $a $k] != [dict get $b $k]} {return 0}
    }
    foreach k {CRPIX1 CRPIX2 CRVAL1 CRVAL2 CD1_1 CD1_2 CD2_1 CD2_2} {
	if {![dict exists $a $k] || ![dict exists $b $k]} {return 0}
	set x [dict get $a $k]; set y [dict get $b $k]
	if {abs($x - $y) > 1e-9 * (abs($x) + 1e-12) + 1e-12} {return 0}
    }
    return 1
}

proc OGFPix2World {d x y} {
    set pi 3.141592653589793
    set r [expr {$pi/180.0}]
    set dx [expr {$x - [dict get $d CRPIX1]}]; set dy [expr {$y - [dict get $d CRPIX2]}]
    set xi  [expr {($r)*([dict get $d CD1_1]*$dx + [dict get $d CD1_2]*$dy)}]
    set eta [expr {($r)*([dict get $d CD2_1]*$dx + [dict get $d CD2_2]*$dy)}]
    set a0 [expr {[dict get $d CRVAL1]*$r}]; set d0 [expr {[dict get $d CRVAL2]*$r}]
    set den [expr {cos($d0) - $eta*sin($d0)}]
    set a [expr {$a0 + atan2($xi, $den)}]
    set dec [expr {atan2((sin($d0) + $eta*cos($d0))*cos($a-$a0), $den)}]
    return [list $a $dec]
}

proc OGFWorld2Pix {d a dec} {
    set pi 3.141592653589793
    set r [expr {$pi/180.0}]
    set a0 [expr {[dict get $d CRVAL1]*$r}]; set d0 [expr {[dict get $d CRVAL2]*$r}]
    set cc [expr {sin($d0)*sin($dec) + cos($d0)*cos($dec)*cos($a-$a0)}]
    set xi  [expr {cos($dec)*sin($a-$a0)/$cc/$r}]
    set eta [expr {(cos($d0)*sin($dec) - sin($d0)*cos($dec)*cos($a-$a0))/$cc/$r}]
    set c11 [dict get $d CD1_1]; set c12 [dict get $d CD1_2]
    set c21 [dict get $d CD2_1]; set c22 [dict get $d CD2_2]
    set det [expr {$c11*$c22 - $c12*$c21}]
    set dx [expr {( $c22*$xi - $c12*$eta)/$det}]
    set dy [expr {(-$c21*$xi + $c11*$eta)/$det}]
    return [list [expr {$dx + [dict get $d CRPIX1]}] [expr {$dy + [dict get $d CRPIX2]}]]
}

# ------------------------------------------------------------ registry
proc OGFBandOfFrame {frame} {
    global ogfband
    foreach b $ogfband(names) {
	if {$ogfband($b,frame) eq $frame} {return $b}
    }
    return {}
}

# Existing frames that are registered bands (detection band first)
proc OGFBandFrames {{all 1}} {
    global ogfband ds9
    set out {}
    set order $ogfband(names)
    if {$ogfband(detect) ne {}} {
	set order [linsert [lsearch -all -inline -not -exact $order $ogfband(detect)] 0 $ogfband(detect)]
    }
    foreach b $order {
	set f $ogfband($b,frame)
	if {[lsearch -exact $ds9(frames) $f] >= 0} {lappend out $f}
    }
    return $out
}

proc OGFBandsActive {} {
    global ogfband
    return [expr {[llength $ogfband(names)] >= 2}]
}

# Are both frames registered bands (-> the catalog is shared between them)?
proc OGFBandsShareCatalog {a b} {
    global ogfband
    if {[llength $ogfband(names)] < 2} {return 0}
    return [expr {[OGFBandOfFrame $a] ne {} && [OGFBandOfFrame $b] ne {}}]
}

proc OGFBandsRecomputeGrid {} {
    global ogfband
    set det $ogfband(detect)
    foreach b $ogfband(names) {
	if {$det eq {} || $b eq $det} {
	    set ogfband($b,same) 1
	} else {
	    set ogfband($b,same) [OGFSameGrid $ogfband($det,info) $ogfband($b,info)]
	}
    }
}

# Map a detection-grid position into the pixel grid of band frame `frame`.
# Returns {x y same}; same=1 when the grids coincide (identity).
proc OGFBandPos {frame num x y} {
    global ogfband
    set b [OGFBandOfFrame $frame]
    if {$b eq {} || $ogfband(detect) eq {} || $b eq $ogfband(detect) || $ogfband($b,same)} {
	return [list $x $y 1]
    }
    set di $ogfband($ogfband(detect),info)
    set bi $ogfband($b,info)
    if {[dict get $di WCSVALID] != 1 || [dict get $bi WCSVALID] != 1} {return [list $x $y 0]}
    lassign [OGFPix2World $di $x $y] a dec
    lassign [OGFWorld2Pix $bi $a $dec] bx by
    return [list $bx $by 0]
}

proc CatalogPanelBandsRegisterFrame {frame name zp fwhm {file {}}} {
    global ogfband
    if {$file eq {}} {
	catch {set file [string trim [$frame get fits file name full] "{}"]}
	regsub {\[.*\]$} $file {} file
    }
    set info [OGFImageInfo $file]
    set pix {}
    if {[dict exists $info PIXSCALE]} {set pix [dict get $info PIXSCALE]}
    if {[lsearch -exact $ogfband(names) $name] >= 0} {
	# re-register under an existing name: replace
	set ogfband($name,frame) $frame
    } else {
	lappend ogfband(names) $name
    }
    set ogfband($name,frame) $frame
    set ogfband($name,file) $file
    set ogfband($name,zp) $zp
    set ogfband($name,fwhm) $fwhm
    set ogfband($name,pscale) $pix
    set ogfband($name,info) $info
    set ogfband($name,pivot) [expr {[dict exists $info PIVOT] ? [dict get $info PIVOT] : {}}]
    catch {OGFSessLog bands.register config {} -tool native -title "Register band $name (ZP $zp)" \
	-payload [dict create name $name zp $zp fwhm $fwhm file [file normalize $file] key $name]}
    if {$ogfband(detect) eq {}} {set ogfband(detect) $name}
    OGFBandsRecomputeGrid
    return $name
}

proc CatalogPanelBandsRegister {} {
    global current ogfband
    set frame $current(frame)
    if {$frame eq {} || ![$frame has fits]} {
	::ogf::cat::set status "Bands: no image in the current frame"
	return
    }
    set file [CatalogPanelGetFITS]
    set info [OGFImageInfo $file]
    set dname [file rootname [file tail $file]]
    if {[dict exists $info FILTER] && [dict get $info FILTER] ne {}} {set dname [dict get $info FILTER]}
    set dzp [::ogf::cat::get param,mag-zeropoint]
    if {[dict exists $info ZP_AB]} {set dzp [dict get $info ZP_AB]}
    set dps {}
    if {[dict exists $info PIXSCALE]} {set dps [dict get $info PIXSCALE]}
    set r [OGFForm "Register Band ($frame)" [list \
	{name "Band / filter name" $dname} \
	{zp "AB zeropoint" $dzp} \
	{fwhm "PSF FWHM (arcsec, optional)" {}}] {} \
	"File: [file tail $file]    pixel scale from header: $dps\"/px"]
    if {$r eq {}} return
    set d $r
    set name [regsub -all {[^A-Za-z0-9_]} [dict get $d name] _]
    if {$name eq {}} return
    set zp [dict get $d zp]
    if {![string is double -strict $zp]} {::ogf::cat::set status "Bands: bad zeropoint"; return}
    CatalogPanelBandsRegisterFrame $frame $name $zp [dict get $d fwhm] $file
    ::ogf::cat::set status "Band $name registered ($frame, ZP=$zp); [llength $ogfband(names)] band(s), detection band = $ogfband(detect)"
}

proc CatalogPanelBandsLoad {} {
    set fn [tk_getOpenFile -title "Load Band (FITS)..." \
	-filetypes {{{FITS files} {.fits .fit .fits.gz}} {{All files} *}}]
    if {$fn eq {}} return
    CatalogPanelBandsLoadFile $fn
}

proc CatalogPanelBandsLoadFile {fn {name {}} {zp {}}} {
    global current
    set had [expr {$current(frame) ne {} && [$current(frame) has fits]}]
    if {$had} {CreateFrame}
    if {[catch {LoadFitsFile $fn {} {}} err]} {
	::ogf::cat::set status "Bands: cannot load $fn: $err"
	return
    }
    if {$name eq {}} {
	CatalogPanelBandsRegister
    } else {
	set info [OGFImageInfo $fn]
	if {$zp eq {} && [dict exists $info ZP_AB]} {set zp [dict get $info ZP_AB]}
	if {$zp eq {}} {set zp [::ogf::cat::get param,mag-zeropoint]}
	CatalogPanelBandsRegisterFrame $current(frame) $name $zp {} $fn
    }
}

proc CatalogPanelBandsSetDetect {name} {
    global ogfband
    set ogfband(detect) $name
    catch {OGFSessLog bands.detect config {} -tool native -title "Detection band = $name" \
	-payload [dict create name $name]}
    OGFBandsRecomputeGrid
    ::ogf::cat::set status "Detection band = $name"
}

proc CatalogPanelBandsRemove {name} {
    global ogfband
    set i [lsearch -exact $ogfband(names) $name]
    if {$i < 0} return
    set ogfband(names) [lreplace $ogfband(names) $i $i]
    catch {OGFSessLog bands.remove config {} -tool native -title "Remove band $name" \
	-payload [dict create name $name]}
    foreach k [array names ogfband "$name,*"] {unset ogfband($k)}
    if {$ogfband(detect) eq $name} {
	set ogfband(detect) [lindex $ogfband(names) 0]
    }
    OGFBandsRecomputeGrid
    ::ogf::cat::set status "Band $name removed ([llength $ogfband(names)] left)"
}

proc OGFBandFrameDeleted {frame} {
    set b [OGFBandOfFrame $frame]
    if {$b ne {}} {CatalogPanelBandsRemove $b}
}

proc CatalogPanelBandsList {} {
    global ogfband
    if {[llength $ogfband(names)] == 0} {
	OGFTextWindow "Bands" "No bands registered.\nUse Bands > Register Current Frame as Band... or Load Band..."
	return
    }
    set t [format "%-3s %-10s %-8s %-7s %-8s %-7s %-6s %s\n" {} Band Frame ZP_AB FWHM\" Pix\" Grid File]
    foreach b $ogfband(names) {
	set m [expr {$b eq $ogfband(detect) ? "*" : " "}]
	append t [format "%-3s %-10s %-8s %-7s %-8s %-7s %-6s %s\n" $m $b $ogfband($b,frame) \
	    $ogfband($b,zp) $ogfband($b,fwhm) $ogfband($b,pscale) \
	    [expr {$ogfband($b,same) ? "same" : "WCS"}] $ogfband($b,file)]
    }
    append t "\n* = detection band.  'Grid WCS' bands are measured via RA/Dec -> pixel."
    OGFTextWindow "Registered Bands" $t
}

proc OGFBandsPostDetect {m} {
    global ogfband
    $m delete 0 end
    foreach b $ogfband(names) {
	$m add radiobutton -label $b -variable ogfband(detect) -value $b \
	    -command [list CatalogPanelBandsSetDetect $b]
    }
    if {[llength $ogfband(names)] == 0} {$m add command -label "(none)" -state disabled}
}

proc OGFBandsPostRemove {m} {
    global ogfband
    $m delete 0 end
    foreach b $ogfband(names) {
	$m add command -label $b -command [list CatalogPanelBandsRemove $b]
    }
    if {[llength $ogfband(names)] == 0} {$m add command -label "(none)" -state disabled}
}

# ------------------------------------------------------------ detection
proc CatalogPanelBandsDetect {} {
    global ogfband
    if {$ogfband(detect) eq {}} {
	::ogf::cat::set status "Bands: register a band first"
	return
    }
    set b $ogfband(detect)
    set fr $ogfband($b,frame)
    CatalogPanelGotoBandFrame $fr
    ::ogf::cat::set param,mag-zeropoint $ogfband($b,zp)
    if {$ogfband($b,pscale) ne {}} {::ogf::cat::set param,pixel-scale $ogfband($b,pscale)}
    CatalogPanelExtract
    set n [expr {[llength [split [::ogf::cat::tsv] \n]] - 1}]
    ::ogf::cat::set status "Detected $n sources in $b (ZP=$ogfband($b,zp)); use Bands > Measure in All Bands"
}

proc CatalogPanelGotoBandFrame {fr} {
    global current
    if {$current(frame) ne $fr} {GotoFrame $fr}
}

# ------------------------------------------------------------ measure
# Order bands by pivot wavelength when known, else registration order.
proc OGFBandsSorted {} {
    global ogfband
    set known 0
    foreach b $ogfband(names) {if {$ogfband($b,pivot) ne {}} {incr known}}
    if {$known != [llength $ogfband(names)]} {return $ogfband(names)}
    set l {}
    foreach b $ogfband(names) {lappend l [list $b $ogfband($b,pivot)]}
    set out {}
    foreach e [lsort -real -index 1 $l] {lappend out [lindex $e 0]}
    return $out
}

proc CatalogPanelBandsMeasure {{snr {}}} {
    global ogfband
    if {$ogfband(detect) eq {}} {
	::ogf::cat::set status "Bands: register bands and detect first"
	return
    }
    if {[::ogf::cat::tsv] eq {}} {
	::ogf::cat::set status "Bands: no catalog - run Detect first"
	return
    }
    if {$snr eq {}} {
	set r [OGFForm "Measure in All Bands" [list \
	    {snr "Non-detection below S/N" $ogfband(snrmin)}] {} \
	    "Forced photometry at the detection-band positions and Kron apertures.\nMAG = 99 when flux < S/N x error."]
	if {$r eq {}} return
	set snr [dict get $r snr]
    }
    if {![string is double -strict $snr]} return
    set ogfband(snrmin) $snr
    set sbin [OGFSextractBin]
    if {![file exists $sbin]} {::ogf::cat::set status "ogfmeas/sextract.py not found"; return}
    set py [OGFPython]
    set det $ogfband(detect)
    set detfile $ogfband($det,file)
    set catfile [CatalogPanelSaveTempCatalog bands]
    catch {OGFSessLog bands.measure auto {} -tool native -requires {catalog bands>=2} \
	-title "Forced photometry in [llength $ogfband(names)] bands (S/N < $snr -> 99)" \
	-payload [dict create snr $snr bands_list [OGFBandsSorted] detect $det \
	    zps_json "\{[join [lmap b [OGFBandsSorted] {format {"%s": "%s"} $b $ogfband($b,zp)}] {, }]\}"]}
    set t0 [clock milliseconds]
    set summary {}
    foreach b [OGFBandsSorted] {
	::ogf::cat::set status "Measuring in $b ..."
	update idletasks
	set args [list $py $sbin $detfile --forced-catalog $catfile \
	    --measure-image $ogfband($b,file) --band $b \
	    --mag-zeropoint $ogfband($b,zp) --snr-min $snr]
	if {[catch {set out [exec {*}$args 2>@stderr]} err]} {
	    ::ogf::cat::set status "Measure error ($b): $err"
	    return
	}
	set rename [OGFBandsRenameCols $out $b]
	CatalogPanelAddColumnsFromTSV $rename [list MAG_$b MAGERR_$b]
	lappend summary "$b:[OGFBandsCountDet $b]"
    }
    OGFBandsAddColors
    set dt [expr {[clock milliseconds] - $t0}]
    set n [expr {[llength [split [::ogf::cat::tsv] \n]] - 1}]
    ::ogf::cat::set status "Measured $n sources in [llength $ogfband(names)] bands ([expr {$dt/1000.0}] s); detected: [join $summary {  }]"
    # keep the current selection visible
    if {[llength [::ogf::cat::selection]]} {catch {OGFApplySelection 0}}
}

# ds9_sextract forced output -> TSV with NUMBER, MAG_<b>, MAGERR_<b>
proc OGFBandsRenameCols {tsv b} {
    set lines [split $tsv \n]
    set h [split [lindex $lines 0] "\t"]
    set ci_n [lsearch -exact $h NUMBER]
    set ci_m [lsearch -exact $h MAG_AUTO_$b]
    set ci_e [lsearch -exact $h MAGERR_AUTO_$b]
    set out "NUMBER\tMAG_$b\tMAGERR_$b"
    foreach line [lrange $lines 1 end] {
	if {[string trim $line] eq {}} continue
	set f [split $line "\t"]
	append out "\n[lindex $f $ci_n]\t[lindex $f $ci_m]\t[lindex $f $ci_e]"
    }
    return $out
}

proc OGFBandsCountDet {b} {
    set lines [split [::ogf::cat::tsv] \n]
    set h [split [lindex $lines 0] "\t"]
    set ci [lsearch -exact $h MAG_$b]
    if {$ci < 0} {return 0}
    set n 0; set tot 0
    foreach line [lrange $lines 1 end] {
	if {[string trim $line] eq {}} continue
	incr tot
	set v [lindex [split $line "\t"] $ci]
	if {[string is double -strict $v] && $v < 90} {incr n}
    }
    return "$n/$tot"
}

# adjacent-band colour columns (wavelength order); 99 = undefined
proc OGFBandsAddColors {} {
    global ogfband
    set order [OGFBandsSorted]
    if {[llength $order] < 2} return
    set lines [split [::ogf::cat::tsv] \n]
    set h [split [lindex $lines 0] "\t"]
    set pairs {}
    for {set i 0} {$i < [llength $order]-1} {incr i} {
	set a [lindex $order $i]; set b [lindex $order [expr {$i+1}]]
	set ia [lsearch -exact $h MAG_$a]; set ib [lsearch -exact $h MAG_$b]
	if {$ia >= 0 && $ib >= 0} {lappend pairs [list $a $b $ia $ib]}
    }
    if {[llength $pairs] == 0} return
    set ni [lsearch -exact $h NUMBER]
    set res "NUMBER"
    set names {}
    foreach p $pairs {append res "\t[lindex $p 0]-[lindex $p 1]"; lappend names "[lindex $p 0]-[lindex $p 1]"}
    foreach line [lrange $lines 1 end] {
	if {[string trim $line] eq {}} continue
	set f [split $line "\t"]
	append res "\n[lindex $f $ni]"
	foreach p $pairs {
	    lassign $p a b ia ib
	    set ma [lindex $f $ia]; set mb [lindex $f $ib]
	    if {[string is double -strict $ma] && [string is double -strict $mb] && $ma < 90 && $mb < 90} {
		append res [format "\t%.3f" [expr {$ma - $mb}]]
	    } else {
		append res "\t99.000"
	    }
	}
    }
    CatalogPanelAddColumnsFromTSV $res $names
}

# ----------------------------------------------------------------- view
proc CatalogPanelBandsTile {} {
    global ogfband current ds9 tile panzoom crosshair scale
    set frs [OGFBandFrames]
    if {[llength $frs] < 2} {
	::ogf::cat::set status "Bands: register at least two bands to tile"
	return
    }
    # hide non-band frames from the tile, show bands
    set cur $current(frame)
    foreach fr $frs {
	catch {$fr show}
    }
    set current(display) tile
    set tile(mode) grid
    DisplayMode
    foreach fr $frs {
	CatalogPanelGotoBandFrame $fr
	set scale(mode) zscale
	ChangeScaleMode
    }
    CatalogPanelGotoBandFrame [lindex $frs 0]
    set panzoom(lock) wcs
    LockFrameCurrent
    set crosshair(lock) wcs
    LockCrosshairCurrent
    ZoomToFit
    if {[llength [::ogf::cat::selection]]} {catch {OGFApplySelection 1}}
    ::ogf::cat::set status "Tiled [llength $frs] bands (WCS-locked pan/zoom and crosshair)"
}

proc CatalogPanelBandsSingleView {} {
    global current
    set current(display) single
    DisplayMode
}
