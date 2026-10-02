# GUI half of the real-data comparison "legacy dialog/proc path vs headless cli block" (scripts/verify_headless_real.py starts this under Xvfb).
# Runs each legacy step on /workspace/fits/m51.fits with the REAL Python drivers and writes the catalog after the step to $OGF_HLREAL_DIR/gui_<name>.tsv
# (plus the deconvolved image).  HOME must be a scratch directory.  Parameters: see scripts/verify_headless_real.py (same values).
global catpanel current
set ::D $::env(OGF_HLREAL_DIR)
set ::fh [open $::D/gui.log w]
proc L {m} {puts $::fh $m; flush $::fh}
proc bgerror {m} {L "BGERROR $m"}
proc tk_messageBox {args} {L "MESSAGEBOX $args"; return ok}
proc wait_idle {ms} {set t [clock milliseconds]; while {[clock milliseconds]-$t < $ms} {update; after 20}}
proc save {name} {L "status: $::catpanel(status)";set fd [open $::D/gui_$name.tsv w]; fconfigure $fd -translation lf; puts -nonewline $fd $::catpanel(alldata); close $fd; L "saved $name rows=[::ogf::cat::nrows]"}
proc restore {} {::ogf::cat::load_tsv $::BASE restored; update; wait_idle 100}
proc run {} {
    global catpanel
    update; wait_idle 500
    set catpanel(param,detect-thresh) 2.5; set catpanel(param,mag-zeropoint) 25.0; set catpanel(param,n-workers) 1
    CatalogPanelExtract
    set t0 [clock milliseconds]; while {$catpanel(alldata) eq {} && [clock milliseconds]-$t0 < 60000} {update; after 100}
    wait_idle 500; save extract; set ::BASE $catpanel(alldata)
    # dual extract (measurement image = a copy of m51 scaled by 0.5)
    after 500 {set ::ed(dual,measure) $::D/m51_half.fits; set ::ed(ok) 1}
    CatalogPanelDualExtract; wait_idle 300; save dual
    restore
    # completeness
    after 500 {set ::ed(comp,ninject) 30; set ::ed(comp,magmin) 20.0; set ::ed(comp,magmax) 24.0; set ::ed(comp,nbins) 4; set ::ed(ok) 1}
    CatalogPanelCompleteness; wait_idle 300; save completeness
    restore
    # multiband: bands = scaled copies of m51 (flux x0.5, x0.25)
    after 500 {.multibandphot.param.tbands insert 1.0 "A:$::D/m51_half.fits\nB:$::D/m51_quarter.fits\n"; set ::ed(ok) 1}
    CatalogPanelMultiBand; wait_idle 300; save multiband
    restore
    # cross-match (live VizieR)
    after 500 {set ::ed(xm,catalog) GAIA_DR3; set ::ed(xm,radius) 2.0; set ::ed(ok) 1}
    CatalogPanelCrossMatch; wait_idle 300; save crossmatch
    # photo-z and SED on the cross-matched catalog (it has MAG_AUTO / MAG_APER)
    CatalogPanelPhotoZ; set d .catphotoz; $d.bands delete 0 end; $d.bands insert 0 g,r; $d.mags delete 0 end; $d.mags insert 0 MAG_AUTO,MAG_APER; CatalogPanelPhotoZRun $d; wait_idle 300; save photoz
    CatalogPanelSEDFit; set d .catsedfit; $d.backend set auto; $d.bands delete 0 end; $d.bands insert 0 g,r; $d.mags delete 0 end; $d.mags insert 0 MAG_AUTO,MAG_APER; $d.pzcol delete 0 end; $d.pzcol insert 0 PHOTO_Z; CatalogPanelSEDFitRun $d; wait_idle 300; save sed
    restore
    # LSBG full pipeline (ensures the shared mask with the LSBG preset first), real ds9_mask.py / ds9_lsbg.py
    set catpanel(lsbg,param,sersic-fit) 0; set catpanel(lsbg,param,multiscale) 0; set catpanel(lsbg,param,pixel-scale) 0.05
    foreach {k v} {mu-eff-min 5 mu-eff-max 40 r-eff-min 0.01 r-eff-max 5000 ellipticity-max 1.0 min-snr 0 sersic-n-filter-min 0 sersic-n-filter-max 50 sersic-chi2-max 1e9} {set catpanel(lsbg,param,$k) $v}
    CatalogPanelLSBGRunAll; wait_idle 500; save lsbg
    foreach {src dst} [list [file join [OGFSessWorkDir] mask_m51_bool.fits] gui_mask_bool.fits [file join [OGFSessWorkDir] lsbg_cleaned_m51.fits] gui_lsbg_cleaned.fits] {catch {file copy -force $src $::D/$dst}}
    for {set i 0} {$i < 4 && [file tail [CatalogPanelGetFITS]] ne "m51.fits"} {incr i} {catch {DeleteCurrentFrame}; wait_idle 200}
    L "image after lsbg: [CatalogPanelGetFITS]"
    # deconvolution (Richardson-Lucy, 6 iterations) with a Gaussian PSF file prepared by the Python half
    set catpanel(psf,file) $::D/psf.fits; set catpanel(psf,has_psf) 1; set catpanel(psf,param,rl-iterations) 6
    file delete -force [file join [OGFSessWorkDir] deconv_result.fits]
    CatalogPanelDeconvolve rl; wait_idle 500
    file copy -force [file join [OGFSessWorkDir] deconv_result.fits] $::D/gui_deconv.fits
    L "image after deconv: [CatalogPanelGetFITS]"
    catch {DeleteCurrentFrame}; wait_idle 300
    L "image after delete: [CatalogPanelGetFITS]"
    # segmentation map
    CatalogPanelSegmentationMap; wait_idle 500
    file copy -force [file join [OGFSessWorkDir] segmap.fits] $::D/gui_segmap.fits
    L "done"; close $::fh; exit
}
after 3000 run
