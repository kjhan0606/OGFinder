# PSF-model plugin (plugins/psfex): Tcl side = remember the outputs, show maps/tables, hand the model to the other steps.
# Catalog keys: psfex,* (registered in ogf_core.tcl) and the existing psf,file (plain PSF image used by the PSF/fit steps)
# plus psf,model (the spatially varying model JSON; the PSF photometry steps add --psf-model when it is set).

proc OGFPsfexFile {name} {return [file join [OGFSessWorkDir] psfex psfex_$name]}

proc OGFPsfexAfter {} {
    foreach {k f} {model_file model.json center_file center.fits stars_file stars.tsv maps_file maps.tsv plot_file diag.png} {
	::ogf::cat::set psfex,$k [OGFPsfexFile $f]
    }
}

# make the model the catalog PSF: psf,file = central stamp (any step that takes a PSF image), psf,model = spatially varying model
proc OGFPsfexUse {} {
    set c [OGFPsfexFile center.fits]
    set m [OGFPsfexFile model.json]
    if {![file exists $c] || ![file exists $m]} {::ogf::status "PSF model: build the model first"; return}
    ::ogf::cat::set psf,file $c
    ::ogf::cat::set psf,model $m
    ::ogf::status "PSF model in use: psf,file = central stamp, psf,model = spatially varying model (PSF photometry, crowded photometry now use it)"
}

proc OGFPsfexText {title name} {
    set f [OGFPsfexFile $name]
    if {![file exists $f]} {::ogf::status "PSF model: $name does not exist yet - build the model first"; return}
    set fd [open $f r]; set txt [read $fd]; close $fd
    OGFTextWindow "$title ($f)" $txt
}
proc OGFPsfexStars {} {OGFPsfexText "PSF stars" stars.tsv}
proc OGFPsfexMaps {} {OGFPsfexText "PSF maps" maps.tsv}

proc OGFPsfexDiag {} {
    set pf [OGFPsfexFile diag.png]
    if {![file exists $pf]} {::ogf::status "PSF model: no plot yet - build the model first"; return}
    set w .ogfpsfexplot
    if {[winfo exists $w]} {destroy $w}
    toplevel $w
    wm title $w "PSF model maps"
    catch {image delete ogfpsfeximg}
    image create photo ogfpsfeximg -file $pf
    label $w.l -image ogfpsfeximg
    pack $w.l -fill both -expand 1
    ttk::button $w.close -text Close -command [list destroy $w]
    pack $w.close -pady 3
    return $w
}
