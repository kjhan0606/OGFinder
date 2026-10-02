# Extended morphology (plugins/morphology, step morph_ext): outputs of ds9_morph_ext.py -> catalog keys morphext,*, growth-curve plot / table windows.

proc OGFMorphExtFile {name} {return [file join [OGFSessWorkDir] morph_ext $name]}

proc OGFMorphExtAfter {} {
    ::ogf::cat::set morphext,growth_file [OGFMorphExtFile morph_ext_growth.tsv]
    ::ogf::cat::set morphext,plot_file [OGFMorphExtFile morph_ext_curves.png]
}

proc OGFMorphExtTable {} {
    set f [OGFMorphExtFile morph_ext_growth.tsv]
    if {![file exists $f]} {::ogf::status "Morphology: no growth curves yet - run Extended Morphology first"; return}
    set fd [open $f r]; set txt [read $fd]; close $fd
    OGFTextWindow "Growth curves ($f)" $txt
}

proc OGFMorphExtCurves {} {
    set png [OGFMorphExtFile morph_ext_curves.png]
    if {![file exists $png]} {::ogf::status "Morphology: no plot yet - run Extended Morphology first"; return}
    set w .ogfmorphextplot
    if {[winfo exists $w]} {destroy $w}
    toplevel $w
    wm title $w "Growth curves and Petrosian eta"
    catch {image delete ogfmorphextimg}
    image create photo ogfmorphextimg -file $png
    label $w.l -image ogfmorphextimg
    pack $w.l -fill both -expand 1
    ttk::button $w.close -text Close -command [list destroy $w]
    pack $w.close -pady 3
    return $w
}
