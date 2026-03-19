#  Copyright (C) 1999-2021
#  Smithsonian Astrophysical Observatory, Cambridge, MA, USA
#  For conditions of distribution and use, see copyright notice in "copyright"

package provide DS9 1.0

proc CreateHeader {} {
    global ds9

    # Panel Frame
    set ds9(header) [ttk::frame $ds9(main).header]
    set ds9(header,sep) [ttk::separator $ds9(main).sheader -orient horizontal]
}

proc CanvasDef {} {
    global canvas
    global ds9

    switch $ds9(wm) {
	x11 {set canvas(width) 738}
	aqua {set canvas(width) 777}
	win32 {set canvas(width) 740}
    }
    set canvas(height) 528
    set canvas(gap) 4

    switch $ds9(wm) {
	x11 {
	    # this is not fool proof. it does not take into account redirecting
	    # the DISPLAY. There must be a better way.
	    global tcl_platform
	    switch -- $tcl_platform(os) {
		Darwin {set canvas(gap,bottom) 14}
		default {set canvas(gap,bottom) 0}
	    }
	}
	aqua  {set canvas(gap,bottom) 14}
	win32 {set canvas(gap,bottom) 0}
    }
}

proc BlinkDef {} {
    global blink
    global iblink
    global pblink

    set iblink(id) {}
    set iblink(index) -1

    set blink(interval) 1000

    array set pblink [array get blink]}

proc FadeDef {} {
    global fade
    global ifade
    global pfade

    set ifade(id) {}
    set ifade(index) -1
    set ifade(alpha) 0

    set fade(blend) screen
    set fade(interval) 2000
    set fade(step) 25

    array set pfade [array get fade]
}

proc TileDef {} {
    global tile
    global itile
    global ptile

    set itile(top) .tile
    set itile(mb) .tilemb

    set tile(mode) grid
    set tile(grid,row) 10
    set tile(grid,col) 10
    set tile(grid,mode) automatic
    set tile(grid,dir) x
    set tile(grid,gap) 4

    array set ptile [array get tile]
}

proc ViewDef {} {
    global view
    global pview

    set view(layout) horizontal
    set view(multi) 1
    set view(info) 1
    set view(panner) 1
    set view(magnifier) 1
    set view(buttons) 1
    set view(icons) 1
    set view(colorbar) 1
    set view(graph,horz) 0
    set view(graph,vert) 0

    set view(info,filename) 1
    set view(info,object) 1
    set view(info,keyvalue) {}
    set view(info,keyword) 0
    set view(info,minmax) 0
    set view(info,lowhigh) 0
    set view(info,bunit) 0
    set view(info,wcs) 1
    foreach l {a b c d e f g h i j k l m n o p q r s t u v w x y z} {
	set "view(info,wcs$l)" 0
    }
    set view(info,detector) 0
    set view(info,amplifier) 0
    set view(info,physical) 1
    set view(info,image) 1
    set view(info,frame) 1

    array set pview [array get view]
}

# canvas

proc CreateCanvas {} {
    global ds9
    global canvas

    set ds9(image) [ttk::frame $ds9(main).f]

    set ds9(canvas) [canvas $ds9(image).c \
			 -width $canvas(width) \
			 -height $canvas(height) \
			 -highlightthickness 0 \
			 -insertofftime 0 \
			 -bg [ThemeTreeBackground] \
			]
    grid rowconfigure $ds9(image) 0 -weight 1
    grid columnconfigure $ds9(image) 0 -weight 1
    grid $ds9(canvas) -row 0 -column 0 -sticky news

    # extra space for window tab
    set ds9(canvas,bottom) {}
    if {$canvas(gap,bottom)>0} {
	set ds9(canvas,bottom) [ttk::frame $ds9(image).b \
				    -width 1 \
				    -height $canvas(gap,bottom) \
				    -style Tree.TFrame \
				   ]
	grid $ds9(canvas,bottom) -row 1 -column 0 -sticky ew
    }

    # needed to realize window so Layout routines will work
    grid $ds9(image)

    switch $ds9(wm) {
	x11 -
	win32 {bind $ds9(canvas) <<ThemeChanged>> {ThemeConfigCanvas %W}}
	aqua {}
    }
}

proc CreateCatalogPanel {} {
    global ds9
    global catpanel

    set f $ds9(catalog_frame)

    # Menubar with dropdown menus
    set catpanel(menubar) [ttk::frame $f.menubar]

    # Flat menubuttons: no border, no indicator arrow
    ttk::style layout CatMenu.TMenubutton {
	Menubutton.focus -sticky nswe -children {
	    Menubutton.padding -sticky we -children {
		Menubutton.label -side left -sticky {}
	    }
	}
    }
    ttk::style configure CatMenu.TMenubutton -relief flat
    ttk::style map CatMenu.TMenubutton -relief {
	pressed flat
	active  flat
    }

    # SExtractor menu
    ttk::menubutton $f.menubar.sextract -text "SExtractor" \
	-menu $f.menubar.sextract.m -style CatMenu.TMenubutton
    menu $f.menubar.sextract.m -tearoff 0
    $f.menubar.sextract.m add command -label "Extract" \
	-command CatalogPanelExtract
    $f.menubar.sextract.m add command -label "Dual-Image Extract..." \
	-command CatalogPanelDualExtract
    $f.menubar.sextract.m add command -label "Settings..." \
	-command CatalogPanelSettingsDialog
    $f.menubar.sextract.m add separator
    $f.menubar.sextract.m add command -label "Trim..." \
	-command CatalogPanelTrimDialog
    $f.menubar.sextract.m add separator
    $f.menubar.sextract.m add command -label "Save Catalog" \
	-command CatalogPanelSaveCatalog
    $f.menubar.sextract.m add command -label "Export Regions (.reg)" \
	-command CatalogPanelExportRegions
    $f.menubar.sextract.m add command -label "Export FITS Table" \
	-command CatalogPanelExportFITS
    $f.menubar.sextract.m add command -label "Load Catalog" \
	-command CatalogPanelLoadCatalog
    $f.menubar.sextract.m add separator
    $f.menubar.sextract.m add command -label "Clear" \
	-command CatalogPanelClear

    # Display menu
    ttk::menubutton $f.menubar.display -text "Display" \
	-menu $f.menubar.display.m -style CatMenu.TMenubutton
    menu $f.menubar.display.m -tearoff 0
    $f.menubar.display.m add command -label "Mark All" \
	-command CatalogPanelMarkAll
    $f.menubar.display.m add command -label "Clear Markers" \
	-command CatalogPanelClearMarkers
    $f.menubar.display.m add separator
    $f.menubar.display.m add checkbutton -label "Show Visible Only" \
	-variable catpanel(visible_mode) -command CatalogPanelShowVisible
    $f.menubar.display.m add separator
    $f.menubar.display.m add checkbutton -label "Add Objects (Click + A)" \
	-variable catpanel(add_objects_mode)
    $f.menubar.display.m add command -label "Delete Selected (Click + D)" \
	-command CatalogPanelDeleteSelected
    $f.menubar.display.m add separator
    $f.menubar.display.m add command \
	-label "Separate Selected (Click + S)" \
	-command CatalogPanelSeparateSelected
    $f.menubar.display.m add command \
	-label "Separate Settings..." \
	-command CatalogPanelSeparateSettings
    $f.menubar.display.m add command \
	-label "Save Separated Catalog..." \
	-command CatalogPanelSeparateSave
    $f.menubar.display.m add command \
	-label "Load Separated Catalog..." \
	-command CatalogPanelSeparateLoad

    # AI Merge menu
    ttk::menubutton $f.menubar.aimerge -text "AI Merge" \
	-menu $f.menubar.aimerge.m -style CatMenu.TMenubutton
    menu $f.menubar.aimerge.m -tearoff 0
    $f.menubar.aimerge.m add command -label "Run AI Merge" \
	-command CatalogPanelAIMerge

    # Galaxy Model menu
    ttk::menubutton $f.menubar.galaxy -text "Galaxy Model" \
	-menu $f.menubar.galaxy.m -style CatMenu.TMenubutton
    menu $f.menubar.galaxy.m -tearoff 0
    $f.menubar.galaxy.m add command -label "AI Morphology Classification" \
	-command CatalogPanelGalaxyMorphology
    $f.menubar.galaxy.m add separator
    $f.menubar.galaxy.m add command -label "Fit Elliptical Model" \
	-command [list CatalogPanelGalaxyFit elliptical]
    $f.menubar.galaxy.m add command -label "Fit Spiral Model" \
	-command [list CatalogPanelGalaxyFit spiral]
    $f.menubar.galaxy.m add separator
    $f.menubar.galaxy.m add command -label "Extract Parameters" \
	-command CatalogPanelGalaxyParams

    # Star(PSF) menu
    ttk::menubutton $f.menubar.starpsf -text "Star(PSF)" \
	-menu $f.menubar.starpsf.m -style CatMenu.TMenubutton
    menu $f.menubar.starpsf.m -tearoff 0

    # Star Finding
    menu $f.menubar.starpsf.m.find -tearoff 0
    $f.menubar.starpsf.m add cascade -label "Find Stars" \
	-menu $f.menubar.starpsf.m.find
    $f.menubar.starpsf.m.find add command \
	-label "Combined" \
	-command [list CatalogPanelFindStars combined]
    $f.menubar.starpsf.m.find add command \
	-label "CLASS_STAR" \
	-command [list CatalogPanelFindStars class_star]
    $f.menubar.starpsf.m.find add command \
	-label "FWHM" \
	-command [list CatalogPanelFindStars fwhm]
    $f.menubar.starpsf.m add separator
    $f.menubar.starpsf.m add command \
	-label "Show Stars" -command CatalogPanelShowStars
    $f.menubar.starpsf.m add command \
	-label "Clear Stars" -command CatalogPanelClearStars
    $f.menubar.starpsf.m add separator

    # PSF Generation
    menu $f.menubar.starpsf.m.psf -tearoff 0
    $f.menubar.starpsf.m add cascade -label "Build PSF" \
	-menu $f.menubar.starpsf.m.psf
    $f.menubar.starpsf.m.psf add command \
	-label "Median Stack" \
	-command [list CatalogPanelBuildPSF median]
    $f.menubar.starpsf.m.psf add command \
	-label "Moffat Fit" \
	-command [list CatalogPanelBuildPSF moffat]
    $f.menubar.starpsf.m.psf add command \
	-label "Gaussian Fit" \
	-command [list CatalogPanelBuildPSF gaussian]
    $f.menubar.starpsf.m.psf add command \
	-label "ePSF" \
	-command [list CatalogPanelBuildPSF epsf]
    $f.menubar.starpsf.m.psf add separator
    $f.menubar.starpsf.m.psf add command \
	-label "Extended PSF..." \
	-command CatalogPanelBuildExtendedPSF
    $f.menubar.starpsf.m.psf add separator
    $f.menubar.starpsf.m.psf add command \
	-label "WebbPSF (JWST)..." \
	-command CatalogPanelSimPSFWebbPSF
    $f.menubar.starpsf.m.psf add command \
	-label "TinyTim (HST)..." \
	-command CatalogPanelSimPSFTinyTim
    $f.menubar.starpsf.m add separator
    $f.menubar.starpsf.m add command \
	-label "View PSF" -command CatalogPanelViewPSF
    $f.menubar.starpsf.m add command \
	-label "Save PSF..." -command CatalogPanelSavePSF
    $f.menubar.starpsf.m add command \
	-label "Load PSF..." -command CatalogPanelLoadPSF
    $f.menubar.starpsf.m add command \
	-label "AI Star Classification" -command CatalogPanelStarFinder
    $f.menubar.starpsf.m add separator
    $f.menubar.starpsf.m add command -label "Settings..." \
	-command CatalogPanelStarPSFSettings

    # Deconvolution menu
    ttk::menubutton $f.menubar.deconv -text "Deconvolution" \
	-menu $f.menubar.deconv.m -style CatMenu.TMenubutton
    menu $f.menubar.deconv.m -tearoff 0
    $f.menubar.deconv.m add command \
	-label "Richardson-Lucy" \
	-command [list CatalogPanelDeconvolve rl]
    $f.menubar.deconv.m add command \
	-label "Richardson-Lucy (Accelerated)" \
	-command [list CatalogPanelDeconvolve rl_accelerated]
    $f.menubar.deconv.m add command \
	-label "Richardson-Lucy (Regularized)" \
	-command [list CatalogPanelDeconvolve rl_tv]
    $f.menubar.deconv.m add separator
    $f.menubar.deconv.m add command \
	-label "Wiener" \
	-command [list CatalogPanelDeconvolve wiener]
    $f.menubar.deconv.m add command \
	-label "Tikhonov" \
	-command [list CatalogPanelDeconvolve tikhonov]
    $f.menubar.deconv.m add separator
    $f.menubar.deconv.m add command \
	-label "CLEAN" \
	-command [list CatalogPanelDeconvolve clean]
    $f.menubar.deconv.m add command \
	-label "Maximum Entropy (MEM)" \
	-command [list CatalogPanelDeconvolve mem]
    $f.menubar.deconv.m add separator
    $f.menubar.deconv.m add command -label "Settings..." \
	-command CatalogPanelDeconvSettings
    $f.menubar.deconv.m add command -label "Quick Deconvolve (RL)" \
	-command CatalogPanelQuickDeconvolve

    # Separate items moved into Display menu above

    # ICL menu
    ttk::menubutton $f.menubar.icl -text "ICL" \
	-menu $f.menubar.icl.m -style CatMenu.TMenubutton
    menu $f.menubar.icl.m -tearoff 0
    $f.menubar.icl.m add command \
	-label "1. Source Masking" \
	-command CatalogPanelICLMask
    $f.menubar.icl.m add command \
	-label "   View Mask" \
	-command CatalogPanelICLViewMask
    $f.menubar.icl.m add separator
    menu $f.menubar.icl.m.bkg -tearoff 0
    $f.menubar.icl.m add cascade -label "2. Background Model" \
	-menu $f.menubar.icl.m.bkg
    $f.menubar.icl.m.bkg add command \
	-label "Polynomial Fit" \
	-command [list CatalogPanelICLBackground polynomial]
    $f.menubar.icl.m.bkg add command \
	-label "Chebyshev Fit" \
	-command [list CatalogPanelICLBackground chebyshev]
    $f.menubar.icl.m.bkg add command \
	-label "SEP Large Mesh" \
	-command [list CatalogPanelICLBackground sep_large]
    $f.menubar.icl.m add command \
	-label "   View Background" \
	-command CatalogPanelICLViewBkg
    $f.menubar.icl.m add separator
    $f.menubar.icl.m add command \
	-label "3. Set BCG Center (Click)" \
	-command CatalogPanelICLSetCenter
    $f.menubar.icl.m add command \
	-label "   Measure Profile" \
	-command CatalogPanelICLProfile
    $f.menubar.icl.m add command \
	-label "   Sector Profile..." \
	-command CatalogPanelICLSectorProfile
    $f.menubar.icl.m add separator
    $f.menubar.icl.m add command \
	-label "4. ICL Measurements" \
	-command CatalogPanelICLMeasure
    $f.menubar.icl.m add command \
	-label "   Multi-Threshold ICL" \
	-command CatalogPanelICLMeasureMulti
    $f.menubar.icl.m add separator
    $f.menubar.icl.m add command \
	-label "5. BCG+ICL Decomposition" \
	-command CatalogPanelICLDecompose
    $f.menubar.icl.m add separator
    $f.menubar.icl.m add command \
	-label "Color Profile..." \
	-command CatalogPanelICLColorProfile
    $f.menubar.icl.m add separator
    $f.menubar.icl.m add command \
	-label "Save Profile..." \
	-command CatalogPanelICLSaveProfile
    $f.menubar.icl.m add command \
	-label "Load Profile..." \
	-command CatalogPanelICLLoadProfile
    $f.menubar.icl.m add separator
    $f.menubar.icl.m add command \
	-label "Settings..." \
	-command CatalogPanelICLSettings

    # LSBG menu
    ttk::menubutton $f.menubar.lsbg -text "LSBG" \
	-menu $f.menubar.lsbg.m -style CatMenu.TMenubutton
    menu $f.menubar.lsbg.m -tearoff 0
    $f.menubar.lsbg.m add command \
	-label "1. Mask Bright Sources" \
	-command CatalogPanelLSBGMask
    $f.menubar.lsbg.m add command \
	-label "   View Masked Image" \
	-command CatalogPanelLSBGViewMask
    $f.menubar.lsbg.m add separator
    menu $f.menubar.lsbg.m.bkg -tearoff 0
    $f.menubar.lsbg.m add cascade -label "2. Background Model" \
	-menu $f.menubar.lsbg.m.bkg
    $f.menubar.lsbg.m.bkg add command \
	-label "SEP Large Mesh" \
	-command [list CatalogPanelLSBGClean sep_large]
    $f.menubar.lsbg.m.bkg add command \
	-label "Polynomial Fit" \
	-command [list CatalogPanelLSBGClean polynomial]
    $f.menubar.lsbg.m.bkg add command \
	-label "Chebyshev Fit" \
	-command [list CatalogPanelLSBGClean chebyshev]
    $f.menubar.lsbg.m add command \
	-label "   View Cleaned Image" \
	-command CatalogPanelLSBGViewClean
    $f.menubar.lsbg.m add separator
    $f.menubar.lsbg.m add command \
	-label "3. Detect LSBG Candidates" \
	-command CatalogPanelLSBGDetect
    $f.menubar.lsbg.m add separator
    $f.menubar.lsbg.m add command \
	-label "4. Photometry" \
	-command CatalogPanelLSBGPhotometry
    $f.menubar.lsbg.m add separator
    $f.menubar.lsbg.m add command \
	-label "5. Sérsic Profile Fit" \
	-command CatalogPanelLSBGSersic
    $f.menubar.lsbg.m add separator
    $f.menubar.lsbg.m add command \
	-label "6. Filter + Grade" \
	-command CatalogPanelLSBGFilter
    $f.menubar.lsbg.m add separator
    $f.menubar.lsbg.m add command \
	-label "Save Catalog..." \
	-command CatalogPanelSaveCatalog
    $f.menubar.lsbg.m add command \
	-label "Load Catalog..." \
	-command CatalogPanelLoadCatalog
    $f.menubar.lsbg.m add separator
    $f.menubar.lsbg.m add command \
	-label "7. Forced Photometry (Multi-Band)" \
	-command CatalogPanelLSBGForcedPhot
    $f.menubar.lsbg.m add separator
    $f.menubar.lsbg.m add command \
	-label "Run Full Pipeline" \
	-command CatalogPanelLSBGRunAll
    $f.menubar.lsbg.m add separator
    $f.menubar.lsbg.m add command \
	-label "Settings..." \
	-command CatalogPanelLSBGSettings

    # Analysis menu
    ttk::menubutton $f.menubar.analysis -text "Analysis" \
	-menu $f.menubar.analysis.m -style CatMenu.TMenubutton
    menu $f.menubar.analysis.m -tearoff 0
    $f.menubar.analysis.m add command \
	-label "Non-Parametric Morphology (CAS/Gini/M20)" \
	-command CatalogPanelMorphometry
    $f.menubar.analysis.m add command \
	-label "Sérsic Fitting" \
	-command CatalogPanelSersicFit
    $f.menubar.analysis.m add separator
    $f.menubar.analysis.m add command \
	-label "PSF Photometry" \
	-command CatalogPanelPSFPhotometry
    $f.menubar.analysis.m add command \
	-label "Multi-Band Photometry..." \
	-command CatalogPanelMultiBand
    $f.menubar.analysis.m add command \
	-label "Crowded Field Photometry" \
	-command CatalogPanelCrowdedPhot
    $f.menubar.analysis.m add separator
    $f.menubar.analysis.m add command \
	-label "Cross-Match (VizieR)..." \
	-command CatalogPanelCrossMatch
    $f.menubar.analysis.m add separator
    $f.menubar.analysis.m add command \
	-label "Segmentation Map" \
	-command CatalogPanelSegmentationMap
    $f.menubar.analysis.m add command \
	-label "Completeness Simulation..." \
	-command CatalogPanelCompleteness

    pack $f.menubar.sextract -side left
    pack $f.menubar.display -side left
    pack $f.menubar.aimerge -side left
    pack $f.menubar.galaxy -side left
    pack $f.menubar.starpsf -side left
    pack $f.menubar.deconv -side left
    pack $f.menubar.icl -side left
    pack $f.menubar.lsbg -side left
    pack $f.menubar.analysis -side left

    # Search/Filter bar
    set catpanel(searchbar) [ttk::frame $f.searchbar]
    ttk::label $f.searchbar.lbl -text "Filter:"
    set catpanel(search_var) {}
    ttk::entry $f.searchbar.entry -textvariable catpanel(search_var) -width 20
    ttk::button $f.searchbar.go -text "Apply" \
	-command CatalogPanelFilter -width 6
    pack $f.searchbar.lbl -side left -padx 4 -pady 2
    pack $f.searchbar.entry -side left -padx 2 -pady 2 -fill x -expand true
    pack $f.searchbar.go -side right -padx 2 -pady 2
    bind $f.searchbar.entry <Return> CatalogPanelFilter

    # Table frame with scrollbars
    set catpanel(tblframe) [ttk::frame $f.tblf]

    set catpanel(tbldb) catpaneltbldb
    global $catpanel(tbldb)

    set catpanel(tbl) [table $f.tblf.t \
			   -state disabled \
			   -usecommand 0 \
			   -variable $catpanel(tbldb) \
			   -colorigin 1 \
			   -roworigin 0 \
			   -cols 19 \
			   -rows 20 \
			   -width -1 \
			   -height -1 \
			   -colwidth 11 \
			   -maxwidth 0 \
			   -maxheight 0 \
			   -titlerows 1 \
			   -resizeborders col \
			   -xscrollcommand [list $f.tblf.xscroll set] \
			   -yscrollcommand [list $f.tblf.yscroll set] \
			   -selecttype row \
			   -selectmode browse \
			   -browsecommand [list CatalogPanelSelectCmd %s %S] \
			   -anchor w \
			   -font [font actual TkDefaultFont] \
			   -fg [ThemeTreeForeground] \
			   -bg [ThemeTreeBackground] \
			  ]

    $catpanel(tbl) tag configure sel \
	-fg [ThemeSelectedForeground] -bg [ThemeSelectedBackground]
    $catpanel(tbl) tag configure title \
	-fg [ThemeForeground] -bg [ThemeBackground]

    ttk::scrollbar $f.tblf.yscroll \
	-command [list $catpanel(tbl) yview] -orient vertical
    ttk::scrollbar $f.tblf.xscroll \
	-command [list $catpanel(tbl) xview] -orient horizontal

    grid $catpanel(tbl) $f.tblf.yscroll -sticky news
    grid $f.tblf.xscroll -sticky news
    grid rowconfigure $f.tblf 0 -weight 1
    grid columnconfigure $f.tblf 0 -weight 1

    # Status bar
    set catpanel(status) {Ready - Load a FITS file to extract sources}
    set catpanel(statusbar) [ttk::frame $f.statusbar]
    ttk::label $f.statusbar.lbl -textvariable catpanel(status) \
	-anchor w -relief sunken
    pack $f.statusbar.lbl -fill x -expand true -padx 2 -pady 2

    # Pack all into catalog frame
    pack $f.menubar -fill x -side top
    pack $f.searchbar -fill x -side top
    pack $f.statusbar -fill x -side bottom
    pack $f.tblf -fill both -expand true -side top

    # Initialize state
    set catpanel(alldata) {}
    set catpanel(filename) {}
    set catpanel(delim) "\t"
    set catpanel(sort,col) {}
    set catpanel(sort,dir) {}

    # Feature B: visible mode
    set catpanel(visible_mode) 0

    # Add Objects mode
    set catpanel(add_objects_mode) 0

    # Mark All cached region string
    set catpanel(markall,on) 0

    # Feature C: merge state
    set catpanel(merge,list) {}
    set catpanel(merge,active) 0

    # Feature D: trim state
    set catpanel(trim,active) 0

    # AI Merge state
    set catpanel(ai,groups) {}
    set catpanel(ai,current) 0
    set catpanel(ai,total) 0
    set catpanel(ai,threshold) 0.7
    set catpanel(ai,active) 0

    # PSF/Deconv state
    set catpanel(psf,stars) {}
    set catpanel(psf,star_indices) {}
    set catpanel(psf,file) [file join [file normalize ~] .ds9 psf_current.fits]
    set catpanel(psf,has_psf) 0
    set catpanel(psf,param,class-star-thresh) 0.8
    set catpanel(psf,param,max-ellipticity) 0.2
    set catpanel(psf,param,fwhm-sigma) 2.0
    set catpanel(psf,param,min-flux-snr) 10.0
    set catpanel(psf,param,psf-size) 51
    set catpanel(psf,param,rl-iterations) 30
    set catpanel(psf,param,wiener-nsr) 0.01
    set catpanel(psf,param,tikhonov-lambda) 0.001
    set catpanel(psf,param,tv-lambda) 0.001
    set catpanel(psf,param,clean-gain) 0.1
    set catpanel(psf,param,clean-niter) 1000
    set catpanel(psf,param,clean-threshold) 0.0
    set catpanel(psf,param,mem-lambda) 0.1
    set catpanel(psf,param,mem-niter) 100

    # Extended PSF params
    set catpanel(psf,param,ext-core-mag-min)     18.0
    set catpanel(psf,param,ext-core-mag-max)     22.0
    set catpanel(psf,param,ext-wing-mag-max)     16.0
    set catpanel(psf,param,ext-core-size)        51
    set catpanel(psf,param,ext-wing-size)        201
    set catpanel(psf,param,ext-blend-inner)      20.0
    set catpanel(psf,param,ext-blend-outer)      30.0
    set catpanel(psf,param,ext-saturation-limit) 60000.0

    # Simulation PSF params
    set catpanel(psf,param,sim-telescope)        auto
    set catpanel(psf,param,sim-instrument)       auto
    set catpanel(psf,param,sim-filter)           auto
    set catpanel(psf,param,sim-psf-size)         201
    set catpanel(psf,param,sim-oversample)       1
    set catpanel(psf,param,sim-jitter-sigma)     0.007
    set catpanel(psf,param,sim-focus-offset)     0.0

    # Simulation availability flags (-1 = unchecked)
    set catpanel(psf,sim_webbpsf_ok) -1
    set catpanel(psf,sim_tinytim_ok) -1

    CatalogPanelPSFParamLoad

    # ICL state
    set catpanel(icl,mask_file)    [file join [file normalize ~] .ds9 icl_mask.fits]
    set catpanel(icl,masked_file)  [file join [file normalize ~] .ds9 icl_masked.fits]
    set catpanel(icl,bkg_file)     [file join [file normalize ~] .ds9 icl_background.fits]
    set catpanel(icl,bgsub_file)   [file join [file normalize ~] .ds9 icl_bgsub.fits]
    set catpanel(icl,profile_file) [file join [file normalize ~] .ds9 icl_profile.tsv]
    set catpanel(icl,has_mask)     0
    set catpanel(icl,has_bkg)      0
    set catpanel(icl,has_profile)  0
    set catpanel(icl,center_x)     {}
    set catpanel(icl,center_y)     {}
    set catpanel(icl,param,expand-factor)          2.5
    set catpanel(icl,param,bright-star-mag-limit)  18.0
    set catpanel(icl,param,bright-star-radius-scale) 10.0
    set catpanel(icl,param,interp-method)          linear
    set catpanel(icl,param,detect-thresh)          1.5
    set catpanel(icl,param,bkg-method)             polynomial
    set catpanel(icl,param,bkg-order)              3
    set catpanel(icl,param,bkg-sigma-clip)         3.0
    set catpanel(icl,param,bkg-sep-mesh)           256
    set catpanel(icl,param,rmin)                   5.0
    set catpanel(icl,param,rmax)                   1000.0
    set catpanel(icl,param,nsteps)                 80
    set catpanel(icl,param,spacing)                log
    set catpanel(icl,param,ellipticity)            0.0
    set catpanel(icl,param,pa)                     0.0
    set catpanel(icl,param,mag-zeropoint)          25.0
    set catpanel(icl,param,pixel-scale)            0.06
    set catpanel(icl,param,mu-threshold)           26.5
    set catpanel(icl,param,mu-levels)              26.0,27.0,28.0
    set catpanel(icl,param,measure-radius)         500.0
    set catpanel(icl,param,bkg-iterative)          0
    set catpanel(icl,param,bkg-n-iterations)       3
    set catpanel(icl,param,bkg-convergence-tol)    0.01
    set catpanel(icl,param,bkg-refine-thresh)      2.0
    CatalogPanelICLParamLoad

    # LSBG state
    set catpanel(lsbg,mask_file)    [file join [file normalize ~] .ds9 lsbg_mask.fits]
    set catpanel(lsbg,masked_file)  [file join [file normalize ~] .ds9 lsbg_masked.fits]
    set catpanel(lsbg,bkg_file)     [file join [file normalize ~] .ds9 lsbg_background.fits]
    set catpanel(lsbg,cleaned_file) [file join [file normalize ~] .ds9 lsbg_cleaned.fits]
    set catpanel(lsbg,segmap_file)  [file join [file normalize ~] .ds9 lsbg_segmap.fits]
    set catpanel(lsbg,catalog_file) [file join [file normalize ~] .ds9 lsbg_catalog.tsv]
    set catpanel(lsbg,has_mask)     0
    set catpanel(lsbg,has_clean)    0
    set catpanel(lsbg,has_detect)   0
    set catpanel(lsbg,has_catalog)  0
    set catpanel(lsbg,detect_data)  {}
    set catpanel(lsbg,param,mask-detect-thresh)         1.5
    set catpanel(lsbg,param,mask-detect-minarea)        5
    set catpanel(lsbg,param,mask-expand-factor)         3.0
    set catpanel(lsbg,param,bright-star-mag-limit)      18.0
    set catpanel(lsbg,param,bright-star-radius-scale)   12.0
    set catpanel(lsbg,param,mask-mag-threshold)         22.0
    set catpanel(lsbg,param,interp-method)              linear
    set catpanel(lsbg,param,lsb-protect)                1
    set catpanel(lsbg,param,lsb-mu-threshold)           24.0
    set catpanel(lsbg,param,bkg-method)                 sep_large
    set catpanel(lsbg,param,bkg-mesh-size)              256
    set catpanel(lsbg,param,bkg-poly-order)             3
    set catpanel(lsbg,param,bkg-sigma-clip)             3.0
    set catpanel(lsbg,param,bkg-n-iterations)           3
    set catpanel(lsbg,param,bkg-refine-thresh)          2.0
    set catpanel(lsbg,param,bkg-rms-quantile)           0.25
    set catpanel(lsbg,param,bkg-convergence-tol)        0.01
    set catpanel(lsbg,param,detect-thresh)              0.8
    set catpanel(lsbg,param,detect-minarea)             50
    set catpanel(lsbg,param,detect-filter-kernel)       gauss5x5
    set catpanel(lsbg,param,deblend-nthresh)            32
    set catpanel(lsbg,param,deblend-mincont)            0.005
    set catpanel(lsbg,param,multiscale)                 1
    set catpanel(lsbg,param,multiscale-factors)         1,2,4
    set catpanel(lsbg,param,sersic-fit)                 1
    set catpanel(lsbg,param,sersic-n-min)               0.2
    set catpanel(lsbg,param,sersic-n-max)               10.0
    set catpanel(lsbg,param,sersic-re-min)              0.5
    set catpanel(lsbg,param,sersic-cutout-scale)        5.0
    set catpanel(lsbg,param,sersic-max-nfev)            500
    set catpanel(lsbg,param,phot-apertures)             5,10,20,40
    set catpanel(lsbg,param,mag-zeropoint)              25.0
    set catpanel(lsbg,param,pixel-scale)                0.06
    set catpanel(lsbg,param,mu-eff-min)                 24.0
    set catpanel(lsbg,param,mu-eff-max)                 30.0
    set catpanel(lsbg,param,r-eff-min)                  2.5
    set catpanel(lsbg,param,r-eff-max)                  60.0
    set catpanel(lsbg,param,ellipticity-max)            0.7
    set catpanel(lsbg,param,min-snr)                    2.0
    set catpanel(lsbg,param,sersic-n-filter-min)        0.3
    set catpanel(lsbg,param,sersic-n-filter-max)        6.0
    set catpanel(lsbg,param,sersic-chi2-max)            10.0
    CatalogPanelLSBGParamLoad

    # Ctrl key tracking (Feature A/C)
    set ::catpanel_ctrl 0
    bind . <KeyPress-Control_L>   {set ::catpanel_ctrl 1}
    bind . <KeyRelease-Control_L> {set ::catpanel_ctrl 0}
    bind . <KeyPress-Control_R>   {set ::catpanel_ctrl 1}
    bind . <KeyRelease-Control_R> {set ::catpanel_ctrl 0}

    # Key bindings (Feature C)
    bind . <Control-Key-m> {CatalogPanelMergeSources}
    bind . <Escape> {+CatalogPanelEscapeKey}

    # Bind table header click for sorting (ButtonRelease to not conflict with tktable)
    bind $catpanel(tbl) <ButtonRelease-1> {+CatalogPanelTableClick %x %y}

    # Mouse wheel scroll for catalog table (natural/macOS direction)
    bind $catpanel(tbl) <Button-4> {
	%W yview scroll 3 units
	break
    }
    bind $catpanel(tbl) <Button-5> {
	%W yview scroll -3 units
	break
    }
    # Horizontal scroll (Shift + wheel, natural/macOS direction)
    bind $catpanel(tbl) <Shift-Button-4> {
	%W xview scroll 3 units
	break
    }
    bind $catpanel(tbl) <Shift-Button-5> {
	%W xview scroll -3 units
	break
    }

    # Initialize extraction parameters
    CatalogPanelParamDef

    # Force ttk widgets to redraw on resize (X11 compositing conflict)
    bind $f <Configure> [list CatalogPanelRedrawTtk $f]
}

proc CatalogPanelRedrawTtk {f} {
    # Debounce: cancel previous scheduled redraw
    catch {after cancel $::catpanel_redraw_id}

    # Schedule redraw after resize settles
    set ::catpanel_redraw_id [after 50 [list CatalogPanelRedrawTtkDo $f]]
}

proc CatalogPanelRedrawTtkDo {f} {
    # Force titlebar and statusbar to re-expose
    catch {
	foreach w [winfo children $f.titlebar] {
	    event generate $w <Expose>
	}
	event generate $f.statusbar.lbl <Expose>
    }
}

# Run source extraction on the currently loaded FITS image
proc CatalogPanelExtract {} {
    global catpanel
    global ds9
    global current
    global loadParam

    # Find the ds9_sextract binary (platform-aware)
    set bindir [file dirname [info nameofexecutable]]
    set os $::tcl_platform(os)

    if {$os eq "Windows NT"} {
	set sextract [file join $bindir ds9_sextract.exe]
    } else {
	set sextract [file join $bindir ds9_sextract]
    }
    if {![file executable $sextract]} {
	set catpanel(status) "ERROR: ds9_sextract not found in $bindir"
	return
    }

    # Get current FITS filename
    set fn {}
    if {$current(frame) != {}} {
	catch {set fn [$current(frame) get fits file name full]}
    }
    if {$fn eq {}} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    # Strip curly braces if present
    set fn [string trim $fn "{}"]

    # Strip FITS HDU extension specifier (e.g. [SCI], [1], [SCI,2])
    regsub {\[.*\]$} $fn {} fn

    if {![file exists $fn]} {
	set catpanel(status) "File not found: $fn"
	return
    }

    # Set log scale with optimized limits
    CatalogPanelSetLogScale

    set catpanel(status) "Extracting sources from [file tail $fn] ..."
    update idletasks

    # Save extraction parameters so AI Merge uses the same values
    foreach pname {detect-thresh detect-minarea deblend-nthresh deblend-mincont \
		   mag-zeropoint back-size back-filtersize} {
	if {[info exists catpanel(param,$pname)]} {
	    set catpanel(extract_param,$pname) $catpanel(param,$pname)
	}
    }

    # Build parameter arguments list
    set paramargs {}
    foreach pname {detect-thresh detect-minarea deblend-nthresh deblend-mincont \
		   phot-aperture mag-zeropoint gain pixel-scale seeing-fwhm \
		   back-size back-filtersize \
		   phot-aperture-2 phot-aperture-3 phot-aperture-5 conv-filter} {
	if {[info exists catpanel(param,$pname)]} {
	    lappend paramargs "--$pname" $catpanel(param,$pname)
	}
    }

    # Platform-specific library path setup and execution
    if {$os eq "Darwin"} {
	# macOS: set DYLD_LIBRARY_PATH
	set libpaths {}
	if {[info exists ::env(CONDA_PREFIX)]} {
	    lappend libpaths "$::env(CONDA_PREFIX)/lib"
	}
	set home_conda [file join [file normalize ~] miniconda3/lib]
	if {[file isdirectory $home_conda]} {
	    lappend libpaths $home_conda
	}
	if {[info exists ::env(DYLD_LIBRARY_PATH)]} {
	    lappend libpaths $::env(DYLD_LIBRARY_PATH)
	}
	if {[llength $libpaths] > 0} {
	    set ::env(DYLD_LIBRARY_PATH) [join $libpaths :]
	}
    } elseif {$os ne "Windows NT"} {
	# Linux/Unix: set LD_LIBRARY_PATH
	set libpaths {}
	if {[info exists ::env(CONDA_PREFIX)]} {
	    lappend libpaths "$::env(CONDA_PREFIX)/lib"
	}
	set home_conda [file join [file normalize ~] miniconda3/lib]
	if {[file isdirectory $home_conda]} {
	    lappend libpaths $home_conda
	}
	if {[info exists ::env(LD_LIBRARY_PATH)]} {
	    lappend libpaths $::env(LD_LIBRARY_PATH)
	}
	if {[llength $libpaths] > 0} {
	    set ::env(LD_LIBRARY_PATH) [join $libpaths :]
	}
    }
    # Windows: DLLs found via PATH automatically

    # Run extraction (cross-platform exec)
    if {[catch {set data [exec $sextract $fn {*}$paramargs 2>@stderr]} err]} {
	set catpanel(status) "Extraction error: $err"
	return
    }

    # Parse TSV output into table
    CatalogPanelLoadTSV $data [file tail $fn]
}

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

    global $catpanel(tbldb)
    $catpanel(tbl) configure -variable {}
    unset -nocomplain $catpanel(tbldb)
    $catpanel(tbl) configure -variable $catpanel(tbldb) \
	-cols 19 -rows 20

    set catpanel(status) {Ready}
    set catpanel(filename) {}
    set catpanel(alldata) {}

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

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
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

    CatalogPanelLoadTSV $data [file tail $fn]
}

proc CatalogPanelFilter {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set pattern $catpanel(search_var)

    global $catpanel(tbldb)

    # Unbind table while modifying
    $catpanel(tbl) configure -variable {}
    unset -nocomplain $catpanel(tbldb)

    set data $catpanel(alldata)
    set lines [split $data \n]

    # Header
    set headers [split [lindex $lines 0] "\t"]
    set ncols [llength $headers]

    for {set c 0} {$c < $ncols} {incr c} {
	set ${catpanel(tbldb)}(0,[expr {$c+1}]) \
	    [string trim [lindex $headers $c]]
    }

    # Filter data rows
    set row 1
    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	if {$pattern ne {} && ![string match -nocase "*${pattern}*" $line]} continue
	set fields [split $line "\t"]
	for {set c 0} {$c < $ncols} {incr c} {
	    set ${catpanel(tbldb)}($row,[expr {$c+1}]) \
		[string trim [lindex $fields $c]]
	}
	incr row
    }

    # Rebind table to trigger full refresh
    $catpanel(tbl) configure -variable $catpanel(tbldb) \
	-cols $ncols -rows $row

    set ndata [expr {$row - 1}]
    if {$pattern eq {}} {
	set catpanel(status) "Showing all $ndata sources"
    } else {
	set catpanel(status) "Filtered: $ndata sources matching '$pattern'"
    }
}

# Hook: automatically extract sources after FITS file is loaded
proc CatalogPanelAutoExtract {} {
    global catpanel
    if {[info exists catpanel(tbl)]} {
	after 500 CatalogPanelExtract
    }
}

# Row selection: navigate to source and mark it on the image
proc CatalogPanelSelectCmd {prev cur} {
    global catpanel

    # cur is "row,col" of current selection
    set row [lindex [split $cur ,] 0]
    if {![string is integer -strict $row] || $row <= 0} return

    after cancel CatalogPanelGotoSource
    after 100 [list CatalogPanelGotoSource $row]
}

proc CatalogPanelGotoSource {row} {
    global catpanel
    global current
    global ds9

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return

    global $catpanel(tbldb)

    # Find column indices from header row
    set ncols [$catpanel(tbl) cget -cols]
    set col_x -1
    set col_y -1
    set col_a -1
    set col_b -1
    set col_theta -1
    set col_ir -1
    set col_reff -1
    for {set c 1} {$c <= $ncols} {incr c} {
	if {[info exists ${catpanel(tbldb)}(0,$c)]} {
	    set hdr [set ${catpanel(tbldb)}(0,$c)]
	    switch -- $hdr {
		X_IMAGE     { set col_x $c }
		Y_IMAGE     { set col_y $c }
		A_IMAGE     { set col_a $c }
		B_IMAGE     { set col_b $c }
		THETA_IMAGE { set col_theta $c }
		ISO_RADIUS  { set col_ir $c }
		R_EFF_PIX   { set col_reff $c }
	    }
	}
    }

    if {$col_x < 0 || $col_y < 0} return

    # Get coordinates from selected row
    if {![info exists ${catpanel(tbldb)}($row,$col_x)]} return
    set x [set ${catpanel(tbldb)}($row,$col_x)]
    set y [set ${catpanel(tbldb)}($row,$col_y)]

    if {![string is double -strict $x] || ![string is double -strict $y]} return

    # Get ellipse parameters (with NaN/Inf safety via catch)
    set iso_radius 10.0
    set a_image 0
    set b_image 0
    set theta 0

    if {$col_ir >= 0 && [info exists ${catpanel(tbldb)}($row,$col_ir)]} {
	set val [set ${catpanel(tbldb)}($row,$col_ir)]
	if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} {
	    set iso_radius $v
	}
    }
    if {$col_a >= 0 && [info exists ${catpanel(tbldb)}($row,$col_a)]} {
	set val [set ${catpanel(tbldb)}($row,$col_a)]
	if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set a_image $v }
    }
    if {$col_b >= 0 && [info exists ${catpanel(tbldb)}($row,$col_b)]} {
	set val [set ${catpanel(tbldb)}($row,$col_b)]
	if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set b_image $v }
    }
    if {$col_theta >= 0 && [info exists ${catpanel(tbldb)}($row,$col_theta)]} {
	set val [set ${catpanel(tbldb)}($row,$col_theta)]
	if {[catch {set v [expr {$val + 0.0}]}] == 0} { set theta $v }
    }

    # Get R_EFF_PIX if available (for LSBG circle display)
    set r_eff_pix 0
    if {$col_reff >= 0 && [info exists ${catpanel(tbldb)}($row,$col_reff)]} {
	set val [set ${catpanel(tbldb)}($row,$col_reff)]
	if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} {
	    set r_eff_pix $v
	}
    }

    # Compute ellipse: ISO_RADIUS as semi-major, scaled by B/A for semi-minor
    set semi_a $iso_radius
    set semi_b $iso_radius
    if {$a_image > 0 && $b_image > 0} {
	set semi_b [expr {$iso_radius * $b_image / $a_image}]
    }

    # Delete previous selection markers
    set frame $current(frame)
    catch {$frame marker catalog sextract_sel delete}

    # Rebuild sextract_all markers from alldata to keep image in sync
    if {[info exists catpanel(markall,on)] && $catpanel(markall,on)} {
	CatalogPanelCreateAllMarkers
    }

    # Use global variable for marker creation (var form requires global access)
    global sextract_sel_reg

    # Create cross point marker (cyan)
    set sextract_sel_reg "image\ncross point($x $y) # color=cyan width=2 point=cross 15 tag={sextract_sel} select=0 edit=0 move=0 rotate=0 delete=1\n"
    catch {$frame marker catalog command ds9 var sextract_sel_reg}

    # If R_EFF_PIX is available, show a circle with effective radius
    if {$r_eff_pix > 0} {
	set sextract_sel_reg "image\ncircle($x $y ${r_eff_pix}i) # color=green width=2 dash=1 tag={sextract_sel} select=0 edit=0 move=0 rotate=0 delete=1\n"
	catch {$frame marker catalog command ds9 var sextract_sel_reg}
    }

    # Create ellipse marker (green, dashed)
    set sextract_sel_reg "image\nellipse($x $y ${semi_a}i ${semi_b}i $theta) # color=green width=2 dash=1 tag={sextract_sel} select=0 edit=0 move=0 rotate=0 delete=1\n"
    catch {$frame marker catalog command ds9 var sextract_sel_reg}

    # Pan to the object
    PanToFrame $current(frame) $x $y image {}

    set catpanel(status) "Source at image ($x, $y)"
}

# --- Source Extractor Parameter Management ---

proc CatalogPanelParamDef {} {
    global catpanel

    set catpanel(param,detect-thresh) 1.5
    set catpanel(param,detect-minarea) 5
    set catpanel(param,deblend-nthresh) 32
    set catpanel(param,deblend-mincont) 0.005
    set catpanel(param,phot-aperture) 5.0
    set catpanel(param,mag-zeropoint) 25.0
    set catpanel(param,gain) 0.0
    set catpanel(param,pixel-scale) 1.0
    set catpanel(param,seeing-fwhm) 3.0
    set catpanel(param,back-size) 64
    set catpanel(param,back-filtersize) 3
    set catpanel(param,phot-aperture-2) 4.0
    set catpanel(param,phot-aperture-3) 6.0
    set catpanel(param,phot-aperture-5) 10.0
    set catpanel(param,conv-filter) default
    set catpanel(param,n-workers) 0

    # Separate (deblend) parameters
    set catpanel(param,sep-deblend-nthresh) 64
    set catpanel(param,sep-deblend-mincont) 0.0001
    set catpanel(param,sep-detect-thresh) 0.8
    set catpanel(param,sep-detect-minarea) 3
    set catpanel(param,sep-radius-factor) 3.0
    set catpanel(param,sep-back-size) 32

    CatalogPanelParamLoad
}

proc CatalogPanelParamLoad {} {
    global catpanel

    set preffile [file join [file normalize ~] .ds9 sextract.prf]
    if {![file exists $preffile]} return
    if {[catch {set fd [open $preffile r]} err]} return
    while {[gets $fd line] >= 0} {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set parts [split $line]
	if {[llength $parts] >= 2} {
	    set key [lindex $parts 0]
	    set val [lindex $parts 1]
	    if {[info exists catpanel(param,$key)]} {
		set catpanel(param,$key) $val
	    }
	}
    }
    close $fd
}

proc CatalogPanelParamSave {} {
    global catpanel

    set prefdir [file join [file normalize ~] .ds9]
    if {![file isdirectory $prefdir]} {
	file mkdir $prefdir
    }
    set preffile [file join $prefdir sextract.prf]
    if {[catch {set fd [open $preffile w]} err]} return
    foreach pname {detect-thresh detect-minarea deblend-nthresh deblend-mincont \
		   phot-aperture mag-zeropoint gain pixel-scale seeing-fwhm \
		   back-size back-filtersize \
		   phot-aperture-2 phot-aperture-3 phot-aperture-5 conv-filter \
		   n-workers \
		   sep-deblend-nthresh sep-deblend-mincont sep-detect-thresh \
		   sep-detect-minarea sep-radius-factor sep-back-size} {
	puts $fd "$pname $catpanel(param,$pname)"
    }
    close $fd
}

proc CatalogPanelParamDefaults {} {
    global ed

    set ed(detect-thresh) 1.5
    set ed(detect-minarea) 5
    set ed(deblend-nthresh) 32
    set ed(deblend-mincont) 0.005
    set ed(phot-aperture) 5.0
    set ed(mag-zeropoint) 25.0
    set ed(gain) 0.0
    set ed(pixel-scale) 1.0
    set ed(seeing-fwhm) 3.0
    set ed(back-size) 64
    set ed(back-filtersize) 3
    set ed(phot-aperture-2) 4.0
    set ed(phot-aperture-3) 6.0
    set ed(phot-aperture-5) 10.0
    set ed(conv-filter) default
    set ed(n-workers) 0
}

proc CatalogPanelSettingsDialog {} {
    global catpanel
    global ed

    set w {.sextractparam}

    set ed(ok) 0

    # Copy current params to ed()
    foreach pname {detect-thresh detect-minarea deblend-nthresh deblend-mincont \
		   phot-aperture mag-zeropoint gain pixel-scale seeing-fwhm \
		   back-size back-filtersize \
		   phot-aperture-2 phot-aperture-3 phot-aperture-5 conv-filter \
		   n-workers} {
	set ed($pname) $catpanel(param,$pname)
    }

    DialogCreate $w {Source Extractor Settings} ed(ok)

    # Param frame
    set f [ttk::frame $w.param]
    set row 0
    foreach {pname plabel} {
	detect-thresh {Detect Threshold}
	detect-minarea {Detect Min Area}
	deblend-nthresh {Deblend NThresh}
	deblend-mincont {Deblend MinCont}
	phot-aperture {Phot Aperture (diam)}
	phot-aperture-2 {Phot Aperture 2}
	phot-aperture-3 {Phot Aperture 3}
	phot-aperture-5 {Phot Aperture 5}
	mag-zeropoint {Mag Zeropoint}
	gain {Gain}
	pixel-scale {Pixel Scale}
	seeing-fwhm {Seeing FWHM}
	back-size {Back Size}
	back-filtersize {Back Filter Size}
    } {
	ttk::label $f.l$row -text "$plabel:" -anchor w
	ttk::entry $f.e$row -textvariable ed($pname) -width 12
	grid $f.l$row $f.e$row -padx 4 -pady 2 -sticky w
	incr row
    }
    # Convolution filter dropdown
    ttk::label $f.l$row -text "Conv Filter:" -anchor w
    ttk::combobox $f.e$row -textvariable ed(conv-filter) -width 12 \
	-values {default gauss5x5 mexhat tophat} -state readonly
    grid $f.l$row $f.e$row -padx 4 -pady 2 -sticky w
    incr row
    # Parallel workers
    ttk::label $f.l$row -text "Parallel Workers (0=auto):" -anchor w
    ttk::entry $f.e$row -textvariable ed(n-workers) -width 12
    grid $f.l$row $f.e$row -padx 4 -pady 2 -sticky w
    incr row

    # Buttons
    set bf [ttk::frame $w.buttons]
    ttk::button $bf.ok -text {OK} -command {set ed(ok) 1} -default active
    ttk::button $bf.cancel -text {Cancel} -command {set ed(ok) 0}
    ttk::button $bf.defaults -text {Defaults} -command CatalogPanelParamDefaults
    ttk::button $bf.save -text {Save} -command {
	foreach pname {detect-thresh detect-minarea deblend-nthresh deblend-mincont \
		       phot-aperture mag-zeropoint gain pixel-scale seeing-fwhm \
		       back-size back-filtersize \
		       phot-aperture-2 phot-aperture-3 phot-aperture-5 conv-filter \
		       n-workers} {
	    set catpanel(param,$pname) $ed($pname)
	}
	CatalogPanelParamSave
    }
    pack $bf.ok $bf.cancel $bf.defaults $bf.save \
	-side left -expand true -padx 2 -pady 4

    bind $w <Return> {set ed(ok) 1}

    # Fini
    ttk::separator $w.sep -orient horizontal
    pack $w.buttons $w.sep -side bottom -fill x
    pack $w.param -side top -fill both -expand true

    DialogWait $w ed(ok) $w.param.e0
    destroy $w

    if {$ed(ok)} {
	foreach pname {detect-thresh detect-minarea deblend-nthresh deblend-mincont \
		       phot-aperture mag-zeropoint gain pixel-scale seeing-fwhm \
		       back-size back-filtersize \
		       phot-aperture-2 phot-aperture-3 phot-aperture-5 conv-filter \
		       n-workers} {
	    set catpanel(param,$pname) $ed($pname)
	}
	CatalogPanelParamSave
    }

    unset ed
}

# --- Mark All Sources ---

# Build region string and create sextract_all markers from catpanel(alldata).
# This is the single source of truth for marker creation.
# Called by: CatalogPanelMarkAll, CatalogPanelMergeSources, AI merge, GotoSource.
proc CatalogPanelCreateAllMarkers {} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set frame $current(frame)

    # Delete previous sextract_all markers
    catch {$frame marker catalog sextract_all delete}

    # Parse directly from alldata (authoritative data source)
    set lines [split $catpanel(alldata) \n]
    if {[llength $lines] < 2} return

    set headers [split [lindex $lines 0] "\t"]
    set ncols [llength $headers]

    # Find column indices (0-based in tab-split fields)
    set col_x -1
    set col_y -1
    set col_a -1
    set col_b -1
    set col_theta -1
    set col_ir -1
    set col_num -1
    for {set c 0} {$c < $ncols} {incr c} {
	set hdr [string trim [lindex $headers $c]]
	switch -- $hdr {
	    NUMBER      { set col_num $c }
	    X_IMAGE     { set col_x $c }
	    Y_IMAGE     { set col_y $c }
	    A_IMAGE     { set col_a $c }
	    B_IMAGE     { set col_b $c }
	    THETA_IMAGE { set col_theta $c }
	    ISO_RADIUS  { set col_ir $c }
	}
    }
    if {$col_x < 0 || $col_y < 0} return

    # Build region strings in batches to avoid DS9 marker command size limits
    set batch_size 500
    set reg "image\n"
    set count 0
    set batch_count 0
    global sextract_all_reg

    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	set fields [split $line "\t"]

	set x [string trim [lindex $fields $col_x]]
	set y [string trim [lindex $fields $col_y]]
	if {![string is double -strict $x] || ![string is double -strict $y]} continue

	# Get source NUMBER for individual tag
	set src_num [expr {$i}]
	if {$col_num >= 0} {
	    set nv [string trim [lindex $fields $col_num]]
	    if {$nv ne {}} { set src_num $nv }
	}

	# Get ellipse parameters (NaN/Inf safe via catch)
	set iso_radius 5.0
	set a_image 0
	set b_image 0
	set theta 0

	if {$col_ir >= 0} {
	    set val [string trim [lindex $fields $col_ir]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} {
		set iso_radius $v
	    }
	}
	if {$col_a >= 0} {
	    set val [string trim [lindex $fields $col_a]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set a_image $v }
	}
	if {$col_b >= 0} {
	    set val [string trim [lindex $fields $col_b]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set b_image $v }
	}
	if {$col_theta >= 0} {
	    set val [string trim [lindex $fields $col_theta]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0} { set theta $v }
	}

	# ISO_RADIUS as semi-major, scaled by B/A for semi-minor
	set semi_a $iso_radius
	set semi_b $iso_radius
	if {$a_image > 0 && $b_image > 0} {
	    set semi_b [expr {$iso_radius * $b_image / $a_image}]
	}

	append reg "ellipse($x $y ${semi_a}i ${semi_b}i $theta) # color=yellow width=1 tag={sextract_all} tag={sextract_src.$src_num} select=0 edit=0 move=0 rotate=0 delete=1 highlite=1 callback=highlite CatalogPanelMarkerCB {$src_num} callback=unhighlite CatalogPanelMarkerUnCB {$src_num}\n"
	incr count
	incr batch_count

	# Flush batch when limit reached
	if {$batch_count >= $batch_size} {
	    set sextract_all_reg $reg
	    catch {$frame marker catalog command ds9 var sextract_all_reg}
	    set reg "image\n"
	    set batch_count 0
	}
    }

    # Flush remaining markers
    if {$batch_count > 0} {
	set sextract_all_reg $reg
	catch {$frame marker catalog command ds9 var sextract_all_reg}
    }

    if {$count == 0} return

    set catpanel(markall,on) 1
    set catpanel(status) "Marked $count sources (yellow ellipses)"
}

proc CatalogPanelMarkAll {} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set frame $current(frame)

    # If already marked, clear first then re-mark
    catch {$frame marker catalog sextract_all delete}
    CatalogPanelCreateAllMarkers
}

proc CatalogPanelClearMarkers {} {
    global catpanel
    global current

    if {$current(frame) == {}} return

    set frame $current(frame)
    catch {$frame marker catalog sextract_all delete}
    set catpanel(markall,on) 0
    set catpanel(status) "Markers cleared"
}

# --- Marker Callbacks (Feature A) ---

proc CatalogPanelMarkerCB {num_str id} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return

    global $catpanel(tbldb)

    # Find NUMBER column index
    set ncols [$catpanel(tbl) cget -cols]
    set col_num -1
    for {set c 1} {$c <= $ncols} {incr c} {
	if {[info exists ${catpanel(tbldb)}(0,$c)]} {
	    set hdr [set ${catpanel(tbldb)}(0,$c)]
	    if {$hdr eq "NUMBER"} {
		set col_num $c
		break
	    }
	}
    }

    # Find table row matching this source NUMBER
    set nrows [$catpanel(tbl) cget -rows]
    set target_row -1

    if {$col_num >= 0} {
	for {set r 1} {$r < $nrows} {incr r} {
	    if {[info exists ${catpanel(tbldb)}($r,$col_num)]} {
		set val [set ${catpanel(tbldb)}($r,$col_num)]
		if {$val eq $num_str} {
		    set target_row $r
		    break
		}
	    }
	}
    }

    if {$target_row < 0} return

    # Select and scroll to row in table
    $catpanel(tbl) selection set $target_row,1
    $catpanel(tbl) see $target_row,1

    # Show selection marker and pan
    CatalogPanelGotoSource $target_row
}

proc CatalogPanelMarkerUnCB {num_str id} {
    # no-op
}

# Click handler called from Button1Frame in none mode
# Find the smallest ellipse source at canvas coordinate (cx, cy).
# Converts canvas→image coords via the frame, then tests all ellipses.
# When ellipses overlap, returns the source NUMBER with the smallest area.
# Returns "" if no ellipse contains the point.
proc CatalogPanelSmallestEllipseAt {frame cx cy} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} { return {} }

    # Convert canvas coords to 1-indexed image coords
    set imgcoord [$frame get coordinates $cx $cy image]
    set imgx [lindex $imgcoord 0]
    set imgy [lindex $imgcoord 1]

    set lines [split $catpanel(alldata) \n]
    if {[llength $lines] < 2} { return {} }

    set headers [split [lindex $lines 0] "\t"]
    set col_x -1; set col_y -1; set col_a -1; set col_b -1
    set col_theta -1; set col_ir -1; set col_num -1
    for {set c 0} {$c < [llength $headers]} {incr c} {
	switch -- [string trim [lindex $headers $c]] {
	    NUMBER      { set col_num $c }
	    X_IMAGE     { set col_x $c }
	    Y_IMAGE     { set col_y $c }
	    A_IMAGE     { set col_a $c }
	    B_IMAGE     { set col_b $c }
	    THETA_IMAGE { set col_theta $c }
	    ISO_RADIUS  { set col_ir $c }
	}
    }
    if {$col_x < 0 || $col_y < 0} { return {} }

    set best_num {}
    set best_area 1e30
    set pi 3.141592653589793

    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	set fields [split $line "\t"]

	# X_IMAGE, Y_IMAGE are 1-indexed image coords
	set sx [string trim [lindex $fields $col_x]]
	set sy [string trim [lindex $fields $col_y]]
	if {![string is double -strict $sx] || ![string is double -strict $sy]} continue

	set src_num $i
	if {$col_num >= 0} {
	    set nv [string trim [lindex $fields $col_num]]
	    if {$nv ne {}} { set src_num $nv }
	}

	# Reconstruct marker ellipse (same logic as CreateAllMarkers)
	set iso_radius 5.0
	set a_image 0; set b_image 0; set theta_deg 0
	if {$col_ir >= 0} {
	    set val [string trim [lindex $fields $col_ir]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set iso_radius $v }
	}
	if {$col_a >= 0} {
	    set val [string trim [lindex $fields $col_a]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set a_image $v }
	}
	if {$col_b >= 0} {
	    set val [string trim [lindex $fields $col_b]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set b_image $v }
	}
	if {$col_theta >= 0} {
	    set val [string trim [lindex $fields $col_theta]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0} { set theta_deg $v }
	}

	set semi_a $iso_radius
	set semi_b $iso_radius
	if {$a_image > 0 && $b_image > 0} {
	    set semi_b [expr {$iso_radius * $b_image / $a_image}]
	}

	# Point-in-ellipse test: rotate (dx,dy) into ellipse frame
	set dx [expr {$imgx - $sx}]
	set dy [expr {$imgy - $sy}]
	set theta_rad [expr {$theta_deg * $pi / 180.0}]
	set cosT [expr {cos($theta_rad)}]
	set sinT [expr {sin($theta_rad)}]
	set rx [expr { $cosT * $dx + $sinT * $dy}]
	set ry [expr {-$sinT * $dx + $cosT * $dy}]

	if {$semi_a <= 0 || $semi_b <= 0} continue
	set t [expr {($rx * $rx) / ($semi_a * $semi_a) + ($ry * $ry) / ($semi_b * $semi_b)}]

	if {$t <= 1.0} {
	    set area [expr {$semi_a * $semi_b}]
	    if {$area < $best_area} {
		set best_area $area
		set best_num $src_num
	    }
	}
    }

    return $best_num
}

proc CatalogPanelMarkerClick {which x y} {
    global catpanel

    if {![info exists catpanel(tbl)]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return
    if {![$which has fits]} return

    # Quick test: is there any marker at this canvas position?
    set id [$which get marker catalog id $x $y]
    if {$id == 0} return

    # Among all overlapping ellipses, pick the smallest one
    set src_num [CatalogPanelSmallestEllipseAt $which $x $y]
    if {$src_num eq {}} {
	# Fallback: use the marker DS9 picked (original behaviour)
	set tags [$which get marker catalog $id tag]
	foreach tag $tags {
	    if {[string match "sextract_src.*" $tag]} {
		set src_num [string range $tag 13 end]
		break
	    }
	}
	if {$src_num eq {}} return
    }

    CatalogPanelMarkerCB $src_num $id
}

# Ctrl+Click handler called from ControlButton1Frame in none mode
proc CatalogPanelMarkerCtrlClick {which x y} {
    global catpanel

    if {![$which has fits]} return

    # Quick test: is there any marker at this canvas position?
    set id [$which get marker catalog id $x $y]
    if {$id == 0} return

    # Check if this is a psf_star marker — if so, remove it
    set tags [$which get marker catalog $id tag]
    foreach tag $tags {
	if {[string match "psf_star.*" $tag]} {
	    set star_num [string range $tag 9 end]
	    # Delete the marker
	    catch {$which marker catalog tag $tag delete}
	    # Remove from star_indices list
	    if {[info exists catpanel(psf,star_indices)]} {
		set idx [lsearch -exact $catpanel(psf,star_indices) $star_num]
		if {$idx >= 0} {
		    set catpanel(psf,star_indices) [lreplace $catpanel(psf,star_indices) $idx $idx]
		}
		set catpanel(status) "Removed star $star_num ([llength $catpanel(psf,star_indices)] stars remaining)"
	    }
	    return
	}
    }

    # Not a star marker — proceed with source merge selection
    if {![info exists catpanel(tbl)]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    # Among all overlapping ellipses, pick the smallest one
    set src_num [CatalogPanelSmallestEllipseAt $which $x $y]
    if {$src_num eq {}} {
	# Fallback: use the marker DS9 picked (original behaviour)
	foreach tag $tags {
	    if {[string match "sextract_src.*" $tag]} {
		set src_num [string range $tag 13 end]
		break
	    }
	}
	if {$src_num eq {}} return
    }

    CatalogPanelCtrlSelect $src_num
}

# --- Visible Filter (Feature B) ---

proc CatalogPanelShowVisible {} {
    global catpanel
    global current
    global ds9

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return
    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return

    # checkbutton already toggled catpanel(visible_mode) before calling us
    if {!$catpanel(visible_mode)} {
	CatalogPanelLoadTSV $catpanel(alldata) "all"
	set catpanel(status) "Showing all sources"
	return
    }

    set frame $current(frame)

    # Get viewport center in image coordinates
    set cursor [$frame get cursor image]
    set cx [lindex $cursor 0]
    set cy [lindex $cursor 1]

    # Get zoom level
    set zoom [$frame get zoom]
    set zx [lindex $zoom 0]
    set zy [lindex $zoom 1]

    # Get canvas size
    set cw [winfo width $ds9(canvas)]
    set ch [winfo height $ds9(canvas)]

    # Compute viewport bounds in image coordinates
    set x_min [expr {$cx - $cw / 2.0 / $zx}]
    set x_max [expr {$cx + $cw / 2.0 / $zx}]
    set y_min [expr {$cy - $ch / 2.0 / $zy}]
    set y_max [expr {$cy + $ch / 2.0 / $zy}]

    # Parse alldata, find X_IMAGE/Y_IMAGE columns
    set lines [split $catpanel(alldata) \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]
    set ncols [llength $headers]

    set idx_x -1
    set idx_y -1
    for {set i 0} {$i < $ncols} {incr i} {
	set h [string trim [lindex $headers $i]]
	if {$h eq "X_IMAGE"} { set idx_x $i }
	if {$h eq "Y_IMAGE"} { set idx_y $i }
    }
    if {$idx_x < 0 || $idx_y < 0} return

    # Filter rows within viewport
    set filtered $header
    set count 0
    set total 0
    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	incr total
	set fields [split $line "\t"]
	set x [string trim [lindex $fields $idx_x]]
	set y [string trim [lindex $fields $idx_y]]
	if {![string is double -strict $x] || ![string is double -strict $y]} continue
	if {$x >= $x_min && $x <= $x_max && $y >= $y_min && $y <= $y_max} {
	    append filtered "\n$line"
	    incr count
	}
    }

    CatalogPanelLoadTSV $filtered "visible"
    set catpanel(status) "Visible: $count of $total sources in current view"
}

# --- Merge Selection (Feature C) ---

proc CatalogPanelCtrlSelect {src_num} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return

    set frame $current(frame)

    # Toggle: if already in list, remove; otherwise add
    set idx [lsearch -exact $catpanel(merge,list) $src_num]
    if {$idx >= 0} {
	# Remove from merge list
	set catpanel(merge,list) [lreplace $catpanel(merge,list) $idx $idx]
	# Delete this source's merge marker
	catch {$frame marker catalog sextract_merge.$src_num delete}
    } else {
	# Add to merge list
	lappend catpanel(merge,list) $src_num

	# Find source position from alldata
	set lines [split $catpanel(alldata) \n]
	set header [lindex $lines 0]
	set headers [split $header "\t"]
	set ncols [llength $headers]

	set idx_num -1
	set idx_x -1
	set idx_y -1
	set idx_a -1
	set idx_b -1
	set idx_theta -1
	set idx_ir -1
	for {set i 0} {$i < $ncols} {incr i} {
	    set h [string trim [lindex $headers $i]]
	    switch -- $h {
		NUMBER      { set idx_num $i }
		X_IMAGE     { set idx_x $i }
		Y_IMAGE     { set idx_y $i }
		A_IMAGE     { set idx_a $i }
		B_IMAGE     { set idx_b $i }
		THETA_IMAGE { set idx_theta $i }
		ISO_RADIUS  { set idx_ir $i }
	    }
	}

	# Find the matching line
	for {set i 1} {$i < [llength $lines]} {incr i} {
	    set line [lindex $lines $i]
	    if {[string trim $line] eq {}} continue
	    set fields [split $line "\t"]
	    set num_val [string trim [lindex $fields $idx_num]]
	    if {$num_val eq $src_num} {
		set x [string trim [lindex $fields $idx_x]]
		set y [string trim [lindex $fields $idx_y]]

		# Get ellipse params
		set iso_radius 5.0
		set a_image 0
		set b_image 0
		set theta 0
		if {$idx_ir >= 0} {
		    set val [string trim [lindex $fields $idx_ir]]
		    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} {
			set iso_radius $v
		    }
		}
		if {$idx_a >= 0} {
		    set val [string trim [lindex $fields $idx_a]]
		    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set a_image $v }
		}
		if {$idx_b >= 0} {
		    set val [string trim [lindex $fields $idx_b]]
		    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set b_image $v }
		}
		if {$idx_theta >= 0} {
		    set val [string trim [lindex $fields $idx_theta]]
		    if {[catch {set v [expr {$val + 0.0}]}] == 0} { set theta $v }
		}

		set semi_a $iso_radius
		set semi_b $iso_radius
		if {$a_image > 0 && $b_image > 0} {
		    set semi_b [expr {$iso_radius * $b_image / $a_image}]
		}

		# Create red thick merge marker
		global sextract_merge_reg
		set sextract_merge_reg "image\nellipse($x $y ${semi_a}i ${semi_b}i $theta) # color=red width=3 tag={sextract_merge} tag={sextract_merge.$src_num} select=0 edit=0 move=0 rotate=0 delete=1\n"
		catch {$frame marker catalog command ds9 var sextract_merge_reg}
		break
	    }
	}
    }

    set catpanel(merge,active) 1
    set n [llength $catpanel(merge,list)]
    if {$n == 0} {
	set catpanel(merge,active) 0
	set catpanel(status) "Merge selection cleared"
    } else {
	set catpanel(status) "Merge: $n sources selected (Ctrl+M to merge, Esc to cancel)"
    }
}

proc CatalogPanelMergeSources {} {
    global catpanel
    global current

    if {!$catpanel(merge,active)} return
    if {[llength $catpanel(merge,list)] < 2} {
	set catpanel(status) "Need at least 2 sources to merge"
	return
    }

    # Parse alldata
    set lines [split $catpanel(alldata) \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]
    set ncols [llength $headers]

    # Find column indices
    set idx_num -1
    set idx_x -1
    set idx_y -1
    set idx_a -1
    set idx_b -1
    set idx_theta -1
    set idx_ir -1
    set idx_flux -1
    set idx_mag -1
    set idx_npix -1
    for {set i 0} {$i < $ncols} {incr i} {
	set h [string trim [lindex $headers $i]]
	switch -- $h {
	    NUMBER      { set idx_num $i }
	    X_IMAGE     { set idx_x $i }
	    Y_IMAGE     { set idx_y $i }
	    A_IMAGE     { set idx_a $i }
	    B_IMAGE     { set idx_b $i }
	    THETA_IMAGE { set idx_theta $i }
	    ISO_RADIUS  { set idx_ir $i }
	    FLUX_AUTO   { set idx_flux $i }
	    MAG_AUTO    { set idx_mag $i }
	    NPIX_ISO    { set idx_npix $i }
	}
    }

    # Collect data for merge sources and find brightest
    set merge_rows {}
    set other_rows {}
    set max_number 0
    set brightest_idx -1
    set brightest_flux -1e30

    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	set fields [split $line "\t"]
	set num_val [string trim [lindex $fields $idx_num]]

	# Track max NUMBER
	if {[string is integer -strict $num_val] && $num_val > $max_number} {
	    set max_number $num_val
	}

	if {[lsearch -exact $catpanel(merge,list) $num_val] >= 0} {
	    lappend merge_rows $fields
	    if {$idx_flux >= 0} {
		set fval [string trim [lindex $fields $idx_flux]]
		if {[string is double -strict $fval] && $fval > $brightest_flux} {
		    set brightest_flux $fval
		    set brightest_idx [expr {[llength $merge_rows] - 1}]
		}
	    }
	} else {
	    lappend other_rows $line
	}
    }

    if {[llength $merge_rows] < 2} {
	set catpanel(status) "Merge error: sources not found in catalog"
	return
    }

    if {$brightest_idx < 0} { set brightest_idx 0 }
    set new_num [expr {$max_number + 1}]

    # Compute merged values
    # Pass 1: flux-weighted centroid for X,Y + total flux/npix
    set total_flux 0.0
    set wx 0.0
    set wy 0.0
    set total_npix 0

    foreach row $merge_rows {
	set flux 1.0
	if {$idx_flux >= 0} {
	    set fv [string trim [lindex $row $idx_flux]]
	    if {[string is double -strict $fv] && $fv > 0} { set flux $fv }
	}
	set x [string trim [lindex $row $idx_x]]
	set y [string trim [lindex $row $idx_y]]
	if {![string is double -strict $x]} { set x 0 }
	if {![string is double -strict $y]} { set y 0 }

	set total_flux [expr {$total_flux + $flux}]
	set wx [expr {$wx + $x * $flux}]
	set wy [expr {$wy + $y * $flux}]

	if {$idx_npix >= 0} {
	    set nv [string trim [lindex $row $idx_npix]]
	    if {[string is integer -strict $nv]} {
		set total_npix [expr {$total_npix + $nv}]
	    }
	}
    }

    if {$total_flux <= 0} { set total_flux 1.0 }

    set new_x [expr {$wx / $total_flux}]
    set new_y [expr {$wy / $total_flux}]

    # Pass 2: second-moment tensor for merged ellipse (A, B, THETA)
    set Ixx 0.0
    set Iyy 0.0
    set Ixy 0.0

    foreach row $merge_rows {
	set flux 1.0
	if {$idx_flux >= 0} {
	    set fv [string trim [lindex $row $idx_flux]]
	    if {[string is double -strict $fv] && $fv > 0} { set flux $fv }
	}
	set x [string trim [lindex $row $idx_x]]
	set y [string trim [lindex $row $idx_y]]
	if {![string is double -strict $x]} { set x 0 }
	if {![string is double -strict $y]} { set y 0 }

	set ak 0.0
	set bk 0.0
	set thetak 0.0
	if {$idx_a >= 0} {
	    set av [string trim [lindex $row $idx_a]]
	    if {[string is double -strict $av]} { set ak $av }
	}
	if {$idx_b >= 0} {
	    set bv [string trim [lindex $row $idx_b]]
	    if {[string is double -strict $bv]} { set bk $bv }
	}
	if {$idx_theta >= 0} {
	    set tv [string trim [lindex $row $idx_theta]]
	    if {[string is double -strict $tv]} { set thetak $tv }
	}

	# Intrinsic second moments of this source's ellipse
	set rad [expr {$thetak * 3.14159265358979 / 180.0}]
	set cosT [expr {cos($rad)}]
	set sinT [expr {sin($rad)}]
	set a2 [expr {$ak * $ak}]
	set b2 [expr {$bk * $bk}]
	set ixx_k [expr {$a2 * $cosT * $cosT + $b2 * $sinT * $sinT}]
	set iyy_k [expr {$a2 * $sinT * $sinT + $b2 * $cosT * $cosT}]
	set ixy_k [expr {($a2 - $b2) * $sinT * $cosT}]

	# Parallel axis theorem: add offset from merged centroid
	set dx [expr {$x - $new_x}]
	set dy [expr {$y - $new_y}]
	set Ixx [expr {$Ixx + $flux * ($ixx_k + $dx * $dx)}]
	set Iyy [expr {$Iyy + $flux * ($iyy_k + $dy * $dy)}]
	set Ixy [expr {$Ixy + $flux * ($ixy_k + $dx * $dy)}]
    }

    # Normalize by total flux
    set Ixx [expr {$Ixx / $total_flux}]
    set Iyy [expr {$Iyy / $total_flux}]
    set Ixy [expr {$Ixy / $total_flux}]

    # Eigenvalue decomposition → A, B, THETA
    set trace [expr {$Ixx + $Iyy}]
    set det [expr {$Ixx * $Iyy - $Ixy * $Ixy}]
    set disc [expr {sqrt(abs(($Ixx - $Iyy) * ($Ixx - $Iyy) + 4.0 * $Ixy * $Ixy))}]
    set lam1 [expr {($trace + $disc) / 2.0}]
    set lam2 [expr {($trace - $disc) / 2.0}]
    if {$lam1 < 0} { set lam1 0.0 }
    if {$lam2 < 0} { set lam2 0.0 }
    set new_a [expr {sqrt($lam1)}]
    set new_b [expr {sqrt($lam2)}]
    set new_theta [expr {0.5 * atan2(2.0 * $Ixy, $Ixx - $Iyy) * 180.0 / 3.14159265358979}]

    # MAG_AUTO from total flux
    set mag_zp 25.0
    if {[info exists catpanel(param,mag-zeropoint)]} {
	set mag_zp $catpanel(param,mag-zeropoint)
    }
    set new_mag [expr {-2.5 * log10($total_flux) + $mag_zp}]

    # ISO_RADIUS from total NPIX
    set new_ir 5.0
    if {$total_npix > 0 && $new_a > 0 && $new_b > 0} {
	set ratio [expr {$new_b / $new_a}]
	if {$ratio <= 0} { set ratio 1.0 }
	set new_ir [expr {sqrt($total_npix / (3.14159265 * $ratio))}]
    }

    # Build merged row: copy from brightest, override computed fields
    set base_row [lindex $merge_rows $brightest_idx]
    set new_fields {}
    for {set c 0} {$c < $ncols} {incr c} {
	set val [string trim [lindex $base_row $c]]
	if {$c == $idx_num} { set val $new_num }
	if {$c == $idx_x} { set val [format "%.4f" $new_x] }
	if {$c == $idx_y} { set val [format "%.4f" $new_y] }
	if {$c == $idx_flux && $idx_flux >= 0} { set val [format "%.6g" $total_flux] }
	if {$c == $idx_mag && $idx_mag >= 0} { set val [format "%.4f" $new_mag] }
	if {$c == $idx_npix && $idx_npix >= 0} { set val $total_npix }
	if {$c == $idx_ir && $idx_ir >= 0} { set val [format "%.4f" $new_ir] }
	if {$c == $idx_a && $idx_a >= 0} { set val [format "%.4f" $new_a] }
	if {$c == $idx_b && $idx_b >= 0} { set val [format "%.4f" $new_b] }
	if {$c == $idx_theta && $idx_theta >= 0} { set val [format "%.4f" $new_theta] }
	lappend new_fields $val
    }
    set new_line [join $new_fields "\t"]

    # Rebuild alldata: header + other rows + merged row
    set newdata $header
    foreach row $other_rows {
	append newdata "\n$row"
    }
    append newdata "\n$new_line"
    set catpanel(alldata) $newdata

    # Clear merge state
    set nmerged [llength $catpanel(merge,list)]
    set catpanel(merge,list) {}
    set catpanel(merge,active) 0

    # Delete merge markers
    if {$current(frame) != {}} {
	catch {$current(frame) marker catalog sextract_merge delete}
    }

    # Reload table and markers
    CatalogPanelLoadTSV $catpanel(alldata) "merged"
    CatalogPanelCreateAllMarkers

    # Find merged source row and auto-select/navigate
    global $catpanel(tbldb)
    set ncols [$catpanel(tbl) cget -cols]
    set nrows [$catpanel(tbl) cget -rows]
    set col_num -1
    for {set c 1} {$c <= $ncols} {incr c} {
	if {[info exists ${catpanel(tbldb)}(0,$c)]} {
	    set hdr [set ${catpanel(tbldb)}(0,$c)]
	    if {$hdr eq "NUMBER"} {
		set col_num $c
		break
	    }
	}
    }
    set merged_row -1
    if {$col_num >= 0} {
	for {set r 1} {$r < $nrows} {incr r} {
	    if {[info exists ${catpanel(tbldb)}($r,$col_num)]} {
		set val [set ${catpanel(tbldb)}($r,$col_num)]
		if {$val eq $new_num} {
		    set merged_row $r
		    break
		}
	    }
	}
    }
    if {$merged_row >= 0} {
	$catpanel(tbl) selection set $merged_row,1
	$catpanel(tbl) see $merged_row,1
	CatalogPanelGotoSource $merged_row
    }

    set catpanel(status) "Merged $nmerged sources into #$new_num: pos=([format %.2f $new_x],[format %.2f $new_y]) mag=[format %.3f $new_mag]"
}

proc CatalogPanelMergeCancel {} {
    global catpanel
    global current

    # Delete all merge markers
    if {$current(frame) != {}} {
	catch {$current(frame) marker catalog sextract_merge delete}
    }

    set catpanel(merge,list) {}
    set catpanel(merge,active) 0
    set catpanel(status) "Merge cancelled"
}

proc CatalogPanelEscapeKey {} {
    global catpanel

    if {$catpanel(ai,active)} {
	CatalogPanelAIDone
	return
    }
    if {$catpanel(merge,active)} {
	CatalogPanelMergeCancel
    }
}

# --- Log Scale ---

proc CatalogPanelSetLogScale {} {
    global current
    global scale

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return

    # Set scale to log
    set scale(type) log
    $current(frame) colorscale log $scale(log)
    $current(frame) colorscale log

    # Get min/max pixel values
    $current(frame) clip mode minmax
    set limits [$current(frame) get clip]
    set pmin [lindex $limits 0]
    set pmax [lindex $limits 1]

    # Guard against NaN/Inf/non-numeric
    if {[catch {expr {$pmin + 0.0}}]} { set pmin 0.001 }
    if {[catch {expr {$pmax + 0.0}}]} { set pmax 1000.0 }

    # Ensure positive values for log
    if {$pmin <= 0} { set pmin 0.001 }
    if {$pmax <= $pmin} { set pmax [expr {$pmin * 1000}] }

    # Compute display max: 0.8*(log(max)-log(min)) + log(min)
    set log_min [expr {log10($pmin)}]
    set log_max [expr {log10($pmax)}]
    set log_disp [expr {0.8 * ($log_max - $log_min) + $log_min}]
    set disp_max [expr {pow(10.0, $log_disp)}]

    # Apply user-defined limits
    set scale(min) $pmin
    set scale(max) $disp_max
    set scale(mode) user
    $current(frame) clip user $pmin $disp_max
    $current(frame) clip mode user
    UpdateScale
}

# --- AI Merge ---

proc CatalogPanelAIMerge {} {
    global catpanel
    global current
    global ds9

    # Need extracted catalog first
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first before AI Merge"
	return
    }

    # Get current FITS filename (same pattern as CatalogPanelExtract)
    set fn {}
    if {$current(frame) != {}} {
	catch {set fn [$current(frame) get fits file name full]}
    }
    set fn [string trim $fn "{}"]
    regsub {\[.*\]$} $fn {} fn
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    # Find ds9_ai_merge.py
    set bindir [file dirname [info nameofexecutable]]
    set script [file join $bindir ds9_ai_merge.py]
    if {![file exists $script]} {
	# Also check library dir
	set libdir [file join [file dirname $bindir] ds9 library]
	set script [file join $libdir ds9_ai_merge.py]
    }
    if {![file exists $script]} {
	set catpanel(status) "ERROR: ds9_ai_merge.py not found"
	return
    }

    # Save catalog to temp TSV for --catalog mode
    set catfile [file join [file normalize ~] .ds9 ai_merge_catalog.tsv]
    catch {file mkdir [file dirname $catfile]}
    if {[catch {
	set fd [open $catfile w]
	puts $fd $catpanel(alldata)
	close $fd
    } err]} {
	set catpanel(status) "AI Merge error: cannot write catalog: $err"
	return
    }

    # Build parameter arguments — use params from Extract time, not current settings
    set paramargs {}
    lappend paramargs "--threshold" $catpanel(ai,threshold)
    lappend paramargs "--catalog" $catfile
    foreach pname {detect-thresh detect-minarea deblend-nthresh deblend-mincont \
		   mag-zeropoint back-size back-filtersize} {
	if {[info exists catpanel(extract_param,$pname)]} {
	    lappend paramargs "--$pname" $catpanel(extract_param,$pname)
	} elseif {[info exists catpanel(param,$pname)]} {
	    lappend paramargs "--$pname" $catpanel(param,$pname)
	}
    }

    # Find checkpoint
    set ckpt [file join [file dirname $bindir] ai_merge data checkpoints mlp_best.pt]
    if {[file exists $ckpt]} {
	lappend paramargs "--checkpoint" $ckpt
    }

    set catpanel(status) "AI Merge: running prediction on [file tail $fn] ..."
    update idletasks

    # Run prediction — capture stderr for error diagnostics
    set errfile [file join [file normalize ~] .ds9 ai_merge_stderr.txt]
    if {[catch {set data [exec python3 $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	catch {file delete $catfile}
	if {$stderr_msg ne ""} {
	    # Show last line of stderr (most relevant error)
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    set catpanel(status) "AI Merge error: $last_err"
	    puts "AI Merge full stderr:\n$stderr_msg"
	} else {
	    set catpanel(status) "AI Merge error: $err"
	}
	return
    }
    catch {file delete $errfile}
    catch {file delete $catfile}

    # Parse output
    set lines [split $data \n]
    set groups {}
    set n_groups 0
    set threshold 0.70

    foreach line $lines {
	if {[string match "#AI_MERGE*" $line]} {
	    # Parse header: #AI_MERGE	N_GROUPS=5	THRESHOLD=0.70
	    foreach field [split $line "\t"] {
		if {[string match "N_GROUPS=*" $field]} {
		    set n_groups [string range $field 9 end]
		}
		if {[string match "THRESHOLD=*" $field]} {
		    set threshold [string range $field 10 end]
		}
	    }
	    continue
	}
	if {[string match "GROUP*" $line] && [string match "*MEMBERS_X*" $line]} {
	    # Skip column header
	    continue
	}
	if {[string trim $line] eq {}} continue

	# Data row: GROUP	N_MEMBERS	CONFIDENCE	MEMBERS_X	MEMBERS_Y	MEMBERS_NUM
	set fields [split $line "\t"]
	if {[llength $fields] < 6} continue
	set g_idx [lindex $fields 0]
	set n_mem [lindex $fields 1]
	set conf  [lindex $fields 2]
	set mem_x [lindex $fields 3]
	set mem_y [lindex $fields 4]
	set mem_num [lindex $fields 5]
	lappend groups [list $g_idx $n_mem $conf $mem_x $mem_y $mem_num]
    }

    if {[llength $groups] == 0} {
	set catpanel(status) "AI Merge: no merge groups found (threshold=$threshold)"
	return
    }

    # Store state
    set catpanel(ai,groups) $groups
    set catpanel(ai,total) [llength $groups]
    set catpanel(ai,current) 0
    set catpanel(ai,active) 1

    # Bind navigation keys
    CatalogPanelAIBindKeys

    # Show first group
    CatalogPanelAIShowGroup 0
}

proc CatalogPanelAIBindKeys {} {
    global ds9

    # Bind on toplevel AND canvas (canvas has focus when mouse is over image)
    foreach w [list . $ds9(canvas)] {
	bind $w <Key-n> {CatalogPanelAINext}
	bind $w <Key-p> {CatalogPanelAIPrev}
	bind $w <Key-a> {CatalogPanelAIAccept}
	bind $w <Key-r> {CatalogPanelAIReject}
	bind $w <Right> {CatalogPanelAINext}
	bind $w <Left>  {CatalogPanelAIPrev}
    }

    # Force focus to canvas so keys work immediately
    focus -force $ds9(canvas)
}

proc CatalogPanelAIUnbindKeys {} {
    global ds9

    foreach w [list . $ds9(canvas)] {
	bind $w <Key-n> {}
	bind $w <Key-p> {}
	bind $w <Key-a> {}
	bind $w <Key-r> {}
	bind $w <Right> {}
	bind $w <Left>  {}
    }
}

proc CatalogPanelAIShowGroup {idx} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    set frame $current(frame)

    # Delete previous ai_merge markers
    catch {$frame marker catalog ai_merge delete}

    set groups $catpanel(ai,groups)
    if {$idx < 0 || $idx >= [llength $groups]} return

    set catpanel(ai,current) $idx
    set group [lindex $groups $idx]

    # Parse group: {g_idx n_mem conf mem_x mem_y mem_num}
    set n_mem [lindex $group 1]
    set conf  [lindex $group 2]
    set xs_str [lindex $group 3]
    set ys_str [lindex $group 4]
    set nums_str [lindex $group 5]
    set xs [split $xs_str ","]
    set ys [split $ys_str ","]
    set nums [split $nums_str ","]

    if {[llength $xs] < 2 || [llength $xs] != [llength $ys]} return

    # Use MEMBERS_NUM directly — no coordinate matching needed
    set matched_nums $nums
    set matched_xs {}
    set matched_ys {}
    for {set m 0} {$m < [llength $xs]} {incr m} {
	set ax [lindex $xs $m]
	set ay [lindex $ys $m]
	if {[catch {expr {$ax + 0.0}}] || [catch {expr {$ay + 0.0}}]} continue
	lappend matched_xs $ax
	lappend matched_ys $ay
    }

    # Draw markers for each member
    set color magenta
    set reg "image\n"
    for {set m 0} {$m < [llength $matched_xs]} {incr m} {
	set mx [lindex $matched_xs $m]
	set my [lindex $matched_ys $m]
	# Draw circle marker for each member
	append reg "circle($mx,$my,8) # color=$color width=3 tag={ai_merge}\n"
    }

    # Draw connecting lines between all pairs
    for {set a 0} {$a < [llength $matched_xs]} {incr a} {
	for {set b [expr {$a + 1}]} {$b < [llength $matched_xs]} {incr b} {
	    set x1 [lindex $matched_xs $a]
	    set y1 [lindex $matched_ys $a]
	    set x2 [lindex $matched_xs $b]
	    set y2 [lindex $matched_ys $b]
	    append reg "line($x1,$y1,$x2,$y2) # color=$color width=1 dash=1 tag={ai_merge}\n"
	}
    }

    # Create markers
    global ai_merge_reg
    set ai_merge_reg $reg
    $frame marker catalog command ds9 var ai_merge_reg

    # Pan to group center
    if {[llength $matched_xs] < 2} return
    set cx 0.0
    set cy 0.0
    foreach mx $matched_xs my $matched_ys {
	catch {set cx [expr {$cx + $mx}]}
	catch {set cy [expr {$cy + $my}]}
    }
    set nm [llength $matched_xs]
    if {$nm > 0} {
	set cx [expr {$cx / $nm}]
	set cy [expr {$cy / $nm}]
	PanTo $cx $cy image {}
    }

    # Status bar
    set g_num [expr {$idx + 1}]
    set total $catpanel(ai,total)
    set num_str [join $matched_nums ","]
    set catpanel(status) "AI Group $g_num/$total (conf=[format %.2f $conf], ${n_mem} sources: $num_str) \[n:Next p:Prev a:Accept r:Reject Esc:Done\]"

    # Ensure canvas has focus so keys work
    global ds9
    focus -force $ds9(canvas)
}

proc CatalogPanelAINext {} {
    global catpanel
    if {!$catpanel(ai,active)} return
    set next [expr {$catpanel(ai,current) + 1}]
    if {$next >= $catpanel(ai,total)} {
	set catpanel(status) "AI Merge: last group reached. Press Esc to finish."
	return
    }
    CatalogPanelAIShowGroup $next
}

proc CatalogPanelAIPrev {} {
    global catpanel
    if {!$catpanel(ai,active)} return
    set prev [expr {$catpanel(ai,current) - 1}]
    if {$prev < 0} {
	set catpanel(status) "AI Merge: already at first group."
	return
    }
    CatalogPanelAIShowGroup $prev
}

proc CatalogPanelAIAccept {} {
    global catpanel
    global current
    if {!$catpanel(ai,active)} return

    set idx $catpanel(ai,current)
    set groups $catpanel(ai,groups)
    set group [lindex $groups $idx]

    # Get member NUMBERs directly from MEMBERS_NUM field
    set nums_str [lindex $group 5]
    set merge_nums [split $nums_str ","]

    if {[llength $merge_nums] < 2} {
	set catpanel(status) "AI Accept: could not match enough sources"
	return
    }

    # Delete AI markers before merge (merge will re-mark)
    catch {$current(frame) marker catalog ai_merge delete}

    # Set up merge and execute
    set catpanel(merge,list) $merge_nums
    set catpanel(merge,active) 1
    CatalogPanelMergeSources

    # Remove accepted group from list
    set catpanel(ai,groups) [lreplace $groups $idx $idx]
    set catpanel(ai,total) [llength $catpanel(ai,groups)]

    # Advance to next (or stay at end)
    if {$catpanel(ai,total) == 0} {
	CatalogPanelAIDone
	return
    }
    if {$idx >= $catpanel(ai,total)} {
	set idx [expr {$catpanel(ai,total) - 1}]
    }
    CatalogPanelAIShowGroup $idx
}

proc CatalogPanelAIReject {} {
    global catpanel
    if {!$catpanel(ai,active)} return

    set idx $catpanel(ai,current)
    # Remove rejected group from list
    set catpanel(ai,groups) [lreplace $catpanel(ai,groups) $idx $idx]
    set catpanel(ai,total) [llength $catpanel(ai,groups)]

    if {$catpanel(ai,total) == 0} {
	CatalogPanelAIDone
	return
    }
    if {$idx >= $catpanel(ai,total)} {
	set idx [expr {$catpanel(ai,total) - 1}]
    }
    CatalogPanelAIShowGroup $idx
}

proc CatalogPanelAIDone {} {
    global catpanel
    global current

    # Delete AI markers
    if {$current(frame) != {}} {
	catch {$current(frame) marker catalog ai_merge delete}
    }

    # Unbind navigation keys
    CatalogPanelAIUnbindKeys

    # Reset state
    set catpanel(ai,groups) {}
    set catpanel(ai,active) 0
    set catpanel(ai,total) 0
    set catpanel(ai,current) 0

    # Re-mark all sources from authoritative data
    CatalogPanelCreateAllMarkers

    set catpanel(status) "AI Merge session ended"
}

# --- Column Header Click Sorting ---

proc CatalogPanelTableClick {x y} {
    global catpanel

    set tbl $catpanel(tbl)
    set idx [$tbl index @$x,$y]
    set row [lindex [split $idx ,] 0]

    # Only handle header row clicks
    if {$row != 0} return

    set col [lindex [split $idx ,] 1]

    global $catpanel(tbldb)
    if {![info exists ${catpanel(tbldb)}(0,$col)]} return
    set colname [set ${catpanel(tbldb)}(0,$col)]

    # Toggle direction if same column clicked again
    if {$catpanel(sort,col) eq $colname} {
	if {$catpanel(sort,dir) eq "ascending"} {
	    set catpanel(sort,dir) descending
	} else {
	    set catpanel(sort,dir) ascending
	}
    } else {
	set catpanel(sort,col) $colname
	set catpanel(sort,dir) ascending
    }

    CatalogPanelSort $colname $catpanel(sort,dir)
}

proc CatalogPanelSort {colname direction} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set lines [split $catpanel(alldata) \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]

    # Find column index
    set colidx -1
    for {set i 0} {$i < [llength $headers]} {incr i} {
	if {[string trim [lindex $headers $i]] eq $colname} {
	    set colidx $i
	    break
	}
    }
    if {$colidx < 0} return

    # Collect data rows (skip header and empty lines)
    set datarows {}
    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	lappend datarows $line
    }

    # Determine sort type: check first non-empty value
    set isnumeric 1
    foreach drow $datarows {
	set val [string trim [lindex [split $drow "\t"] $colidx]]
	if {$val ne {}} {
	    if {![string is double -strict $val]} {
		set isnumeric 0
	    }
	    break
	}
    }

    # Sort
    if {$isnumeric} {
	set cmd [list CatalogPanelSortCmpNum $colidx]
    } else {
	set cmd [list CatalogPanelSortCmpStr $colidx]
    }
    if {$direction eq "descending"} {
	set sortedrows [lsort -decreasing -command $cmd $datarows]
    } else {
	set sortedrows [lsort -command $cmd $datarows]
    }

    # Rebuild alldata with sorted rows
    set newdata $header
    foreach drow $sortedrows {
	append newdata "\n$drow"
    }
    set catpanel(alldata) $newdata

    # Reload table
    CatalogPanelLoadTSV $catpanel(alldata) "sorted"

    set catpanel(status) "Sorted by $colname $direction"
}

proc CatalogPanelSortCmpNum {colidx a b} {
    set va [string trim [lindex [split $a "\t"] $colidx]]
    set vb [string trim [lindex [split $b "\t"] $colidx]]
    if {![string is double -strict $va]} { set va 0 }
    if {![string is double -strict $vb]} { set vb 0 }
    if {$va < $vb} { return -1 }
    if {$va > $vb} { return 1 }
    return 0
}

proc CatalogPanelSortCmpStr {colidx a b} {
    set va [string trim [lindex [split $a "\t"] $colidx]]
    set vb [string trim [lindex [split $b "\t"] $colidx]]
    return [string compare $va $vb]
}

# --- Trim Filter (Feature D) ---

proc CatalogPanelTrimDialog {} {
    global catpanel
    global ed

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "No catalog data to trim"
	return
    }

    set w {.sextracttrim}

    set ed(ok) 0

    # Get column names from alldata header
    set lines [split $catpanel(alldata) \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]
    set ncols [llength $headers]

    set ed(trim,cols) {}
    for {set i 0} {$i < $ncols} {incr i} {
	set colname [string trim [lindex $headers $i]]
	lappend ed(trim,cols) $colname
	# Initialize from existing trim values or empty
	if {[info exists catpanel(trim,$colname,min)]} {
	    set ed(trim,$colname,min) $catpanel(trim,$colname,min)
	} else {
	    set ed(trim,$colname,min) {}
	}
	if {[info exists catpanel(trim,$colname,max)]} {
	    set ed(trim,$colname,max) $catpanel(trim,$colname,max)
	} else {
	    set ed(trim,$colname,max) {}
	}
    }

    DialogCreate $w {Trim - Column Filter} ed(ok)

    # Scrollable frame for columns
    set sf [ttk::frame $w.param]
    set canvas_w [canvas $sf.c -width 400 -height 300 \
		      -yscrollcommand [list $sf.vs set]]
    ttk::scrollbar $sf.vs -orient vertical -command [list $canvas_w yview]
    set inner [ttk::frame $canvas_w.inner]
    $canvas_w create window 0 0 -anchor nw -window $inner

    # Header labels
    ttk::label $inner.hcol -text "Column" -font {Helvetica 10 bold} -width 15
    ttk::label $inner.hmin -text "Min" -font {Helvetica 10 bold} -width 12
    ttk::label $inner.htilde -text "" -width 2
    ttk::label $inner.hmax -text "Max" -font {Helvetica 10 bold} -width 12
    grid $inner.hcol $inner.hmin $inner.htilde $inner.hmax \
	-padx 2 -pady 2 -sticky w

    set row 1
    foreach colname $ed(trim,cols) {
	ttk::label $inner.l$row -text "$colname:" -anchor w -width 15
	ttk::entry $inner.emin$row -textvariable ed(trim,$colname,min) -width 12
	ttk::label $inner.tilde$row -text "~" -width 2
	ttk::entry $inner.emax$row -textvariable ed(trim,$colname,max) -width 12
	grid $inner.l$row $inner.emin$row $inner.tilde$row $inner.emax$row \
	    -padx 2 -pady 1 -sticky w
	incr row
    }

    # Update scroll region after layout
    bind $inner <Configure> [list $canvas_w configure -scrollregion \
				 [$canvas_w bbox all]]

    pack $canvas_w -side left -fill both -expand true
    pack $sf.vs -side right -fill y

    # Buttons
    set bf [ttk::frame $w.buttons]
    ttk::button $bf.ok -text {Apply} -command {set ed(ok) 1} -default active
    ttk::button $bf.cancel -text {Cancel} -command {set ed(ok) 0}
    ttk::button $bf.reset -text {Reset} -command {
	foreach col $ed(trim,cols) {
	    set ed(trim,$col,min) {}
	    set ed(trim,$col,max) {}
	}
    }
    ttk::button $bf.save -text {Save} -command {
	CatalogPanelTrimSaveFromEd
    }
    ttk::button $bf.load -text {Load} -command {
	CatalogPanelTrimLoadToEd
    }
    pack $bf.ok $bf.cancel $bf.reset $bf.save $bf.load \
	-side left -expand true -padx 2 -pady 4

    bind $w <Return> {set ed(ok) 1}

    # Fini
    ttk::separator $w.sep -orient horizontal
    pack $w.buttons $w.sep -side bottom -fill x
    pack $w.param -side top -fill both -expand true

    DialogWait $w ed(ok)
    destroy $w

    if {$ed(ok)} {
	# Copy trim values from ed to catpanel
	foreach colname $ed(trim,cols) {
	    set catpanel(trim,$colname,min) $ed(trim,$colname,min)
	    set catpanel(trim,$colname,max) $ed(trim,$colname,max)
	}
	CatalogPanelTrimApply
    }

    unset ed
}

proc CatalogPanelTrimSaveFromEd {} {
    global ed

    set prefdir [file join [file normalize ~] .ds9]
    if {![file isdirectory $prefdir]} {
	file mkdir $prefdir
    }
    set preffile [file join $prefdir sextract_trim.prf]
    if {[catch {set fd [open $preffile w]} err]} return

    foreach colname $ed(trim,cols) {
	puts $fd "$colname\t$ed(trim,$colname,min)\t$ed(trim,$colname,max)"
    }
    close $fd
}

proc CatalogPanelTrimLoadToEd {} {
    global ed

    set preffile [file join [file normalize ~] .ds9 sextract_trim.prf]
    if {![file exists $preffile]} return
    if {[catch {set fd [open $preffile r]} err]} return

    while {[gets $fd line] >= 0} {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set parts [split $line "\t"]
	if {[llength $parts] >= 3} {
	    set colname [lindex $parts 0]
	    set minval [lindex $parts 1]
	    set maxval [lindex $parts 2]
	    if {[lsearch -exact $ed(trim,cols) $colname] >= 0} {
		set ed(trim,$colname,min) $minval
		set ed(trim,$colname,max) $maxval
	    }
	}
    }
    close $fd
}

proc CatalogPanelTrimApply {} {
    global catpanel
    global current

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set lines [split $catpanel(alldata) \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]
    set ncols [llength $headers]

    # Build list of active trim conditions
    set conditions {}
    for {set i 0} {$i < $ncols} {incr i} {
	set colname [string trim [lindex $headers $i]]
	set has_min 0
	set has_max 0
	set minval 0
	set maxval 0
	if {[info exists catpanel(trim,$colname,min)] && $catpanel(trim,$colname,min) ne {}} {
	    if {[string is double -strict $catpanel(trim,$colname,min)]} {
		set has_min 1
		set minval $catpanel(trim,$colname,min)
	    }
	}
	if {[info exists catpanel(trim,$colname,max)] && $catpanel(trim,$colname,max) ne {}} {
	    if {[string is double -strict $catpanel(trim,$colname,max)]} {
		set has_max 1
		set maxval $catpanel(trim,$colname,max)
	    }
	}
	if {$has_min || $has_max} {
	    lappend conditions [list $i $has_min $minval $has_max $maxval]
	}
    }

    # If no conditions, show all
    if {[llength $conditions] == 0} {
	set catpanel(trim,active) 0
	CatalogPanelLoadTSV $catpanel(alldata) "all"
	set catpanel(status) "Trim cleared - showing all sources"
	return
    }

    # Filter rows
    set filtered $header
    set count 0
    set total 0
    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	incr total
	set fields [split $line "\t"]
	set pass 1

	foreach cond $conditions {
	    set cidx [lindex $cond 0]
	    set has_min [lindex $cond 1]
	    set minval [lindex $cond 2]
	    set has_max [lindex $cond 3]
	    set maxval [lindex $cond 4]

	    set val [string trim [lindex $fields $cidx]]
	    if {![string is double -strict $val]} {
		set pass 0
		break
	    }
	    if {$has_min && $val < $minval} {
		set pass 0
		break
	    }
	    if {$has_max && $val > $maxval} {
		set pass 0
		break
	    }
	}

	if {$pass} {
	    append filtered "\n$line"
	    incr count
	}
    }

    set catpanel(trim,active) 1
    CatalogPanelLoadTSV $filtered "trimmed"

    # Re-mark from authoritative data
    CatalogPanelCreateAllMarkers

    set catpanel(status) "Trimmed: $count of $total sources match conditions"
}

proc CatalogPanelTrimSave {} {
    global catpanel

    set prefdir [file join [file normalize ~] .ds9]
    if {![file isdirectory $prefdir]} {
	file mkdir $prefdir
    }
    set preffile [file join $prefdir sextract_trim.prf]
    if {[catch {set fd [open $preffile w]} err]} return

    set lines [split $catpanel(alldata) \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]
    foreach h $headers {
	set colname [string trim $h]
	set minval {}
	set maxval {}
	if {[info exists catpanel(trim,$colname,min)]} {
	    set minval $catpanel(trim,$colname,min)
	}
	if {[info exists catpanel(trim,$colname,max)]} {
	    set maxval $catpanel(trim,$colname,max)
	}
	puts $fd "$colname\t$minval\t$maxval"
    }
    close $fd
}

proc CatalogPanelTrimLoad {} {
    global catpanel

    set preffile [file join [file normalize ~] .ds9 sextract_trim.prf]
    if {![file exists $preffile]} return
    if {[catch {set fd [open $preffile r]} err]} return

    while {[gets $fd line] >= 0} {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set parts [split $line "\t"]
	if {[llength $parts] >= 3} {
	    set colname [lindex $parts 0]
	    set catpanel(trim,$colname,min) [lindex $parts 1]
	    set catpanel(trim,$colname,max) [lindex $parts 2]
	}
    }
    close $fd
}

proc ThemeConfigCanvas {w} {
    global ds9

    $w configure -bg [ThemeTreeBackground]

    $w itemconfigure colorbar -fg [ThemeTreeForeground]
    $w itemconfigure colorbar -bg [ThemeTreeBackground]

    foreach ff $ds9(frames) {
	$w itemconfigure $ff -fg [ThemeTreeForeground]
	$w itemconfigure $ff -bg [ThemeTreeBackground]

	$w itemconfigure ${ff}cb -fg [ThemeTreeForeground]
	$w itemconfigure ${ff}cb -bg [ThemeTreeBackground]

	# since graphs are created, but maybe not realized
	# must update manually
	set varname ${ff}gr
	global $varname
	ThemeConfigGraph [subst $${varname}(horz)]
	ThemeConfigGraph [subst $${varname}(vert)]
    }
}

proc InitCanvas {} {
    global ds9

    # must wait until now
    bind $ds9(canvas) <Configure> [list LayoutView]
    BindEventsCanvas
}

proc BindEventsCanvas {} {
    global ds9

    # Bindings
    bind $ds9(canvas) <Tab> [list NextFrame]
    bind $ds9(canvas) <Shift-Tab> [list PrevFrame]
    switch $ds9(wm) {
	x11 {bind $ds9(canvas) <ISO_Left_Tab> [list PrevFrame]}
	aqua -
	win32 {}
    }

    # iis
    bind $ds9(canvas) <Key> {}
    # freeze
    bind $ds9(canvas) <f> {ToggleFreeze}

    # keyboard focus
    switch $ds9(wm) {
	x11 -
	aqua {
	    bind $ds9(canvas) <Enter> [list focus $ds9(canvas)]
	    bind $ds9(canvas) <Leave> [list focus {}]
	}
	win32 {}
    }
    switch $ds9(wm) {
	x11 {}
	aqua -
	win32 {bind $ds9(canvas) <MouseWheel> [list MouseWheelFrame %x %y %D]}
    }

    # backward compatible bindings
    switch $ds9(wm) {
	x11 -
	win32 {
	    bind $ds9(canvas) <Button-3> {Button3Canvas %x %y}
	    bind $ds9(canvas) <B3-Motion> {Motion3Canvas %x %y}
	    bind $ds9(canvas) <ButtonRelease-3> {Release3Canvas %x %y}
	}
	aqua {
	    # swap button-2 and button-3 on the mighty mouse
	    bind $ds9(canvas) <Button-2> {Button3Canvas %x %y}
	    bind $ds9(canvas) <B2-Motion> {Motion3Canvas %x %y}
	    bind $ds9(canvas) <ButtonRelease-2> {Release3Canvas %x %y}

	    # x11 command key emulation
	    bind $ds9(canvas) <Command-Button-1> {Button3Canvas %x %y}
	    bind $ds9(canvas) <Command-B1-Motion> {Motion3Canvas %x %y}
	    bind $ds9(canvas) <Command-ButtonRelease-1> {Release3Canvas %x %y}
	}
    }
}

proc UnBindEventsCanvas {} {
    global ds9

    # Bindings
    bind $ds9(canvas) <Tab> {}
    bind $ds9(canvas) <Shift-Tab> {}
    switch $ds9(wm) {
	x11 {bind $ds9(canvas) <ISO_Left_Tab> {}}
	aqua -
	win32 {}
    }

    # iis
    bind $ds9(canvas) <Key> {}
    # freeze
    bind $ds9(canvas) <f> {}

    # keyboard focus
    switch $ds9(wm) {
	x11 -
	aqua {
	    bind $ds9(canvas) <Enter> {}
	    bind $ds9(canvas) <Leave> {}
	}
	win32 {}
    }
    switch $ds9(wm) {
	x11 {}
	aqua -
	win32 {bind $ds9(canvas) <MouseWheel> {}}
    }

    # backward compatible bindings
    switch $ds9(wm) {
	x11 -
	win32 {
	    bind $ds9(canvas) <Button-3> {}
	    bind $ds9(canvas) <B3-Motion> {}
	    bind $ds9(canvas) <ButtonRelease-3> {}
	}
	aqua {
	    # swap button-2 and button-3 on the mighty mouse
	    bind $ds9(canvas) <Button-2> {}
	    bind $ds9(canvas) <B2-Motion> {}
	    bind $ds9(canvas) <ButtonRelease-2> {}

	    # x11 command key emulation
	    bind $ds9(canvas) <Command-Button-1> {}
	    bind $ds9(canvas) <Command-B1-Motion> {}
	    bind $ds9(canvas) <Command-ButtonRelease-1> {}
	}
    }
}

proc Button3Canvas {x y} {
    global ds9
    global current

    global debug
    if {$debug(tcl,events)} {
	puts stderr "Button3Canvas"
    }

    set ds9(b3) 1
    if {$current(frame) != {}} {
	ColorbarButton3 $current(frame) $x $y
    }
}

proc Motion3Canvas {x y} {
    global ds9
    global current

    global debug
    if {$debug(tcl,events)} {
	puts stderr "Motion3Canvas"
    }

    if {$current(frame) != {}} {
	ColorbarMotion3 $current(frame) $x $y
    }
}

proc Release3Canvas {x y} {
    global ds9
    global current

    global debug
    if {$debug(tcl,events)} {
	puts stderr "Release3Canvas"
    }

    set ds9(b3) 0
    if {$current(frame) != {}} {
	ColorbarRelease3 $current(frame) $x $y
    }
}

proc UnBindEventsCanvasItems {} {
    global ds9

    foreach ff $ds9(active) {
	UnBindEventsFrame $ff
	UnBindEventsColorbar ${ff}cb
	UnBindEventsGraph $ff
    }
}

proc BindEventsCanvasItems {} {
    global ds9

    foreach ff $ds9(active) {
	BindEventsFrame $ff
	BindEventsColorbar ${ff}cb
	BindEventsGraph $ff
    }
}

proc LayoutRaise {id} {
    global ds9

    set ll [$ds9(canvas) find withtag {graphic}]
    if {$ll != {}} {
	$ds9(canvas) lower $id [lindex $ll 0]
    } else {
	$ds9(canvas) raise $id
    }
}

proc LayoutView {} {
    global view

    global debug
    if {$debug(tcl,layout)} {
	puts stderr "LayoutView"
    }

    LayoutViewInit
    switch $view(layout) {
	horizontal {LayoutViewHorz}
	vertical {LayoutViewVert}
	basic {LayoutViewBasic}
	advanced {LayoutViewAdvanced}
    }

    LayoutInfoPanel
    LayoutButtons
    LayoutFrames

    UpdateViewMenu
}

proc LayoutViewInit {} {
    global ds9

    # reset weights
    grid rowconfigure $ds9(main) 0 -weight 0
    grid columnconfigure $ds9(main) 0 -weight 0
    grid rowconfigure $ds9(main) 2 -weight 0
    grid columnconfigure $ds9(main) 2 -weight 0
    grid rowconfigure $ds9(main) 4 -weight 0
    grid columnconfigure $ds9(main) 4 -weight 0

    grid forget $ds9(image)
    grid forget $ds9(header)
    grid forget $ds9(header,sep)
    grid forget $ds9(buttons,frame)
    grid forget $ds9(buttons,sep)
    grid forget $ds9(icons,top)
    grid forget $ds9(icons,top,sep)
    grid forget $ds9(icons,left)
    grid forget $ds9(icons,left,sep)
    grid forget $ds9(icons,bottom)
    grid forget $ds9(icons,bottom,sep)

    pack forget $ds9(panner)
    pack forget $ds9(panner,align)
    pack forget $ds9(panner,center)
    pack forget $ds9(magnifier)
    pack forget $ds9(magnifier,plus)
    pack forget $ds9(magnifier,minus)
    pack forget $ds9(info)
}

proc LayoutViewHorz {} {
    global ds9
    global view

    # ds9(main) weight
    grid rowconfigure $ds9(main) 4 -weight 1
    grid columnconfigure $ds9(main) 0 -weight 1

    # info panel
    if {$view(info) || $view(magnifier) || $view(panner)} {
	grid $ds9(header) -row 0 -column 0 -sticky ew
	$ds9(header,sep) configure -orient horizontal
	grid $ds9(header,sep) -row 1 -column 0 -sticky ew
    }

    if {$view(info)} {
	pack $ds9(info) -side left -anchor nw -padx 2 -pady 2 \
	    -fill x -expand true
    }

    if {$view(panner)} {
	pack $ds9(panner) -side right -padx 2 -pady 2
    }

    if {$view(magnifier)} {
	pack $ds9(magnifier) -side right -padx 2 -pady 2
	if {$view(panner)} {
	    pack $ds9(magnifier) -before $ds9(panner)
	}
    }

    # buttons
    if {$view(buttons)} {
	grid $ds9(buttons,frame) -row 2 -sticky ew -columnspan 3
	$ds9(buttons,sep) configure -orient horizontal
	grid $ds9(buttons,sep) -row 3 -column 0 -sticky ew -columnspan 3
    }

    # image
    grid $ds9(image) -row 4 -column 0 -sticky news
}

proc LayoutViewVert {} {
    global ds9
    global view

    # ds9(main) weight
    grid rowconfigure $ds9(main) 0 -weight 1
    grid columnconfigure $ds9(main) 4 -weight 1

    # info panel
    if {$view(info) || $view(magnifier) || $view(panner)} {
	grid $ds9(header) -row 0 -column 0 -sticky ns
	$ds9(header,sep) configure -orient vertical
	grid $ds9(header,sep) -row 0 -column 1 -sticky ns
    }

    if {$view(magnifier)} {
	pack $ds9(magnifier) -side top -padx 2 -pady 2
    }

    if {$view(info)} {
	pack $ds9(info) -side top -padx 2 -pady 2 -fill y -expand true
	if {$view(magnifier)} {
	    pack $ds9(info) -after $ds9(magnifier)
	}
    }

    if {$view(panner)} {
	pack $ds9(panner) -side bottom -padx 2 -pady 2
    }

    # buttons
    if {$view(buttons)} {
	grid $ds9(buttons,frame) -row 0 -column 2 -sticky ns
	$ds9(buttons,sep) configure -orient vertical
	grid $ds9(buttons,sep) -row 0 -column 3 -sticky ns
    }

    # image
    grid $ds9(image) -row 0 -column 4 -sticky news
}

proc LayoutViewBasic {} {
    global ds9
    global view

    # ds9(main) weight
    grid rowconfigure $ds9(main) 0 -weight 1
    grid columnconfigure $ds9(main) 0 -weight 1

    # image
    grid $ds9(image) -row 0 -column 0 -sticky news
}

proc LayoutViewAdvanced {} {
    global ds9
    global view

    # ds9(main) weight
    grid rowconfigure $ds9(main) 2 -weight 1
    grid columnconfigure $ds9(main) 2 -weight 1

    # info panel
    if {$view(info) || $view(magnifier) || $view(panner)} {
	$ds9(header,sep) configure -orient vertical
	grid $ds9(header,sep) -row 2 -column 3 -sticky ns
	grid $ds9(header) -row 2 -column 4 -sticky ns
    }

    if {$view(panner)} {
	pack $ds9(panner) -side top -padx 2 -pady 2
	if {$view(icons)} {
	    pack $ds9(panner,align) -side left
	    pack $ds9(panner,center) -side left
	}
    }

    if {$view(magnifier)} {
	pack $ds9(magnifier) -side top -padx 2 -pady 2
	if {$view(icons)} {
	    pack $ds9(magnifier,minus) -side left
	    pack $ds9(magnifier,plus) -side left
	}
    }

    if {$view(info)} {
	pack $ds9(info) -side bottom -padx 2 -pady 2 -fill y -expand true
	if {$view(magnifier)} {
	    pack $ds9(info) -after $ds9(magnifier)
	}
    }

    # buttons
    if {$view(buttons)} {
	$ds9(buttons,sep) configure -orient vertical
	grid $ds9(buttons,sep) -row 2 -column 5 -sticky ns
	grid $ds9(buttons,frame) -row 2 -column 6 -sticky ns
    }

    # icons
    if {$view(icons)} {
	grid $ds9(icons,top) -row 0 -column 0 -sticky ew -columnspan 7
	grid $ds9(icons,top,sep) -row 1 -column 0 -sticky ew -columnspan 7
	grid $ds9(icons,left) -row 2 -column 0 -sticky ns
	grid $ds9(icons,left,sep) -row 2 -column 1 -sticky ns
	grid $ds9(icons,bottom,sep) -row 3 -column 0 -sticky ew -columnspan 7
	grid $ds9(icons,bottom) -row 4 -column 0 -sticky ew -columnspan 7
    }

    # image
    grid $ds9(image) -row 2 -column 2 -sticky news
}

proc LayoutFrames {} {
    global ds9
    global current
    global tile
    global view
    global colorbar

    # turn off default colorbar
    colorbar hide

    # turn off default graphs
    GraphHide graph horz
    GraphHide graph vert

    # all frames turn everything off
    foreach ff $ds9(frames) {
	$ff hide
	$ff highlite off
	$ff panner off
	$ff magnifier off

	# colorbar
	${ff}cb hide

	# graphs
	GraphHide $ff horz
	GraphHide $ff vert
    }

    # be sure colorbar/graph sizes are correct
    LayoutColorbarAdjust
    LayoutGraphsAdjust

    if {[llength $ds9(active)] > 0} {
	LayoutFramesOneOrMore
    } else {
	LayoutFramesNone
    }

    # after all layed out, update data cut for graphs if needed
    #  one problem- if single mode, non-current graphs are incorrectly updated
    switch -- $current(mode) {
	crosshair {
	    if {$view(graph,horz) || $view(graph,vert)} {
		update idletasks
		foreach ff $ds9(active) {
		    set vv [$ff get crosshair canvas]
		    UpdateGraphsData $ff [lindex $vv 0] [lindex $vv 1] canvas
		}
	    }
	}
    }
}

proc LayoutFramesNone {} {
    global ds9
    global current
    global colorbar
    global view

    catch {CatalogPanelSaveFrameState $current(frame)}
    set current(frame) {}
    set current(colorbar) colorbar

    set colorbar(map) [colorbar get name]
    set colorbar(invert) [colorbar get invert]

    # panner
    if {$view(panner)} {
	panner clear
    }

    # magnifier
    if {$view(magnifier)} {
	magnifier clear
    }

    # colorbar
    if {$view(colorbar)} {
	if {[LayoutColorbar colorbar 0 0 [winfo width $ds9(canvas)] [winfo height $ds9(canvas)]]} {
	    colorbar show
	    LayoutRaise colorbar
#	    $ds9(canvas) raise colorbar
	}
    }

    # graphs
    if {$view(graph,horz)} {
	LayoutGraphHorz graph 0 0 \
	    [winfo width $ds9(canvas)] [winfo height $ds9(canvas)]
	GraphShow graph horz
    }
    if {$view(graph,vert)} {
	LayoutGraphVert graph 0 0 \
	    [winfo width $ds9(canvas)] [winfo height $ds9(canvas)]
	GraphShow graph vert
    }

    # update menus/dialogs
    UpdateDS9
}

proc LayoutFramesOneOrMore {} {
    global ds9
    global view

    switch -- $ds9(display) {
	fade -
	blink -
	single {LayoutFrameOne}
	tile {
	    if {[llength $ds9(active)] > 1} {
		if {$view(multi)} {
		    LayoutFrame
		} else {
		    LayoutFrameNone
		}
	    } else {
		LayoutFrameOne
	    }
	}
    }
}

proc LayoutFrameOne {} {
    global ds9
    global view
    global current
    global colorbar

    set ww [winfo width $ds9(canvas)]
    set hh [winfo height $ds9(canvas)]

    foreach ff $ds9(active) {
	set fw $ww
	set fh $hh

	# frame
	LayoutFrameAdjust fw fh
	$ff configure -x 0 -y 0 -width $fw -height $fh -anchor nw

	# colorbar
	if {$view(colorbar)} {
	    LayoutColorbar ${ff}cb 0 0 $ww $hh
	}

	# graphs
	if {$view(graph,horz)} {
	    LayoutGraphHorz $ff 0 0 $ww $hh
	    UpdateGraphAxis $ff horz
	}
	if {$view(graph,vert)} {
	    LayoutGraphVert $ff 0 0 $ww $hh
	    UpdateGraphAxis $ff vert
    	}
    }

    # frame
    $current(frame) show
    LayoutRaise $current(frame)
#    $ds9(canvas) raise $current(frame)

    # colorbar
    if {$view(colorbar)} {
	$current(colorbar) show
	LayoutRaise $current(colorbar)
#	$ds9(canvas) raise $current(colorbar)
    }

    # graphs
    if {$view(graph,horz)} {
	GraphShow $current(frame) horz
    }
    if {$view(graph,vert)} {
	GraphShow $current(frame) vert
    }

    FrameToFront
}

proc LayoutFrame {} {
    global ds9
    global tile

    set num [llength $ds9(active)]
    switch -- $tile(mode) {
	row {
	    TileRect 1 $num
	}
	column {
	    TileRect $num 1
	}
	grid {
	    switch -- $tile(grid,mode) {
		automatic {
		    TileRect \
			[expr int(sqrt($num-1))+1] [expr int(sqrt($num)+.5)]
		}
		manual {
		    set cnt [expr $tile(grid,col)*$tile(grid,row)]
		    if {[llength $ds9(active)] > $cnt} {
			Error "Too many Frames to display manual, using automatic"
			TileRect \
			    [expr int(sqrt($num-1))+1] [expr int(sqrt($num)+.5)]
		    } else {
			TileRect $tile(grid,col) $tile(grid,row)
		    }
		}
	    }
	}
    }
}

proc LayoutFrameNone {} {
    global ds9
    global tile

    set num [llength $ds9(active)]
    switch -- $tile(mode) {
	row {
	    TileRectNone 1 $num
	}
	column {
	    TileRectNone $num 1
	}
	grid {
	    switch -- $tile(grid,mode) {
		automatic {
		    TileRectNone \
			[expr int(sqrt($num-1))+1] [expr int(sqrt($num)+.5)]
		}
		manual {
		    set cnt [expr $tile(grid,col)*$tile(grid,row)]
		    if {[llength $ds9(active)] > $cnt} {
			Error "Too many Frames to display manual, using automatic"
			TileRectNone \
			    [expr int(sqrt($num-1))+1] [expr int(sqrt($num)+.5)]
		    } else {
			TileRectNone $tile(grid,col) $tile(grid,row)
		    }
		}
	    }
	}
    }
}

proc TileRect {numx numy} {
    global ds9
    global tile
    global current
    global view
    global colorbar

    set ww [expr int(([winfo width  $ds9(canvas)]-($tile(grid,gap)*($numx-1)))/$numx)]
    set hh [expr int(([winfo height $ds9(canvas)]-($tile(grid,gap)*($numy-1)))/$numy)]

    switch $tile(grid,dir) {
	x {
	    for {set jj 0} {$jj<$numy} {incr jj} {
		for {set ii 0} {$ii<$numx} {incr ii} {
		    set nn [expr $jj*$numx + $ii]
		    set xx($nn) [expr ($ww+$tile(grid,gap))*$ii]
		    set yy($nn) [expr ($hh+$tile(grid,gap))*$jj]
		}
	    }
	}
	y {
	    for {set ii 0} {$ii<$numx} {incr ii} {
		for {set jj 0} {$jj<$numy} {incr jj} {
		    set nn [expr $ii*$numy + $jj]
		    set xx($nn) [expr ($ww+$tile(grid,gap))*$ii]
		    set yy($nn) [expr ($hh+$tile(grid,gap))*$jj]
		}
	    }
	}
    }

    set ii 0
    foreach ff $ds9(active) {
	set fw $ww
	set fh $hh

	# frame
	LayoutFrameAdjust fw fh
	$ff configure -x $xx($ii) -y $yy($ii) -width $fw -height $fh -anchor nw
	$ff show
	LayoutRaise $ff
#	$ds9(canvas) raise $ff

	# colorbar
	if {$view(colorbar)} {
	    LayoutColorbar ${ff}cb $xx($ii) $yy($ii) $ww $hh
	    ${ff}cb show
	    LayoutRaise ${ff}cb
#	    $ds9(canvas) raise ${ff}cb
	}

	# graphs
	if {$view(graph,horz)} {
	    LayoutGraphHorz $ff $xx($ii) $yy($ii) $ww $hh
	    UpdateGraphAxis $ff horz
	    GraphShow $ff horz
	}
	if {$view(graph,vert)} {
	    LayoutGraphVert $ff $xx($ii) $yy($ii) $ww $hh
	    UpdateGraphAxis $ff vert
	    GraphShow $ff vert
	}

	incr ii
    }

    FrameToFront
}

proc TileRectNone {numx numy} {
    global ds9
    global tile
    global current
    global view
    global colorbar

    set fw [winfo width $ds9(canvas)]
    set fh [winfo height $ds9(canvas)]
    LayoutFrameAdjust fw fh

    set ww [expr int(($fw-($tile(grid,gap)*($numx-1)))/$numx)]
    set hh [expr int(($fh-($tile(grid,gap)*($numy-1)))/$numy)]

    switch $tile(grid,dir) {
	x {
	    for {set jj 0} {$jj<$numy} {incr jj} {
		for {set ii 0} {$ii<$numx} {incr ii} {
		    set nn [expr $jj*$numx + $ii]
		    set xx($nn) [expr ($ww+$tile(grid,gap))*$ii]
		    set yy($nn) [expr ($hh+$tile(grid,gap))*$jj]
		}
	    }
	}
	y {
	    for {set ii 0} {$ii<$numx} {incr ii} {
		for {set jj 0} {$jj<$numy} {incr jj} {
		    set nn [expr $ii*$numy + $jj]
		    set xx($nn) [expr ($ww+$tile(grid,gap))*$ii]
		    set yy($nn) [expr ($hh+$tile(grid,gap))*$jj]
		}
	    }
	}
    }

    # frames
    set ii 0
    set cnt [expr $numx*$numy]
    foreach ff $ds9(active) {
	# sanity check
	if {$xx($ii)>=0 && $yy($ii)>=0 && $ww>=0 && $hh>=0} {
	    $ff configure -x $xx($ii) -y $yy($ii) \
		-width $ww -height $hh -anchor nw
	    $ff show
	    LayoutRaise $ff
#	    $ds9(canvas) raise $ff
	}

	if {$view(colorbar)} {
	    LayoutColorbar ${ff}cb 0 0 \
		[winfo width $ds9(canvas)] [winfo height $ds9(canvas)]
	}

	if {$view(graph,horz)} {
	    LayoutGraphHorz $ff 0 0 \
		[winfo width $ds9(canvas)] [winfo height $ds9(canvas)]
	    UpdateGraphAxis $ff horz
	}

	if {$view(graph,vert)} {
	    LayoutGraphVert $ff 0 0 \
		[winfo width $ds9(canvas)] [winfo height $ds9(canvas)]
	    UpdateGraphAxis $ff vert
	}

	incr ii
	if {$ii>=$cnt} {
	    break
	}
    }

    # set colorbar/graph for current frame
    set ff $current(frame)

    # colorbar
    if {$view(colorbar)} {
	${ff}cb show
	LayoutRaise ${ff}cb
#	$ds9(canvas) raise ${ff}cb
    }

    # graphs
    if {$view(graph,horz)} {
	GraphShow $ff horz
    }
    if {$view(graph,vert)} {
	GraphShow $ff vert
    }

    FrameToFront
}

proc LayoutFrameAdjust {wvar hvar} {
    global canvas
    global view
    global colorbar
    global igraph
    global dgraph
    global graph

    upvar $wvar ww
    upvar $hvar hh

    set cbh [expr $view(colorbar) && !$colorbar(orientation)]
    set cbv [expr $view(colorbar) &&  $colorbar(orientation)]
    set grh $view(graph,horz)
    set grv $view(graph,vert)

    # cbh
    if {$cbh && !$cbv && !$grh && !$grv} {
	incr hh -$colorbar(horizontal,height)
	incr hh -$canvas(gap)
    }
    # cbv
    if {!$cbh && $cbv && !$grh && !$grv} {
	incr ww -$colorbar(vertical,width)
	incr ww -$canvas(gap)
    }

    # cbhgrh
    if {$cbh && !$cbv && $grh && !$grv} {
	incr hh -$colorbar(horizontal,height)
	incr hh -$canvas(gap)
	incr hh -$graph(size)
	incr ww -$dgraph(horz,offset)
    }
    # cbhgrv
    if {$cbh && !$cbv && !$grh && $grv} {
	incr hh -$colorbar(horizontal,height)
	incr hh -$canvas(gap)
	incr ww -$graph(size)
    }
    # cbhgrhgrv
    if {$cbh && !$cbv && $grh && $grv} {
	incr hh -$colorbar(horizontal,height)
	incr hh -$canvas(gap)
	incr hh -$graph(size)
	incr ww -$graph(size)
    }

    # cbvgrh
    if {!$cbh && $cbv && $grh && !$grv} {
	incr ww -$colorbar(vertical,width)
	incr ww -$canvas(gap)
	incr hh -$graph(size)
    }
    # cbvgrv
    if {!$cbh && $cbv && !$grh && $grv} {
	incr ww -$colorbar(vertical,width)
	incr ww -$canvas(gap)
	incr ww -$graph(size)
	incr hh -$dgraph(vert,offset)
    }
    # cbvgrhgrv
    if {!$cbh && $cbv && $grh && $grv} {
	incr ww -$colorbar(vertical,width)
	incr ww -$canvas(gap)
	incr ww -$graph(size)
	incr hh -$graph(size)
    }

    # grh
    if {!$cbh && !$cbv && $grh && !$grv} {
	incr hh -$graph(size)
	incr hh -$canvas(gap)
	incr ww -$dgraph(horz,offset)
    }
    # grv
    if {!$cbh && !$cbv && !$grh && $grv} {
	incr ww -$graph(size)
	incr ww -$canvas(gap)
	incr hh -$dgraph(vert,offset)
    }
    # grhgrv
    if {!$cbh && !$cbv && $grh && $grv} {
	incr ww -$graph(size)
	incr ww -$canvas(gap)
	incr hh -$graph(size)
	incr hh -$canvas(gap)
    }

    # sanity check
    if {$ww<0} {
	set ww 1
    }
    if {$hh<0} {
	set hh 1
    }
}

proc LayoutChangeWidth {ww} {
    global ds9

    set cw [winfo width $ds9(canvas)]
    set tw [winfo width $ds9(top)]
    set th [winfo height $ds9(top)]
    set dw $ww-$cw

    # change window size
    wm geometry $ds9(top) "[expr $tw+$dw]x${th}"
    LayoutView
}

proc LayoutChangeHeight {hh} {
    global ds9

    set ch [winfo height $ds9(canvas)]
    set tw [winfo width $ds9(top)]
    set th [winfo height $ds9(top)]
    set dh $hh-$ch

    # change window size
    wm geometry $ds9(top) "${tw}x[expr $th+$dh]"
    LayoutView
}

proc LayoutChangeSize {ww hh} {
    global ds9

    set cw [winfo width $ds9(canvas)]
    set ch [winfo height $ds9(canvas)]
    set tw [winfo width $ds9(top)]
    set th [winfo height $ds9(top)]
    set dw $ww-$cw
    set dh $hh-$ch

    # change window size
    wm geometry $ds9(top) "[expr $tw+$dw]x[expr $th+$dh]"
    LayoutView
}

proc DisplayDefaultDialog {} {
    global ed
    global ds9

    set w {.defdpy}

    set ed(ok) 0
    set ed(x) [winfo width $ds9(canvas)]
    set ed(y) [winfo height $ds9(canvas)]

    DialogCreate $w [msgcat::mc {Display Size}] ed(ok)

    # Param
    set f [ttk::frame $w.param]

    ttk::label $f.xTitle -text {X}
    ttk::label $f.yTitle -text {Y}
    ttk::entry $f.x -textvariable ed(x) -width 10
    ttk::entry $f.y -textvariable ed(y) -width 10
    ttk::label $f.xunit -text [msgcat::mc {Pixels}]
    ttk::label $f.yunit -text [msgcat::mc {Pixels}]

    grid $f.xTitle $f.x $f.xunit -padx 2 -pady 2 -sticky w
    grid $f.yTitle $f.y $f.yunit -padx 2 -pady 2 -sticky w

    # Buttons
    set f [ttk::frame $w.buttons]
    ttk::button $f.ok -text [msgcat::mc {OK}] -command {set ed(ok) 1} \
	-default active
    ttk::button $f.cancel -text [msgcat::mc {Cancel}] -command {set ed(ok) 0}
    pack $f.ok $f.cancel -side left -expand true -padx 2 -pady 4

    bind $w <Return> {set ed(ok) 1}

    # Fini
    ttk::separator $w.sep -orient horizontal
    pack $w.buttons $w.sep -side bottom -fill x
    pack $w.param -side top -fill both -expand true

    $w.param.x select range 0 end
    DialogWait $w ed(ok) $w.param.x
    destroy $w

    if {$ed(ok)} {
	LayoutChangeSize $ed(x) $ed(y)
    }

    set rr $ed(ok)
    unset ed
    return $rr
}

# Process Cmds

proc ProcessHeightCmd {varname iname} {
    upvar $varname var
    upvar $iname i

    # we need to be realized
    # can't use ProcessRealize
    RealizeDS9

    height::YY_FLUSH_BUFFER
    height::yy_scan_string [lrange $var $i end]
    height::yyparse
    incr i [expr $height::yycnt-1]
}

proc ProcessSendHeightCmd {proc id param {sock {}} {fn {}}} {
    global ds9
    $proc $id "[winfo height $ds9(canvas)]\n"
}

proc ProcessWidthCmd {varname iname} {
    upvar $varname var
    upvar $iname i

    # we need to be realized
    # can't use ProcessRealize
    RealizeDS9

    width::YY_FLUSH_BUFFER
    width::yy_scan_string [lrange $var $i end]
    width::yyparse
    incr i [expr $width::yycnt-1]
}

proc ProcessSendWidthCmd {proc id param {sock {}} {fn {}}} {
    global ds9
    $proc $id "[winfo width $ds9(canvas)]\n"
}

proc ProcessViewCmd {varname iname} {
    upvar $varname var
    upvar $iname i

    view::YY_FLUSH_BUFFER
    view::yy_scan_string [lrange $var $i end]
    view::yyparse
    incr i [expr $view::yycnt-1]
}

proc ProcessSendViewCmd {proc id param {sock {}} {fn {}}} {
    global parse
    set parse(proc) $proc
    set parse(id) $id

    viewsend::YY_FLUSH_BUFFER
    viewsend::yy_scan_string $param
    viewsend::yyparse
}

# --- Galaxy Model Fitting ---

proc CatalogPanelGalaxyFit {model} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "No sources — run Extract first"
	return
    }

    set catpanel(status) "Galaxy $model fitting — not yet implemented"
}

proc CatalogPanelGalaxyParams {} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "No sources — run Extract first"
	return
    }

    set catpanel(status) "Galaxy parameter extraction — not yet implemented"
}

proc CatalogPanelGalaxyMorphology {} {
    global catpanel
    global current
    global ds9

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "No sources — run Extract first"
	return
    }

    # Get current FITS filename
    set fn {}
    if {$current(frame) != {}} {
	catch {set fn [$current(frame) get fits file name full]}
    }
    set fn [string trim $fn "{}"]
    regsub {\[.*\]$} $fn {} fn
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    # Find ds9_galaxy_morph.py
    set bindir [file dirname [info nameofexecutable]]
    set script [file join $bindir ds9_galaxy_morph.py]
    if {![file exists $script]} {
	set libdir [file join [file dirname $bindir] ds9 library]
	set script [file join $libdir ds9_galaxy_morph.py]
    }
    if {![file exists $script]} {
	set catpanel(status) "ERROR: ds9_galaxy_morph.py not found"
	return
    }

    # Save catalog to temp TSV
    set catfile [file join [file normalize ~] .ds9 morph_catalog.tsv]
    catch {file mkdir [file dirname $catfile]}
    if {[catch {
	set fd [open $catfile w]
	puts $fd $catpanel(alldata)
	close $fd
    } err]} {
	set catpanel(status) "Morphology error: cannot write catalog: $err"
	return
    }

    # Build arguments
    set paramargs {}
    lappend paramargs "--catalog" $catfile

    # Find checkpoint
    set ckpt [file join [file dirname $bindir] galaxy_morph data checkpoints cnn_morph_best.pt]
    if {[file exists $ckpt]} {
	lappend paramargs "--checkpoint" $ckpt
    }

    set catpanel(status) "Morphology: classifying sources on [file tail $fn] ..."
    update idletasks

    # Run classification
    set errfile [file join [file normalize ~] .ds9 morph_stderr.txt]
    if {[catch {set data [exec python3 $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	catch {file delete $catfile}
	if {$stderr_msg ne ""} {
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    set catpanel(status) "Morphology error: $last_err"
	    puts "Morphology full stderr:\n$stderr_msg"
	} else {
	    set catpanel(status) "Morphology error: $err"
	}
	return
    }
    catch {file delete $errfile}
    catch {file delete $catfile}

    # Parse results
    CatalogPanelMorphParseResults $data
}

proc CatalogPanelMorphParseResults {data} {
    global catpanel
    global current

    set lines [split $data \n]
    set n_classified 0

    # Parse header: #GALAXY_MORPH	N_CLASSIFIED=245	N_SOURCES=300
    foreach line $lines {
	if {[string match "#GALAXY_MORPH*" $line]} {
	    foreach field [split $line "\t"] {
		if {[string match "N_CLASSIFIED=*" $field]} {
		    set n_classified [string range $field 13 end]
		}
	    }
	    continue
	}
    }

    # Build morph_map: NUMBER -> {morph_type morph_conf color}
    array unset catpanel morph,*
    set catpanel(morph,map) {}

    foreach line $lines {
	if {[string match "#*" $line]} continue
	if {[string match "NUMBER*" $line]} continue
	if {[string trim $line] eq {}} continue

	# NUMBER MORPH_TYPE MORPH_DESC MORPH_CONF TOP1_CLASS TOP1_PROB ... MORPH_COLOR
	set fields [split $line "\t"]
	if {[llength $fields] < 11} continue

	set src_num [lindex $fields 0]
	set morph_type [lindex $fields 1]
	set morph_desc [lindex $fields 2]
	set morph_conf [lindex $fields 3]
	set color [lindex $fields 10]

	set catpanel(morph,$src_num) [list $morph_type $morph_desc $morph_conf $color]
	lappend catpanel(morph,map) $src_num
    }

    if {[llength $catpanel(morph,map)] == 0} {
	set catpanel(status) "Morphology: no galaxies classified"
	return
    }

    # Add MORPH_TYPE and MORPH_CONF columns to alldata
    CatalogPanelMorphAddColumns

    # Recolor markers by morphology
    CatalogPanelMorphColorMarkers

    set catpanel(status) "Morphology: $n_classified galaxies classified"
}

proc CatalogPanelMorphAddColumns {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set lines [split $catpanel(alldata) \n]
    if {[llength $lines] < 2} return

    # Parse header
    set headers [split [lindex $lines 0] "\t"]
    set ncols [llength $headers]

    # Find NUMBER column
    set col_num -1
    for {set c 0} {$c < $ncols} {incr c} {
	if {[string trim [lindex $headers $c]] eq "NUMBER"} {
	    set col_num $c
	    break
	}
    }
    if {$col_num < 0} return

    # Check if columns already exist (avoid duplicates)
    set has_morph 0
    for {set c 0} {$c < $ncols} {incr c} {
	if {[string trim [lindex $headers $c]] eq "MORPH_TYPE"} {
	    set has_morph 1
	    break
	}
    }

    # Build new data with added columns
    set newdata {}

    if {$has_morph} {
	# Find existing MORPH_TYPE, MORPH_DESC, MORPH_CONF column indices
	set col_mt -1
	set col_md -1
	set col_mc -1
	for {set c 0} {$c < $ncols} {incr c} {
	    set h [string trim [lindex $headers $c]]
	    if {$h eq "MORPH_TYPE"} { set col_mt $c }
	    if {$h eq "MORPH_DESC"} { set col_md $c }
	    if {$h eq "MORPH_CONF"} { set col_mc $c }
	}

	# Update existing columns
	append newdata [lindex $lines 0]
	for {set i 1} {$i < [llength $lines]} {incr i} {
	    set line [lindex $lines $i]
	    if {[string trim $line] eq {}} continue
	    set fields [split $line "\t"]
	    set src_num [string trim [lindex $fields $col_num]]

	    if {[info exists catpanel(morph,$src_num)]} {
		set info $catpanel(morph,$src_num)
		if {$col_mt >= 0} {
		    lset fields $col_mt [lindex $info 0]
		}
		if {$col_md >= 0} {
		    lset fields $col_md [lindex $info 1]
		}
		if {$col_mc >= 0} {
		    lset fields $col_mc [lindex $info 2]
		}
	    }
	    append newdata "\n" [join $fields "\t"]
	}
    } else {
	# Add new columns
	append newdata [lindex $lines 0] "\tMORPH_TYPE\tMORPH_DESC\tMORPH_CONF"
	for {set i 1} {$i < [llength $lines]} {incr i} {
	    set line [lindex $lines $i]
	    if {[string trim $line] eq {}} continue
	    set fields [split $line "\t"]
	    set src_num [string trim [lindex $fields $col_num]]

	    set mt ""
	    set md ""
	    set mc ""
	    if {[info exists catpanel(morph,$src_num)]} {
		set info $catpanel(morph,$src_num)
		set mt [lindex $info 0]
		set md [lindex $info 1]
		set mc [lindex $info 2]
	    }
	    append newdata "\n" $line "\t" $mt "\t" $md "\t" $mc
	}
    }

    set catpanel(alldata) $newdata
    CatalogPanelLoadTSV $catpanel(alldata) "morphology"
}

proc CatalogPanelMorphColorMarkers {} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    set frame $current(frame)

    # Delete previous sextract_all markers and recreate with morph colors
    catch {$frame marker catalog sextract_all delete}

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set lines [split $catpanel(alldata) \n]
    if {[llength $lines] < 2} return

    set headers [split [lindex $lines 0] "\t"]
    set ncols [llength $headers]

    # Find needed columns
    set col_x -1
    set col_y -1
    set col_a -1
    set col_b -1
    set col_theta -1
    set col_ir -1
    set col_num -1
    for {set c 0} {$c < $ncols} {incr c} {
	set hdr [string trim [lindex $headers $c]]
	switch -- $hdr {
	    NUMBER      { set col_num $c }
	    X_IMAGE     { set col_x $c }
	    Y_IMAGE     { set col_y $c }
	    A_IMAGE     { set col_a $c }
	    B_IMAGE     { set col_b $c }
	    THETA_IMAGE { set col_theta $c }
	    ISO_RADIUS  { set col_ir $c }
	}
    }
    if {$col_x < 0 || $col_y < 0} return

    # Build region strings with morph-specific colors
    set batch_size 500
    set reg "image\n"
    set count 0
    set batch_count 0
    global sextract_all_reg

    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	set fields [split $line "\t"]

	set x [string trim [lindex $fields $col_x]]
	set y [string trim [lindex $fields $col_y]]
	if {![string is double -strict $x] || ![string is double -strict $y]} continue

	set src_num [expr {$i}]
	if {$col_num >= 0} {
	    set nv [string trim [lindex $fields $col_num]]
	    if {$nv ne {}} { set src_num $nv }
	}

	# Ellipse parameters
	set iso_radius 5.0
	set a_image 0
	set b_image 0
	set theta 0

	if {$col_ir >= 0} {
	    set val [string trim [lindex $fields $col_ir]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} {
		set iso_radius $v
	    }
	}
	if {$col_a >= 0} {
	    set val [string trim [lindex $fields $col_a]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set a_image $v }
	}
	if {$col_b >= 0} {
	    set val [string trim [lindex $fields $col_b]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set b_image $v }
	}
	if {$col_theta >= 0} {
	    set val [string trim [lindex $fields $col_theta]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0} { set theta $v }
	}

	set semi_a $iso_radius
	set semi_b $iso_radius
	if {$a_image > 0 && $b_image > 0} {
	    set semi_b [expr {$iso_radius * $b_image / $a_image}]
	}

	# Color: use morph color if classified, else yellow
	set color yellow
	if {[info exists catpanel(morph,$src_num)]} {
	    set color [lindex $catpanel(morph,$src_num) 3]
	}

	append reg "ellipse($x $y ${semi_a}i ${semi_b}i $theta) # color=$color width=1 tag={sextract_all} tag={sextract_src.$src_num} select=0 edit=0 move=0 rotate=0 delete=1 highlite=1 callback=highlite CatalogPanelMarkerCB {$src_num} callback=unhighlite CatalogPanelMarkerUnCB {$src_num}\n"
	incr count
	incr batch_count

	if {$batch_count >= $batch_size} {
	    set sextract_all_reg $reg
	    catch {$frame marker catalog command ds9 var sextract_all_reg}
	    set reg "image\n"
	    set batch_count 0
	}
    }

    if {$batch_count > 0} {
	set sextract_all_reg $reg
	catch {$frame marker catalog command ds9 var sextract_all_reg}
    }

    set catpanel(markall,on) 1
}

# ============================================================================
# AI Star/Galaxy Classification
# ============================================================================

proc CatalogPanelStarFinder {} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "No sources — run Extract first"
	return
    }

    # Get FITS filename
    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    # Find script
    set script [CatalogPanelGetScript ds9_star_finder.py]
    if {![file exists $script]} {
	set catpanel(status) "ERROR: ds9_star_finder.py not found"
	return
    }

    # Save catalog to temp TSV
    set catfile [CatalogPanelSaveTempCatalog star]
    if {$catfile eq {}} {
	set catpanel(status) "Star Finder: cannot write catalog"
	return
    }

    # Build arguments
    set paramargs {}
    lappend paramargs "--catalog" $catfile

    # Pass PSF file if available
    if {[info exists catpanel(psf,file)] && $catpanel(psf,file) ne {} &&
	[file exists $catpanel(psf,file)]} {
	lappend paramargs "--psf" $catpanel(psf,file)
    }

    set catpanel(status) "AI Star Classification: classifying sources on [file tail $fn] ..."
    update idletasks

    # Run classification
    set errfile [file join [file normalize ~] .ds9 star_finder_stderr.txt]
    if {[catch {set data [exec python3 $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	catch {file delete $catfile}
	if {$stderr_msg ne ""} {
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    set catpanel(status) "Star Finder error: $last_err"
	    puts "Star Finder full stderr:\n$stderr_msg"
	} else {
	    set catpanel(status) "Star Finder error: $err"
	}
	return
    }
    catch {file delete $errfile}
    catch {file delete $catfile}

    # Parse results and add columns
    CatalogPanelStarFinderParse $data
}

proc CatalogPanelStarFinderParse {data} {
    global catpanel

    set lines [split $data \n]
    set n_classified 0

    # Parse header: #STAR_FINDER	N_CLASSIFIED=245	N_SOURCES=300
    foreach line $lines {
	if {[string match "#STAR_FINDER*" $line]} {
	    foreach field [split $line "\t"] {
		if {[string match "N_CLASSIFIED=*" $field]} {
		    set n_classified [string range $field 13 end]
		}
	    }
	    continue
	}
    }

    # Collect result lines (skip header/comments)
    set result_lines {}
    set header_line ""
    foreach line $lines {
	if {[string match "#*" $line]} continue
	if {[string match "NUMBER*" $line]} {
	    set header_line $line
	    continue
	}
	if {[string trim $line] eq {}} continue
	lappend result_lines $line
    }

    if {[llength $result_lines] == 0} {
	set catpanel(status) "Star Finder: no sources classified"
	return
    }

    # Rebuild TSV result for AddColumnsFromTSV
    set result_data $header_line
    foreach line $result_lines {
	append result_data "\n" $line
    }

    # Add AI_STAR and AI_STAR_CONF columns
    CatalogPanelAddColumnsFromTSV $result_data {AI_STAR AI_STAR_CONF}

    # Recolor markers
    CatalogPanelStarFinderColorMarkers $result_lines

    set catpanel(status) "AI Star Classification: $n_classified sources classified"
}

proc CatalogPanelStarFinderColorMarkers {result_lines} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    set frame $current(frame)

    # Build color lookup: NUMBER -> color
    array set star_color {}
    foreach line $result_lines {
	set fields [split $line "\t"]
	if {[llength $fields] < 4} continue
	set src_num [lindex $fields 0]
	set color [lindex $fields 3]
	set star_color($src_num) $color
    }

    # Delete existing markers and recreate with star/galaxy colors
    catch {$frame marker catalog sextract_all delete}

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set lines [split $catpanel(alldata) \n]
    if {[llength $lines] < 2} return

    set headers [split [lindex $lines 0] "\t"]
    set ncols [llength $headers]

    # Find needed columns
    set col_x -1
    set col_y -1
    set col_a -1
    set col_b -1
    set col_theta -1
    set col_ir -1
    set col_num -1
    for {set c 0} {$c < $ncols} {incr c} {
	set hdr [string trim [lindex $headers $c]]
	switch -- $hdr {
	    NUMBER      { set col_num $c }
	    X_IMAGE     { set col_x $c }
	    Y_IMAGE     { set col_y $c }
	    A_IMAGE     { set col_a $c }
	    B_IMAGE     { set col_b $c }
	    THETA_IMAGE { set col_theta $c }
	    ISO_RADIUS  { set col_ir $c }
	}
    }
    if {$col_x < 0 || $col_y < 0} return

    # Build region strings with star/galaxy colors
    set batch_size 500
    set reg "image\n"
    set count 0
    set batch_count 0
    global sextract_all_reg

    for {set i 1} {$i < [llength $lines]} {incr i} {
	set line [lindex $lines $i]
	if {[string trim $line] eq {}} continue
	set fields [split $line "\t"]

	set x [string trim [lindex $fields $col_x]]
	set y [string trim [lindex $fields $col_y]]
	if {![string is double -strict $x] || ![string is double -strict $y]} continue

	set src_num [expr {$i}]
	if {$col_num >= 0} {
	    set nv [string trim [lindex $fields $col_num]]
	    if {$nv ne {}} { set src_num $nv }
	}

	# Ellipse parameters
	set iso_radius 5.0
	set a_image 0
	set b_image 0
	set theta 0

	if {$col_ir >= 0} {
	    set val [string trim [lindex $fields $col_ir]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} {
		set iso_radius $v
	    }
	}
	if {$col_a >= 0} {
	    set val [string trim [lindex $fields $col_a]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set a_image $v }
	}
	if {$col_b >= 0} {
	    set val [string trim [lindex $fields $col_b]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set b_image $v }
	}
	if {$col_theta >= 0} {
	    set val [string trim [lindex $fields $col_theta]]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0} { set theta $v }
	}

	set semi_a $iso_radius
	set semi_b $iso_radius
	if {$a_image > 0 && $b_image > 0} {
	    set semi_b [expr {$iso_radius * $b_image / $a_image}]
	}

	# Color: star=cyan, galaxy=yellow, unclassified=yellow
	set color yellow
	if {[info exists star_color($src_num)]} {
	    set color $star_color($src_num)
	}

	append reg "ellipse($x $y ${semi_a}i ${semi_b}i $theta) # color=$color width=1 tag={sextract_all} tag={sextract_src.$src_num} select=0 edit=0 move=0 rotate=0 delete=1 highlite=1 callback=highlite CatalogPanelMarkerCB {$src_num} callback=unhighlite CatalogPanelMarkerUnCB {$src_num}\n"
	incr count
	incr batch_count

	if {$batch_count >= $batch_size} {
	    set sextract_all_reg $reg
	    catch {$frame marker catalog command ds9 var sextract_all_reg}
	    set reg "image\n"
	    set batch_count 0
	}
    }

    if {$batch_count > 0} {
	set sextract_all_reg $reg
	catch {$frame marker catalog command ds9 var sextract_all_reg}
    }

    set catpanel(markall,on) 1
}

# ============================================================================
# PSF/Deconv procedures
# ============================================================================

proc CatalogPanelPSFParamLoad {} {
    global catpanel

    set preffile [file join [file normalize ~] .ds9 psf_deconv.prf]
    if {![file exists $preffile]} return
    if {[catch {set fd [open $preffile r]} err]} return
    while {[gets $fd line] >= 0} {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set parts [split $line]
	if {[llength $parts] >= 2} {
	    set key [lindex $parts 0]
	    set val [lindex $parts 1]
	    if {[info exists catpanel(psf,param,$key)]} {
		set catpanel(psf,param,$key) $val
	    }
	}
    }
    close $fd
}

proc CatalogPanelPSFParamSave {} {
    global catpanel

    set prefdir [file join [file normalize ~] .ds9]
    if {![file isdirectory $prefdir]} {
	file mkdir $prefdir
    }
    set preffile [file join $prefdir psf_deconv.prf]
    if {[catch {set fd [open $preffile w]} err]} return
    foreach pname {class-star-thresh max-ellipticity fwhm-sigma min-flux-snr \
		   psf-size rl-iterations wiener-nsr tikhonov-lambda tv-lambda \
		   clean-gain clean-niter clean-threshold mem-lambda mem-niter \
		   ext-core-mag-min ext-core-mag-max ext-wing-mag-max \
		   ext-core-size ext-wing-size ext-blend-inner ext-blend-outer \
		   ext-saturation-limit \
		   sim-telescope sim-instrument sim-filter sim-psf-size \
		   sim-oversample sim-jitter-sigma sim-focus-offset} {
	puts $fd "$pname $catpanel(psf,param,$pname)"
    }
    close $fd
}

proc CatalogPanelPSFGetScript {} {
    set bindir [file dirname [info nameofexecutable]]
    set script [file join $bindir ds9_psf_deconv.py]
    if {![file exists $script]} {
	set libdir [file join [file dirname $bindir] ds9 library]
	set script [file join $libdir ds9_psf_deconv.py]
    }
    return $script
}

proc CatalogPanelPSFGetFITS {} {
    global current

    set fn {}
    if {$current(frame) != {}} {
	catch {set fn [$current(frame) get fits file name full]}
    }
    set fn [string trim $fn "{}"]
    regsub {\[.*\]$} $fn {} fn
    return $fn
}

# --- Star Finding ---

proc CatalogPanelFindStars {method} {
    global catpanel
    global current

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first before finding stars"
	return
    }

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	set catpanel(status) "ERROR: ds9_psf_deconv.py not found"
	return
    }

    # Save catalog to temp TSV
    set catfile [file join [file normalize ~] .ds9 psf_catalog.tsv]
    catch {file mkdir [file dirname $catfile]}
    if {[catch {
	set fd [open $catfile w]
	puts $fd $catpanel(alldata)
	close $fd
    } err]} {
	set catpanel(status) "Star finding error: cannot write catalog: $err"
	return
    }

    # Build arguments
    set paramargs {}
    lappend paramargs "--mode" "find_stars"
    lappend paramargs "--catalog" $catfile
    lappend paramargs "--method" $method
    lappend paramargs "--class-star-thresh" $catpanel(psf,param,class-star-thresh)
    lappend paramargs "--max-ellipticity" $catpanel(psf,param,max-ellipticity)
    lappend paramargs "--fwhm-sigma" $catpanel(psf,param,fwhm-sigma)
    lappend paramargs "--min-flux-snr" $catpanel(psf,param,min-flux-snr)

    set catpanel(status) "Finding stars ($method) ..."
    update idletasks

    set errfile [file join [file normalize ~] .ds9 psf_stderr.txt]
    if {[catch {set data [exec python3 $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	catch {file delete $catfile}
	if {$stderr_msg ne ""} {
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    set catpanel(status) "Star finding error: $last_err"
	    puts "Star finding stderr:\n$stderr_msg"
	} else {
	    set catpanel(status) "Star finding error: $err"
	}
	return
    }
    catch {file delete $errfile}
    catch {file delete $catfile}

    # Parse output
    set lines [split $data \n]
    set star_indices {}
    set n_stars 0

    foreach line $lines {
	if {[string match "#PSF_STARS*" $line]} {
	    foreach field [split $line "\t"] {
		if {[string match "N_STARS=*" $field]} {
		    set n_stars [string range $field 8 end]
		}
	    }
	    continue
	}
	# Skip column header
	if {[string match "NUMBER*" $line]} continue
	if {$line eq {}} continue

	set fields [split $line "\t"]
	if {[llength $fields] >= 3} {
	    set num [lindex $fields 0]
	    lappend star_indices $num
	}
    }

    set catpanel(psf,star_indices) $star_indices
    set catpanel(status) "Found $n_stars stars ($method)"

    # Show star markers
    CatalogPanelShowStars
}

proc CatalogPanelShowStars {} {
    global catpanel
    global current

    if {$current(frame) eq {}} return
    set frame $current(frame)

    # Clear existing star markers
    catch {$frame marker catalog tag psf_star delete}

    if {![info exists catpanel(psf,star_indices)] || $catpanel(psf,star_indices) eq {}} {
	set catpanel(status) "No stars found — run Star Finding first"
	return
    }

    # Parse alldata to find star positions
    set lines [split $catpanel(alldata) \n]
    if {[llength $lines] < 2} return

    set header [lindex $lines 0]
    set cols [split $header "\t"]
    set num_idx -1
    set x_idx -1
    set y_idx -1
    set a_idx -1
    set b_idx -1
    for {set i 0} {$i < [llength $cols]} {incr i} {
	set col [string trim [lindex $cols $i]]
	switch $col {
	    NUMBER {set num_idx $i}
	    X_IMAGE {set x_idx $i}
	    Y_IMAGE {set y_idx $i}
	    A_IMAGE {set a_idx $i}
	    B_IMAGE {set b_idx $i}
	}
    }
    if {$num_idx < 0 || $x_idx < 0 || $y_idx < 0} return

    set reg "image\n"
    set count 0
    foreach line [lrange $lines 1 end] {
	set fields [split $line "\t"]
	if {[llength $fields] <= $num_idx} continue
	set num [string trim [lindex $fields $num_idx]]
	if {[lsearch -exact $catpanel(psf,star_indices) $num] < 0} continue

	set x [string trim [lindex $fields $x_idx]]
	set y [string trim [lindex $fields $y_idx]]
	set a 5.0
	set b 5.0
	if {$a_idx >= 0} {set a [expr {max(3.0, [string trim [lindex $fields $a_idx]] * 2)}]}
	if {$b_idx >= 0} {set b [expr {max(3.0, [string trim [lindex $fields $b_idx]] * 2)}]}

	append reg "circle($x $y ${a}i) # color=purple width=2 dash=1 tag={psf_star} tag={psf_star.$num} select=0 edit=0 move=0 rotate=0 delete=1\n"
	incr count
    }

    if {$count > 0} {
	set psf_star_reg $reg
	catch {$frame marker catalog command ds9 var psf_star_reg}
	set catpanel(status) "Showing $count star markers"
    }
}

proc CatalogPanelClearStars {} {
    global catpanel
    global current

    set catpanel(psf,star_indices) {}
    if {$current(frame) ne {}} {
	catch {$current(frame) marker catalog tag psf_star delete}
    }
    set catpanel(status) "Star markers cleared"
}

# --- PSF Generation ---

proc CatalogPanelBuildPSF {method} {
    global catpanel
    global current

    if {![info exists catpanel(psf,star_indices)] || $catpanel(psf,star_indices) eq {}} {
	set catpanel(status) "Find stars first before building PSF"
	return
    }

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	set catpanel(status) "ERROR: ds9_psf_deconv.py not found"
	return
    }

    # Save catalog
    set catfile [file join [file normalize ~] .ds9 psf_catalog.tsv]
    catch {file mkdir [file dirname $catfile]}
    if {[catch {
	set fd [open $catfile w]
	puts $fd $catpanel(alldata)
	close $fd
    } err]} {
	set catpanel(status) "PSF build error: cannot write catalog: $err"
	return
    }

    set star_list [join $catpanel(psf,star_indices) ","]

    set paramargs {}
    lappend paramargs "--mode" "build_psf"
    lappend paramargs "--catalog" $catfile
    lappend paramargs "--star-indices" $star_list
    lappend paramargs "--psf-method" $method
    lappend paramargs "--psf-size" $catpanel(psf,param,psf-size)
    lappend paramargs "--psf-output" $catpanel(psf,file)

    set catpanel(status) "Building PSF ($method) from [llength $catpanel(psf,star_indices)] stars ..."
    update idletasks

    set errfile [file join [file normalize ~] .ds9 psf_stderr.txt]
    if {[catch {set data [exec python3 $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	catch {file delete $catfile}
	if {$stderr_msg ne ""} {
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    set catpanel(status) "PSF build error: $last_err"
	    puts "PSF build stderr:\n$stderr_msg"
	} else {
	    set catpanel(status) "PSF build error: $err"
	}
	return
    }
    catch {file delete $errfile}
    catch {file delete $catfile}

    set catpanel(psf,has_psf) 1

    # Parse info from output
    foreach line [split $data \n] {
	if {[string match "#PSF_BUILT*" $line]} {
	    set info_str [string range $line 10 end]
	    set catpanel(status) "PSF built: $info_str"
	    break
	}
    }
}

# --- Extended PSF ---

proc CatalogPanelBuildExtendedPSF {} {
    global catpanel
    global ed

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first before building extended PSF"
	return
    }

    set w .extpsfdlg
    if {[winfo exists $w]} {
	raise $w
	return
    }

    toplevel $w
    wm title $w "Extended PSF"
    wm geometry $w 360x420

    # Copy current values
    foreach pname {ext-core-mag-min ext-core-mag-max ext-wing-mag-max \
		   ext-core-size ext-wing-size ext-blend-inner ext-blend-outer \
		   ext-saturation-limit} {
	set ed(psf,$pname) $catpanel(psf,param,$pname)
    }

    set f [ttk::frame $w.content]
    pack $f -fill both -expand true -padx 8 -pady 8

    set r 0

    # Core Stars section
    ttk::label $f.hcore -text "Core Stars" -font TkHeadingFont
    grid $f.hcore -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {4 2}
    incr r

    ttk::label $f.lmagmin -text "Magnitude range:"
    ttk::entry $f.emagmin -textvariable ed(psf,ext-core-mag-min) -width 8
    ttk::label $f.ltilde -text "~"
    ttk::entry $f.emagmax -textvariable ed(psf,ext-core-mag-max) -width 8
    grid $f.lmagmin -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.emagmin -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.ltilde  -row $r -column 2 -padx 2 -pady 2
    grid $f.emagmax -row $r -column 3 -sticky w -padx 2 -pady 2
    incr r

    ttk::label $f.lcsize -text "Cutout size:"
    ttk::entry $f.ecsize -textvariable ed(psf,ext-core-size) -width 8
    ttk::label $f.lcpx -text "px"
    grid $f.lcsize -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ecsize -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lcpx   -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    # Wing Stars section
    ttk::label $f.hwing -text "Wing Stars" -font TkHeadingFont
    grid $f.hwing -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}
    incr r

    ttk::label $f.lwmag -text "Max magnitude:"
    ttk::entry $f.ewmag -textvariable ed(psf,ext-wing-mag-max) -width 8
    grid $f.lwmag -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ewmag -row $r -column 1 -sticky w -padx 2 -pady 2
    incr r

    ttk::label $f.lwsize -text "Cutout size:"
    ttk::entry $f.ewsize -textvariable ed(psf,ext-wing-size) -width 8
    ttk::label $f.lwpx -text "px"
    grid $f.lwsize -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ewsize -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lwpx   -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    ttk::label $f.lsat -text "Saturation:"
    ttk::entry $f.esat -textvariable ed(psf,ext-saturation-limit) -width 8
    ttk::label $f.ladu -text "ADU"
    grid $f.lsat -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.esat -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.ladu -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    # Blending section
    ttk::label $f.hblend -text "Blending" -font TkHeadingFont
    grid $f.hblend -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}
    incr r

    ttk::label $f.lbin -text "Inner radius:"
    ttk::entry $f.ebin -textvariable ed(psf,ext-blend-inner) -width 8
    ttk::label $f.lbinpx -text "px"
    grid $f.lbin   -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ebin   -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lbinpx -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    ttk::label $f.lbout -text "Outer radius:"
    ttk::entry $f.ebout -textvariable ed(psf,ext-blend-outer) -width 8
    ttk::label $f.lbopx -text "px"
    grid $f.lbout  -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ebout  -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lbopx  -row $r -column 2 -sticky w -padx 2 -pady 2

    # Buttons
    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8

    ttk::button $bf.build -text "Build" \
	-command [list CatalogPanelBuildExtendedPSFExec $w]
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.build -side right -padx 4
}

proc CatalogPanelBuildExtendedPSFExec {w} {
    global catpanel
    global ed

    # Apply params
    foreach pname {ext-core-mag-min ext-core-mag-max ext-wing-mag-max \
		   ext-core-size ext-wing-size ext-blend-inner ext-blend-outer \
		   ext-saturation-limit} {
	set catpanel(psf,param,$pname) $ed(psf,$pname)
    }
    CatalogPanelPSFParamSave

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	set catpanel(status) "ERROR: ds9_psf_deconv.py not found"
	return
    }

    # Save catalog
    set catfile [file join [file normalize ~] .ds9 psf_catalog.tsv]
    catch {file mkdir [file dirname $catfile]}
    if {[catch {
	set fd [open $catfile w]
	puts $fd $catpanel(alldata)
	close $fd
    } err]} {
	set catpanel(status) "Extended PSF error: cannot write catalog: $err"
	return
    }

    set paramargs {}
    lappend paramargs "--mode" "build_psf_extended"
    lappend paramargs "--catalog" $catfile
    lappend paramargs "--ext-core-mag-min" $catpanel(psf,param,ext-core-mag-min)
    lappend paramargs "--ext-core-mag-max" $catpanel(psf,param,ext-core-mag-max)
    lappend paramargs "--ext-wing-mag-max" $catpanel(psf,param,ext-wing-mag-max)
    lappend paramargs "--ext-core-size" $catpanel(psf,param,ext-core-size)
    lappend paramargs "--ext-wing-size" $catpanel(psf,param,ext-wing-size)
    lappend paramargs "--ext-blend-inner" $catpanel(psf,param,ext-blend-inner)
    lappend paramargs "--ext-blend-outer" $catpanel(psf,param,ext-blend-outer)
    lappend paramargs "--ext-saturation-limit" $catpanel(psf,param,ext-saturation-limit)
    lappend paramargs "--psf-output" $catpanel(psf,file)

    set catpanel(status) "Building extended PSF ..."
    update idletasks

    set errfile [file join [file normalize ~] .ds9 psf_stderr.txt]
    if {[catch {set data [exec python3 $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	catch {file delete $catfile}
	if {$stderr_msg ne ""} {
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    set catpanel(status) "Extended PSF error: $last_err"
	    puts "Extended PSF stderr:\n$stderr_msg"
	} else {
	    set catpanel(status) "Extended PSF error: $err"
	}
	return
    }
    catch {file delete $errfile}
    catch {file delete $catfile}

    set catpanel(psf,has_psf) 1

    foreach line [split $data \n] {
	if {[string match "#PSF_EXTENDED*" $line]} {
	    set catpanel(status) "Extended PSF: [string range $line 15 end]"
	    break
	}
    }

    CatalogPanelViewPSF
}

# --- Simulation PSF: Check Availability ---

proc CatalogPanelCheckSimAvail {} {
    global catpanel

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	set catpanel(psf,sim_webbpsf_ok) 0
	set catpanel(psf,sim_tinytim_ok) 0
	return
    }

    # Use a dummy fits arg for check_sim mode
    if {[catch {set data [exec python3 $script dummy.fits --mode check_sim 2>/dev/null]} err]} {
	set catpanel(psf,sim_webbpsf_ok) 0
	set catpanel(psf,sim_tinytim_ok) 0
	return
    }

    foreach line [split $data \n] {
	if {[string match "#SIM_STATUS*" $line]} {
	    foreach part [split $line \t] {
		if {[string match "WEBBPSF=*" $part]} {
		    set catpanel(psf,sim_webbpsf_ok) [string range $part 8 end]
		}
		if {[string match "TINYTIM=*" $part]} {
		    set catpanel(psf,sim_tinytim_ok) [string range $part 8 end]
		}
	    }
	}
    }
}

# --- Simulation PSF: WebbPSF (JWST) ---

proc CatalogPanelSimPSFWebbPSF {} {
    global catpanel
    global ed

    set w .webbpsfdlg
    if {[winfo exists $w]} {
	raise $w
	return
    }

    # Check availability if not yet done
    if {$catpanel(psf,sim_webbpsf_ok) == -1} {
	CatalogPanelCheckSimAvail
    }

    toplevel $w
    wm title $w "WebbPSF (JWST)"
    wm geometry $w 360x380

    set ed(psf,sim-instrument) $catpanel(psf,param,sim-instrument)
    set ed(psf,sim-filter) $catpanel(psf,param,sim-filter)
    set ed(psf,sim-psf-size) $catpanel(psf,param,sim-psf-size)
    set ed(psf,sim-oversample) $catpanel(psf,param,sim-oversample)
    set ed(psf,sim-jitter-sigma) $catpanel(psf,param,sim-jitter-sigma)
    set ed(psf,sim-focus-offset) $catpanel(psf,param,sim-focus-offset)

    set f [ttk::frame $w.content]
    pack $f -fill both -expand true -padx 8 -pady 8

    set r 0

    # Availability status
    if {$catpanel(psf,sim_webbpsf_ok) == 1} {
	set statxt "WebbPSF: Available"
    } else {
	set statxt "WebbPSF: Not found (pip install webbpsf)"
    }
    ttk::label $f.status -text $statxt -foreground \
	[expr {$catpanel(psf,sim_webbpsf_ok) == 1 ? "green" : "red"}]
    grid $f.status -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {4 8}
    incr r

    # Instrument
    ttk::label $f.linst -text "Instrument:"
    ttk::combobox $f.cinst -textvariable ed(psf,sim-instrument) -width 14 \
	-values {NIRCAM MIRI NIRISS NIRSPEC FGS} -state readonly
    grid $f.linst -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.cinst -row $r -column 1 -sticky w -padx 2 -pady 2
    incr r

    # Filter
    ttk::label $f.lfilt -text "Filter:"
    ttk::combobox $f.cfilt -textvariable ed(psf,sim-filter) -width 14
    grid $f.lfilt -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.cfilt -row $r -column 1 -sticky w -padx 2 -pady 2
    incr r

    # Dynamic filter update
    bind $f.cinst <<ComboboxSelected>> [list CatalogPanelSimPSFUpdateFilters $f.cfilt jwst]
    # Initialize filter list
    CatalogPanelSimPSFUpdateFilters $f.cfilt jwst

    # PSF size
    ttk::label $f.lsz -text "PSF size:"
    ttk::entry $f.esz -textvariable ed(psf,sim-psf-size) -width 8
    ttk::label $f.lpx -text "px"
    grid $f.lsz -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.esz -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lpx -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    # Oversample
    ttk::label $f.lover -text "Oversample:"
    ttk::entry $f.eover -textvariable ed(psf,sim-oversample) -width 8
    grid $f.lover -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.eover -row $r -column 1 -sticky w -padx 2 -pady 2
    incr r

    # Jitter
    ttk::label $f.ljit -text "Jitter sigma:"
    ttk::entry $f.ejit -textvariable ed(psf,sim-jitter-sigma) -width 8
    ttk::label $f.ljitas -text "arcsec"
    grid $f.ljit   -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ejit   -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.ljitas -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    # Focus
    ttk::label $f.lfoc -text "Focus offset:"
    ttk::entry $f.efoc -textvariable ed(psf,sim-focus-offset) -width 8
    ttk::label $f.lfocw -text "waves"
    grid $f.lfoc  -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.efoc  -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lfocw -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    # Auto-detect button
    ttk::button $f.autodet -text "Auto-detect from FITS" \
	-command [list CatalogPanelSimPSFAutoDetect $f.cinst $f.cfilt jwst]
    grid $f.autodet -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}

    # Buttons
    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8

    ttk::button $bf.gen -text "Generate" \
	-command [list CatalogPanelSimPSFWebbPSFExec $w]
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.gen -side right -padx 4
}

proc CatalogPanelSimPSFWebbPSFExec {w} {
    global catpanel
    global ed

    foreach pname {sim-instrument sim-filter sim-psf-size sim-oversample \
		   sim-jitter-sigma sim-focus-offset} {
	set catpanel(psf,param,$pname) $ed(psf,$pname)
    }
    set catpanel(psf,param,sim-telescope) jwst
    CatalogPanelPSFParamSave

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	set catpanel(status) "ERROR: ds9_psf_deconv.py not found"
	return
    }

    set paramargs {}
    lappend paramargs "--mode" "sim_psf"
    lappend paramargs "--sim-telescope" "jwst"
    lappend paramargs "--sim-instrument" $catpanel(psf,param,sim-instrument)
    lappend paramargs "--sim-filter" $catpanel(psf,param,sim-filter)
    lappend paramargs "--sim-psf-size" $catpanel(psf,param,sim-psf-size)
    lappend paramargs "--sim-oversample" $catpanel(psf,param,sim-oversample)
    lappend paramargs "--sim-jitter-sigma" $catpanel(psf,param,sim-jitter-sigma)
    lappend paramargs "--sim-focus-offset" $catpanel(psf,param,sim-focus-offset)
    lappend paramargs "--psf-output" $catpanel(psf,file)

    set catpanel(status) "Generating WebbPSF ($catpanel(psf,param,sim-instrument) / $catpanel(psf,param,sim-filter)) ..."
    update idletasks

    set errfile [file join [file normalize ~] .ds9 psf_stderr.txt]
    if {[catch {set data [exec python3 $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	if {$stderr_msg ne ""} {
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    set catpanel(status) "WebbPSF error: $last_err"
	    puts "WebbPSF stderr:\n$stderr_msg"
	} else {
	    set catpanel(status) "WebbPSF error: $err"
	}
	return
    }
    catch {file delete $errfile}

    set catpanel(psf,has_psf) 1

    foreach line [split $data \n] {
	if {[string match "#PSF_SIM*" $line]} {
	    set catpanel(status) "Sim PSF: [string range $line 9 end]"
	    break
	}
    }

    CatalogPanelViewPSF
}

# --- Simulation PSF: TinyTim (HST) ---

proc CatalogPanelSimPSFTinyTim {} {
    global catpanel
    global ed

    set w .tinytimdlg
    if {[winfo exists $w]} {
	raise $w
	return
    }

    if {$catpanel(psf,sim_tinytim_ok) == -1} {
	CatalogPanelCheckSimAvail
    }

    toplevel $w
    wm title $w "TinyTim (HST)"
    wm geometry $w 360x340

    set ed(psf,sim-instrument) $catpanel(psf,param,sim-instrument)
    set ed(psf,sim-filter) $catpanel(psf,param,sim-filter)
    set ed(psf,sim-psf-size) $catpanel(psf,param,sim-psf-size)
    set ed(psf,sim-oversample) $catpanel(psf,param,sim-oversample)
    set ed(psf,sim-focus-offset) $catpanel(psf,param,sim-focus-offset)

    set f [ttk::frame $w.content]
    pack $f -fill both -expand true -padx 8 -pady 8

    set r 0

    # Availability status
    if {$catpanel(psf,sim_tinytim_ok) == 1} {
	set statxt "TinyTim: Available"
    } else {
	set statxt "TinyTim: Not found (tiny1/tiny2/tiny3 not on PATH)"
    }
    ttk::label $f.status -text $statxt -foreground \
	[expr {$catpanel(psf,sim_tinytim_ok) == 1 ? "green" : "red"}]
    grid $f.status -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {4 8}
    incr r

    # Instrument
    ttk::label $f.linst -text "Instrument:"
    ttk::combobox $f.cinst -textvariable ed(psf,sim-instrument) -width 14 \
	-values {ACS_WFC ACS_HRC WFC3_UVIS WFC3_IR WFPC2} -state readonly
    grid $f.linst -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.cinst -row $r -column 1 -sticky w -padx 2 -pady 2
    incr r

    # Filter
    ttk::label $f.lfilt -text "Filter:"
    ttk::combobox $f.cfilt -textvariable ed(psf,sim-filter) -width 14
    grid $f.lfilt -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.cfilt -row $r -column 1 -sticky w -padx 2 -pady 2
    incr r

    bind $f.cinst <<ComboboxSelected>> [list CatalogPanelSimPSFUpdateFilters $f.cfilt hst]
    CatalogPanelSimPSFUpdateFilters $f.cfilt hst

    # PSF size
    ttk::label $f.lsz -text "PSF size:"
    ttk::entry $f.esz -textvariable ed(psf,sim-psf-size) -width 8
    ttk::label $f.lpx -text "px"
    grid $f.lsz -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.esz -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lpx -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    # Oversample
    ttk::label $f.lover -text "Oversample:"
    ttk::entry $f.eover -textvariable ed(psf,sim-oversample) -width 8
    grid $f.lover -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.eover -row $r -column 1 -sticky w -padx 2 -pady 2
    incr r

    # Focus
    ttk::label $f.lfoc -text "Focus offset:"
    ttk::entry $f.efoc -textvariable ed(psf,sim-focus-offset) -width 8
    ttk::label $f.lfocum -text "um"
    grid $f.lfoc   -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.efoc   -row $r -column 1 -sticky w -padx 2 -pady 2
    grid $f.lfocum -row $r -column 2 -sticky w -padx 2 -pady 2
    incr r

    # Auto-detect
    ttk::button $f.autodet -text "Auto-detect from FITS" \
	-command [list CatalogPanelSimPSFAutoDetect $f.cinst $f.cfilt hst]
    grid $f.autodet -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}

    # Buttons
    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8

    ttk::button $bf.gen -text "Generate" \
	-command [list CatalogPanelSimPSFTinyTimExec $w]
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.gen -side right -padx 4
}

proc CatalogPanelSimPSFTinyTimExec {w} {
    global catpanel
    global ed

    foreach pname {sim-instrument sim-filter sim-psf-size sim-oversample \
		   sim-focus-offset} {
	set catpanel(psf,param,$pname) $ed(psf,$pname)
    }
    set catpanel(psf,param,sim-telescope) hst
    CatalogPanelPSFParamSave

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	set catpanel(status) "ERROR: ds9_psf_deconv.py not found"
	return
    }

    set paramargs {}
    lappend paramargs "--mode" "sim_psf"
    lappend paramargs "--sim-telescope" "hst"
    lappend paramargs "--sim-instrument" $catpanel(psf,param,sim-instrument)
    lappend paramargs "--sim-filter" $catpanel(psf,param,sim-filter)
    lappend paramargs "--sim-psf-size" $catpanel(psf,param,sim-psf-size)
    lappend paramargs "--sim-oversample" $catpanel(psf,param,sim-oversample)
    lappend paramargs "--sim-focus-offset" $catpanel(psf,param,sim-focus-offset)
    lappend paramargs "--psf-output" $catpanel(psf,file)

    set catpanel(status) "Generating TinyTim PSF ($catpanel(psf,param,sim-instrument) / $catpanel(psf,param,sim-filter)) ..."
    update idletasks

    set errfile [file join [file normalize ~] .ds9 psf_stderr.txt]
    if {[catch {set data [exec python3 $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	if {$stderr_msg ne ""} {
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    set catpanel(status) "TinyTim error: $last_err"
	    puts "TinyTim stderr:\n$stderr_msg"
	} else {
	    set catpanel(status) "TinyTim error: $err"
	}
	return
    }
    catch {file delete $errfile}

    set catpanel(psf,has_psf) 1

    foreach line [split $data \n] {
	if {[string match "#PSF_SIM*" $line]} {
	    set catpanel(status) "Sim PSF: [string range $line 9 end]"
	    break
	}
    }

    CatalogPanelViewPSF
}

# --- Simulation PSF: Shared Helpers ---

proc CatalogPanelSimPSFUpdateFilters {cfilt telescope} {
    global ed

    # Filter lists per instrument
    array set jwst_filters {
	NIRCAM  {F070W F090W F115W F140M F150W F162M F200W F210M F250M F277W F300M F322W2 F335M F356W F360M F410M F430M F444W F460M F480M}
	MIRI    {F560W F770W F1000W F1065C F1130W F1140C F1280W F1500W F1550C F1800W F2100W F2300C F2550W}
	NIRISS  {F090W F115W F140M F150W F158M F200W F277W F356W F380M F430M F444W F480M}
	NIRSPEC {F070LP F100LP F170LP F290LP CLEAR}
	FGS     {FGS}
    }
    array set hst_filters {
	ACS_WFC   {F435W F475W F502N F550M F555W F606W F625W F658N F775W F814W F850LP}
	ACS_HRC   {F220W F250W F330W F435W F475W F555W F606W F625W F775W F814W F850LP}
	WFC3_UVIS {F218W F225W F275W F336W F390W F438W F475W F555W F606W F625W F775W F814W F850LP}
	WFC3_IR   {F098M F105W F110W F125W F140W F160W}
	WFPC2     {F300W F336W F439W F450W F555W F606W F675W F702W F791W F814W}
    }

    set inst $ed(psf,sim-instrument)
    set inst [string toupper $inst]

    if {$telescope eq "jwst" && [info exists jwst_filters($inst)]} {
	$cfilt configure -values $jwst_filters($inst)
    } elseif {$telescope eq "hst" && [info exists hst_filters($inst)]} {
	$cfilt configure -values $hst_filters($inst)
    } else {
	$cfilt configure -values {}
    }
}

proc CatalogPanelSimPSFAutoDetect {cinst cfilt telescope} {
    global catpanel
    global ed

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	set catpanel(status) "ERROR: ds9_psf_deconv.py not found"
	return
    }

    # Run sim_psf with all auto to just get detection output
    set paramargs {}
    lappend paramargs "--mode" "sim_psf"
    lappend paramargs "--sim-telescope" "auto"
    lappend paramargs "--sim-instrument" "auto"
    lappend paramargs "--sim-filter" "auto"
    lappend paramargs "--psf-output" "/dev/null"

    # We expect this to potentially fail (no webbpsf/tinytim), but stderr has detection info
    set errfile [file join [file normalize ~] .ds9 psf_stderr.txt]
    catch {exec python3 $script $fn {*}$paramargs 2>$errfile}

    set stderr_msg ""
    catch {
	set fd [open $errfile r]
	set stderr_msg [read $fd]
	close $fd
    }

    # Parse "Auto-detected: {'telescope': ..., 'instrument': ..., 'filter': ...}"
    if {[regexp {instrument.*?:\s*'([^']+)'} $stderr_msg -> inst_val]} {
	set ed(psf,sim-instrument) [string toupper $inst_val]
	$cinst set [string toupper $inst_val]
    }
    if {[regexp {filter.*?:\s*'([^']+)'} $stderr_msg -> filt_val]} {
	set ed(psf,sim-filter) [string toupper $filt_val]
	$cfilt set [string toupper $filt_val]
    }

    # Update filter list for the detected instrument
    CatalogPanelSimPSFUpdateFilters $cfilt $telescope

    set catpanel(status) "Auto-detected: $ed(psf,sim-instrument) / $ed(psf,sim-filter)"
}

proc CatalogPanelViewPSF {} {
    global catpanel

    if {!$catpanel(psf,has_psf) || ![file exists $catpanel(psf,file)]} {
	set catpanel(status) "No PSF available — build PSF first"
	return
    }

    set w .psfviewer
    if {[winfo exists $w]} {
	raise $w
	CatalogPanelViewPSFRender $w
	return
    }

    toplevel $w
    wm title $w "PSF Viewer"

    # Image display
    ttk::label $w.img -anchor center
    pack $w.img -padx 8 -pady 8 -fill both -expand true

    # Info label
    ttk::label $w.info -text "" -anchor center
    pack $w.info -padx 8 -pady {0 4}

    # Buttons: Save / Load / Close
    ttk::frame $w.btn
    ttk::button $w.btn.save -text "Save PSF..." \
	-command CatalogPanelSavePSF
    ttk::button $w.btn.load -text "Load PSF..." \
	-command [list CatalogPanelViewPSFLoad $w]
    ttk::button $w.btn.close -text "Close" \
	-command [list destroy $w]
    pack $w.btn.save -side left -padx 4
    pack $w.btn.load -side left -padx 4
    pack $w.btn.close -side right -padx 4
    pack $w.btn -fill x -padx 8 -pady 8

    CatalogPanelViewPSFRender $w
}

proc CatalogPanelViewPSFRender {w} {
    global catpanel

    set psffile $catpanel(psf,file)
    set tmpimg [file join [file normalize ~] .ds9 psf_view.ppm]

    # Write render script to temp file
    set tmpscript [file join [file normalize ~] .ds9 psf_render.py]
    set fd [open $tmpscript w]
    puts $fd {import numpy as np
from astropy.io import fits
from scipy.ndimage import zoom
import sys

psffile = sys.argv[1]
outfile = sys.argv[2]

with fits.open(psffile) as hdul:
    data = hdul[0].data.astype(np.float64)

h0, w0 = data.shape
vmin, vmax = float(data.min()), float(data.max())

# Asinh stretch
if vmax > vmin:
    norm = (data - vmin) / (vmax - vmin)
    stretched = np.arcsinh(norm * 10) / np.arcsinh(10)
    img = (stretched * 255).clip(0, 255).astype(np.uint8)
else:
    img = np.zeros_like(data, dtype=np.uint8)

# Scale up to at least 256x256
scale = max(1, 256 // max(img.shape))
if scale > 1:
    img = zoom(img, scale, order=0)

h, w = img.shape

# Write PPM (P6 RGB) — Tk reads this natively
with open(outfile, 'wb') as f:
    f.write(f'P6\n{w} {h}\n255\n'.encode())
    rgb = np.stack([img, img, img], axis=-1)
    f.write(rgb.tobytes())

print(f'{w0}x{h0}  peak={vmax:.4g}')
}
    close $fd

    if {[catch {set info [exec python3 $tmpscript $psffile $tmpimg]} err]} {
	catch {$w.info configure -text "Render error: $err"}
	catch {file delete $tmpscript}
	return
    }
    catch {file delete $tmpscript}

    # Load into Tk photo image
    catch {image delete psfviewimg}
    image create photo psfviewimg -file $tmpimg
    $w.img configure -image psfviewimg
    $w.info configure -text "PSF: $info"

    # Resize window to fit image + buttons
    set iw [image width psfviewimg]
    set ih [image height psfviewimg]
    set ww [expr {max($iw + 16, 280)}]
    set wh [expr {$ih + 90}]
    wm geometry $w ${ww}x${wh}
}

proc CatalogPanelViewPSFLoad {w} {
    global catpanel

    set types {
	{{FITS Files} {.fits .fit .fts}}
	{{All Files} *}
    }
    set infile [tk_getOpenFile -filetypes $types \
		    -title "Load PSF FITS"]
    if {$infile eq {}} return

    file copy -force $infile $catpanel(psf,file)
    set catpanel(psf,has_psf) 1
    set catpanel(status) "PSF loaded from $infile"

    # Refresh the viewer
    CatalogPanelViewPSFRender $w
}

proc CatalogPanelSavePSF {} {
    global catpanel

    if {!$catpanel(psf,has_psf) || ![file exists $catpanel(psf,file)]} {
	set catpanel(status) "No PSF available — build PSF first"
	return
    }

    set types {
	{{FITS Files} {.fits .fit .fts}}
	{{All Files} *}
    }
    set outfile [tk_getSaveFile -filetypes $types \
		     -title "Save PSF FITS" \
		     -initialfile "psf.fits"]
    if {$outfile eq {}} return

    file copy -force $catpanel(psf,file) $outfile
    set catpanel(status) "PSF saved to $outfile"
}

proc CatalogPanelLoadPSF {} {
    global catpanel

    set types {
	{{FITS Files} {.fits .fit .fts}}
	{{All Files} *}
    }
    set infile [tk_getOpenFile -filetypes $types \
		    -title "Load PSF FITS"]
    if {$infile eq {}} return

    file copy -force $infile $catpanel(psf,file)
    set catpanel(psf,has_psf) 1
    set catpanel(status) "PSF loaded from $infile"

    # Refresh viewer if open
    if {[winfo exists .psfviewer]} {
	CatalogPanelViewPSFRender .psfviewer
    }
}

# --- Deconvolution ---

proc CatalogPanelDeconvolve {algorithm} {
    global catpanel
    global current

    if {!$catpanel(psf,has_psf) || ![file exists $catpanel(psf,file)]} {
	set catpanel(status) "No PSF available — build or load PSF first"
	return
    }

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	set catpanel(status) "ERROR: ds9_psf_deconv.py not found"
	return
    }

    set outfile [file join [file normalize ~] .ds9 deconv_result.fits]

    set paramargs {}
    lappend paramargs "--mode" "deconvolve"
    lappend paramargs "--psf" $catpanel(psf,file)
    lappend paramargs "--algorithm" $algorithm
    lappend paramargs "--output" $outfile

    # Algorithm-specific parameters
    switch $algorithm {
	rl - rl_accelerated {
	    lappend paramargs "--iterations" $catpanel(psf,param,rl-iterations)
	}
	rl_tv {
	    lappend paramargs "--iterations" $catpanel(psf,param,rl-iterations)
	    lappend paramargs "--tv-lambda" $catpanel(psf,param,tv-lambda)
	}
	wiener {
	    lappend paramargs "--wiener-nsr" $catpanel(psf,param,wiener-nsr)
	}
	tikhonov {
	    lappend paramargs "--tikhonov-lambda" $catpanel(psf,param,tikhonov-lambda)
	}
	clean {
	    lappend paramargs "--clean-gain" $catpanel(psf,param,clean-gain)
	    lappend paramargs "--clean-niter" $catpanel(psf,param,clean-niter)
	    lappend paramargs "--clean-threshold" $catpanel(psf,param,clean-threshold)
	}
	mem {
	    lappend paramargs "--mem-lambda" $catpanel(psf,param,mem-lambda)
	    lappend paramargs "--mem-niter" $catpanel(psf,param,mem-niter)
	}
    }

    set catpanel(status) "Deconvolving ($algorithm) ..."
    update idletasks

    set errfile [file join [file normalize ~] .ds9 psf_stderr.txt]
    if {[catch {set data [exec python3 $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	if {$stderr_msg ne ""} {
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    set catpanel(status) "Deconvolution error: $last_err"
	    puts "Deconvolution stderr:\n$stderr_msg"
	} else {
	    set catpanel(status) "Deconvolution error: $err"
	}
	return
    }
    catch {file delete $errfile}

    # Load result in a new frame (preserve original)
    if {[file exists $outfile]} {
	CreateFrame
	if {[catch {LoadFitsFile $outfile {} {}} loaderr]} {
	    set catpanel(status) "Deconvolution error: cannot load result: $loaderr"
	    return
	}
	# Apply zscale via DS9's standard scale API
	global scale
	set scale(mode) zscale
	ChangeScaleMode
	set catpanel(status) "Deconvolution complete ($algorithm) — result in new frame"
    } else {
	set catpanel(status) "Deconvolution complete but output file not found"
    }
}

# --- Quick Deconvolve (one-click) ---

proc CatalogPanelQuickDeconvolve {} {
    global catpanel
    global current

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first"
	return
    }

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set catpanel(status) "Quick Deconvolve: Step 1/3 — Finding stars ..."
    update idletasks

    # Step 1: Find stars
    CatalogPanelFindStars combined

    if {![info exists catpanel(psf,star_indices)] || $catpanel(psf,star_indices) eq {}} {
	set catpanel(status) "Quick Deconvolve failed: no stars found"
	return
    }

    set catpanel(status) "Quick Deconvolve: Step 2/3 — Building PSF ..."
    update idletasks

    # Step 2: Build PSF
    CatalogPanelBuildPSF median

    if {!$catpanel(psf,has_psf)} {
	set catpanel(status) "Quick Deconvolve failed: PSF build failed"
	return
    }

    set catpanel(status) "Quick Deconvolve: Step 3/3 — Richardson-Lucy deconvolution ..."
    update idletasks

    # Step 3: Deconvolve
    CatalogPanelDeconvolve rl
}

# --- Settings Dialog ---

proc CatalogPanelStarPSFSettings {} {
    global catpanel
    global ed

    set w .starpsfsettings
    if {[winfo exists $w]} {
	raise $w
	return
    }

    toplevel $w
    wm title $w "Star / PSF Settings"
    wm geometry $w 400x460

    # Copy current values to edit vars
    foreach pname {class-star-thresh max-ellipticity fwhm-sigma min-flux-snr psf-size \
		   ext-core-mag-min ext-core-mag-max ext-wing-mag-max \
		   ext-core-size ext-wing-size ext-blend-inner ext-blend-outer \
		   ext-saturation-limit \
		   sim-telescope sim-instrument sim-filter sim-psf-size \
		   sim-oversample sim-jitter-sigma sim-focus-offset} {
	set ed(psf,$pname) $catpanel(psf,param,$pname)
    }

    ttk::notebook $w.nb
    pack $w.nb -fill both -expand true -padx 8 -pady 8

    # --- Tab 1: Star Finding ---
    set t1 [ttk::frame $w.nb.stars]
    $w.nb add $t1 -text "Star Finding"

    set r 0
    ttk::label $t1.lcs -text "CLASS_STAR threshold:"
    ttk::entry $t1.ecs -textvariable ed(psf,class-star-thresh) -width 10
    grid $t1.lcs -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.ecs -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lell -text "Max ellipticity:"
    ttk::entry $t1.eell -textvariable ed(psf,max-ellipticity) -width 10
    grid $t1.lell -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.eell -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lfwhm -text "FWHM sigma:"
    ttk::entry $t1.efwhm -textvariable ed(psf,fwhm-sigma) -width 10
    grid $t1.lfwhm -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.efwhm -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lsnr -text "Min flux S/N:"
    ttk::entry $t1.esnr -textvariable ed(psf,min-flux-snr) -width 10
    grid $t1.lsnr -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.esnr -row $r -column 1 -sticky w -padx 4 -pady 4

    # --- Tab 2: PSF ---
    set t2 [ttk::frame $w.nb.psf]
    $w.nb add $t2 -text "PSF"

    set r 0
    ttk::label $t2.lsz -text "PSF size (pixels):"
    ttk::entry $t2.esz -textvariable ed(psf,psf-size) -width 10
    grid $t2.lsz -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.esz -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    # Extended PSF sub-section
    ttk::label $t2.hext -text "Extended PSF:" -font TkHeadingFont
    grid $t2.hext -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}
    incr r

    ttk::label $t2.lcmin -text "Core mag min:"
    ttk::entry $t2.ecmin -textvariable ed(psf,ext-core-mag-min) -width 10
    grid $t2.lcmin -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $t2.ecmin -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $t2.lcmax -text "Core mag max:"
    ttk::entry $t2.ecmax -textvariable ed(psf,ext-core-mag-max) -width 10
    grid $t2.lcmax -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $t2.ecmax -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $t2.lwmag -text "Wing mag max:"
    ttk::entry $t2.ewmag -textvariable ed(psf,ext-wing-mag-max) -width 10
    grid $t2.lwmag -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $t2.ewmag -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $t2.lcsz -text "Core cutout size:"
    ttk::entry $t2.ecsz -textvariable ed(psf,ext-core-size) -width 10
    grid $t2.lcsz -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $t2.ecsz -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $t2.lwsz -text "Wing cutout size:"
    ttk::entry $t2.ewsz -textvariable ed(psf,ext-wing-size) -width 10
    grid $t2.lwsz -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $t2.ewsz -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $t2.lbi -text "Blend inner (px):"
    ttk::entry $t2.ebi -textvariable ed(psf,ext-blend-inner) -width 10
    grid $t2.lbi -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $t2.ebi -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $t2.lbo -text "Blend outer (px):"
    ttk::entry $t2.ebo -textvariable ed(psf,ext-blend-outer) -width 10
    grid $t2.lbo -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $t2.ebo -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $t2.lsat -text "Saturation (ADU):"
    ttk::entry $t2.esat -textvariable ed(psf,ext-saturation-limit) -width 10
    grid $t2.lsat -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $t2.esat -row $r -column 1 -sticky w -padx 4 -pady 2

    # --- Tab 3: Simulation ---
    set t3 [ttk::frame $w.nb.sim]
    $w.nb add $t3 -text "Simulation"

    set r 0
    ttk::label $t3.ltel -text "Telescope:"
    ttk::combobox $t3.etel -textvariable ed(psf,sim-telescope) -width 10 \
	-values {auto jwst hst}
    grid $t3.ltel -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.etel -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.linst -text "Instrument:"
    ttk::entry $t3.einst -textvariable ed(psf,sim-instrument) -width 14
    grid $t3.linst -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.einst -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lfilt -text "Filter:"
    ttk::entry $t3.efilt -textvariable ed(psf,sim-filter) -width 14
    grid $t3.lfilt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.efilt -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lsz -text "PSF size (px):"
    ttk::entry $t3.esz -textvariable ed(psf,sim-psf-size) -width 10
    grid $t3.lsz -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.esz -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lover -text "Oversample:"
    ttk::entry $t3.eover -textvariable ed(psf,sim-oversample) -width 10
    grid $t3.lover -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.eover -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.ljit -text "Jitter sigma (arcsec):"
    ttk::entry $t3.ejit -textvariable ed(psf,sim-jitter-sigma) -width 10
    grid $t3.ljit -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.ejit -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lfoc -text "Focus offset:"
    ttk::entry $t3.efoc -textvariable ed(psf,sim-focus-offset) -width 10
    grid $t3.lfoc -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.efoc -row $r -column 1 -sticky w -padx 4 -pady 4

    # Buttons
    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8

    ttk::button $bf.apply -text "Apply" -command [list CatalogPanelStarPSFSettingsApply $w]
    ttk::button $bf.defaults -text "Defaults" -command CatalogPanelStarPSFSettingsDefaults
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.defaults -side right -padx 4
    pack $bf.apply -side right -padx 4
}

proc CatalogPanelStarPSFSettingsApply {w} {
    global catpanel
    global ed

    foreach pname {class-star-thresh max-ellipticity fwhm-sigma min-flux-snr psf-size \
		   ext-core-mag-min ext-core-mag-max ext-wing-mag-max \
		   ext-core-size ext-wing-size ext-blend-inner ext-blend-outer \
		   ext-saturation-limit \
		   sim-telescope sim-instrument sim-filter sim-psf-size \
		   sim-oversample sim-jitter-sigma sim-focus-offset} {
	set catpanel(psf,param,$pname) $ed(psf,$pname)
    }
    CatalogPanelPSFParamSave
    set catpanel(status) "Star/PSF settings applied and saved"
}

proc CatalogPanelStarPSFSettingsDefaults {} {
    global ed

    set ed(psf,class-star-thresh) 0.8
    set ed(psf,max-ellipticity) 0.2
    set ed(psf,fwhm-sigma) 2.0
    set ed(psf,min-flux-snr) 10.0
    set ed(psf,psf-size) 51

    set ed(psf,ext-core-mag-min) 18.0
    set ed(psf,ext-core-mag-max) 22.0
    set ed(psf,ext-wing-mag-max) 16.0
    set ed(psf,ext-core-size) 51
    set ed(psf,ext-wing-size) 201
    set ed(psf,ext-blend-inner) 20.0
    set ed(psf,ext-blend-outer) 30.0
    set ed(psf,ext-saturation-limit) 60000.0

    set ed(psf,sim-telescope) auto
    set ed(psf,sim-instrument) auto
    set ed(psf,sim-filter) auto
    set ed(psf,sim-psf-size) 201
    set ed(psf,sim-oversample) 1
    set ed(psf,sim-jitter-sigma) 0.007
    set ed(psf,sim-focus-offset) 0.0
}

proc CatalogPanelDeconvSettings {} {
    global catpanel
    global ed

    set w .deconvsettings
    if {[winfo exists $w]} {
	raise $w
	return
    }

    toplevel $w
    wm title $w "Deconvolution Settings"
    wm geometry $w 380x480

    # Copy current values to edit vars
    foreach pname {rl-iterations wiener-nsr tikhonov-lambda tv-lambda \
		   clean-gain clean-niter clean-threshold mem-lambda mem-niter} {
	set ed(psf,$pname) $catpanel(psf,param,$pname)
    }

    set f [ttk::frame $w.content]
    pack $f -fill both -expand true -padx 8 -pady 8

    set r 0
    ttk::label $f.h1 -text "Richardson-Lucy:" -font TkHeadingFont
    grid $f.h1 -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}
    incr r

    ttk::label $f.liter -text "  Iterations:"
    ttk::entry $f.eiter -textvariable ed(psf,rl-iterations) -width 10
    grid $f.liter -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.eiter -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $f.ltv -text "  TV lambda:"
    ttk::entry $f.etv -textvariable ed(psf,tv-lambda) -width 10
    grid $f.ltv -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.etv -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $f.h2 -text "Wiener:" -font TkHeadingFont
    grid $f.h2 -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}
    incr r

    ttk::label $f.lnsr -text "  NSR:"
    ttk::entry $f.ensr -textvariable ed(psf,wiener-nsr) -width 10
    grid $f.lnsr -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ensr -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $f.h3 -text "Tikhonov:" -font TkHeadingFont
    grid $f.h3 -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}
    incr r

    ttk::label $f.ltik -text "  Lambda:"
    ttk::entry $f.etik -textvariable ed(psf,tikhonov-lambda) -width 10
    grid $f.ltik -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.etik -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $f.h4 -text "CLEAN:" -font TkHeadingFont
    grid $f.h4 -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}
    incr r

    ttk::label $f.lcg -text "  Gain:"
    ttk::entry $f.ecg -textvariable ed(psf,clean-gain) -width 10
    grid $f.lcg -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ecg -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $f.lcn -text "  Iterations:"
    ttk::entry $f.ecn -textvariable ed(psf,clean-niter) -width 10
    grid $f.lcn -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ecn -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $f.lct -text "  Threshold:"
    ttk::entry $f.ect -textvariable ed(psf,clean-threshold) -width 10
    grid $f.lct -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ect -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $f.h5 -text "MEM:" -font TkHeadingFont
    grid $f.h5 -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}
    incr r

    ttk::label $f.lml -text "  Lambda:"
    ttk::entry $f.eml -textvariable ed(psf,mem-lambda) -width 10
    grid $f.lml -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.eml -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $f.lmn -text "  Iterations:"
    ttk::entry $f.emn -textvariable ed(psf,mem-niter) -width 10
    grid $f.lmn -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.emn -row $r -column 1 -sticky w -padx 4 -pady 2

    # Buttons
    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8

    ttk::button $bf.apply -text "Apply" -command [list CatalogPanelDeconvSettingsApply $w]
    ttk::button $bf.defaults -text "Defaults" -command CatalogPanelDeconvSettingsDefaults
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.defaults -side right -padx 4
    pack $bf.apply -side right -padx 4
}

proc CatalogPanelDeconvSettingsApply {w} {
    global catpanel
    global ed

    foreach pname {rl-iterations wiener-nsr tikhonov-lambda tv-lambda \
		   clean-gain clean-niter clean-threshold mem-lambda mem-niter} {
	set catpanel(psf,param,$pname) $ed(psf,$pname)
    }
    CatalogPanelPSFParamSave
    set catpanel(status) "Deconvolution settings applied and saved"
}

proc CatalogPanelDeconvSettingsDefaults {} {
    global ed

    set ed(psf,rl-iterations) 30
    set ed(psf,wiener-nsr) 0.01
    set ed(psf,tikhonov-lambda) 0.001
    set ed(psf,tv-lambda) 0.001
    set ed(psf,clean-gain) 0.1
    set ed(psf,clean-niter) 1000
    set ed(psf,clean-threshold) 0.0
    set ed(psf,mem-lambda) 0.1
    set ed(psf,mem-niter) 100
}

# Keep backward compatibility — old proc name routes to star/PSF settings
proc CatalogPanelPSFDeconvSettings {} {
    CatalogPanelStarPSFSettings
}
proc CatalogPanelPSFSettingsApply {w} {
    CatalogPanelStarPSFSettingsApply $w
}
proc CatalogPanelPSFSettingsDefaults {} {
    CatalogPanelStarPSFSettingsDefaults
}

# ============================================================================
# Separate — Source Deblending / Splitting
# ============================================================================

proc CatalogPanelSeparateGetScript {} {
    set bindir [file dirname [info nameofexecutable]]
    set script [file join $bindir ds9_separate.py]
    if {![file exists $script]} {
	set libdir [file join [file dirname $bindir] ds9 library]
	set script [file join $libdir ds9_separate.py]
    }
    return $script
}

# Get the currently selected source NUMBER from the table
proc CatalogPanelGetSelectedSource {} {
    global catpanel

    if {![info exists catpanel(tbl)]} { return {} }
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} { return {} }

    # Get selected row
    set sel [$catpanel(tbl) curselection]
    if {$sel eq {}} { return {} }
    # sel is "row,col" — extract row
    set row [lindex [split [lindex $sel 0] ","] 0]
    if {$row <= 0} { return {} }

    global $catpanel(tbldb)
    set ncols [$catpanel(tbl) cget -cols]

    # Find column indices
    set col_num -1; set col_x -1; set col_y -1
    set col_a -1; set col_b -1; set col_theta -1; set col_ir -1
    for {set c 1} {$c <= $ncols} {incr c} {
	if {[info exists ${catpanel(tbldb)}(0,$c)]} {
	    switch -- [set ${catpanel(tbldb)}(0,$c)] {
		NUMBER      { set col_num $c }
		X_IMAGE     { set col_x $c }
		Y_IMAGE     { set col_y $c }
		A_IMAGE     { set col_a $c }
		B_IMAGE     { set col_b $c }
		THETA_IMAGE { set col_theta $c }
		ISO_RADIUS  { set col_ir $c }
	    }
	}
    }
    if {$col_x < 0 || $col_y < 0} { return {} }

    set result [dict create]
    dict set result row $row

    if {$col_num >= 0 && [info exists ${catpanel(tbldb)}($row,$col_num)]} {
	dict set result number [set ${catpanel(tbldb)}($row,$col_num)]
    } else {
	dict set result number $row
    }
    dict set result x [set ${catpanel(tbldb)}($row,$col_x)]
    dict set result y [set ${catpanel(tbldb)}($row,$col_y)]

    set a 10.0; set b 10.0; set theta 0.0; set ir 10.0
    if {$col_a >= 0 && [info exists ${catpanel(tbldb)}($row,$col_a)]} {
	set val [set ${catpanel(tbldb)}($row,$col_a)]
	if {[string is double -strict $val] && $val > 0} { set a $val }
    }
    if {$col_b >= 0 && [info exists ${catpanel(tbldb)}($row,$col_b)]} {
	set val [set ${catpanel(tbldb)}($row,$col_b)]
	if {[string is double -strict $val] && $val > 0} { set b $val }
    }
    if {$col_theta >= 0 && [info exists ${catpanel(tbldb)}($row,$col_theta)]} {
	set val [set ${catpanel(tbldb)}($row,$col_theta)]
	if {[string is double -strict $val]} { set theta $val }
    }
    if {$col_ir >= 0 && [info exists ${catpanel(tbldb)}($row,$col_ir)]} {
	set val [set ${catpanel(tbldb)}($row,$col_ir)]
	if {[string is double -strict $val] && $val > 0} { set ir $val }
    }
    dict set result a $a
    dict set result b $b
    dict set result theta $theta
    dict set result iso_radius $ir

    return $result
}

proc CatalogPanelDeleteSelected {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "No catalog data"
	return
    }

    # Get selected source info
    set src [CatalogPanelGetSelectedSource]
    if {$src eq {}} {
	set catpanel(status) "No source selected — click a source first"
	return
    }

    set src_num [dict get $src number]

    # Remove the row from alldata
    set lines [split $catpanel(alldata) \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]

    # Find NUMBER column index
    set num_idx -1
    for {set i 0} {$i < [llength $headers]} {incr i} {
	if {[string trim [lindex $headers $i]] eq "NUMBER"} {
	    set num_idx $i
	    break
	}
    }

    set new_lines [list $header]
    set deleted 0
    foreach line [lrange $lines 1 end] {
	if {$line eq {}} continue
	set fields [split $line "\t"]

	set this_num ""
	if {$num_idx >= 0 && [llength $fields] > $num_idx} {
	    set this_num [string trim [lindex $fields $num_idx]]
	}

	if {$this_num eq [string trim $src_num]} {
	    set deleted 1
	} else {
	    lappend new_lines $line
	}
    }

    if {!$deleted} {
	set catpanel(status) "Source $src_num not found in catalog"
	return
    }

    set catpanel(alldata) [join $new_lines \n]

    # Delete the source marker
    global current
    if {$current(frame) ne {}} {
	catch {$current(frame) marker catalog sextract_src.$src_num delete}
	catch {$current(frame) marker catalog sextract_sel delete}
    }

    # Reload table and refresh markers
    CatalogPanelLoadTSV $catpanel(alldata) "deleted"
    if {[info exists catpanel(markall,on)] && $catpanel(markall,on)} {
	CatalogPanelCreateAllMarkers
    }

    set catpanel(status) "Deleted source $src_num"
}

proc CatalogPanelSeparateSelected {} {
    global catpanel
    global current

    if {$current(frame) eq {}} return

    set src [CatalogPanelGetSelectedSource]
    if {$src eq {}} {
	set catpanel(status) "Select a source first (click on an ellipse)"
	return
    }

    # Get FITS filename
    set fn {}
    if {[$current(frame) has fits]} {
	set fn [$current(frame) get fits file name full]
    }
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelSeparateGetScript]
    if {![file exists $script]} {
	set catpanel(status) "ERROR: ds9_separate.py not found"
	return
    }

    set src_num [dict get $src number]
    set src_x [dict get $src x]
    set src_y [dict get $src y]
    set src_a [dict get $src a]
    set src_b [dict get $src b]
    set src_theta [dict get $src theta]
    set src_ir [dict get $src iso_radius]

    set catpanel(status) "Separating source $src_num ..."
    update idletasks

    set paramargs {}
    lappend paramargs "--x" $src_x
    lappend paramargs "--y" $src_y
    lappend paramargs "--a" $src_a
    lappend paramargs "--b" $src_b
    lappend paramargs "--theta" $src_theta
    lappend paramargs "--iso-radius" $src_ir
    lappend paramargs "--parent-number" $src_num
    lappend paramargs "--deblend-nthresh" $catpanel(param,sep-deblend-nthresh)
    lappend paramargs "--deblend-mincont" $catpanel(param,sep-deblend-mincont)
    lappend paramargs "--detect-thresh" $catpanel(param,sep-detect-thresh)
    lappend paramargs "--detect-minarea" $catpanel(param,sep-detect-minarea)
    lappend paramargs "--radius-factor" $catpanel(param,sep-radius-factor)
    lappend paramargs "--back-size" $catpanel(param,sep-back-size)

    set errfile [file join [file normalize ~] .ds9 separate_stderr.txt]
    if {[catch {set data [exec python3 $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	if {$stderr_msg ne ""} {
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    set catpanel(status) "Separate error: $last_err"
	} else {
	    set catpanel(status) "Separate error: $err"
	}
	return
    }
    catch {file delete $errfile}

    # Parse output
    set lines [split $data \n]
    set n_sub 0
    set sub_lines {}

    foreach line $lines {
	if {[string match "#SEPARATE*" $line]} {
	    foreach field [split $line "\t"] {
		if {[string match "N_SUB=*" $field]} {
		    set n_sub [string range $field 6 end]
		}
	    }
	    continue
	}
	if {[string match "NUMBER*" $line]} continue
	if {$line eq {}} continue
	lappend sub_lines $line
    }

    if {$n_sub == 0 || [llength $sub_lines] < 2} {
	set catpanel(status) "No sub-components found for source $src_num"
	return
    }

    # Replace the parent source in alldata with sub-sources
    CatalogPanelSeparateReplace $src_num $sub_lines

    set catpanel(status) "Source $src_num separated into $n_sub sub-components"
}

proc CatalogPanelSeparateReplace {parent_num sub_lines} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set all_lines [split $catpanel(alldata) \n]
    set header [lindex $all_lines 0]
    set headers [split $header "\t"]

    # Find NUMBER column index
    set num_idx -1
    for {set i 0} {$i < [llength $headers]} {incr i} {
	if {[string trim [lindex $headers $i]] eq "NUMBER"} {
	    set num_idx $i
	    break
	}
    }

    # Sub-source output has fixed columns; map to parent catalog columns
    set sub_cols {NUMBER X_IMAGE Y_IMAGE A_IMAGE B_IMAGE THETA_IMAGE ISO_RADIUS FLUX_AUTO MAG_AUTO FLAGS}

    # Build column index map: sub_col_name → parent_col_index
    set col_map {}
    for {set si 0} {$si < [llength $sub_cols]} {incr si} {
	set scol [lindex $sub_cols $si]
	for {set pi 0} {$pi < [llength $headers]} {incr pi} {
	    if {[string trim [lindex $headers $pi]] eq $scol} {
		lappend col_map [list $si $pi]
		break
	    }
	}
    }

    # Build new alldata lines
    set new_lines [list $header]
    set ncols [llength $headers]

    foreach line [lrange $all_lines 1 end] {
	if {$line eq {}} continue
	set fields [split $line "\t"]

	set this_num ""
	if {$num_idx >= 0 && [llength $fields] > $num_idx} {
	    set this_num [string trim [lindex $fields $num_idx]]
	}

	if {$this_num eq [string trim $parent_num]} {
	    # Replace parent with sub-sources
	    foreach sub_line $sub_lines {
		set sub_fields [split $sub_line "\t"]
		# Start with empty row matching parent column count
		set new_row {}
		for {set c 0} {$c < $ncols} {incr c} {
		    lappend new_row ""
		}
		# Fill in mapped columns from sub-source
		foreach mapping $col_map {
		    set si [lindex $mapping 0]
		    set pi [lindex $mapping 1]
		    if {$si < [llength $sub_fields]} {
			lset new_row $pi [lindex $sub_fields $si]
		    }
		}
		lappend new_lines [join $new_row "\t"]
	    }
	} else {
	    lappend new_lines $line
	}
    }

    set catpanel(alldata) [join $new_lines \n]

    # Reload table and markers
    CatalogPanelLoadTSV $catpanel(alldata) "separated"
    if {[info exists catpanel(markall,on)] && $catpanel(markall,on)} {
	CatalogPanelCreateAllMarkers
    }
}

# --- Separate Settings Dialog ---

proc CatalogPanelSeparateSettings {} {
    global catpanel
    global ed

    set w .separatesettings
    if {[winfo exists $w]} {
	raise $w
	return
    }

    toplevel $w
    wm title $w "Separate Settings"
    wm geometry $w 380x300

    foreach pname {sep-deblend-nthresh sep-deblend-mincont sep-detect-thresh \
		   sep-detect-minarea sep-radius-factor sep-back-size} {
	set ed($pname) $catpanel(param,$pname)
    }

    set f [ttk::frame $w.content]
    pack $f -fill both -expand true -padx 8 -pady 8

    set r 0
    ttk::label $f.h1 -text "Deblending:" -font TkHeadingFont
    grid $f.h1 -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}
    incr r

    ttk::label $f.lnt -text "  Sub-thresholds:"
    ttk::entry $f.ent -textvariable ed(sep-deblend-nthresh) -width 10
    grid $f.lnt -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ent -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $f.lmc -text "  Min contrast:"
    ttk::entry $f.emc -textvariable ed(sep-deblend-mincont) -width 10
    grid $f.lmc -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.emc -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $f.h2 -text "Detection:" -font TkHeadingFont
    grid $f.h2 -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}
    incr r

    ttk::label $f.ldt -text "  Threshold (sigma):"
    ttk::entry $f.edt -textvariable ed(sep-detect-thresh) -width 10
    grid $f.ldt -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.edt -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $f.lma -text "  Min area (pixels):"
    ttk::entry $f.ema -textvariable ed(sep-detect-minarea) -width 10
    grid $f.lma -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ema -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $f.h3 -text "Cutout:" -font TkHeadingFont
    grid $f.h3 -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady {8 2}
    incr r

    ttk::label $f.lrf -text "  Radius factor:"
    ttk::entry $f.erf -textvariable ed(sep-radius-factor) -width 10
    grid $f.lrf -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.erf -row $r -column 1 -sticky w -padx 4 -pady 2
    incr r

    ttk::label $f.lbs -text "  Background size:"
    ttk::entry $f.ebs -textvariable ed(sep-back-size) -width 10
    grid $f.lbs -row $r -column 0 -sticky w -padx 8 -pady 2
    grid $f.ebs -row $r -column 1 -sticky w -padx 4 -pady 2

    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8

    ttk::button $bf.apply -text "Apply" -command [list CatalogPanelSeparateSettingsApply $w]
    ttk::button $bf.defaults -text "Defaults" -command CatalogPanelSeparateSettingsDefaults
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.defaults -side right -padx 4
    pack $bf.apply -side right -padx 4
}

proc CatalogPanelSeparateSettingsApply {w} {
    global catpanel
    global ed

    foreach pname {sep-deblend-nthresh sep-deblend-mincont sep-detect-thresh \
		   sep-detect-minarea sep-radius-factor sep-back-size} {
	set catpanel(param,$pname) $ed($pname)
    }
    CatalogPanelParamSave
    set catpanel(status) "Separate settings applied and saved"
}

proc CatalogPanelSeparateSettingsDefaults {} {
    global ed

    set ed(sep-deblend-nthresh) 64
    set ed(sep-deblend-mincont) 0.0001
    set ed(sep-detect-thresh) 0.8
    set ed(sep-detect-minarea) 3
    set ed(sep-radius-factor) 3.0
    set ed(sep-back-size) 32
}

# --- Separate Save/Load Catalog ---

proc CatalogPanelSeparateSave {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "No catalog to save"
	return
    }

    set types {
	{{TSV Files} {.tsv .txt}}
	{{All Files} *}
    }
    set outfile [tk_getSaveFile -filetypes $types \
		     -title "Save Catalog" \
		     -initialfile "catalog_separated.tsv"]
    if {$outfile eq {}} return

    if {[catch {
	set fd [open $outfile w]
	puts -nonewline $fd $catpanel(alldata)
	close $fd
    } err]} {
	set catpanel(status) "Save error: $err"
	return
    }
    set catpanel(status) "Catalog saved to $outfile"
}

proc CatalogPanelSeparateLoad {} {
    global catpanel

    set types {
	{{TSV Files} {.tsv .txt}}
	{{All Files} *}
    }
    set infile [tk_getOpenFile -filetypes $types \
		    -title "Load Catalog"]
    if {$infile eq {}} return

    if {[catch {
	set fd [open $infile r]
	set data [read $fd]
	close $fd
    } err]} {
	set catpanel(status) "Load error: $err"
	return
    }

    set catpanel(alldata) [string trim $data]
    CatalogPanelLoadTSV $catpanel(alldata) "loaded"
    if {[info exists catpanel(markall,on)] && $catpanel(markall,on)} {
	CatalogPanelCreateAllMarkers
    }
    set catpanel(status) "Catalog loaded from $infile"
}

# ============================================================================
# Add Objects — Detect and add source at cursor position (Ctrl+A)
# ============================================================================

proc CatalogPanelAddObjectGetScript {} {
    set bindir [file dirname [info nameofexecutable]]
    set script [file join $bindir ds9_add_source.py]
    if {![file exists $script]} {
	set libdir [file join [file dirname $bindir] ds9 library]
	set script [file join $libdir ds9_add_source.py]
    }
    return $script
}

proc CatalogPanelAddObjectAtPosition {which imgx imgy} {
    global catpanel
    global current

    # Check if add objects mode is enabled
    if {![info exists catpanel(add_objects_mode)] || !$catpanel(add_objects_mode)} {
	set catpanel(status) "Enable Add Objects in Display menu first"
	return
    }

    if {$which ne $current(frame) || $which eq {}} return
    if {![$which has fits]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first, then enable Add Objects"
	return
    }

    # Get FITS filename
    set fn [$which get fits file name full]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelAddObjectGetScript]
    if {![file exists $script]} {
	set catpanel(status) "ERROR: ds9_add_source.py not found"
	return
    }

    # Find max NUMBER in current catalog
    set max_num 0
    set lines [split $catpanel(alldata) \n]
    set header [lindex $lines 0]
    set headers [split $header "\t"]
    set num_idx -1
    for {set i 0} {$i < [llength $headers]} {incr i} {
	if {[string trim [lindex $headers $i]] eq "NUMBER"} {
	    set num_idx $i
	    break
	}
    }
    if {$num_idx >= 0} {
	foreach line [lrange $lines 1 end] {
	    if {$line eq {}} continue
	    set fields [split $line "\t"]
	    set nv [string trim [lindex $fields $num_idx]]
	    # Handle fractional numbers like "42.1" from separate
	    set int_part [lindex [split $nv "."] 0]
	    if {[string is integer -strict $int_part] && $int_part > $max_num} {
		set max_num $int_part
	    }
	}
    }
    set new_num [expr {$max_num + 1}]

    set catpanel(status) "Detecting source at ($imgx, $imgy) ..."
    update idletasks

    set paramargs {}
    lappend paramargs "--x" $imgx
    lappend paramargs "--y" $imgy
    lappend paramargs "--number" $new_num
    lappend paramargs "--detect-thresh" $catpanel(param,detect-thresh)
    lappend paramargs "--detect-minarea" $catpanel(param,detect-minarea)
    lappend paramargs "--deblend-nthresh" $catpanel(param,deblend-nthresh)
    lappend paramargs "--deblend-mincont" $catpanel(param,deblend-mincont)
    lappend paramargs "--back-size" $catpanel(param,back-size)
    lappend paramargs "--mag-zeropoint" $catpanel(param,mag-zeropoint)
    lappend paramargs "--phot-aperture" $catpanel(param,phot-aperture)

    set errfile [file join [file normalize ~] .ds9 add_source_stderr.txt]
    if {[catch {set data [exec python3 $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	if {$stderr_msg ne ""} {
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    set catpanel(status) "Add source error: $last_err"
	} else {
	    set catpanel(status) "Add source error: $err"
	}
	return
    }
    catch {file delete $errfile}

    # Parse output
    set result_lines [split $data \n]
    set found 0
    set source_line {}

    foreach line $result_lines {
	if {[string match "#ADD_SOURCE*" $line]} {
	    foreach field [split $line "\t"] {
		if {[string match "FOUND=*" $field]} {
		    set found [string range $field 6 end]
		}
	    }
	    continue
	}
	# Skip column header
	if {[string match "NUMBER*" $line]} continue
	if {$line eq {}} continue
	set source_line $line
    }

    if {$found == 0 || $source_line eq {}} {
	set catpanel(status) "No source detected at ($imgx, $imgy)"
	return
    }

    # Add the new source to alldata
    # Map output columns to existing catalog columns
    set out_cols {NUMBER X_IMAGE Y_IMAGE A_IMAGE B_IMAGE THETA_IMAGE ELLIPTICITY KRON_RADIUS ISO_RADIUS FLUX_AUTO FLUX_APER MAG_AUTO MAG_APER PEAK CLASS_STAR FLAGS FWHM_IMAGE}
    set src_fields [split $source_line "\t"]

    set ncols [llength $headers]
    set new_row {}
    for {set c 0} {$c < $ncols} {incr c} {
	lappend new_row ""
    }

    # Map each output column to the catalog column
    for {set si 0} {$si < [llength $out_cols]} {incr si} {
	set scol [lindex $out_cols $si]
	for {set pi 0} {$pi < $ncols} {incr pi} {
	    if {[string trim [lindex $headers $pi]] eq $scol} {
		if {$si < [llength $src_fields]} {
		    lset new_row $pi [lindex $src_fields $si]
		}
		break
	    }
	}
    }

    # Append to alldata
    append catpanel(alldata) "\n" [join $new_row "\t"]

    # Reload table and markers
    CatalogPanelLoadTSV $catpanel(alldata) "added"
    if {[info exists catpanel(markall,on)] && $catpanel(markall,on)} {
	CatalogPanelCreateAllMarkers
    }

    set catpanel(status) "Added source #$new_num at ($imgx, $imgy)"
}

# ============================================================================
# Analysis — Advanced Source Extraction Features (Phase A-D)
# ============================================================================

# --- Helper: Locate a Python script (check bindir, then library dir) ---

proc CatalogPanelGetScript {scriptname} {
    set bindir [file dirname [info nameofexecutable]]
    set script [file join $bindir $scriptname]
    if {![file exists $script]} {
	set libdir [file join [file dirname $bindir] ds9 library]
	set script [file join $libdir $scriptname]
    }
    return $script
}

# --- Helper: Get FITS filename from current frame ---

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

# --- Helper: Save current catalog to temp TSV ---

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

# --- Helper: Add columns from result TSV to alldata ---

proc CatalogPanelAddColumnsFromTSV {result_data col_names} {
    global catpanel

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

# --- A3: Export Regions (.reg) ---

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

# --- B9: Export FITS Table ---

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

    if {[catch {
	set data [exec python3 $script --input $tmpfile --output $fn 2>@stderr]
    } err]} {
	set catpanel(status) "FITS export error: $err"
	return
    }

    set catpanel(status) "Exported FITS table to [file tail $fn]"
}

# --- B8: Segmentation Map ---

proc CatalogPanelSegmentationMap {} {
    global catpanel current ds9

    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_segmap.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_segmap.py"
	return
    }

    set catpanel(status) "Generating segmentation map..."
    update idletasks

    set args [list python3 $script $fn]
    if {[info exists catpanel(param,detect-thresh)]} {
	lappend args --detect-thresh $catpanel(param,detect-thresh)
    }
    if {[info exists catpanel(param,detect-minarea)]} {
	lappend args --detect-minarea $catpanel(param,detect-minarea)
    }

    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "Segmentation map error: $err"
	return
    }

    # Parse output: "OK N_SOURCES OUTPUT_PATH"
    set parts [split $data]
    if {[llength $parts] >= 3 && [lindex $parts 0] eq "OK"} {
	set nsrc [lindex $parts 1]
	set outpath [lindex $parts 2]

	# Load in new frame
	if {[catch {
	    CreateFrame
	    LoadFitsFile $outpath {} {}
	} err2]} {
	    set catpanel(status) "Error loading segmap: $err2"
	    return
	}

	set catpanel(status) "Segmentation map: $nsrc sources"
    } else {
	set catpanel(status) "Segmentation map: unexpected output"
    }
}

# --- B10: Non-Parametric Morphology (CAS/Gini/M20) ---

proc CatalogPanelMorphometry {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_morphometry.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_morphometry.py"
	return
    }

    # Save temp catalog
    set tmpcat [CatalogPanelSaveTempCatalog "morphometry"]
    if {$tmpcat eq {}} {
	set catpanel(status) "Failed to save temp catalog"
	return
    }

    set catpanel(status) "Measuring morphometry (CAS/Gini/M20)..."
    update idletasks

    if {[catch {
	set data [exec python3 $script $fn --catalog $tmpcat \
	    --n-workers $catpanel(param,n-workers) 2>@stderr]
    } err]} {
	set catpanel(status) "Morphometry error: $err"
	return
    }

    if {[string trim $data] eq {}} {
	set catpanel(status) "Morphometry: no output"
	return
    }

    # Add columns to alldata
    CatalogPanelAddColumnsFromTSV $data {CONC ASYM GINI M20 R_PETRO}
    set catpanel(status) "Morphometry complete (CAS/Gini/M20/Petrosian)"
}

# --- D17: Sérsic Fitting ---

proc CatalogPanelSersicFit {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_sersic.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_sersic.py"
	return
    }

    set tmpcat [CatalogPanelSaveTempCatalog "sersic"]
    if {$tmpcat eq {}} {
	set catpanel(status) "Failed to save temp catalog"
	return
    }

    set catpanel(status) "Sérsic profile fitting..."
    update idletasks

    set args [list python3 $script $fn --catalog $tmpcat]
    if {[info exists catpanel(param,mag-zeropoint)]} {
	lappend args --mag-zeropoint $catpanel(param,mag-zeropoint)
    }
    lappend args --n-workers $catpanel(param,n-workers)

    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "Sérsic fit error: $err"
	return
    }

    CatalogPanelAddColumnsFromTSV $data \
	{SERSIC_N SERSIC_RE SERSIC_IE SERSIC_ELLIP SERSIC_THETA SERSIC_CHI2}
    set catpanel(status) "Sérsic fitting complete"
}

# --- C13: PSF Photometry ---

proc CatalogPanelPSFPhotometry {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    # Check for PSF file
    if {![info exists catpanel(psf,file)] || ![file exists $catpanel(psf,file)]} {
	set catpanel(status) "No PSF file — build PSF first (Reconstruction > PSF Generation)"
	return
    }

    set script [CatalogPanelGetScript ds9_psf_phot.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_psf_phot.py"
	return
    }

    set tmpcat [CatalogPanelSaveTempCatalog "psfphot"]
    if {$tmpcat eq {}} {
	set catpanel(status) "Failed to save temp catalog"
	return
    }

    set catpanel(status) "PSF photometry..."
    update idletasks

    set args [list python3 $script $fn --catalog $tmpcat --psf $catpanel(psf,file)]
    if {[info exists catpanel(param,mag-zeropoint)]} {
	lappend args --mag-zeropoint $catpanel(param,mag-zeropoint)
    }
    lappend args --n-workers $catpanel(param,n-workers)

    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "PSF photometry error: $err"
	return
    }

    CatalogPanelAddColumnsFromTSV $data \
	{FLUX_PSF FLUXERR_PSF MAG_PSF MAGERR_PSF CHI2_PSF X_PSF Y_PSF}
    set catpanel(status) "PSF photometry complete"
}

# --- D16: Multi-Band Photometry ---

proc CatalogPanelMultiBand {} {
    global catpanel ds9

    set w {.multibandphot}
    set ed(ok) 0

    DialogCreate $w {Multi-Band Photometry} ed(ok)

    set f [ttk::frame $w.param]

    # Detection image
    ttk::label $f.ldet -text "Detection Image:" -anchor w
    ttk::entry $f.edet -textvariable ed(mb,detect) -width 40
    ttk::button $f.bdet -text "Browse..." -command {
	set ff [tk_getOpenFile -title "Detection Image" \
	    -filetypes {{{FITS Files} {.fits .fit .fts}} {{All Files} {*}}}]
	if {$ff ne {}} { set ed(mb,detect) $ff }
    }
    grid $f.ldet $f.edet $f.bdet -padx 4 -pady 2 -sticky w

    # Current frame as default
    set ed(mb,detect) [CatalogPanelGetFITS]

    # Band list
    ttk::label $f.lbands -text "Bands (name:file, one per line):" -anchor w
    text $f.tbands -width 50 -height 6
    grid $f.lbands -padx 4 -pady 2 -sticky w -columnspan 3
    grid $f.tbands -padx 4 -pady 2 -sticky we -columnspan 3

    # Buttons
    set bf [ttk::frame $w.buttons]
    ttk::button $bf.ok -text {Run} -command {set ed(ok) 1} -default active
    ttk::button $bf.cancel -text {Cancel} -command {set ed(ok) 0}
    pack $bf.ok $bf.cancel -side left -expand true -padx 2 -pady 4

    bind $w <Return> {set ed(ok) 1}
    ttk::separator $w.sep -orient horizontal
    pack $w.buttons $w.sep -side bottom -fill x
    pack $w.param -side top -fill both -expand true

    DialogWait $w ed(ok) $f.edet
    set detect_img $ed(mb,detect)
    set bands_text [$f.tbands get 1.0 end]
    destroy $w

    if {!$ed(ok)} { unset ed; return }
    unset ed

    if {$detect_img eq {} || ![file exists $detect_img]} {
	set catpanel(status) "Detection image not found"
	return
    }

    # Parse band lines
    set bands {}
    foreach line [split $bands_text \n] {
	set line [string trim $line]
	if {$line eq {}} continue
	lappend bands $line
    }
    if {[llength $bands] == 0} {
	set catpanel(status) "No bands specified"
	return
    }

    set script [CatalogPanelGetScript ds9_multiband.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_multiband.py"
	return
    }

    set catpanel(status) "Multi-band photometry ([llength $bands] bands)..."
    update idletasks

    set args [list python3 $script --detect-image $detect_img \
	--bands [join $bands ","]]

    # Use existing catalog if available
    if {[info exists catpanel(alldata)] && $catpanel(alldata) ne {}} {
	set tmpcat [CatalogPanelSaveTempCatalog "multiband"]
	if {$tmpcat ne {}} {
	    lappend args --catalog $tmpcat
	}
    }
    if {[info exists catpanel(param,mag-zeropoint)]} {
	lappend args --mag-zeropoint $catpanel(param,mag-zeropoint)
    }
    lappend args --n-workers $catpanel(param,n-workers)

    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "Multi-band error: $err"
	return
    }

    # Load as new catalog (replaces current)
    CatalogPanelLoadTSV $data "multi-band"
    set catpanel(status) "Multi-band photometry complete"
}

# --- D19: Crowded Field Photometry ---

proc CatalogPanelCrowdedPhot {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    if {![info exists catpanel(psf,file)] || ![file exists $catpanel(psf,file)]} {
	set catpanel(status) "No PSF file — build PSF first"
	return
    }

    set script [CatalogPanelGetScript ds9_crowded_phot.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_crowded_phot.py"
	return
    }

    set tmpcat [CatalogPanelSaveTempCatalog "crowded"]
    if {$tmpcat eq {}} {
	set catpanel(status) "Failed to save temp catalog"
	return
    }

    set catpanel(status) "Crowded field photometry..."
    update idletasks

    set args [list python3 $script $fn --catalog $tmpcat \
	--psf $catpanel(psf,file)]
    if {[info exists catpanel(param,mag-zeropoint)]} {
	lappend args --mag-zeropoint $catpanel(param,mag-zeropoint)
    }
    lappend args --n-workers $catpanel(param,n-workers)

    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "Crowded phot error: $err"
	return
    }

    CatalogPanelAddColumnsFromTSV $data \
	{FLUX_CROWD FLUXERR_CROWD MAG_CROWD X_CROWD Y_CROWD N_NEIGHBORS}
    set catpanel(status) "Crowded field photometry complete"
}

# --- C14: Cross-Match (VizieR) ---

proc CatalogPanelCrossMatch {} {
    global catpanel ed

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first"
	return
    }

    set w {.crossmatch}
    set ed(ok) 0

    DialogCreate $w {Cross-Match (VizieR)} ed(ok)

    set f [ttk::frame $w.param]

    # Catalog selection
    ttk::label $f.lcat -text "VizieR Catalog:" -anchor w
    ttk::combobox $f.ecat -textvariable ed(xm,catalog) -width 25 \
	-values {GAIA_DR3 2MASS PanSTARRS_DR1 SDSS_DR17 ALLWISE} -state readonly
    set ed(xm,catalog) GAIA_DR3
    grid $f.lcat $f.ecat -padx 4 -pady 2 -sticky w

    # Match radius
    ttk::label $f.lrad -text "Match Radius (arcsec):" -anchor w
    ttk::entry $f.erad -textvariable ed(xm,radius) -width 12
    set ed(xm,radius) 2.0
    grid $f.lrad $f.erad -padx 4 -pady 2 -sticky w

    # Buttons
    set bf [ttk::frame $w.buttons]
    ttk::button $bf.ok -text {Run} -command {set ed(ok) 1} -default active
    ttk::button $bf.cancel -text {Cancel} -command {set ed(ok) 0}
    pack $bf.ok $bf.cancel -side left -expand true -padx 2 -pady 4

    bind $w <Return> {set ed(ok) 1}
    ttk::separator $w.sep -orient horizontal
    pack $w.buttons $w.sep -side bottom -fill x
    pack $w.param -side top -fill both -expand true

    DialogWait $w ed(ok) $f.ecat
    set vizcat $ed(xm,catalog)
    set matchrad $ed(xm,radius)
    destroy $w

    if {!$ed(ok)} { unset ed; return }
    unset ed

    set script [CatalogPanelGetScript ds9_crossmatch.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_crossmatch.py"
	return
    }

    set tmpcat [CatalogPanelSaveTempCatalog "crossmatch"]
    if {$tmpcat eq {}} {
	set catpanel(status) "Failed to save temp catalog"
	return
    }

    set catpanel(status) "Cross-matching with $vizcat (r=${matchrad}\")..."
    update idletasks

    if {[catch {
	set data [exec python3 $script --catalog $tmpcat \
	    --vizier-cat $vizcat --radius $matchrad 2>@stderr]
    } err]} {
	set catpanel(status) "Cross-match error: $err"
	return
    }

    CatalogPanelAddColumnsFromTSV $data {MATCH_DIST MATCH_ID}
    set catpanel(status) "Cross-match complete ($vizcat)"
}

# --- C12: Dual-Image Extract ---

proc CatalogPanelDualExtract {} {
    global catpanel ds9 current ed

    set w {.dualextract}
    set ed(ok) 0

    DialogCreate $w {Dual-Image Extract} ed(ok)

    set f [ttk::frame $w.param]

    # Detection image
    ttk::label $f.ldet -text "Detection Image:" -anchor w
    ttk::entry $f.edet -textvariable ed(dual,detect) -width 40
    ttk::button $f.bdet -text "Browse..." -command {
	set ff [tk_getOpenFile -title "Detection Image" \
	    -filetypes {{{FITS Files} {.fits .fit .fts}} {{All Files} {*}}}]
	if {$ff ne {}} { set ed(dual,detect) $ff }
    }
    grid $f.ldet $f.edet $f.bdet -padx 4 -pady 2 -sticky w

    # Measurement image
    ttk::label $f.lmeas -text "Measurement Image:" -anchor w
    ttk::entry $f.emeas -textvariable ed(dual,measure) -width 40
    ttk::button $f.bmeas -text "Browse..." -command {
	set ff [tk_getOpenFile -title "Measurement Image" \
	    -filetypes {{{FITS Files} {.fits .fit .fts}} {{All Files} {*}}}]
	if {$ff ne {}} { set ed(dual,measure) $ff }
    }
    grid $f.lmeas $f.emeas $f.bmeas -padx 4 -pady 2 -sticky w

    # Set defaults from current frame
    set ed(dual,detect) [CatalogPanelGetFITS]
    set ed(dual,measure) {}

    # Buttons
    set bf [ttk::frame $w.buttons]
    ttk::button $bf.ok -text {Extract} -command {set ed(ok) 1} -default active
    ttk::button $bf.cancel -text {Cancel} -command {set ed(ok) 0}
    pack $bf.ok $bf.cancel -side left -expand true -padx 2 -pady 4

    bind $w <Return> {set ed(ok) 1}
    ttk::separator $w.sep -orient horizontal
    pack $w.buttons $w.sep -side bottom -fill x
    pack $w.param -side top -fill both -expand true

    DialogWait $w ed(ok) $f.edet
    set detect_img $ed(dual,detect)
    set measure_img $ed(dual,measure)
    destroy $w

    if {!$ed(ok)} { unset ed; return }
    unset ed

    if {$detect_img eq {} || ![file exists $detect_img]} {
	set catpanel(status) "Detection image not found"
	return
    }
    if {$measure_img eq {} || ![file exists $measure_img]} {
	set catpanel(status) "Measurement image not found"
	return
    }

    set script [CatalogPanelGetScript ds9_dual_extract.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_dual_extract.py"
	return
    }

    set catpanel(status) "Dual-image extraction..."
    update idletasks

    set args [list python3 $script \
	--detect-image $detect_img --measure-image $measure_img]
    if {[info exists catpanel(param,detect-thresh)]} {
	lappend args --detect-thresh $catpanel(param,detect-thresh)
    }
    if {[info exists catpanel(param,detect-minarea)]} {
	lappend args --detect-minarea $catpanel(param,detect-minarea)
    }
    if {[info exists catpanel(param,phot-aperture)]} {
	lappend args --phot-aperture $catpanel(param,phot-aperture)
    }
    if {[info exists catpanel(param,mag-zeropoint)]} {
	lappend args --mag-zeropoint $catpanel(param,mag-zeropoint)
    }

    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "Dual extract error: $err"
	return
    }

    CatalogPanelLoadTSV $data "dual-image"
    set catpanel(status) "Dual-image extraction complete"
}

# --- D18: Completeness Simulation ---

proc CatalogPanelCompleteness {} {
    global catpanel ed

    set fn [CatalogPanelGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	set catpanel(status) "No FITS image loaded"
	return
    }

    set w {.completeness}
    set ed(ok) 0

    DialogCreate $w {Completeness Simulation} ed(ok)

    set f [ttk::frame $w.param]

    set ed(comp,ninject) 500
    set ed(comp,magmin) 20.0
    set ed(comp,magmax) 28.0
    set ed(comp,nbins) 16

    set row 0
    foreach {vname vlabel} {
	comp,ninject {N Inject (per bin)}
	comp,magmin {Mag Min}
	comp,magmax {Mag Max}
	comp,nbins {N Bins}
    } {
	ttk::label $f.l$row -text "$vlabel:" -anchor w
	ttk::entry $f.e$row -textvariable ed($vname) -width 12
	grid $f.l$row $f.e$row -padx 4 -pady 2 -sticky w
	incr row
    }

    set bf [ttk::frame $w.buttons]
    ttk::button $bf.ok -text {Run} -command {set ed(ok) 1} -default active
    ttk::button $bf.cancel -text {Cancel} -command {set ed(ok) 0}
    pack $bf.ok $bf.cancel -side left -expand true -padx 2 -pady 4

    bind $w <Return> {set ed(ok) 1}
    ttk::separator $w.sep -orient horizontal
    pack $w.buttons $w.sep -side bottom -fill x
    pack $w.param -side top -fill both -expand true

    DialogWait $w ed(ok) $f.e0
    set ninject $ed(comp,ninject)
    set magmin $ed(comp,magmin)
    set magmax $ed(comp,magmax)
    set nbins $ed(comp,nbins)
    destroy $w

    if {!$ed(ok)} { unset ed; return }
    unset ed

    set script [CatalogPanelGetScript ds9_completeness.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_completeness.py"
	return
    }

    set catpanel(status) "Completeness simulation ($nbins bins)..."
    update idletasks

    set args [list python3 $script $fn \
	--n-inject $ninject --mag-min $magmin --mag-max $magmax --n-bins $nbins]
    if {[info exists catpanel(param,detect-thresh)]} {
	lappend args --detect-thresh $catpanel(param,detect-thresh)
    }
    if {[info exists catpanel(param,mag-zeropoint)]} {
	lappend args --mag-zeropoint $catpanel(param,mag-zeropoint)
    }
    lappend args --n-workers $catpanel(param,n-workers)

    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "Completeness error: $err"
	return
    }

    # Load completeness results as a new catalog view
    CatalogPanelLoadTSV $data "completeness"
    set catpanel(status) "Completeness simulation complete"
}

# ============================================================================
# ICL (Intra-Cluster Light) Detection Pipeline
# ============================================================================

proc CatalogPanelICLParamLoad {} {
    global catpanel

    set preffile [file join [file normalize ~] .ds9 icl.prf]
    if {![file exists $preffile]} return
    if {[catch {set fd [open $preffile r]} err]} return
    while {[gets $fd line] >= 0} {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set parts [split $line]
	if {[llength $parts] >= 2} {
	    set key [lindex $parts 0]
	    set val [lindex $parts 1]
	    if {[info exists catpanel(icl,param,$key)]} {
		set catpanel(icl,param,$key) $val
	    }
	}
    }
    close $fd
}

proc CatalogPanelICLParamSave {} {
    global catpanel

    set prefdir [file join [file normalize ~] .ds9]
    if {![file isdirectory $prefdir]} {
	file mkdir $prefdir
    }
    set preffile [file join $prefdir icl.prf]
    if {[catch {set fd [open $preffile w]} err]} return
    foreach pname {expand-factor bright-star-mag-limit bright-star-radius-scale \
		   interp-method detect-thresh \
		   bkg-method bkg-order bkg-sigma-clip bkg-sep-mesh \
		   bkg-iterative bkg-n-iterations bkg-convergence-tol bkg-refine-thresh \
		   rmin rmax nsteps spacing ellipticity pa \
		   mag-zeropoint pixel-scale \
		   mu-threshold mu-levels measure-radius} {
	puts $fd "$pname $catpanel(icl,param,$pname)"
    }
    close $fd
}

# --- 1. Source Masking ---

proc CatalogPanelICLMask {} {
    global catpanel current

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "ICL: No FITS file loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Creating source mask..."
    update idletasks

    set args [list python3 $script $fn --mode mask \
	--expand-factor $catpanel(icl,param,expand-factor) \
	--bright-star-mag-limit $catpanel(icl,param,bright-star-mag-limit) \
	--bright-star-radius-scale $catpanel(icl,param,bright-star-radius-scale) \
	--interp-method $catpanel(icl,param,interp-method) \
	--detect-thresh $catpanel(icl,param,detect-thresh) \
	--mask-output $catpanel(icl,mask_file) \
	--masked-output $catpanel(icl,masked_file)]

    # Add catalog if available
    if {[info exists catpanel(alldata)] && $catpanel(alldata) ne {}} {
	set catfile [CatalogPanelSaveTempCatalog icl]
	if {$catfile ne {}} {
	    lappend args --catalog $catfile
	}
    }

    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL mask error: $err"
	return
    }

    set catpanel(icl,has_mask) 1

    # Auto-display masked image in new frame
    if {[file exists $catpanel(icl,masked_file)]} {
	CreateFrame
	if {![catch {LoadFitsFile $catpanel(icl,masked_file) {} {}}]} {
	    global scale
	    set scale(mode) zscale
	    ChangeScaleMode
	}
    }

    set catpanel(status) "ICL: Source mask created"
}

proc CatalogPanelICLViewMask {} {
    global catpanel

    if {!$catpanel(icl,has_mask) || ![file exists $catpanel(icl,masked_file)]} {
	set catpanel(status) "ICL: No mask available — run Source Masking first"
	return
    }

    CreateFrame
    if {[catch {LoadFitsFile $catpanel(icl,masked_file) {} {}} err]} {
	set catpanel(status) "ICL: Error loading masked image: $err"
	return
    }
    global scale
    set scale(mode) zscale
    ChangeScaleMode
    set catpanel(status) "ICL: Masked image loaded in new frame"
}

# --- 2. Background Model ---

proc CatalogPanelICLBackground {method} {
    global catpanel

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "ICL: No FITS file loaded"
	return
    }

    # Use masked image if available, otherwise raw
    set input $fn
    if {$catpanel(icl,has_mask) && [file exists $catpanel(icl,masked_file)]} {
	set input $catpanel(icl,masked_file)
    }

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Fitting background ($method)..."
    update idletasks

    set args [list python3 $script $input --mode background \
	--bkg-method $method \
	--bkg-order $catpanel(icl,param,bkg-order) \
	--bkg-sigma-clip $catpanel(icl,param,bkg-sigma-clip) \
	--bkg-sep-mesh $catpanel(icl,param,bkg-sep-mesh) \
	--bkg-output $catpanel(icl,bkg_file) \
	--bgsub-output $catpanel(icl,bgsub_file)]

    if {$catpanel(icl,has_mask) && [file exists $catpanel(icl,mask_file)]} {
	lappend args --mask $catpanel(icl,mask_file)
    }

    # Iterative background refinement
    if {$catpanel(icl,param,bkg-iterative)} {
	lappend args --iterative \
	    --interp-method $catpanel(icl,param,interp-method) \
	    --bkg-n-iterations $catpanel(icl,param,bkg-n-iterations) \
	    --bkg-convergence-tol $catpanel(icl,param,bkg-convergence-tol) \
	    --bkg-refine-thresh $catpanel(icl,param,bkg-refine-thresh) \
	    --mask-output $catpanel(icl,mask_file)
    }

    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL background error: $err"
	return
    }

    set catpanel(icl,has_bkg) 1

    # Auto-display bgsub image in new frame
    if {[file exists $catpanel(icl,bgsub_file)]} {
	CreateFrame
	if {![catch {LoadFitsFile $catpanel(icl,bgsub_file) {} {}}]} {
	    global scale
	    set scale(mode) zscale
	    ChangeScaleMode
	}
    }

    set catpanel(status) "ICL: Background model ($method) complete"
}

proc CatalogPanelICLViewBkg {} {
    global catpanel

    if {!$catpanel(icl,has_bkg) || ![file exists $catpanel(icl,bgsub_file)]} {
	set catpanel(status) "ICL: No background model — run Background Model first"
	return
    }

    CreateFrame
    if {[catch {LoadFitsFile $catpanel(icl,bgsub_file) {} {}} err]} {
	set catpanel(status) "ICL: Error loading background-subtracted image: $err"
	return
    }
    global scale
    set scale(mode) zscale
    ChangeScaleMode
    set catpanel(status) "ICL: Background-subtracted image loaded in new frame"
}

# --- 3. BCG Center + Profile ---

proc CatalogPanelICLSetCenter {} {
    global catpanel

    # Use currently selected source as BCG center
    if {![info exists catpanel(selected_row)] || $catpanel(selected_row) < 0} {
	set catpanel(status) "ICL: Select a source (click a row) to set as BCG center"
	return
    }

    set lines [split $catpanel(alldata) \n]
    set headers [split [lindex $lines 0] "\t"]

    # Find X_IMAGE, Y_IMAGE columns
    set xcol -1
    set ycol -1
    for {set c 0} {$c < [llength $headers]} {incr c} {
	set h [string trim [lindex $headers $c]]
	if {$h eq "X_IMAGE"} { set xcol $c }
	if {$h eq "Y_IMAGE"} { set ycol $c }
    }
    if {$xcol < 0 || $ycol < 0} {
	set catpanel(status) "ICL: Cannot find X_IMAGE/Y_IMAGE columns"
	return
    }

    set row_idx [expr {$catpanel(selected_row) + 1}]
    set row_data [split [lindex $lines $row_idx] "\t"]
    if {[llength $row_data] <= $xcol || [llength $row_data] <= $ycol} {
	set catpanel(status) "ICL: Cannot read coordinates from selected row"
	return
    }

    # Store as 0-indexed
    set catpanel(icl,center_x) [expr {[lindex $row_data $xcol] - 1.0}]
    set catpanel(icl,center_y) [expr {[lindex $row_data $ycol] - 1.0}]

    set catpanel(status) "ICL: BCG center set to ([format %.1f [expr {$catpanel(icl,center_x)+1}]], [format %.1f [expr {$catpanel(icl,center_y)+1}]])"
}

proc CatalogPanelICLProfile {} {
    global catpanel

    if {$catpanel(icl,center_x) eq {} || $catpanel(icl,center_y) eq {}} {
	set catpanel(status) "ICL: Set BCG center first"
	return
    }

    # Use bgsub image if available, otherwise raw
    set fn [CatalogPanelGetFITS]
    if {$catpanel(icl,has_bkg) && [file exists $catpanel(icl,bgsub_file)]} {
	set fn $catpanel(icl,bgsub_file)
    }
    if {$fn eq {}} {
	set catpanel(status) "ICL: No image available"
	return
    }

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Measuring SB profile..."
    update idletasks

    set center "$catpanel(icl,center_x),$catpanel(icl,center_y)"

    set args [list python3 $script $fn --mode profile \
	--center $center \
	--rmin $catpanel(icl,param,rmin) \
	--rmax $catpanel(icl,param,rmax) \
	--nsteps $catpanel(icl,param,nsteps) \
	--spacing $catpanel(icl,param,spacing) \
	--ellipticity $catpanel(icl,param,ellipticity) \
	--pa $catpanel(icl,param,pa) \
	--mag-zeropoint $catpanel(icl,param,mag-zeropoint) \
	--pixel-scale $catpanel(icl,param,pixel-scale) \
	--profile-output $catpanel(icl,profile_file)]

    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL profile error: $err"
	return
    }

    set catpanel(icl,has_profile) 1

    # Display profile in catalog table
    CatalogPanelLoadTSV $data "icl_profile"

    # Draw annulus markers
    CatalogPanelICLDrawAnnuli

    set catpanel(status) "ICL: SB profile measured"
}

proc CatalogPanelICLDrawAnnuli {} {
    global catpanel current

    set frame $current(frame)
    if {$frame eq {}} return

    # Delete existing annulus markers
    catch {$frame marker catalog icl_annulus delete}

    # BCG center (0-indexed → 1-indexed for ds9 markers)
    set cx [expr {$catpanel(icl,center_x) + 1.0}]
    set cy [expr {$catpanel(icl,center_y) + 1.0}]

    # Draw BCG center cross
    set MARKER "physical; cross point $cx $cy 20 # color=cyan tag={icl_annulus}"
    set marker_var _icl_ann_center
    global $marker_var
    set $marker_var $MARKER
    $frame marker catalog command ds9 var $marker_var

    # Draw annulus rings at 25%, 50%, 75%, 100% of rmax
    set rmin $catpanel(icl,param,rmin)
    set rmax $catpanel(icl,param,rmax)
    foreach frac {0.25 0.50 0.75 1.00} {
	set r [expr {$rmin + ($rmax - $rmin) * $frac}]
	set ri [expr {int($r)}]
	set MARKER "physical; circle $cx $cy $r # color=green dash=1 width=1 tag={icl_annulus} text={r=${ri}}"
	set marker_var _icl_ann_ring_$ri
	global $marker_var
	set $marker_var $MARKER
	$frame marker catalog command ds9 var $marker_var
    }
}

proc CatalogPanelICLSectorProfile {} {
    global catpanel ed

    if {$catpanel(icl,center_x) eq {} || $catpanel(icl,center_y) eq {}} {
	set catpanel(status) "ICL: Set BCG center first"
	return
    }

    set w .iclsector
    if {[winfo exists $w]} {
	raise $w
	return
    }

    toplevel $w
    wm title $w "ICL Sector Profile"
    wm geometry $w 300x180

    set ed(icl,sector-pa) 0.0
    set ed(icl,sector-width) 90.0

    set f [ttk::frame $w.content]
    pack $f -fill both -expand true -padx 8 -pady 8

    set r 0
    ttk::label $f.lpa -text "Sector PA (deg):"
    ttk::entry $f.epa -textvariable ed(icl,sector-pa) -width 10
    grid $f.lpa -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $f.epa -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $f.lw -text "Sector Width (deg):"
    ttk::entry $f.ew -textvariable ed(icl,sector-width) -width 10
    grid $f.lw -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $f.ew -row $r -column 1 -sticky w -padx 4 -pady 4

    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8

    ttk::button $bf.run -text "Measure" -command [list CatalogPanelICLSectorProfileRun $w]
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.run -side right -padx 4
}

proc CatalogPanelICLSectorProfileRun {w} {
    global catpanel ed

    set fn [CatalogPanelGetFITS]
    if {$catpanel(icl,has_bkg) && [file exists $catpanel(icl,bgsub_file)]} {
	set fn $catpanel(icl,bgsub_file)
    }
    if {$fn eq {}} {
	set catpanel(status) "ICL: No image available"
	return
    }

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Measuring sector profile..."
    update idletasks

    set center "$catpanel(icl,center_x),$catpanel(icl,center_y)"

    set args [list python3 $script $fn --mode profile \
	--center $center \
	--rmin $catpanel(icl,param,rmin) \
	--rmax $catpanel(icl,param,rmax) \
	--nsteps $catpanel(icl,param,nsteps) \
	--spacing $catpanel(icl,param,spacing) \
	--mag-zeropoint $catpanel(icl,param,mag-zeropoint) \
	--pixel-scale $catpanel(icl,param,pixel-scale) \
	--sector-pa $ed(icl,sector-pa) \
	--sector-width $ed(icl,sector-width) \
	--profile-output $catpanel(icl,profile_file)]

    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL sector profile error: $err"
	return
    }

    destroy $w
    set catpanel(icl,has_profile) 1
    CatalogPanelLoadTSV $data "icl_sector_profile"
    set catpanel(status) "ICL: Sector profile measured (PA=$ed(icl,sector-pa), width=$ed(icl,sector-width))"
}

# --- 4. ICL Measurements ---

proc CatalogPanelICLMeasure {} {
    global catpanel

    if {!$catpanel(icl,has_profile) || ![file exists $catpanel(icl,profile_file)]} {
	set catpanel(status) "ICL: No profile available — run Measure Profile first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} { set fn "dummy.fits" }

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Computing ICL measurements..."
    update idletasks

    set args [list python3 $script $fn --mode measure \
	--profile-file $catpanel(icl,profile_file) \
	--mu-threshold $catpanel(icl,param,mu-threshold) \
	--mu-levels $catpanel(icl,param,mu-levels) \
	--pixel-scale $catpanel(icl,param,pixel-scale)]

    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL measure error: $err"
	return
    }

    CatalogPanelLoadTSV $data "icl_measurements"

    # Draw isophotal radius markers
    CatalogPanelICLDrawIsophotes $data

    set catpanel(status) "ICL: Measurements complete"
}

proc CatalogPanelICLDrawIsophotes {data} {
    global catpanel current

    set frame $current(frame)
    if {$frame eq {}} return
    if {$catpanel(icl,center_x) eq {} || $catpanel(icl,center_y) eq {}} return

    # Delete existing annulus markers
    catch {$frame marker catalog icl_annulus delete}

    set cx [expr {$catpanel(icl,center_x) + 1.0}]
    set cy [expr {$catpanel(icl,center_y) + 1.0}]

    # BCG center cross
    set MARKER "physical; cross point $cx $cy 20 # color=cyan tag={icl_annulus}"
    set marker_var _icl_iso_center
    global $marker_var
    set $marker_var $MARKER
    $frame marker catalog command ds9 var $marker_var

    # Parse TSV for R_MU* columns
    set lines [split $data \n]
    if {[llength $lines] < 2} return
    set headers [split [lindex $lines 0] "\t"]
    set vals [split [lindex $lines 1] "\t"]

    for {set c 0} {$c < [llength $headers]} {incr c} {
	set h [string trim [lindex $headers $c]]
	if {[string match "R_MU*" $h] && ![string match "*ARCSEC" $h]} {
	    set v [string trim [lindex $vals $c]]
	    if {$v ne "NaN" && [string is double $v]} {
		set r [expr {double($v)}]
		if {$r > 0 && $r < 1e6} {
		    # Extract mu level from column name (R_MU26 → 26)
		    set mu_label [string range $h 4 end]
		    set MARKER "physical; circle $cx $cy $r # color=magenta dash=1 width=1 tag={icl_annulus} text={mu=$mu_label}"
		    set marker_var _icl_iso_$c
		    global $marker_var
		    set $marker_var $MARKER
		    $frame marker catalog command ds9 var $marker_var
		}
	    }
	}
    }
}

# --- Multi-Threshold ICL ---

proc CatalogPanelICLMeasureMulti {} {
    global catpanel

    if {!$catpanel(icl,has_profile) || ![file exists $catpanel(icl,profile_file)]} {
	set catpanel(status) "ICL: No profile available — run Measure Profile first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} { set fn "dummy.fits" }

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Computing multi-threshold ICL..."
    update idletasks

    set args [list python3 $script $fn --mode measure-multi \
	--profile-file $catpanel(icl,profile_file) \
	--pixel-scale $catpanel(icl,param,pixel-scale)]

    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL multi-threshold error: $err"
	return
    }

    CatalogPanelLoadTSV $data "icl_multi_threshold"
    set catpanel(status) "ICL: Multi-threshold measurements complete"
}

# --- BCG+ICL Decomposition ---

proc CatalogPanelICLDecompose {} {
    global catpanel

    if {!$catpanel(icl,has_profile) || ![file exists $catpanel(icl,profile_file)]} {
	set catpanel(status) "ICL: No profile available — run Measure Profile first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} { set fn "dummy.fits" }

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: BCG+ICL decomposition..."
    update idletasks

    set args [list python3 $script $fn --mode decompose \
	--profile-file $catpanel(icl,profile_file) \
	--pixel-scale $catpanel(icl,param,pixel-scale) \
	--mag-zeropoint $catpanel(icl,param,mag-zeropoint)]

    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL decompose error: $err"
	return
    }

    CatalogPanelLoadTSV $data "icl_decomposition"
    set catpanel(status) "ICL: BCG+ICL decomposition complete"
}

# --- Color Profile Dialog ---

proc CatalogPanelICLColorProfile {} {
    global catpanel ed

    set w .iclcolor
    if {[winfo exists $w]} {
	raise $w
	return
    }

    toplevel $w
    wm title $w "ICL Color Profile"
    wm geometry $w 450x280

    set f [ttk::frame $w.content]
    pack $f -fill both -expand true -padx 8 -pady 8

    ttk::label $f.linfo -text "Specify 2-4 bands for color profile measurement:"
    grid $f.linfo -row 0 -column 0 -columnspan 3 -sticky w -padx 4 -pady 4

    for {set i 1} {$i <= 4} {incr i} {
	set ed(icl,color_name_$i) {}
	set ed(icl,color_file_$i) {}

	ttk::label $f.ln$i -text "Band $i name:"
	ttk::entry $f.en$i -textvariable ed(icl,color_name_$i) -width 8
	ttk::entry $f.ef$i -textvariable ed(icl,color_file_$i) -width 25
	ttk::button $f.bb$i -text "Browse" -command [list CatalogPanelICLColorBrowse $i]

	grid $f.ln$i -row $i -column 0 -sticky w -padx 4 -pady 2
	grid $f.en$i -row $i -column 1 -sticky w -padx 2 -pady 2
	grid $f.ef$i -row $i -column 2 -sticky ew -padx 2 -pady 2
	grid $f.bb$i -row $i -column 3 -sticky w -padx 2 -pady 2
    }

    grid columnconfigure $f 2 -weight 1

    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8
    ttk::button $bf.run -text "Run" -command [list CatalogPanelICLColorProfileRun $w]
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.run -side right -padx 4
}

proc CatalogPanelICLColorBrowse {idx} {
    global ed

    set fname [tk_getOpenFile -filetypes {{{FITS} {.fits .fit .fts}} {{All} *}}]
    if {$fname ne {}} {
	set ed(icl,color_file_$idx) $fname
    }
}

proc CatalogPanelICLColorProfileRun {w} {
    global catpanel ed

    if {$catpanel(icl,center_x) eq {} || $catpanel(icl,center_y) eq {}} {
	set catpanel(status) "ICL: Set BCG center first"
	return
    }

    # Build band spec
    set bands {}
    for {set i 1} {$i <= 4} {incr i} {
	set name [string trim $ed(icl,color_name_$i)]
	set file [string trim $ed(icl,color_file_$i)]
	if {$name ne {} && $file ne {} && [file exists $file]} {
	    if {$bands ne {}} { append bands , }
	    append bands "$name:$file"
	}
    }

    if {$bands eq {}} {
	set catpanel(status) "ICL: Need at least 2 bands with name and file"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} { set fn "dummy.fits" }

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Measuring color profile..."
    update idletasks

    set center "$catpanel(icl,center_x),$catpanel(icl,center_y)"

    set args [list python3 $script $fn --mode color \
	--center $center \
	--bands $bands \
	--rmin $catpanel(icl,param,rmin) \
	--rmax $catpanel(icl,param,rmax) \
	--nsteps $catpanel(icl,param,nsteps) \
	--spacing $catpanel(icl,param,spacing) \
	--mag-zeropoint $catpanel(icl,param,mag-zeropoint) \
	--pixel-scale $catpanel(icl,param,pixel-scale)]

    if {$catpanel(icl,has_mask) && [file exists $catpanel(icl,mask_file)]} {
	lappend args --mask $catpanel(icl,mask_file)
    }

    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL color profile error: $err"
	return
    }

    destroy $w
    CatalogPanelLoadTSV $data "icl_color_profile"
    set catpanel(status) "ICL: Color profile complete"
}

# --- Save/Load Profile ---

proc CatalogPanelICLSaveProfile {} {
    global catpanel

    if {!$catpanel(icl,has_profile) || ![file exists $catpanel(icl,profile_file)]} {
	set catpanel(status) "ICL: No profile available"
	return
    }

    set fname [tk_getSaveFile -defaultextension .tsv \
	-filetypes {{{TSV} {.tsv}} {{All} *}} \
	-initialfile icl_profile.tsv]
    if {$fname eq {}} return

    if {[catch {file copy -force $catpanel(icl,profile_file) $fname} err]} {
	set catpanel(status) "ICL: Save error: $err"
	return
    }
    set catpanel(status) "ICL: Profile saved to $fname"
}

proc CatalogPanelICLLoadProfile {} {
    global catpanel

    set fname [tk_getOpenFile -filetypes {{{TSV} {.tsv}} {{All} *}}]
    if {$fname eq {} || ![file exists $fname]} return

    if {[catch {set fd [open $fname r]} err]} {
	set catpanel(status) "ICL: Load error: $err"
	return
    }
    set data [read $fd]
    close $fd

    # Copy to standard location
    if {[catch {file copy -force $fname $catpanel(icl,profile_file)} err]} {
	# Non-fatal
    }

    set catpanel(icl,has_profile) 1
    CatalogPanelLoadTSV [string trim $data] "icl_profile"
    set catpanel(status) "ICL: Profile loaded from $fname"
}

# --- Settings Dialog ---

proc CatalogPanelICLSettings {} {
    global catpanel
    global ed

    set w .iclsettings
    if {[winfo exists $w]} {
	raise $w
	return
    }

    toplevel $w
    wm title $w "ICL Detection Settings"
    wm geometry $w 420x520

    # Copy current values
    foreach pname {expand-factor bright-star-mag-limit bright-star-radius-scale \
		   interp-method detect-thresh \
		   bkg-order bkg-sigma-clip bkg-sep-mesh \
		   bkg-iterative bkg-n-iterations bkg-convergence-tol bkg-refine-thresh \
		   rmin rmax nsteps spacing ellipticity pa \
		   mag-zeropoint pixel-scale \
		   mu-threshold mu-levels measure-radius} {
	set ed(icl,$pname) $catpanel(icl,param,$pname)
    }

    ttk::notebook $w.nb
    pack $w.nb -fill both -expand true -padx 8 -pady 8

    # --- Tab 1: Masking ---
    set t1 [ttk::frame $w.nb.mask]
    $w.nb add $t1 -text "Masking"

    set r 0
    ttk::label $t1.lef -text "Expand factor:"
    ttk::entry $t1.eef -textvariable ed(icl,expand-factor) -width 10
    grid $t1.lef -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.eef -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lml -text "Bright star mag limit:"
    ttk::entry $t1.eml -textvariable ed(icl,bright-star-mag-limit) -width 10
    grid $t1.lml -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.eml -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lrs -text "Bright star radius scale:"
    ttk::entry $t1.ers -textvariable ed(icl,bright-star-radius-scale) -width 10
    grid $t1.lrs -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.ers -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lim -text "Interpolation method:"
    ttk::combobox $t1.eim -textvariable ed(icl,interp-method) -width 10 \
	-values {linear cubic nearest}
    grid $t1.lim -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.eim -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.ldt -text "Detect threshold (sigma):"
    ttk::entry $t1.edt -textvariable ed(icl,detect-thresh) -width 10
    grid $t1.ldt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.edt -row $r -column 1 -sticky w -padx 4 -pady 4

    # --- Tab 2: Background ---
    set t2 [ttk::frame $w.nb.bkg]
    $w.nb add $t2 -text "Background"

    set r 0
    ttk::label $t2.lbo -text "Polynomial order:"
    ttk::entry $t2.ebo -textvariable ed(icl,bkg-order) -width 10
    grid $t2.lbo -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.ebo -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lsc -text "Sigma clip:"
    ttk::entry $t2.esc -textvariable ed(icl,bkg-sigma-clip) -width 10
    grid $t2.lsc -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.esc -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lsm -text "SEP mesh size:"
    ttk::entry $t2.esm -textvariable ed(icl,bkg-sep-mesh) -width 10
    grid $t2.lsm -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.esm -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::separator $t2.sep1 -orient horizontal
    grid $t2.sep1 -row $r -column 0 -columnspan 2 -sticky ew -padx 8 -pady 6
    incr r

    ttk::checkbutton $t2.cbi -text "Iterative refinement" \
	-variable ed(icl,bkg-iterative)
    grid $t2.cbi -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady 4
    incr r

    ttk::label $t2.lni -text "N iterations:"
    ttk::entry $t2.eni -textvariable ed(icl,bkg-n-iterations) -width 10
    grid $t2.lni -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.eni -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lct -text "Convergence tol:"
    ttk::entry $t2.ect -textvariable ed(icl,bkg-convergence-tol) -width 10
    grid $t2.lct -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.ect -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lrt -text "Refine threshold (sigma):"
    ttk::entry $t2.ert -textvariable ed(icl,bkg-refine-thresh) -width 10
    grid $t2.lrt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.ert -row $r -column 1 -sticky w -padx 4 -pady 4

    # --- Tab 3: Profile ---
    set t3 [ttk::frame $w.nb.prof]
    $w.nb add $t3 -text "Profile"

    set r 0
    ttk::label $t3.lrn -text "R min (px):"
    ttk::entry $t3.ern -textvariable ed(icl,rmin) -width 10
    grid $t3.lrn -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.ern -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lrx -text "R max (px):"
    ttk::entry $t3.erx -textvariable ed(icl,rmax) -width 10
    grid $t3.lrx -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.erx -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lns -text "N steps:"
    ttk::entry $t3.ens -textvariable ed(icl,nsteps) -width 10
    grid $t3.lns -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.ens -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lsp -text "Spacing:"
    ttk::combobox $t3.esp -textvariable ed(icl,spacing) -width 10 \
	-values {log linear}
    grid $t3.lsp -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.esp -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lel -text "Ellipticity:"
    ttk::entry $t3.eel -textvariable ed(icl,ellipticity) -width 10
    grid $t3.lel -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.eel -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lpa -text "PA (deg):"
    ttk::entry $t3.epa -textvariable ed(icl,pa) -width 10
    grid $t3.lpa -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.epa -row $r -column 1 -sticky w -padx 4 -pady 4

    # --- Tab 4: Calibration & ICL ---
    set t4 [ttk::frame $w.nb.cal]
    $w.nb add $t4 -text "Calibration"

    set r 0
    ttk::label $t4.lzp -text "Mag zeropoint:"
    ttk::entry $t4.ezp -textvariable ed(icl,mag-zeropoint) -width 10
    grid $t4.lzp -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.ezp -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lps -text "Pixel scale (arcsec/px):"
    ttk::entry $t4.eps -textvariable ed(icl,pixel-scale) -width 10
    grid $t4.lps -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.eps -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lmt -text "ICL mu threshold (mag/arcsec2):"
    ttk::entry $t4.emt -textvariable ed(icl,mu-threshold) -width 10
    grid $t4.lmt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.emt -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lml -text "Isophotal mu levels:"
    ttk::entry $t4.eml -textvariable ed(icl,mu-levels) -width 18
    grid $t4.lml -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.eml -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lmr -text "Measure radius (px):"
    ttk::entry $t4.emr -textvariable ed(icl,measure-radius) -width 10
    grid $t4.lmr -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.emr -row $r -column 1 -sticky w -padx 4 -pady 4

    # Buttons
    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8

    ttk::button $bf.apply -text "Apply" -command [list CatalogPanelICLSettingsApply $w]
    ttk::button $bf.defaults -text "Defaults" -command CatalogPanelICLSettingsDefaults
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.defaults -side right -padx 4
    pack $bf.apply -side right -padx 4
}

proc CatalogPanelICLSettingsApply {w} {
    global catpanel
    global ed

    foreach pname {expand-factor bright-star-mag-limit bright-star-radius-scale \
		   interp-method detect-thresh \
		   bkg-order bkg-sigma-clip bkg-sep-mesh \
		   bkg-iterative bkg-n-iterations bkg-convergence-tol bkg-refine-thresh \
		   rmin rmax nsteps spacing ellipticity pa \
		   mag-zeropoint pixel-scale \
		   mu-threshold mu-levels measure-radius} {
	set catpanel(icl,param,$pname) $ed(icl,$pname)
    }
    CatalogPanelICLParamSave
    set catpanel(status) "ICL settings applied and saved"
}

proc CatalogPanelICLSettingsDefaults {} {
    global ed

    set ed(icl,expand-factor) 2.5
    set ed(icl,bright-star-mag-limit) 18.0
    set ed(icl,bright-star-radius-scale) 10.0
    set ed(icl,interp-method) linear
    set ed(icl,detect-thresh) 1.5
    set ed(icl,bkg-order) 3
    set ed(icl,bkg-sigma-clip) 3.0
    set ed(icl,bkg-sep-mesh) 256
    set ed(icl,bkg-iterative) 0
    set ed(icl,bkg-n-iterations) 3
    set ed(icl,bkg-convergence-tol) 0.01
    set ed(icl,bkg-refine-thresh) 2.0
    set ed(icl,rmin) 5.0
    set ed(icl,rmax) 1000.0
    set ed(icl,nsteps) 80
    set ed(icl,spacing) log
    set ed(icl,ellipticity) 0.0
    set ed(icl,pa) 0.0
    set ed(icl,mag-zeropoint) 25.0
    set ed(icl,pixel-scale) 0.06
    set ed(icl,mu-threshold) 26.5
    set ed(icl,mu-levels) 26.0,27.0,28.0
    set ed(icl,measure-radius) 500.0
}

# ============================================================
# LSBG (Low Surface Brightness Galaxy) Detection Pipeline
# ============================================================

# --- LSBG Param Load/Save ---

proc CatalogPanelLSBGParamLoad {} {
    global catpanel

    set preffile [file join [file normalize ~] .ds9 lsbg.prf]
    if {![file exists $preffile]} return
    if {[catch {set fd [open $preffile r]} err]} return
    while {[gets $fd line] >= 0} {
	set line [string trim $line]
	if {$line eq {} || [string index $line 0] eq "#"} continue
	set parts [split $line]
	if {[llength $parts] >= 2} {
	    set key [lindex $parts 0]
	    set val [lindex $parts 1]
	    if {[info exists catpanel(lsbg,param,$key)]} {
		set catpanel(lsbg,param,$key) $val
	    }
	}
    }
    close $fd
}

proc CatalogPanelLSBGParamSave {} {
    global catpanel

    set prefdir [file join [file normalize ~] .ds9]
    if {![file isdirectory $prefdir]} {
	file mkdir $prefdir
    }
    set preffile [file join $prefdir lsbg.prf]
    if {[catch {set fd [open $preffile w]} err]} return
    foreach pname {mask-detect-thresh mask-detect-minarea mask-expand-factor \
		   bright-star-mag-limit bright-star-radius-scale \
		   mask-mag-threshold interp-method \
		   lsb-protect lsb-mu-threshold \
		   bkg-method bkg-mesh-size bkg-poly-order \
		   bkg-sigma-clip bkg-n-iterations bkg-refine-thresh \
		   bkg-rms-quantile bkg-convergence-tol \
		   detect-thresh detect-minarea detect-filter-kernel \
		   deblend-nthresh deblend-mincont \
		   multiscale multiscale-factors \
		   sersic-fit sersic-n-min sersic-n-max sersic-re-min \
		   sersic-cutout-scale sersic-max-nfev \
		   phot-apertures mag-zeropoint pixel-scale \
		   mu-eff-min mu-eff-max r-eff-min r-eff-max \
		   ellipticity-max min-snr \
		   sersic-n-filter-min sersic-n-filter-max sersic-chi2-max} {
	puts $fd "$pname $catpanel(lsbg,param,$pname)"
    }
    close $fd
}

# --- 1. Mask Bright Sources ---

proc CatalogPanelLSBGMask {} {
    global catpanel current

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Masking bright sources..."
    update idletasks

    set args [list python3 $script $fn --mode mask \
	--mask-detect-thresh $catpanel(lsbg,param,mask-detect-thresh) \
	--mask-detect-minarea $catpanel(lsbg,param,mask-detect-minarea) \
	--mask-expand-factor $catpanel(lsbg,param,mask-expand-factor) \
	--bright-star-mag-limit $catpanel(lsbg,param,bright-star-mag-limit) \
	--bright-star-radius-scale $catpanel(lsbg,param,bright-star-radius-scale) \
	--mask-mag-threshold $catpanel(lsbg,param,mask-mag-threshold) \
	--interp-method $catpanel(lsbg,param,interp-method) \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--lsb-mu-threshold $catpanel(lsbg,param,lsb-mu-threshold) \
	--mask-output $catpanel(lsbg,mask_file) \
	--masked-output $catpanel(lsbg,masked_file) \
	--n-workers $catpanel(param,n-workers)]
    if {$catpanel(lsbg,param,lsb-protect)} {
	lappend args --lsb-protect
    } else {
	lappend args --no-lsb-protect
    }

    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG mask error: $err"
	return
    }

    set catpanel(lsbg,has_mask) 1

    # Auto-display masked image in new frame
    CreateFrame
    if {[catch {LoadFitsFile $catpanel(lsbg,masked_file) {} {}} err]} {
	set catpanel(status) "LSBG: Mask done but could not display: $err"
    } else {
	global scale
	set scale(mode) zscale
	ChangeScaleMode
    }

    set catpanel(status) "LSBG: Bright source mask created (new frame)"
}

proc CatalogPanelLSBGViewMask {} {
    global catpanel

    if {!$catpanel(lsbg,has_mask) || ![file exists $catpanel(lsbg,masked_file)]} {
	set catpanel(status) "LSBG: No mask available — run Mask Bright Sources first"
	return
    }

    CreateFrame
    if {[catch {LoadFitsFile $catpanel(lsbg,masked_file) {} {}} err]} {
	set catpanel(status) "LSBG: Error loading masked image: $err"
	return
    }
    global scale
    set scale(mode) zscale
    ChangeScaleMode
    set catpanel(status) "LSBG: Masked image loaded in new frame"
}

# --- 2. Background Model (Iterative Cleaning) ---

proc CatalogPanelLSBGClean {method} {
    global catpanel

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }

    if {!$catpanel(lsbg,has_mask) || ![file exists $catpanel(lsbg,mask_file)]} {
	set catpanel(status) "LSBG: No mask — run Mask Bright Sources first"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Iterative background ($method)..."
    update idletasks

    set args [list python3 $script $fn --mode clean \
	--mask $catpanel(lsbg,mask_file) \
	--bkg-method $method \
	--bkg-mesh-size $catpanel(lsbg,param,bkg-mesh-size) \
	--bkg-poly-order $catpanel(lsbg,param,bkg-poly-order) \
	--bkg-sigma-clip $catpanel(lsbg,param,bkg-sigma-clip) \
	--bkg-n-iterations $catpanel(lsbg,param,bkg-n-iterations) \
	--bkg-refine-thresh $catpanel(lsbg,param,bkg-refine-thresh) \
	--bkg-rms-quantile $catpanel(lsbg,param,bkg-rms-quantile) \
	--bkg-convergence-tol $catpanel(lsbg,param,bkg-convergence-tol) \
	--interp-method $catpanel(lsbg,param,interp-method) \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--bkg-output $catpanel(lsbg,bkg_file) \
	--cleaned-output $catpanel(lsbg,cleaned_file) \
	--n-workers $catpanel(param,n-workers)]

    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG clean error: $err"
	return
    }

    set catpanel(lsbg,has_clean) 1

    # Auto-display cleaned image in new frame
    CreateFrame
    if {[catch {LoadFitsFile $catpanel(lsbg,cleaned_file) {} {}} err]} {
	set catpanel(status) "LSBG: Clean done but could not display: $err"
    } else {
	global scale
	set scale(mode) zscale
	ChangeScaleMode
    }

    set catpanel(status) "LSBG: Background cleaned ($method) (new frame)"
}

proc CatalogPanelLSBGViewClean {} {
    global catpanel

    if {!$catpanel(lsbg,has_clean) || ![file exists $catpanel(lsbg,cleaned_file)]} {
	set catpanel(status) "LSBG: No cleaned image — run Background Model first"
	return
    }

    CreateFrame
    if {[catch {LoadFitsFile $catpanel(lsbg,cleaned_file) {} {}} err]} {
	set catpanel(status) "LSBG: Error loading cleaned image: $err"
	return
    }
    global scale
    set scale(mode) zscale
    ChangeScaleMode
    set catpanel(status) "LSBG: Cleaned image loaded in new frame"
}

# --- 3. Detect LSBG Candidates ---

proc CatalogPanelLSBGDetect {} {
    global catpanel

    if {!$catpanel(lsbg,has_clean) || ![file exists $catpanel(lsbg,cleaned_file)]} {
	set catpanel(status) "LSBG: No cleaned image — run Background Model first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Detecting candidates..."
    update idletasks

    set args [list python3 $script $fn --mode detect \
	--cleaned $catpanel(lsbg,cleaned_file) \
	--detect-thresh $catpanel(lsbg,param,detect-thresh) \
	--detect-minarea $catpanel(lsbg,param,detect-minarea) \
	--detect-filter-kernel $catpanel(lsbg,param,detect-filter-kernel) \
	--deblend-nthresh $catpanel(lsbg,param,deblend-nthresh) \
	--deblend-mincont $catpanel(lsbg,param,deblend-mincont) \
	--multiscale-factors $catpanel(lsbg,param,multiscale-factors) \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--segmap-output $catpanel(lsbg,segmap_file) \
	--n-workers $catpanel(param,n-workers)]
    if {$catpanel(lsbg,param,multiscale)} {
	lappend args --multiscale
    } else {
	lappend args --no-multiscale
    }

    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG detect error: $err"
	return
    }

    # Parse detection results into table
    set catpanel(lsbg,detect_data) $result
    set catpanel(lsbg,has_detect) 1

    # Auto-display segmentation map in new frame
    if {[file exists $catpanel(lsbg,segmap_file)]} {
	CreateFrame
	if {[catch {LoadFitsFile $catpanel(lsbg,segmap_file) {} {}} err]} {
	    set catpanel(status) "LSBG: Detect done but could not display segmap: $err"
	} else {
	    global scale
	    set scale(mode) zscale
	    ChangeScaleMode
	}
    }

    # Load into panel table
    set catpanel(alldata) $result
    CatalogPanelDisplayTable

    # Count detections
    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]
    set catpanel(status) "LSBG: $nsrc candidates detected (segmap in new frame)"
}

# --- 4. Photometry ---

proc CatalogPanelLSBGPhotometry {} {
    global catpanel

    if {!$catpanel(lsbg,has_clean) || ![file exists $catpanel(lsbg,cleaned_file)]} {
	set catpanel(status) "LSBG: No cleaned image — run Background Model first"
	return
    }

    if {!$catpanel(lsbg,has_detect) || ![file exists $catpanel(lsbg,segmap_file)]} {
	set catpanel(status) "LSBG: No detections — run Detect first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Measuring photometry..."
    update idletasks

    set args [list python3 $script $fn --mode photometry \
	--cleaned $catpanel(lsbg,cleaned_file) \
	--segmap $catpanel(lsbg,segmap_file) \
	--detect-thresh $catpanel(lsbg,param,detect-thresh) \
	--detect-minarea $catpanel(lsbg,param,detect-minarea) \
	--detect-filter-kernel $catpanel(lsbg,param,detect-filter-kernel) \
	--deblend-nthresh $catpanel(lsbg,param,deblend-nthresh) \
	--deblend-mincont $catpanel(lsbg,param,deblend-mincont) \
	--phot-apertures $catpanel(lsbg,param,phot-apertures) \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--n-workers $catpanel(param,n-workers)]

    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG photometry error: $err"
	return
    }

    set catpanel(alldata) $result
    set catpanel(lsbg,has_catalog) 1
    CatalogPanelDisplayTable

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    CatalogPanelMarkAll
    set catpanel(status) "LSBG: $nsrc sources — photometry done"
}

# --- 5. Sérsic Profile Fit ---

proc CatalogPanelLSBGSersic {} {
    global catpanel

    if {!$catpanel(lsbg,has_clean) || ![file exists $catpanel(lsbg,cleaned_file)]} {
	set catpanel(status) "LSBG: No cleaned image — run Background Model first"
	return
    }

    if {!$catpanel(lsbg,has_detect) || ![file exists $catpanel(lsbg,segmap_file)]} {
	set catpanel(status) "LSBG: No detections — run Detect first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Fitting Sérsic profiles..."
    update idletasks

    set args [list python3 $script $fn --mode sersic \
	--cleaned $catpanel(lsbg,cleaned_file) \
	--segmap $catpanel(lsbg,segmap_file) \
	--detect-thresh $catpanel(lsbg,param,detect-thresh) \
	--detect-minarea $catpanel(lsbg,param,detect-minarea) \
	--detect-filter-kernel $catpanel(lsbg,param,detect-filter-kernel) \
	--deblend-nthresh $catpanel(lsbg,param,deblend-nthresh) \
	--deblend-mincont $catpanel(lsbg,param,deblend-mincont) \
	--phot-apertures $catpanel(lsbg,param,phot-apertures) \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--sersic-n-min $catpanel(lsbg,param,sersic-n-min) \
	--sersic-n-max $catpanel(lsbg,param,sersic-n-max) \
	--sersic-re-min $catpanel(lsbg,param,sersic-re-min) \
	--sersic-cutout-scale $catpanel(lsbg,param,sersic-cutout-scale) \
	--sersic-max-nfev $catpanel(lsbg,param,sersic-max-nfev) \
	--n-workers $catpanel(param,n-workers)]

    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG Sérsic fit error: $err"
	return
    }

    set catpanel(alldata) $result
    set catpanel(lsbg,has_catalog) 1
    CatalogPanelDisplayTable

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    CatalogPanelMarkAll
    set catpanel(status) "LSBG: $nsrc sources — Sérsic fit done"
}

# --- 6. Filter + Grade ---

proc CatalogPanelLSBGFilter {} {
    global catpanel

    if {!$catpanel(lsbg,has_clean) || ![file exists $catpanel(lsbg,cleaned_file)]} {
	set catpanel(status) "LSBG: No cleaned image — run Background Model first"
	return
    }

    if {!$catpanel(lsbg,has_detect) || ![file exists $catpanel(lsbg,segmap_file)]} {
	set catpanel(status) "LSBG: No detections — run Detect first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Filtering + grading candidates..."
    update idletasks

    set args [list python3 $script $fn --mode filter \
	--cleaned $catpanel(lsbg,cleaned_file) \
	--segmap $catpanel(lsbg,segmap_file) \
	--detect-thresh $catpanel(lsbg,param,detect-thresh) \
	--detect-minarea $catpanel(lsbg,param,detect-minarea) \
	--detect-filter-kernel $catpanel(lsbg,param,detect-filter-kernel) \
	--deblend-nthresh $catpanel(lsbg,param,deblend-nthresh) \
	--deblend-mincont $catpanel(lsbg,param,deblend-mincont) \
	--phot-apertures $catpanel(lsbg,param,phot-apertures) \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--mu-eff-min $catpanel(lsbg,param,mu-eff-min) \
	--mu-eff-max $catpanel(lsbg,param,mu-eff-max) \
	--r-eff-min $catpanel(lsbg,param,r-eff-min) \
	--r-eff-max $catpanel(lsbg,param,r-eff-max) \
	--ellipticity-max $catpanel(lsbg,param,ellipticity-max) \
	--min-snr $catpanel(lsbg,param,min-snr) \
	--sersic-n-min $catpanel(lsbg,param,sersic-n-min) \
	--sersic-n-max $catpanel(lsbg,param,sersic-n-max) \
	--sersic-re-min $catpanel(lsbg,param,sersic-re-min) \
	--sersic-cutout-scale $catpanel(lsbg,param,sersic-cutout-scale) \
	--sersic-max-nfev $catpanel(lsbg,param,sersic-max-nfev) \
	--sersic-n-filter-min $catpanel(lsbg,param,sersic-n-filter-min) \
	--sersic-n-filter-max $catpanel(lsbg,param,sersic-n-filter-max) \
	--sersic-chi2-max $catpanel(lsbg,param,sersic-chi2-max) \
	--n-workers $catpanel(param,n-workers)]
    if {$catpanel(lsbg,param,sersic-fit)} {
	lappend args --sersic-fit
    } else {
	lappend args --no-sersic-fit
    }

    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG filter error: $err"
	return
    }

    set catpanel(alldata) $result
    set catpanel(lsbg,has_catalog) 1
    CatalogPanelDisplayTable

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    CatalogPanelMarkAll
    set catpanel(status) "LSBG: $nsrc candidates passed filtering"
}

# --- Run Full Pipeline ---

proc CatalogPanelLSBGRunAll {} {
    global catpanel

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Running full pipeline..."
    update idletasks

    set args [list python3 $script $fn --mode run \
	--mask-detect-thresh $catpanel(lsbg,param,mask-detect-thresh) \
	--mask-detect-minarea $catpanel(lsbg,param,mask-detect-minarea) \
	--mask-expand-factor $catpanel(lsbg,param,mask-expand-factor) \
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
	--mask-output $catpanel(lsbg,mask_file) \
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

    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG pipeline error: $err"
	return
    }

    # Update state
    set catpanel(lsbg,has_mask) 1
    set catpanel(lsbg,has_clean) 1
    set catpanel(lsbg,has_detect) 1
    set catpanel(lsbg,has_catalog) 1

    # Load cleaned image in new frame
    CreateFrame
    if {[catch {LoadFitsFile $catpanel(lsbg,cleaned_file) {} {}} err]} {
	set catpanel(status) "LSBG: Warning — could not load cleaned image"
    } else {
	global scale
	set scale(mode) zscale
	ChangeScaleMode
    }

    # Load catalog
    set catpanel(alldata) $result
    CatalogPanelDisplayTable

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    # Mark all LSBG candidates
    CatalogPanelMarkAll
    set catpanel(status) "LSBG: Pipeline complete — $nsrc candidates"
}

# --- 7. Forced Photometry (Multi-Band) ---

proc CatalogPanelLSBGForcedPhot {} {
    global catpanel current

    # Check that we have a catalog
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "LSBG: No catalog — run pipeline first"
	return
    }

    # Select the other band FITS file
    set band_fits [tk_getOpenFile \
	-title "Select Other Band FITS Image" \
	-filetypes {
	    {{FITS Files} {.fits .fit .fts .fits.gz .fit.gz}}
	    {{All Files} {*}}
	}]
    if {$band_fits eq {}} return

    # Ask for band name
    set w .lsbg_bandname
    catch {destroy $w}
    toplevel $w
    wm title $w "Band Name"
    wm geometry $w 300x100
    wm transient $w .
    ttk::label $w.l -text "Enter band name (e.g. F606W):"
    ttk::entry $w.e -textvariable catpanel(lsbg,tmp_bandname)
    set catpanel(lsbg,tmp_bandname) ""
    ttk::frame $w.btns
    ttk::button $w.btns.ok -text "OK" -command [list set catpanel(lsbg,band_dialog_done) 1]
    ttk::button $w.btns.cancel -text "Cancel" -command [list set catpanel(lsbg,band_dialog_done) 0]
    pack $w.l -padx 10 -pady 5
    pack $w.e -padx 10 -fill x
    pack $w.btns -pady 5
    pack $w.btns.ok $w.btns.cancel -side left -padx 10
    focus $w.e
    bind $w.e <Return> [list set catpanel(lsbg,band_dialog_done) 1]
    tkwait variable catpanel(lsbg,band_dialog_done)
    set band_name $catpanel(lsbg,tmp_bandname)
    catch {destroy $w}
    if {!$catpanel(lsbg,band_dialog_done) || $band_name eq {}} return

    # Save current catalog to temp file
    set tmpcat [CatalogPanelSaveTempCatalog lsbg_forced]
    if {$tmpcat eq {}} {
	set catpanel(status) "LSBG: Failed to save temp catalog"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Forced photometry ($band_name)..."
    update idletasks

    set args [list python3 $script $band_fits --mode forced \
	--catalog $tmpcat \
	--band-name $band_name \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--n-workers $catpanel(param,n-workers)]

    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG forced phot error: $err"
	return
    }

    # Merge result columns into current catalog
    set col_names [list FLUX_$band_name FLUXERR_$band_name MAG_$band_name MAGERR_$band_name]
    CatalogPanelAddColumnsFromTSV $result $col_names

    # Load the other band image in a new frame
    CreateFrame
    if {[catch {LoadFitsFile $band_fits {} {}} err]} {
	set catpanel(status) "LSBG: Forced phot done but could not display band image"
    } else {
	global scale
	set scale(mode) zscale
	ChangeScaleMode
    }

    set catpanel(status) "LSBG: Forced photometry ($band_name) complete — columns added"
}

# --- LSBG Settings Dialog ---

proc CatalogPanelLSBGSettings {} {
    global catpanel
    global ed

    set w .lsbgsettings
    if {[winfo exists $w]} {
	raise $w
	return
    }

    toplevel $w
    wm title $w "LSBG Detection Settings"
    wm geometry $w 500x620

    # Copy current values
    foreach pname {mask-detect-thresh mask-detect-minarea mask-expand-factor \
		   bright-star-mag-limit bright-star-radius-scale \
		   mask-mag-threshold interp-method \
		   lsb-protect lsb-mu-threshold \
		   bkg-method bkg-mesh-size bkg-poly-order \
		   bkg-sigma-clip bkg-n-iterations bkg-refine-thresh \
		   bkg-rms-quantile bkg-convergence-tol \
		   detect-thresh detect-minarea detect-filter-kernel \
		   deblend-nthresh deblend-mincont \
		   multiscale multiscale-factors \
		   sersic-fit sersic-n-min sersic-n-max sersic-re-min \
		   sersic-cutout-scale sersic-max-nfev \
		   phot-apertures mag-zeropoint pixel-scale \
		   mu-eff-min mu-eff-max r-eff-min r-eff-max \
		   ellipticity-max min-snr \
		   sersic-n-filter-min sersic-n-filter-max sersic-chi2-max} {
	set ed(lsbg,$pname) $catpanel(lsbg,param,$pname)
    }

    ttk::notebook $w.nb
    pack $w.nb -fill both -expand true -padx 8 -pady 8

    # --- Tab 1: Masking ---
    set t1 [ttk::frame $w.nb.mask]
    $w.nb add $t1 -text "Masking"

    set r 0
    ttk::label $t1.lmt -text "Mask mag threshold:"
    ttk::entry $t1.emt -textvariable ed(lsbg,mask-mag-threshold) -width 10
    grid $t1.lmt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.emt -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lef -text "Expand factor:"
    ttk::entry $t1.eef -textvariable ed(lsbg,mask-expand-factor) -width 10
    grid $t1.lef -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.eef -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lml -text "Bright star mag limit:"
    ttk::entry $t1.eml -textvariable ed(lsbg,bright-star-mag-limit) -width 10
    grid $t1.lml -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.eml -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lrs -text "Bright star radius scale:"
    ttk::entry $t1.ers -textvariable ed(lsbg,bright-star-radius-scale) -width 10
    grid $t1.lrs -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.ers -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.lim -text "Interpolation method:"
    ttk::combobox $t1.eim -textvariable ed(lsbg,interp-method) -width 10 \
	-values {linear cubic nearest}
    grid $t1.lim -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.eim -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t1.ldt -text "Mask detect threshold:"
    ttk::entry $t1.edt -textvariable ed(lsbg,mask-detect-thresh) -width 10
    grid $t1.ldt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.edt -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::checkbutton $t1.clp -text "LSB structure protection" \
	-variable ed(lsbg,lsb-protect)
    grid $t1.clp -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady 4
    incr r

    ttk::label $t1.llm -text "LSB mu threshold:"
    ttk::entry $t1.elm -textvariable ed(lsbg,lsb-mu-threshold) -width 10
    grid $t1.llm -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.elm -row $r -column 1 -sticky w -padx 4 -pady 4

    # --- Tab 2: Background ---
    set t2 [ttk::frame $w.nb.bkg]
    $w.nb add $t2 -text "Background"

    set r 0
    ttk::label $t2.lbm -text "Method:"
    ttk::combobox $t2.ebm -textvariable ed(lsbg,bkg-method) -width 12 \
	-values {sep_large polynomial chebyshev}
    grid $t2.lbm -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.ebm -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lms -text "Mesh size (px):"
    ttk::entry $t2.ems -textvariable ed(lsbg,bkg-mesh-size) -width 10
    grid $t2.lms -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.ems -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lbo -text "Polynomial order:"
    ttk::entry $t2.ebo -textvariable ed(lsbg,bkg-poly-order) -width 10
    grid $t2.lbo -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.ebo -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lsc -text "Sigma clip:"
    ttk::entry $t2.esc -textvariable ed(lsbg,bkg-sigma-clip) -width 10
    grid $t2.lsc -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.esc -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lni -text "Iterations:"
    ttk::entry $t2.eni -textvariable ed(lsbg,bkg-n-iterations) -width 10
    grid $t2.lni -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.eni -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lrt -text "Refine threshold (sigma):"
    ttk::entry $t2.ert -textvariable ed(lsbg,bkg-refine-thresh) -width 10
    grid $t2.lrt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.ert -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lrq -text "RMS quantile:"
    ttk::entry $t2.erq -textvariable ed(lsbg,bkg-rms-quantile) -width 10
    grid $t2.lrq -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.erq -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t2.lct -text "Convergence tolerance:"
    ttk::entry $t2.ect -textvariable ed(lsbg,bkg-convergence-tol) -width 10
    grid $t2.lct -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t2.ect -row $r -column 1 -sticky w -padx 4 -pady 4

    # --- Tab 3: Detection ---
    set t3 [ttk::frame $w.nb.det]
    $w.nb add $t3 -text "Detection"

    set r 0
    ttk::label $t3.ldt -text "Detect threshold (sigma):"
    ttk::entry $t3.edt -textvariable ed(lsbg,detect-thresh) -width 10
    grid $t3.ldt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.edt -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lma -text "Min area (px):"
    ttk::entry $t3.ema -textvariable ed(lsbg,detect-minarea) -width 10
    grid $t3.lma -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.ema -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lfk -text "Filter kernel:"
    ttk::combobox $t3.efk -textvariable ed(lsbg,detect-filter-kernel) -width 12 \
	-values {none gauss3x3 gauss5x5 gauss7x7 gauss9x9 tophat5 tophat7 mexhat}
    grid $t3.lfk -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.efk -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lnt -text "Deblend N thresholds:"
    ttk::entry $t3.ent -textvariable ed(lsbg,deblend-nthresh) -width 10
    grid $t3.lnt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.ent -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t3.lmc -text "Deblend min contrast:"
    ttk::entry $t3.emc -textvariable ed(lsbg,deblend-mincont) -width 10
    grid $t3.lmc -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.emc -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::checkbutton $t3.cms -text "Multi-scale detection" \
	-variable ed(lsbg,multiscale)
    grid $t3.cms -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady 4
    incr r

    ttk::label $t3.lmf -text "Scale factors:"
    ttk::entry $t3.emf -textvariable ed(lsbg,multiscale-factors) -width 10
    grid $t3.lmf -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t3.emf -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::checkbutton $t3.csf -text "Sérsic profile fitting" \
	-variable ed(lsbg,sersic-fit)
    grid $t3.csf -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady 4

    # --- Tab 4: Calibration & Filter ---
    set t4 [ttk::frame $w.nb.cal]
    $w.nb add $t4 -text "Calibration & Filter"

    set r 0
    ttk::label $t4.lzp -text "Mag zeropoint:"
    ttk::entry $t4.ezp -textvariable ed(lsbg,mag-zeropoint) -width 10
    grid $t4.lzp -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.ezp -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lps -text "Pixel scale (arcsec/px):"
    ttk::entry $t4.eps -textvariable ed(lsbg,pixel-scale) -width 10
    grid $t4.lps -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.eps -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lmn -text "mu_eff min (mag/arcsec2):"
    ttk::entry $t4.emn -textvariable ed(lsbg,mu-eff-min) -width 10
    grid $t4.lmn -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.emn -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lmx -text "mu_eff max (mag/arcsec2):"
    ttk::entry $t4.emx -textvariable ed(lsbg,mu-eff-max) -width 10
    grid $t4.lmx -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.emx -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lrn -text "R_eff min (arcsec):"
    ttk::entry $t4.ern -textvariable ed(lsbg,r-eff-min) -width 10
    grid $t4.lrn -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.ern -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lrx -text "R_eff max (arcsec):"
    ttk::entry $t4.erx -textvariable ed(lsbg,r-eff-max) -width 10
    grid $t4.lrx -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.erx -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lex -text "Ellipticity max:"
    ttk::entry $t4.eex -textvariable ed(lsbg,ellipticity-max) -width 10
    grid $t4.lex -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.eex -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lsn -text "Min SNR:"
    ttk::entry $t4.esn -textvariable ed(lsbg,min-snr) -width 10
    grid $t4.lsn -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.esn -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lap -text "Phot apertures (px):"
    ttk::entry $t4.eap -textvariable ed(lsbg,phot-apertures) -width 18
    grid $t4.lap -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.eap -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::separator $t4.sep -orient horizontal
    grid $t4.sep -row $r -column 0 -columnspan 2 -sticky ew -padx 8 -pady 6
    incr r

    ttk::label $t4.lsnm -text "Sérsic n min:"
    ttk::entry $t4.esnm -textvariable ed(lsbg,sersic-n-filter-min) -width 10
    grid $t4.lsnm -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.esnm -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lsnx -text "Sérsic n max:"
    ttk::entry $t4.esnx -textvariable ed(lsbg,sersic-n-filter-max) -width 10
    grid $t4.lsnx -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.esnx -row $r -column 1 -sticky w -padx 4 -pady 4
    incr r

    ttk::label $t4.lsc2 -text "Sérsic chi2 max:"
    ttk::entry $t4.esc2 -textvariable ed(lsbg,sersic-chi2-max) -width 10
    grid $t4.lsc2 -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.esc2 -row $r -column 1 -sticky w -padx 4 -pady 4

    # Buttons
    set bf [ttk::frame $w.buttons]
    pack $bf -fill x -padx 8 -pady 8

    ttk::button $bf.apply -text "Apply" -command [list CatalogPanelLSBGSettingsApply $w]
    ttk::button $bf.defaults -text "Defaults" -command CatalogPanelLSBGSettingsDefaults
    ttk::button $bf.close -text "Close" -command [list destroy $w]
    pack $bf.close -side right -padx 4
    pack $bf.defaults -side right -padx 4
    pack $bf.apply -side right -padx 4
}

proc CatalogPanelLSBGSettingsApply {w} {
    global catpanel
    global ed

    foreach pname {mask-detect-thresh mask-detect-minarea mask-expand-factor \
		   bright-star-mag-limit bright-star-radius-scale \
		   mask-mag-threshold interp-method \
		   lsb-protect lsb-mu-threshold \
		   bkg-method bkg-mesh-size bkg-poly-order \
		   bkg-sigma-clip bkg-n-iterations bkg-refine-thresh \
		   bkg-rms-quantile bkg-convergence-tol \
		   detect-thresh detect-minarea detect-filter-kernel \
		   deblend-nthresh deblend-mincont \
		   multiscale multiscale-factors \
		   sersic-fit sersic-n-min sersic-n-max sersic-re-min \
		   sersic-cutout-scale sersic-max-nfev \
		   phot-apertures mag-zeropoint pixel-scale \
		   mu-eff-min mu-eff-max r-eff-min r-eff-max \
		   ellipticity-max min-snr \
		   sersic-n-filter-min sersic-n-filter-max sersic-chi2-max} {
	set catpanel(lsbg,param,$pname) $ed(lsbg,$pname)
    }
    CatalogPanelLSBGParamSave
    set catpanel(status) "LSBG settings applied and saved"
}

proc CatalogPanelLSBGSettingsDefaults {} {
    global ed

    set ed(lsbg,mask-detect-thresh) 1.5
    set ed(lsbg,mask-detect-minarea) 5
    set ed(lsbg,mask-expand-factor) 3.0
    set ed(lsbg,bright-star-mag-limit) 18.0
    set ed(lsbg,bright-star-radius-scale) 12.0
    set ed(lsbg,mask-mag-threshold) 22.0
    set ed(lsbg,interp-method) linear
    set ed(lsbg,lsb-protect) 1
    set ed(lsbg,lsb-mu-threshold) 24.0
    set ed(lsbg,bkg-method) sep_large
    set ed(lsbg,bkg-mesh-size) 256
    set ed(lsbg,bkg-poly-order) 3
    set ed(lsbg,bkg-sigma-clip) 3.0
    set ed(lsbg,bkg-n-iterations) 3
    set ed(lsbg,bkg-refine-thresh) 2.0
    set ed(lsbg,bkg-rms-quantile) 0.25
    set ed(lsbg,bkg-convergence-tol) 0.01
    set ed(lsbg,detect-thresh) 0.8
    set ed(lsbg,detect-minarea) 50
    set ed(lsbg,detect-filter-kernel) gauss5x5
    set ed(lsbg,deblend-nthresh) 32
    set ed(lsbg,deblend-mincont) 0.005
    set ed(lsbg,multiscale) 1
    set ed(lsbg,multiscale-factors) 1,2,4
    set ed(lsbg,sersic-fit) 1
    set ed(lsbg,sersic-n-min) 0.2
    set ed(lsbg,sersic-n-max) 10.0
    set ed(lsbg,sersic-re-min) 0.5
    set ed(lsbg,sersic-cutout-scale) 5.0
    set ed(lsbg,sersic-max-nfev) 500
    set ed(lsbg,phot-apertures) 5,10,20,40
    set ed(lsbg,mag-zeropoint) 25.0
    set ed(lsbg,pixel-scale) 0.06
    set ed(lsbg,mu-eff-min) 24.0
    set ed(lsbg,mu-eff-max) 30.0
    set ed(lsbg,r-eff-min) 2.5
    set ed(lsbg,r-eff-max) 60.0
    set ed(lsbg,ellipticity-max) 0.7
    set ed(lsbg,min-snr) 2.0
    set ed(lsbg,sersic-n-filter-min) 0.3
    set ed(lsbg,sersic-n-filter-max) 6.0
    set ed(lsbg,sersic-chi2-max) 10.0
}

# ===== Per-Frame Catalog Panel State Management =====

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
	set catpanel(icl,has_mask) 0
	set catpanel(icl,has_bkg) 0
	set catpanel(icl,has_profile) 0
	set catpanel(icl,center_x) {}
	set catpanel(icl,center_y) {}
	set catpanel(icl,mask_file) [file join [file normalize ~] .ds9 icl_mask.fits]
	set catpanel(icl,masked_file) [file join [file normalize ~] .ds9 icl_masked.fits]
	set catpanel(icl,bkg_file) [file join [file normalize ~] .ds9 icl_background.fits]
	set catpanel(icl,bgsub_file) [file join [file normalize ~] .ds9 icl_bgsub.fits]
	set catpanel(icl,profile_file) [file join [file normalize ~] .ds9 icl_profile.tsv]
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
    if {$old_frame ne {} && $old_frame ne $new_frame} {
	CatalogPanelSaveFrameState $old_frame
    }
    CatalogPanelRestoreFrameState $new_frame
}
