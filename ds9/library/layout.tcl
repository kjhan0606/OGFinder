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
    # buttons bar duplicates the menubar, so hide it by default
    set view(buttons) 0
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
    # the image tab strip (ogf_tile.tcl, 24 px) is carved out of the canvas request: default window size unchanged
    $ds9(canvas) configure -height [expr {$canvas(height)-24}]

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

    # image tab strip (frames as tabs + Single / Tile), takes its height from the canvas so the window
    # size does not change (ogf_tile.tcl)
    catch {OGFTileBuildStrip $ds9(image)}

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

    # OGFinder core: plugin registry (plugins/*/plugin.json) and the workflow UI built from it
    OGFCoreStart

    # Info area: same height as the left pane header so the catalog
    # table lines up with the image display
    set catpanel(detached) 0
    set catpanel(infoarea) [ttk::frame $f.info -height 154]
    pack propagate $f.info 0
    OGFUIBuild $f
    ttk::separator $f.infosep -orient horizontal
    set catpanel(hdrw) [expr {[info exists ds9(header)] ? $ds9(header) : {}}]

    # Search/Filter bar
    set catpanel(searchbar) [ttk::frame $f.searchbar]
    ttk::label $f.searchbar.lbl -text "Filter:"
    set catpanel(search_var) {}
    ttk::entry $f.searchbar.entry -textvariable catpanel(search_var) -width 20
    # Compact button so the Filter row has the same height (21 px) as
    # the File row in the left pane
    ttk::style configure CatApply.TButton -padding {2 0}
    ttk::button $f.searchbar.go -text "Apply" \
	-command CatalogPanelFilter -width 6 -style CatApply.TButton
    pack $f.searchbar.lbl -side left -padx 2 -pady 0
    pack $f.searchbar.entry -side left -padx 2 -pady 0 -fill x -expand true
    pack $f.searchbar.go -side right -padx 2 -pady 0 -fill y
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
    pack $f.statusbar.lbl -fill x -expand true -padx 2 -pady 0

    # Pack all into catalog frame
    # Selected source summary (Extract / Mark All / Clear moved to the workflow tabs)
    set catpanel(sel,text) {No source selected}
    ttk::label $f.selinfo -textvariable catpanel(sel,text) \
	-anchor nw -justify left -relief groove -padding 4

    pack $f.menubar -fill x -side top
    pack $f.info -fill x -side top
    pack $f.infosep -fill x -side top
    pack $f.info.tabs -in $f.info -fill x -side top
    pack $f.info.content -in $f.info -fill x -side top -pady 1
    pack $f.searchbar -in $f.info -fill x -side top -pady 2
    pack $f.statusbar -in $f.info -fill x -side top
    pack $f.selinfo -in $f.info -fill both -expand true -side top \
	-padx 2 -pady 2
    pack $f.tblf -fill both -expand true -side top

    # Keep the info area height equal to the left header height
    if {$catpanel(hdrw) ne {}} {
	bind $catpanel(hdrw) <Configure> {+CatalogPanelSyncInfoHeight}
    }

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

    # Ensure ~/.ds9 directory exists
    set ds9dir [file join [file normalize ~] .ds9]
    if {![file isdirectory $ds9dir]} {
	file mkdir $ds9dir
    }

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

    # ICL state (default paths; updated per-FITS by CatalogPanelICLUpdateFiles)
    set catpanel(icl,fits_base)    {}
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
    set catpanel(icl,click_mode)   0
    set catpanel(icl,cmdlog)       {}
    set catpanel(icl,param,expand-factor)          1.5
    set catpanel(icl,param,bright-star-mag-limit)  18.0
    set catpanel(icl,param,bright-star-radius-scale) 10.0
    set catpanel(icl,param,interp-method)          linear
    set catpanel(icl,param,detect-thresh)          5.0
    set catpanel(icl,param,max-dilate-radius)      20
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

    # LSBG state (default paths; updated per-FITS by CatalogPanelLSBGUpdateFiles)
    set catpanel(lsbg,fits_base)    {}
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
    set catpanel(lsbg,cmdlog)       {}
    set catpanel(lsbg,param,mask-detect-thresh)         1.5
    set catpanel(lsbg,param,mask-detect-minarea)        5
    set catpanel(lsbg,param,mask-expand-factor)         1.5
    set catpanel(lsbg,param,max-dilate-radius)          30
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
    set catpanel(lsbg,param,svm-classify)               0
    set catpanel(lsbg,param,svm-threshold)              0.3
    set catpanel(lsbg,param,svm-checkpoint)             {}
    CatalogPanelLSBGParamLoad

    # Interactive Plot state
    set catpanel(plot,counter) 0

    # Photo-z state
    set catpanel(photoz,param,bands)       {g,r,i,z}
    set catpanel(photoz,param,mag-columns) {}
    set catpanel(photoz,param,checkpoint)  {}
    CatalogPanelPhotoZParamLoad

    # SED Fitting state
    set catpanel(sed,param,bands)       {g,r,i,z}
    set catpanel(sed,param,mag-columns) {}
    set catpanel(sed,param,photoz-column) PHOTO_Z
    set catpanel(sed,param,checkpoint-emulator) {}
    set catpanel(sed,param,checkpoint-inverse)  {}
    set catpanel(sed,param,backend) auto
    CatalogPanelSEDParamLoad

    # Bulge+Disk state
    set catpanel(bd,param,max-sources)   100
    set catpanel(bd,param,free-bulge-n)  0
    set catpanel(bd,param,mag-zeropoint) 25.0
    set catpanel(bd,param,pixel-scale)   0.263
    CatalogPanelBDParamLoad

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

    # OGFinder extension modules (link / bands / mask)
    OGFBandsInit
    OGFMaskInit
    OGFLinkInit
    OGFTDInit
    OGFSessInit
    OGFCoreReady

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

proc CatalogPanelFilter {} {
    global catpanel

    # time-domain views filter their own rows
    if {[catch {::ogf::td::active} _tda]} {set _tda 0}
    if {$_tda} {::ogf::td::fill; return}
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

    catch {OGFTDAppendKindColumn $ncols $row}
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
# Detach the catalog panel into its own window, or attach it back to the
# main window.  catpanel(detached) holds the requested state.
proc CatalogPanelToggleDetach {} {
    global ds9
    global catpanel

    set f $ds9(catalog_frame)
    set tl [winfo toplevel $f]
    set is_detached [expr {$tl eq $f}]

    if {$catpanel(detached) && !$is_detached} {
	# Detach
	set cw [winfo width $f]
	set ch [winfo height $f]
	set mw [winfo width .]
	set mh [winfo height .]
	set catpanel(detach,cw) $cw
	set catpanel(detach,ch) $ch
	$ds9(toppw) forget $f
	wm manage $f
	wm title $f "Catalog - [wm title .]"
	wm protocol $f WM_DELETE_WINDOW {
	    set catpanel(detached) 0
	    CatalogPanelToggleDetach
	}
	wm geometry $f ${cw}x${ch}
	# give the width back to the main window
	set nw [expr {max($mw - $cw - 5, 400)}]
	wm geometry . ${nw}x${mh}
	update idletasks
	CatalogPanelSyncInfoHeight
    } elseif {!$catpanel(detached) && $is_detached} {
	# Attach
	set cw [winfo width $f]
	set mw [winfo width .]
	set mh [winfo height .]
	wm forget $f
	$ds9(toppw) add $f -weight 2
	wm geometry . [expr {$mw + $cw + 5}]x${mh}
	update idletasks
	CatalogPanelSyncInfoHeight
    }
}

proc CatalogPanelSyncInfoHeight {} {
    global catpanel
    global ds9
    if {![info exists catpanel(hdrw)] || ![winfo exists $catpanel(hdrw)]} return
    set h [winfo height $catpanel(hdrw)]
    if {$h > 1} {
	catch {$ds9(catalog_frame).info configure -height $h}
    }
}

# Show key columns of the selected catalog row in the info area
proc CatalogPanelUpdateSelInfo {row} {
    global catpanel
    if {![catch {::ogf::td::selinfo $row} _tds] && $_tds} return
    global $catpanel(tbldb)

    set want {NUMBER X_IMAGE Y_IMAGE ALPHA_J2000 DELTA_J2000 MAG_AUTO
	FWHM_IMAGE ELLIPTICITY CLASS_STAR}
    set idx {}
    set ncols [$catpanel(tbl) cget -cols]
    for {set c 1} {$c <= $ncols} {incr c} {
	if {[info exists ${catpanel(tbldb)}(0,$c)]} {
	    dict set idx [set ${catpanel(tbldb)}(0,$c)] $c
	}
    }
    set v {}
    foreach k $want {
	set val {-}
	if {[dict exists $idx $k]} {
	    set c [dict get $idx $k]
	    if {[info exists ${catpanel(tbldb)}($row,$c)]} {
		set val [set ${catpanel(tbldb)}($row,$c)]
	    }
	}
	if {[string is double -strict $val] && $k ne {NUMBER}} {
	    set val [format %.4g $val]
	}
	dict set v $k $val
    }
    set catpanel(sel,text) [format \
	"Source #%s   x,y = %s, %s\nRA,Dec = %s, %s\nMAG_AUTO = %s   FWHM = %s   e = %s   Star = %s" \
	[dict get $v NUMBER] [dict get $v X_IMAGE] [dict get $v Y_IMAGE] \
	[dict get $v ALPHA_J2000] [dict get $v DELTA_J2000] \
	[dict get $v MAG_AUTO] [dict get $v FWHM_IMAGE] \
	[dict get $v ELLIPTICITY] [dict get $v CLASS_STAR]]
}

proc CatalogPanelSelectCmd {prev cur} {
    global catpanel

    # cur is "row,col" of current selection
    set row [lindex [split $cur ,] 0]
    if {![string is integer -strict $row] || $row <= 0} return

    catch {CatalogPanelUpdateSelInfo $row}
    # keep the link state in step with clicks in the table itself
    set num [OGFNumberOfRow $row]
    if {$num ne {}} {
	set catpanel(sel,nums) [list $num]
	catch {$catpanel(tbl) tag delete msel}
	catch {$catpanel(tbl) tag configure msel -bg #9cc7f5 -fg black}
	OGFSetSelBase $catpanel(sel,text)
    }
    after cancel CatalogPanelGotoSource
    after 100 [list CatalogPanelGotoSource $row]
}

proc CatalogPanelGotoSource {row {pan 1}} {
    global catpanel
    global current
    global ds9

    if {![catch {::ogf::td::goto $row $pan} _tdg] && $_tdg} return

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

    # Delete previous selection markers (in every registered band frame)
    set frame $current(frame)
    set frames [OGFSelFrames]
    foreach fr $frames {
	catch {$fr marker catalog sextract_sel delete}
    }

    # Rebuild sextract_all markers from alldata to keep image in sync
    if {[info exists catpanel(markall,on)] && $catpanel(markall,on)} {
	CatalogPanelCreateAllMarkers
    }

    # Use global variable for marker creation (var form requires global access)
    global sextract_sel_reg

    # Marker in each band frame (positions mapped through the WCS when the
    # band has a different pixel grid; identity on a shared grid)
    set num [OGFNumberOfRow $row]
    foreach fr $frames {
	if {![$fr has fits]} continue
	lassign [OGFBandPos $fr $num $x $y] bx by same
	set isdet [expr {$same}]
	set sextract_sel_reg "image\ncross point($bx $by) # color=cyan width=2 point=cross 15 tag={sextract_sel} select=0 edit=0 move=0 rotate=0 delete=1\n"
	catch {$fr marker catalog command ds9 var sextract_sel_reg}
	if {$fr ne $frame && !$isdet} continue
	if {$r_eff_pix > 0} {
	    set sextract_sel_reg "image\ncircle($bx $by ${r_eff_pix}i) # color=green width=2 dash=1 tag={sextract_sel} select=0 edit=0 move=0 rotate=0 delete=1\n"
	    catch {$fr marker catalog command ds9 var sextract_sel_reg}
	}
	set sextract_sel_reg "image\nellipse($bx $by ${semi_a}i ${semi_b}i $theta) # color=green width=2 dash=1 tag={sextract_sel} select=0 edit=0 move=0 rotate=0 delete=1\n"
	catch {$fr marker catalog command ds9 var sextract_sel_reg}
    }

    # Pan to the object (locked band frames follow through the WCS lock)
    if {$pan} {
	PanToFrame $current(frame) $x $y image {}
    }

    set catpanel(status) "Source at image ($x, $y)"
}

# --- Source Extractor Parameter Management ---

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
    OGFMarkDelete $frame

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
	    OGFMarkSend $frame
	    set reg "image\n"
	    set batch_count 0
	}
    }

    # Flush remaining markers
    if {$batch_count > 0} {
	set sextract_all_reg $reg
	OGFMarkSend $frame
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
    OGFMarkDelete $frame
    CatalogPanelCreateAllMarkers
}

proc CatalogPanelClearMarkers {} {
    global catpanel
    global current

    if {$current(frame) == {}} return

    set frame $current(frame)
    OGFMarkDelete $frame
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

    # Select, scroll, summary, markers (shared link path)
    CatalogPanelLinkSelect $num_str replace 1
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
    # time-domain markers (moving / transient / detection): select the table row of that object
    if {![catch {$which get marker catalog id $x $y} _mid] && $_mid != 0} {
	set _tdk [OGFTDKeyFromTags [$which get marker catalog $_mid tag]]
	if {$_tdk ne {}} {CatalogPanelLinkSelect $_tdk replace 1; return}
    }
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

# --- Log Scale ---

# --- AI Merge ---

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

    if {[catch {::ogf::td::active} _tda]} {set _tda 0}
    if {$_tda} {::ogf::td::sort $colname $direction; return}
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

    OGFSessLog catalog.sort auto {} -tool native -title "Sort catalog by $colname $direction" \
	-payload [dict create col $colname dir $direction] -requires catalog
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

# ============================================================================
# AI Star/Galaxy Classification
# ============================================================================

# ============================================================================
# PSF/Deconv procedures
# ============================================================================

# --- Star Finding ---

# --- PSF Generation ---

# --- Extended PSF ---

# --- Simulation PSF: Check Availability ---

# --- Simulation PSF: WebbPSF (JWST) ---

# --- Simulation PSF: TinyTim (HST) ---

# --- Simulation PSF: Shared Helpers ---

# --- Deconvolution ---

# --- Quick Deconvolve (one-click) ---

# --- Settings Dialog ---

# ============================================================================
# Separate — Source Deblending / Splitting
# ============================================================================

# --- Separate Settings Dialog ---

# --- Separate Save/Load Catalog ---

# ============================================================================
# Add Objects — Detect and add source at cursor position (Ctrl+A)
# ============================================================================

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

# --- Helper: Extract FITS base name (without path and extensions) ---

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

# --- Helper: Update ICL intermediate file paths based on FITS name ---

# --- Helper: Update LSBG intermediate file paths based on FITS name ---

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

    OGFSessLog catalog.export_fits auto [list [OGFPython] $script --input $tmpfile --output $fn] -title "Export FITS table [file tail $fn]" -useroutputs [list $fn]
    if {[catch {
	set data [exec [OGFPython] $script --input $tmpfile --output $fn 2>@stderr]
    } err]} {
	set catpanel(status) "FITS export error: $err"
	return
    }

    set catpanel(status) "Exported FITS table to [file tail $fn]"
}

# --- B8: Segmentation Map ---

# --- B10: Non-Parametric Morphology (CAS/Gini/M20) ---

# --- D17: Sérsic Fitting ---

# --- C13: PSF Photometry ---

# --- D16: Multi-Band Photometry ---

# --- D19: Crowded Field Photometry ---

# --- C14: Cross-Match (VizieR) ---

# --- C12: Dual-Image Extract ---

# --- D18: Completeness Simulation ---

# ============================================================================
# ICL (Intra-Cluster Light) Detection Pipeline
# ============================================================================

# ===========================================================
# CLI Script Export / Import Engine (ICL & LSBG)
# ===========================================================

# --- ICL fallback: generate from current settings ---

# --- LSBG fallback: generate --mode run from current settings ---

# --- Export entry points ---

# --- Import Script Engine ---
# Handles both Export-generated scripts ($INPUT/$OUTDIR/$SCRIPT_DIR)
# and standalone reproduce scripts (arbitrary shell variables, line continuations).

# --- 1. Source Masking ---

# --- 2. Background Model ---

# --- 3. BCG Center + Profile ---

# --- 4. ICL Measurements ---

# --- Multi-Threshold ICL ---

# --- BCG+ICL Decomposition ---

# --- Color Profile Dialog ---

# --- Save/Load Profile ---

# --- Settings Dialog ---

# ============================================================
# LSBG (Low Surface Brightness Galaxy) Detection Pipeline
# ============================================================

# --- LSBG Param Load/Save ---

# --- 1. Mask Bright Sources ---

# --- 2. Background Model (Iterative Cleaning) ---

# --- 3. Detect LSBG Candidates ---

# --- 4. Photometry ---

# --- 5. Sérsic Profile Fit ---

# --- 6. Filter + Grade ---

# --- SVM Classify ---

# --- Run Full Pipeline ---

# --- 7. Forced Photometry (Multi-Band) ---

# --- LSBG Settings Dialog ---

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

# ============================================================
# Feature: Interactive Plot (BLT Graph)
# ============================================================

# ============================================================
# Feature: Photo-z (AI Photometric Redshift)
# ============================================================

# ============================================================
# Feature: SED Fitting (AI SPS Emulator)
# ============================================================

# ============================================================
# Feature: Analysis Viewer (standalone GUI)
# ============================================================

# ============================================================
# Feature: Bulge+Disk Decomposition
# ============================================================
