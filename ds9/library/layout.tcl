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

proc CatalogPanelICLUpdateFiles {fn} {
    global catpanel
    if {$fn eq {}} return
    set base [CatalogPanelFitsBaseName $fn]
    if {$base eq {}} return
    # Skip if already set for this base
    if {[info exists catpanel(icl,fits_base)] &&
	$catpanel(icl,fits_base) eq $base} return
    set catpanel(icl,fits_base) $base
    set ds9dir [file join [file normalize ~] .ds9]
    # shared mask (ds9_mask.py): boolean view + on-demand interpolated image
    set catpanel(icl,mask_file)    [file join $ds9dir "mask_${base}_bool.fits"]
    set catpanel(icl,masked_file)  [file join $ds9dir "mask_${base}_masked.fits"]
    set catpanel(icl,bkg_file)     [file join $ds9dir "icl_background_${base}.fits"]
    set catpanel(icl,bgsub_file)   [file join $ds9dir "icl_bgsub_${base}.fits"]
    set catpanel(icl,profile_file) [file join $ds9dir "icl_profile_${base}.tsv"]
    # Detect if previous results exist for this FITS
    set catpanel(icl,has_mask)    [file exists $catpanel(icl,mask_file)]
    set catpanel(icl,has_bkg)     [file exists $catpanel(icl,bkg_file)]
    set catpanel(icl,has_profile) [file exists $catpanel(icl,profile_file)]
}

# --- Helper: Update LSBG intermediate file paths based on FITS name ---

proc CatalogPanelLSBGUpdateFiles {fn} {
    global catpanel
    if {$fn eq {}} return
    set base [CatalogPanelFitsBaseName $fn]
    if {$base eq {}} return
    # Skip if already set for this base
    if {[info exists catpanel(lsbg,fits_base)] &&
	$catpanel(lsbg,fits_base) eq $base} return
    set catpanel(lsbg,fits_base) $base
    set ds9dir [file join [file normalize ~] .ds9]
    set catpanel(lsbg,mask_file)    [file join $ds9dir "mask_${base}_bool.fits"]
    set catpanel(lsbg,masked_file)  [file join $ds9dir "mask_${base}_masked.fits"]
    set catpanel(lsbg,bkg_file)     [file join $ds9dir "lsbg_background_${base}.fits"]
    set catpanel(lsbg,cleaned_file) [file join $ds9dir "lsbg_cleaned_${base}.fits"]
    set catpanel(lsbg,segmap_file)  [file join $ds9dir "lsbg_segmap_${base}.fits"]
    set catpanel(lsbg,catalog_file) [file join $ds9dir "lsbg_catalog_${base}.tsv"]
    # Detect if previous results exist for this FITS
    set catpanel(lsbg,has_mask)    [file exists $catpanel(lsbg,mask_file)]
    set catpanel(lsbg,has_clean)   [file exists $catpanel(lsbg,cleaned_file)]
    set catpanel(lsbg,has_detect)  [file exists $catpanel(lsbg,segmap_file)]
    set catpanel(lsbg,has_catalog) [file exists $catpanel(lsbg,catalog_file)]
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
    foreach pname {expand-factor max-dilate-radius \
		   bright-star-mag-limit bright-star-radius-scale \
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

# ===========================================================
# CLI Script Export / Import Engine (ICL & LSBG)
# ===========================================================

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

# --- ICL fallback: generate from current settings ---

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

# --- LSBG fallback: generate --mode run from current settings ---

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

# --- Export entry points ---

proc CatalogPanelICLExportScript {} {
    CatalogPanelExportCLIScript icl
}

proc CatalogPanelLSBGExportScript {} {
    CatalogPanelExportCLIScript lsbg
}

# --- Import Script Engine ---
# Handles both Export-generated scripts ($INPUT/$OUTDIR/$SCRIPT_DIR)
# and standalone reproduce scripts (arbitrary shell variables, line continuations).

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

# --- 1. Source Masking ---

proc CatalogPanelICLMask {} {
    OGFMaskPipelineMask icl
}

proc CatalogPanelICLViewMask {} {
    global ogfmask
    set ogfmask(overlay) 1
    CatalogPanelMaskToggleOverlay
}

proc CatalogPanelICLSaveMask {} {
    CatalogPanelMaskSaveAs
}

proc CatalogPanelICLImportMask {} {
    CatalogPanelMaskImport
}

# --- 2. Background Model ---

proc CatalogPanelICLBackground {method} {
    global catpanel

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "ICL: No FITS file loaded"
	return
    }
    CatalogPanelICLUpdateFiles $fn

    # Shared mask -> interpolated image produced now (on demand); raw if no mask
    set input $fn
    if {[OGFMaskExists]} {
	set input [OGFMaskMaskedFor icl $fn]
    }

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Fitting background ($method)..."
    update idletasks

    set args [list [OGFPython] $script $input --mode background \
	--bkg-method $method \
	--bkg-order $catpanel(icl,param,bkg-order) \
	--bkg-sigma-clip $catpanel(icl,param,bkg-sigma-clip) \
	--bkg-sep-mesh $catpanel(icl,param,bkg-sep-mesh) \
	--bkg-output $catpanel(icl,bkg_file) \
	--bgsub-output $catpanel(icl,bgsub_file)]

    if {[OGFMaskExists]} {
	lappend args --mask [OGFMaskBoolPath]
    }

    # Iterative background refinement
    if {$catpanel(icl,param,bkg-iterative)} {
	lappend args --iterative \
	    --interp-method $catpanel(icl,param,interp-method) \
	    --bkg-n-iterations $catpanel(icl,param,bkg-n-iterations) \
	    --bkg-convergence-tol $catpanel(icl,param,bkg-convergence-tol) \
	    --bkg-refine-thresh $catpanel(icl,param,bkg-refine-thresh) \
	    --mask-output [OGFMaskRefinedPath icl]
    }

    CatalogPanelCmdLog icl $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL background error: $err"
	return
    }

    set catpanel(icl,has_bkg) 1
    if {$catpanel(icl,param,bkg-iterative) && [file exists [OGFMaskRefinedPath icl]]} {
	# iterative refinement result feeds the profile step; shared mask untouched
	set catpanel(icl,mask_file) [OGFMaskRefinedPath icl]
    }

    # Auto-display bgsub image in new frame.  (CreateFrame resets the
    # per-frame panel state, so grab the path first.)
    set _bgsub $catpanel(icl,bgsub_file)
    if {[file exists $_bgsub]} {
	CreateFrame
	if {![catch {LoadFitsFile $_bgsub {} {}}]} {
	    global scale
	    set scale(mode) zscale
	    ChangeScaleMode
	}
    }

    set catpanel(status) "ICL: Background model ($method) complete"
}

proc CatalogPanelICLViewBkg {} {
    global catpanel

    if {![file exists $catpanel(icl,bgsub_file)]} {
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
    set catpanel(icl,has_bkg) 1
    set catpanel(status) "ICL: Background-subtracted image loaded in new frame"
}

# --- 3. BCG Center + Profile ---

proc CatalogPanelICLSetCenter {} {
    global catpanel catpanel_fdata ds9 current

    # Use the first frame's catalog (original image SExtractor result)
    set catalog_data {}
    set first_frame [lindex $ds9(frames) 0]
    if {$first_frame ne {} &&
	[info exists catpanel_fdata($first_frame,alldata)] &&
	$catpanel_fdata($first_frame,alldata) ne {}} {
	set catalog_data $catpanel_fdata($first_frame,alldata)
    } elseif {[info exists catpanel(alldata)] && $catpanel(alldata) ne {}} {
	set catalog_data $catpanel(alldata)
    }

    if {$catalog_data eq {}} {
	set catpanel(status) "ICL: No SExtractor catalog found — run SExtract first"
	return
    }

    # Parse headers to find NUMBER, X_IMAGE, Y_IMAGE, MAG_AUTO
    set lines [split $catalog_data \n]
    set headers [split [lindex $lines 0] "\t"]
    set numcol -1
    set xcol -1
    set ycol -1
    set magcol -1
    for {set c 0} {$c < [llength $headers]} {incr c} {
	set h [string trim [lindex $headers $c]]
	if {$h eq "NUMBER"} { set numcol $c }
	if {$h eq "X_IMAGE"} { set xcol $c }
	if {$h eq "Y_IMAGE"} { set ycol $c }
	if {$h eq "MAG_AUTO"} { set magcol $c }
    }
    if {$numcol < 0 || $xcol < 0 || $ycol < 0} {
	set catpanel(status) "ICL: Catalog missing NUMBER/X_IMAGE/Y_IMAGE"
	return
    }

    # Find brightest source as default suggestion
    set best_id {}
    set best_mag 99.0
    for {set i 1} {$i < [llength $lines]} {incr i} {
	set row [split [lindex $lines $i] "\t"]
	if {[llength $row] <= $xcol} continue
	if {$magcol >= 0 && [llength $row] > $magcol} {
	    set m [string trim [lindex $row $magcol]]
	    if {[string is double $m] && $m < $best_mag} {
		set best_mag $m
		set best_id [string trim [lindex $row $numcol]]
	    }
	}
    }

    # Show dialog
    set w .icl_bcg_id
    if {[winfo exists $w]} { destroy $w }
    toplevel $w
    wm title $w "Set BCG Center"
    wm geometry $w 320x120
    wm resizable $w 0 0

    set default_id $best_id
    if {$default_id eq {}} { set default_id 1 }

    ttk::label $w.lbl -text "Enter source NUMBER from SExtractor catalog:"
    pack $w.lbl -padx 10 -pady {10 2} -anchor w

    ttk::frame $w.ef
    pack $w.ef -padx 10 -pady 2 -fill x
    ttk::label $w.ef.idlbl -text "Source ID:"
    ttk::entry $w.ef.entry -width 10 -textvariable icl_bcg_entry_id
    set ::icl_bcg_entry_id $default_id
    pack $w.ef.idlbl -side left -padx {0 5}
    pack $w.ef.entry -side left

    if {$best_id ne {} && $magcol >= 0} {
	ttk::label $w.ef.hint -text "(brightest: #${best_id}, mag=[format %.1f $best_mag])" \
	    -foreground gray
	pack $w.ef.hint -side left -padx 5
    }

    ttk::frame $w.bf
    pack $w.bf -padx 10 -pady 8 -fill x
    ttk::button $w.bf.ok -text "OK" \
	-command [list CatalogPanelICLSetCenterByID $w [list $catalog_data] $numcol $xcol $ycol]
    ttk::button $w.bf.click -text "Click on Image" \
	-command [list CatalogPanelICLSetCenterClickMode $w]
    ttk::button $w.bf.cancel -text "Cancel" -command [list destroy $w]
    pack $w.bf.cancel -side right -padx 2
    pack $w.bf.click -side right -padx 2
    pack $w.bf.ok -side right -padx 2

    bind $w.ef.entry <Return> [list $w.bf.ok invoke]
    focus $w.ef.entry
    $w.ef.entry selection range 0 end
}

proc CatalogPanelICLSetCenterByID {w catalog_data numcol xcol ycol} {
    global catpanel current

    set src_id [string trim $::icl_bcg_entry_id]
    if {$src_id eq {} || ![string is integer $src_id]} {
	set catpanel(status) "ICL: Enter a valid source NUMBER"
	return
    }

    # Search catalog for this NUMBER
    set lines [split $catalog_data \n]
    for {set i 1} {$i < [llength $lines]} {incr i} {
	set row [split [lindex $lines $i] "\t"]
	if {[llength $row] <= $xcol || [llength $row] <= $ycol} continue
	set n [string trim [lindex $row $numcol]]
	if {$n eq $src_id} {
	    set ix [string trim [lindex $row $xcol]]
	    set iy [string trim [lindex $row $ycol]]

	    # Store as 0-indexed
	    set catpanel(icl,center_x) [expr {$ix - 1.0}]
	    set catpanel(icl,center_y) [expr {$iy - 1.0}]

	    # Draw cyan cross
	    set frame $current(frame)
	    if {$frame ne {}} {
		catch {$frame marker catalog icl_bcg delete}
		global icl_bcg_reg
		set icl_bcg_reg "image\ncross point([format %.1f $ix] [format %.1f $iy]) # color=cyan width=2 point=cross 20 tag={icl_bcg} select=0 edit=0 move=0 rotate=0 delete=1\n"
		catch {$frame marker catalog command ds9 var icl_bcg_reg}
	    }

	    set catpanel(status) "ICL: BCG center set to source #$src_id ([format %.1f $ix], [format %.1f $iy])"
	    destroy $w
	    return
	}
    }

    set catpanel(status) "ICL: Source #$src_id not found in catalog"
}

proc CatalogPanelICLSetCenterClickMode {w} {
    global catpanel
    destroy $w
    set catpanel(icl,click_mode) 1
    set catpanel(status) "ICL: Click on image to set BCG center..."
}

proc CatalogPanelICLClickSetCenter {frame x y} {
    global catpanel

    # Convert canvas coords to image coords (1-based)
    set imgc [$frame get coordinates $x $y image]
    set ix [lindex $imgc 0]
    set iy [lindex $imgc 1]

    # Store as 0-indexed (Python convention)
    set catpanel(icl,center_x) [expr {$ix - 1.0}]
    set catpanel(icl,center_y) [expr {$iy - 1.0}]

    # Draw a cyan cross at BCG center
    catch {$frame marker catalog icl_bcg delete}

    global icl_bcg_reg
    set icl_bcg_reg "image\ncross point([format %.1f $ix] [format %.1f $iy]) # color=cyan width=2 point=cross 20 tag={icl_bcg} select=0 edit=0 move=0 rotate=0 delete=1\n"
    catch {$frame marker catalog command ds9 var icl_bcg_reg}

    set catpanel(icl,click_mode) 0
    set catpanel(status) "ICL: BCG center set to ([format %.1f $ix], [format %.1f $iy])"
}

proc CatalogPanelICLProfile {} {
    global catpanel

    if {$catpanel(icl,center_x) eq {} || $catpanel(icl,center_y) eq {}} {
	set catpanel(status) "ICL: Set BCG center first"
	return
    }

    # Use bgsub image if available, otherwise raw
    set fn [CatalogPanelGetFITS]
    if {[file exists $catpanel(icl,bgsub_file)]} {
	set fn $catpanel(icl,bgsub_file)
    }
    if {$fn eq {}} {
	set catpanel(status) "ICL: No image available"
	return
    }
    CatalogPanelICLUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Measuring SB profile..."
    update idletasks

    set center "$catpanel(icl,center_x),$catpanel(icl,center_y)"

    set args [list [OGFPython] $script $fn --mode profile \
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

    CatalogPanelCmdLog icl $args
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
    global icl_ann_reg
    set icl_ann_reg "image\ncross point($cx $cy) # color=cyan width=2 point=cross 20 tag={icl_annulus} select=0 edit=0 move=0 rotate=0 delete=1\n"
    catch {$frame marker catalog command ds9 var icl_ann_reg}

    # Draw annulus rings at 25%, 50%, 75%, 100% of rmax
    set rmin $catpanel(icl,param,rmin)
    set rmax $catpanel(icl,param,rmax)
    foreach frac {0.25 0.50 0.75 1.00} {
	set r [expr {$rmin + ($rmax - $rmin) * $frac}]
	set ri [expr {int($r)}]
	set icl_ann_reg "image\ncircle($cx $cy ${r}i) # color=green dash=1 width=1 tag={icl_annulus} select=0 edit=0 move=0 rotate=0 delete=1 text={r=${ri}}\n"
	catch {$frame marker catalog command ds9 var icl_ann_reg}
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
    if {[file exists $catpanel(icl,bgsub_file)]} {
	set fn $catpanel(icl,bgsub_file)
    }
    if {$fn eq {}} {
	set catpanel(status) "ICL: No image available"
	return
    }
    CatalogPanelICLUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Measuring sector profile..."
    update idletasks

    set center "$catpanel(icl,center_x),$catpanel(icl,center_y)"

    set args [list [OGFPython] $script $fn --mode profile \
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

    CatalogPanelCmdLog icl $args
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
    CatalogPanelICLUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Computing ICL measurements..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode measure \
	--profile-file $catpanel(icl,profile_file) \
	--mu-threshold $catpanel(icl,param,mu-threshold) \
	--mu-levels $catpanel(icl,param,mu-levels) \
	--pixel-scale $catpanel(icl,param,pixel-scale)]

    CatalogPanelCmdLog icl $args
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
    global icl_ann_reg
    set icl_ann_reg "image\ncross point($cx $cy) # color=cyan width=2 point=cross 20 tag={icl_annulus} select=0 edit=0 move=0 rotate=0 delete=1\n"
    catch {$frame marker catalog command ds9 var icl_ann_reg}

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
		    set icl_ann_reg "image\ncircle($cx $cy ${r}i) # color=magenta dash=1 width=1 tag={icl_annulus} select=0 edit=0 move=0 rotate=0 delete=1 text={mu=$mu_label}\n"
		    catch {$frame marker catalog command ds9 var icl_ann_reg}
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
    CatalogPanelICLUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Computing multi-threshold ICL..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode measure-multi \
	--profile-file $catpanel(icl,profile_file) \
	--pixel-scale $catpanel(icl,param,pixel-scale)]

    CatalogPanelCmdLog icl $args
    if {[catch {set data [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "ICL multi-threshold error: $err"
	return
    }

    CatalogPanelLoadTSV $data "icl_multi_threshold"
    set catpanel(status) "ICL: Multi-threshold measurements complete"
}

# --- BCG+ICL Decomposition ---

proc CatalogPanelICLDecompose {} {
    global catpanel ds9 catpanel_fdata

    if {$catpanel(icl,center_x) eq {} || $catpanel(icl,center_y) eq {}} {
	set catpanel(status) "ICL: Set BCG center first"
	return
    }

    # Use the ORIGINAL image (first frame) for decompose —
    # the masked image has BCG removed, making Sérsic fitting impossible
    set fn {}
    set first_frame [lindex $ds9(frames) 0]
    if {$first_frame ne {}} {
	if {[info exists catpanel_fdata($first_frame,filename)] &&
	    $catpanel_fdata($first_frame,filename) ne {}} {
	    set fn $catpanel_fdata($first_frame,filename)
	}
    }
    if {$fn eq {}} {
	set fn [CatalogPanelGetFITS]
    }
    if {$fn eq {}} {
	set catpanel(status) "ICL: No FITS image available"
	return
    }
    CatalogPanelICLUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: BCG+ICL decomposition (original image)..."
    update idletasks

    # Measure profile on original image and decompose in one step
    set center "$catpanel(icl,center_x),$catpanel(icl,center_y)"

    # First: measure profile on original (unmasked) image
    set prof_args [list [OGFPython] $script $fn --mode profile \
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

    CatalogPanelCmdLog icl $prof_args
    if {[catch {exec {*}$prof_args 2>@stderr} err]} {
	set catpanel(status) "ICL decompose: profile error: $err"
	return
    }

    # Then: decompose using that profile
    set args [list [OGFPython] $script $fn --mode decompose \
	--profile-file $catpanel(icl,profile_file) \
	--pixel-scale $catpanel(icl,param,pixel-scale) \
	--mag-zeropoint $catpanel(icl,param,mag-zeropoint)]

    CatalogPanelCmdLog icl $args
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
    CatalogPanelICLUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_icl.py]
    if {![file exists $script]} {
	set catpanel(status) "ICL: ds9_icl.py not found"
	return
    }

    set catpanel(status) "ICL: Measuring color profile..."
    update idletasks

    set center "$catpanel(icl,center_x),$catpanel(icl,center_y)"

    set args [list [OGFPython] $script $fn --mode color \
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

    catch {OGFSessFromCmdLog icl $args}
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
	-initialfile [file tail $catpanel(icl,profile_file)]]
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
    foreach pname {expand-factor max-dilate-radius \
		   bright-star-mag-limit bright-star-radius-scale \
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

    ttk::label $t1.lmdr -text "Max dilate radius (px):"
    ttk::entry $t1.emdr -textvariable ed(icl,max-dilate-radius) -width 10
    grid $t1.lmdr -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.emdr -row $r -column 1 -sticky w -padx 4 -pady 4
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

    foreach pname {expand-factor max-dilate-radius \
		   bright-star-mag-limit bright-star-radius-scale \
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
    set ed(icl,max-dilate-radius) 50
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
		   max-dilate-radius \
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
		   sersic-n-filter-min sersic-n-filter-max sersic-chi2-max \
		   svm-classify svm-threshold svm-checkpoint} {
	puts $fd "$pname $catpanel(lsbg,param,$pname)"
    }
    close $fd
}

# --- 1. Mask Bright Sources ---

proc CatalogPanelLSBGMask {} {
    OGFMaskPipelineMask lsbg
}

proc CatalogPanelLSBGViewMask {} {
    global ogfmask
    set ogfmask(overlay) 1
    CatalogPanelMaskToggleOverlay
}

proc CatalogPanelLSBGSaveMask {} {
    CatalogPanelMaskSaveAs
}

proc CatalogPanelLSBGImportMask {} {
    CatalogPanelMaskImport
}

# --- 2. Background Model (Iterative Cleaning) ---

proc CatalogPanelLSBGClean {method} {
    global catpanel

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    if {![OGFMaskEnsure lsbg]} {
	set catpanel(status) "LSBG: could not obtain a mask - run Mask > Auto Mask first"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Iterative background ($method)..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode clean \
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

    CatalogPanelCmdLog lsbg $args
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

    if {![file exists $catpanel(lsbg,cleaned_file)]} {
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

    if {![file exists $catpanel(lsbg,cleaned_file)]} {
	set catpanel(status) "LSBG: No cleaned image — run Background Model first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Detecting candidates..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode detect \
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

    CatalogPanelCmdLog lsbg $args
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
    CatalogPanelLoadTSV $catpanel(alldata) "lsbg"

    # Count detections
    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]
    set catpanel(status) "LSBG: $nsrc candidates detected (segmap in new frame)"
}

# --- 4. Photometry ---

proc CatalogPanelLSBGPhotometry {} {
    global catpanel

    if {![file exists $catpanel(lsbg,cleaned_file)]} {
	set catpanel(status) "LSBG: No cleaned image — run Background Model first"
	return
    }

    if {![file exists $catpanel(lsbg,segmap_file)]} {
	set catpanel(status) "LSBG: No detections — run Detect first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Measuring photometry..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode photometry \
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

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG photometry error: $err"
	return
    }

    set catpanel(alldata) $result
    set catpanel(lsbg,has_catalog) 1
    CatalogPanelLoadTSV $catpanel(alldata) "lsbg"

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    CatalogPanelMarkAll
    set catpanel(status) "LSBG: $nsrc sources — photometry done"
}

# --- 5. Sérsic Profile Fit ---

proc CatalogPanelLSBGSersic {} {
    global catpanel

    if {![file exists $catpanel(lsbg,cleaned_file)]} {
	set catpanel(status) "LSBG: No cleaned image — run Background Model first"
	return
    }

    if {![file exists $catpanel(lsbg,segmap_file)]} {
	set catpanel(status) "LSBG: No detections — run Detect first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Fitting Sérsic profiles..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode sersic \
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

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG Sérsic fit error: $err"
	return
    }

    set catpanel(alldata) $result
    set catpanel(lsbg,has_catalog) 1
    CatalogPanelLoadTSV $catpanel(alldata) "lsbg"

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    CatalogPanelMarkAll
    set catpanel(status) "LSBG: $nsrc sources — Sérsic fit done"
}

# --- 6. Filter + Grade ---

proc CatalogPanelLSBGFilter {} {
    global catpanel

    if {![file exists $catpanel(lsbg,cleaned_file)]} {
	set catpanel(status) "LSBG: No cleaned image — run Background Model first"
	return
    }

    if {![file exists $catpanel(lsbg,segmap_file)]} {
	set catpanel(status) "LSBG: No detections — run Detect first"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    set catpanel(status) "LSBG: Filtering + grading candidates..."
    update idletasks

    set args [list [OGFPython] $script $fn --mode filter \
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

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG filter error: $err"
	return
    }

    set catpanel(alldata) $result
    set catpanel(lsbg,has_catalog) 1
    CatalogPanelLoadTSV $catpanel(alldata) "lsbg"

    set nlines [llength [split $result \n]]
    set nsrc [expr {$nlines - 1}]

    CatalogPanelMarkAll
    set catpanel(status) "LSBG: $nsrc candidates passed filtering"
}

# --- SVM Classify ---

proc CatalogPanelLSBGSVMClassify {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "LSBG SVM: No catalog loaded"
	return
    }

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG SVM: ds9_lsbg.py not found"
	return
    }

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG SVM: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set catpanel(status) "LSBG SVM: Classifying candidates..."
    update idletasks

    # Save current catalog to temp file
    set tmpcat [file join [file normalize ~] .ds9 lsbg_svm_tmp.tsv]
    if {[catch {
	set fd [open $tmpcat w]
	puts $fd $catpanel(alldata)
	close $fd
    } err]} {
	set catpanel(status) "LSBG SVM: Cannot write temp catalog"
	return
    }

    set args [list [OGFPython] $script $fn --mode svm-classify \
	--catalog $tmpcat \
	--svm-threshold $catpanel(lsbg,param,svm-threshold)]
    if {$catpanel(lsbg,param,svm-checkpoint) ne {}} {
	lappend args --svm-checkpoint $catpanel(lsbg,param,svm-checkpoint)
    }

    CatalogPanelCmdLog lsbg $args
    if {[catch {set result [exec {*}$args 2>@stderr]} err]} {
	set catpanel(status) "LSBG SVM error: $err"
	catch {file delete $tmpcat}
	return
    }
    catch {file delete $tmpcat}

    # Add SVM columns to catalog
    CatalogPanelAddColumnsFromTSV $result

    # Color markers by SVM result
    CatalogPanelLSBGSVMColorMarkers $result

    # Count results
    set lines [split $result \n]
    set n_lsbg 0
    set n_total 0
    foreach line [lrange $lines 1 end] {
	if {[string trim $line] eq {}} continue
	incr n_total
	set fields [split $line "\t"]
	if {[llength $fields] >= 2 && [lindex $fields 1] eq "1"} {
	    incr n_lsbg
	}
    }
    set catpanel(status) "LSBG SVM: $n_lsbg LSBG / $n_total total"
}

proc CatalogPanelLSBGSVMColorMarkers {result_tsv} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    set frame $current(frame)

    # Parse SVM results: NUMBER -> SVM_LSBG
    array set svm_class {}
    set lines [split $result_tsv \n]
    foreach line [lrange $lines 1 end] {
	set fields [split $line "\t"]
	if {[llength $fields] < 2} continue
	set src_num [lindex $fields 0]
	set svm_lsbg [lindex $fields 1]
	set svm_class($src_num) $svm_lsbg
    }

    # Delete existing markers and recreate with SVM colors
    OGFMarkDelete $frame

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set alllines [split $catpanel(alldata) \n]
    if {[llength $alllines] < 2} return

    set headers [split [lindex $alllines 0] "\t"]
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

    # Build region strings with SVM colors
    set batch_size 500
    set reg "image\n"
    set count 0
    set batch_count 0
    global sextract_all_reg

    for {set i 1} {$i < [llength $alllines]} {incr i} {
	set line [lindex $alllines $i]
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

	# Color: LSBG=green, contaminant=red, unclassified=yellow
	set color yellow
	if {[info exists svm_class($src_num)]} {
	    if {$svm_class($src_num) eq "1"} {
		set color green
	    } else {
		set color red
	    }
	}

	append reg "ellipse($x $y ${semi_a}i ${semi_b}i $theta) # color=$color width=1 tag={sextract_all} tag={sextract_src.$src_num} select=0 edit=0 move=0 rotate=0 delete=1 highlite=1 callback=highlite CatalogPanelMarkerCB {$src_num} callback=unhighlite CatalogPanelMarkerUnCB {$src_num}\n"

	incr count
	if {$count >= $batch_size} {
	    set sextract_all_reg $reg
	    OGFMarkSend $frame
	    set reg "image\n"
	    set count 0
	    incr batch_count
	}
    }

    # Flush remaining
    if {$count > 0} {
	set sextract_all_reg $reg
	OGFMarkSend $frame
    }
}

# --- Run Full Pipeline ---

proc CatalogPanelLSBGRunAll {} {
    global catpanel

    set fn [CatalogPanelGetFITS]
    if {$fn eq {}} {
	set catpanel(status) "LSBG: No FITS file loaded"
	return
    }
    CatalogPanelLSBGUpdateFiles $fn

    set script [CatalogPanelGetScript ds9_lsbg.py]
    if {![file exists $script]} {
	set catpanel(status) "LSBG: ds9_lsbg.py not found"
	return
    }

    OGFMaskEnsure lsbg
    set catpanel(status) "LSBG: Running full pipeline..."
    update idletasks

    # Reset cmdlog for new session
    set catpanel(lsbg,cmdlog) {}

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

    CatalogPanelCmdLog lsbg $args
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
    CatalogPanelLoadTSV $catpanel(alldata) "lsbg"

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

    set args [list [OGFPython] $script $band_fits --mode forced \
	--catalog $tmpcat \
	--band-name $band_name \
	--mag-zeropoint $catpanel(lsbg,param,mag-zeropoint) \
	--pixel-scale $catpanel(lsbg,param,pixel-scale) \
	--n-workers $catpanel(param,n-workers)]

    CatalogPanelCmdLog lsbg $args
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
		   sersic-n-filter-min sersic-n-filter-max sersic-chi2-max \
		   svm-classify svm-threshold svm-checkpoint} {
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

    ttk::label $t1.lmdr -text "Max dilate radius (px):"
    ttk::entry $t1.emdr -textvariable ed(lsbg,max-dilate-radius) -width 10
    grid $t1.lmdr -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t1.emdr -row $r -column 1 -sticky w -padx 4 -pady 4
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
    incr r

    ttk::separator $t4.sep2 -orient horizontal
    grid $t4.sep2 -row $r -column 0 -columnspan 2 -sticky ew -padx 8 -pady 6
    incr r

    ttk::checkbutton $t4.csvm -text "SVM Classification" \
	-variable ed(lsbg,svm-classify)
    grid $t4.csvm -row $r -column 0 -columnspan 2 -sticky w -padx 8 -pady 4
    incr r

    ttk::label $t4.lsvmt -text "SVM threshold:"
    ttk::entry $t4.esvmt -textvariable ed(lsbg,svm-threshold) -width 10
    grid $t4.lsvmt -row $r -column 0 -sticky w -padx 8 -pady 4
    grid $t4.esvmt -row $r -column 1 -sticky w -padx 4 -pady 4

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
		   max-dilate-radius \
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
		   sersic-n-filter-min sersic-n-filter-max sersic-chi2-max \
		   svm-classify svm-threshold svm-checkpoint} {
	set catpanel(lsbg,param,$pname) $ed(lsbg,$pname)
    }
    CatalogPanelLSBGParamSave
    set catpanel(status) "LSBG settings applied and saved"
}

proc CatalogPanelLSBGSettingsDefaults {} {
    global ed

    set ed(lsbg,mask-detect-thresh) 1.5
    set ed(lsbg,mask-detect-minarea) 5
    set ed(lsbg,mask-expand-factor) 1.5
    set ed(lsbg,max-dilate-radius) 30
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
    set ed(lsbg,svm-classify) 0
    set ed(lsbg,svm-threshold) 0.5
    set ed(lsbg,svm-checkpoint) {}
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
