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

    # Title bar
    set catpanel(titlebar) [ttk::frame $f.titlebar]
    ttk::label $f.titlebar.title -text "Source Extractor" \
	-font {Helvetica 11 bold} -anchor w
    ttk::button $f.titlebar.extract -text "Extract" \
	-command CatalogPanelExtract -width 8
    ttk::button $f.titlebar.settings -text "Settings..." \
	-command CatalogPanelSettingsDialog -width 10
    ttk::button $f.titlebar.markall -text "Mark All" \
	-command CatalogPanelMarkAll -width 8
    ttk::button $f.titlebar.visible -text "Visible" \
	-command CatalogPanelShowVisible -width 8
    ttk::button $f.titlebar.trim -text "Trim..." \
	-command CatalogPanelTrimDialog -width 6
    ttk::button $f.titlebar.save -text "Save" \
	-command CatalogPanelSaveCatalog -width 6
    ttk::button $f.titlebar.load -text "Load" \
	-command CatalogPanelLoadCatalog -width 6
    ttk::button $f.titlebar.clear -text "Clear" \
	-command CatalogPanelClear -width 6
    pack $f.titlebar.title -side left -padx 4 -pady 2
    pack $f.titlebar.clear -side right -padx 2 -pady 2
    pack $f.titlebar.load -side right -padx 2 -pady 2
    pack $f.titlebar.save -side right -padx 2 -pady 2
    pack $f.titlebar.trim -side right -padx 2 -pady 2
    pack $f.titlebar.visible -side right -padx 2 -pady 2
    pack $f.titlebar.markall -side right -padx 2 -pady 2
    pack $f.titlebar.settings -side right -padx 2 -pady 2
    pack $f.titlebar.extract -side right -padx 2 -pady 2

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
    pack $f.titlebar -fill x -side top
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

    # Feature C: merge state
    set catpanel(merge,list) {}
    set catpanel(merge,active) 0

    # Feature D: trim state
    set catpanel(trim,active) 0

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

    # Initialize extraction parameters
    CatalogPanelParamDef
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

    if {![file exists $fn]} {
	set catpanel(status) "File not found: $fn"
	return
    }

    set catpanel(status) "Extracting sources from [file tail $fn] ..."
    update idletasks

    # Build parameter arguments list
    set paramargs {}
    foreach pname {detect-thresh detect-minarea deblend-nthresh deblend-mincont \
		   phot-aperture mag-zeropoint gain pixel-scale seeing-fwhm \
		   back-size back-filtersize} {
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

    # Reset visible mode
    set catpanel(visible_mode) 0

    # Reset trim state
    set catpanel(trim,active) 0
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

    # Compute ellipse: ISO_RADIUS as semi-major, scaled by B/A for semi-minor
    set semi_a $iso_radius
    set semi_b $iso_radius
    if {$a_image > 0 && $b_image > 0} {
	set semi_b [expr {$iso_radius * $b_image / $a_image}]
    }

    # Delete previous selection markers
    set frame $current(frame)
    catch {$frame marker catalog sextract_sel delete}

    # Use global variable for marker creation (var form requires global access)
    global sextract_sel_reg

    # Create cross point marker (cyan)
    set sextract_sel_reg "image\ncross point($x $y) # color=cyan width=2 point=cross 15 tag={sextract_sel} select=0 edit=0 move=0 rotate=0 delete=1\n"
    catch {$frame marker catalog command ds9 var sextract_sel_reg}

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
		   back-size back-filtersize} {
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
}

proc CatalogPanelSettingsDialog {} {
    global catpanel
    global ed

    set w {.sextractparam}

    set ed(ok) 0

    # Copy current params to ed()
    foreach pname {detect-thresh detect-minarea deblend-nthresh deblend-mincont \
		   phot-aperture mag-zeropoint gain pixel-scale seeing-fwhm \
		   back-size back-filtersize} {
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
	phot-aperture {Phot Aperture}
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

    # Buttons
    set bf [ttk::frame $w.buttons]
    ttk::button $bf.ok -text {OK} -command {set ed(ok) 1} -default active
    ttk::button $bf.cancel -text {Cancel} -command {set ed(ok) 0}
    ttk::button $bf.defaults -text {Defaults} -command CatalogPanelParamDefaults
    ttk::button $bf.save -text {Save} -command {
	foreach pname {detect-thresh detect-minarea deblend-nthresh deblend-mincont \
		       phot-aperture mag-zeropoint gain pixel-scale seeing-fwhm \
		       back-size back-filtersize} {
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
		       back-size back-filtersize} {
	    set catpanel(param,$pname) $ed($pname)
	}
	CatalogPanelParamSave
    }

    unset ed
}

# --- Mark All Sources ---

proc CatalogPanelMarkAll {} {
    global catpanel
    global current

    if {$current(frame) == {}} return
    if {![$current(frame) has fits]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return

    set frame $current(frame)

    # Delete previous "mark all" markers
    catch {$frame marker catalog sextract_all delete}

    global $catpanel(tbldb)

    # Find column indices
    set ncols [$catpanel(tbl) cget -cols]
    set col_x -1
    set col_y -1
    set col_a -1
    set col_b -1
    set col_theta -1
    set col_ir -1
    set col_num -1
    for {set c 1} {$c <= $ncols} {incr c} {
	if {[info exists ${catpanel(tbldb)}(0,$c)]} {
	    set hdr [set ${catpanel(tbldb)}(0,$c)]
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
    }
    if {$col_x < 0 || $col_y < 0} return

    set catpanel(status) "Marking all sources..."
    update idletasks

    # Build one big region string for all sources
    set reg "image\n"
    set nrows [$catpanel(tbl) cget -rows]
    set count 0

    for {set r 1} {$r < $nrows} {incr r} {
	if {![info exists ${catpanel(tbldb)}($r,$col_x)]} continue
	set x [set ${catpanel(tbldb)}($r,$col_x)]
	set y [set ${catpanel(tbldb)}($r,$col_y)]
	if {![string is double -strict $x] || ![string is double -strict $y]} continue

	# Get source NUMBER for individual tag
	set src_num $r
	if {$col_num >= 0 && [info exists ${catpanel(tbldb)}($r,$col_num)]} {
	    set src_num [set ${catpanel(tbldb)}($r,$col_num)]
	}

	# Get ellipse parameters (NaN/Inf safe via catch)
	set iso_radius 5.0
	set a_image 0
	set b_image 0
	set theta 0

	if {$col_ir >= 0 && [info exists ${catpanel(tbldb)}($r,$col_ir)]} {
	    set val [set ${catpanel(tbldb)}($r,$col_ir)]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} {
		set iso_radius $v
	    }
	}
	if {$col_a >= 0 && [info exists ${catpanel(tbldb)}($r,$col_a)]} {
	    set val [set ${catpanel(tbldb)}($r,$col_a)]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set a_image $v }
	}
	if {$col_b >= 0 && [info exists ${catpanel(tbldb)}($r,$col_b)]} {
	    set val [set ${catpanel(tbldb)}($r,$col_b)]
	    if {[catch {set v [expr {$val + 0.0}]}] == 0 && $v > 0} { set b_image $v }
	}
	if {$col_theta >= 0 && [info exists ${catpanel(tbldb)}($r,$col_theta)]} {
	    set val [set ${catpanel(tbldb)}($r,$col_theta)]
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
    }

    if {$count == 0} return

    # Create all markers in one call using global variable
    global sextract_all_reg
    set sextract_all_reg $reg
    catch {$frame marker catalog command ds9 var sextract_all_reg}

    set catpanel(status) "Marked $count sources (yellow ellipses)"
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
proc CatalogPanelMarkerClick {which x y} {
    global catpanel

    if {![info exists catpanel(tbl)]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return
    if {![$which has fits]} return

    # Check if click is on a catalog marker
    set id [$which get marker catalog id $x $y]
    if {$id == 0} return

    # Get tags of this marker
    set tags [$which get marker catalog $id tag]

    # Look for sextract_src.NUMBER tag
    set src_num {}
    foreach tag $tags {
	if {[string match "sextract_src.*" $tag]} {
	    set src_num [string range $tag 13 end]
	    break
	}
    }
    if {$src_num eq {}} return

    # Found a sextract marker — navigate to it
    CatalogPanelMarkerCB $src_num $id
}

# Ctrl+Click handler called from ControlButton1Frame in none mode
proc CatalogPanelMarkerCtrlClick {which x y} {
    global catpanel

    if {![info exists catpanel(tbl)]} return
    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} return
    if {![$which has fits]} return

    # Check if click is on a catalog marker
    set id [$which get marker catalog id $x $y]
    if {$id == 0} return

    # Get tags of this marker
    set tags [$which get marker catalog $id tag]

    # Look for sextract_src.NUMBER tag
    set src_num {}
    foreach tag $tags {
	if {[string match "sextract_src.*" $tag]} {
	    set src_num [string range $tag 13 end]
	    break
	}
    }
    if {$src_num eq {}} return

    # Ctrl+Click: merge selection
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

    # Toggle mode
    if {$catpanel(visible_mode)} {
	set catpanel(visible_mode) 0
	CatalogPanelLoadTSV $catpanel(alldata) "all"
	set catpanel(status) "Showing all sources"
	return
    }

    set catpanel(visible_mode) 1

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
    # Flux-weighted centroid for X,Y
    set total_flux 0.0
    set wx 0.0
    set wy 0.0
    set total_npix 0
    set wa 0.0
    set wb 0.0
    set wtheta 0.0

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
	if {$idx_a >= 0} {
	    set av [string trim [lindex $row $idx_a]]
	    if {[string is double -strict $av]} {
		set wa [expr {$wa + $av * $flux}]
	    }
	}
	if {$idx_b >= 0} {
	    set bv [string trim [lindex $row $idx_b]]
	    if {[string is double -strict $bv]} {
		set wb [expr {$wb + $bv * $flux}]
	    }
	}
	if {$idx_theta >= 0} {
	    set tv [string trim [lindex $row $idx_theta]]
	    if {[string is double -strict $tv]} {
		set wtheta [expr {$wtheta + $tv * $flux}]
	    }
	}
    }

    if {$total_flux <= 0} { set total_flux 1.0 }

    set new_x [expr {$wx / $total_flux}]
    set new_y [expr {$wy / $total_flux}]
    set new_a [expr {$wa / $total_flux}]
    set new_b [expr {$wb / $total_flux}]
    set new_theta [expr {$wtheta / $total_flux}]

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
    CatalogPanelMarkAll

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

    if {$catpanel(merge,active)} {
	CatalogPanelMergeCancel
    }
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

    # Re-mark if markers were present
    if {$current(frame) != {}} {
	catch {$current(frame) marker catalog sextract_all delete}
    }
    CatalogPanelMarkAll

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
