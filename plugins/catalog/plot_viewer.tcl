# catalog/plot_viewer.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/catalog/plugin.json.

proc CatalogPanelPlotDialog {} {
    global catpanel

    if {![info exists catpanel(alldata)] || $catpanel(alldata) eq {}} {
	set catpanel(status) "Extract sources first"
	return
    }

    # Get column names from header
    set lines [split $catpanel(alldata) "\n"]
    set header [lindex $lines 0]
    set colnames [split $header "\t"]

    set w .catplotdlg
    if {[winfo exists $w]} {
	raise $w
	return
    }
    toplevel $w
    wm title $w "Interactive Plot"
    wm geometry $w 340x260

    # Plot type
    ttk::label $w.lttype -text "Plot Type:"
    ttk::combobox $w.ptype -values {Scatter Histogram} \
	-state readonly -width 14
    $w.ptype set "Scatter"

    # X axis column
    ttk::label $w.ltx -text "X Column:"
    ttk::combobox $w.xcol -values $colnames -state readonly -width 18
    if {[llength $colnames] > 0} { $w.xcol set [lindex $colnames 0] }

    # Y axis column
    ttk::label $w.lty -text "Y Column:"
    ttk::combobox $w.ycol -values $colnames -state readonly -width 18
    if {[llength $colnames] > 1} { $w.ycol set [lindex $colnames 1] }

    # Log scale
    set catpanel(plot,logx) 0
    set catpanel(plot,logy) 0
    ttk::checkbutton $w.logx -text "Log X" -variable catpanel(plot,logx)
    ttk::checkbutton $w.logy -text "Log Y" -variable catpanel(plot,logy)

    # Histogram bins
    ttk::label $w.ltbins -text "Hist Bins:"
    ttk::spinbox $w.bins -from 10 -to 200 -increment 10 -width 6
    $w.bins set 50

    # Buttons
    ttk::frame $w.btns
    ttk::button $w.btns.plot -text "Plot" -command [list CatalogPanelPlotRun $w]
    ttk::button $w.btns.close -text "Close" -command [list destroy $w]
    pack $w.btns.plot $w.btns.close -side left -padx 5

    grid $w.lttype $w.ptype     -padx 5 -pady 3 -sticky w
    grid $w.ltx    $w.xcol      -padx 5 -pady 3 -sticky w
    grid $w.lty    $w.ycol      -padx 5 -pady 3 -sticky w
    grid $w.logx   $w.logy      -padx 5 -pady 3 -sticky w
    grid $w.ltbins $w.bins      -padx 5 -pady 3 -sticky w
    grid $w.btns   -            -padx 5 -pady 10
}

proc CatalogPanelPlotRun {dlg} {
    global catpanel

    set ptype [$dlg.ptype get]
    set xcol  [$dlg.xcol get]
    set ycol  [$dlg.ycol get]
    set logx  $catpanel(plot,logx)
    set logy  $catpanel(plot,logy)
    set nbins [$dlg.bins get]

    if {$xcol eq {}} {
	set catpanel(status) "Select X column"
	return
    }

    # Parse data
    set lines [split $catpanel(alldata) "\n"]
    set header [split [lindex $lines 0] "\t"]

    set xi -1
    set yi -1
    set ni -1
    for {set i 0} {$i < [llength $header]} {incr i} {
	set h [lindex $header $i]
	if {$h eq $xcol} { set xi $i }
	if {$h eq $ycol} { set yi $i }
	if {$h eq "NUMBER"} { set ni $i }
    }
    if {$xi < 0} {
	set catpanel(status) "Column $xcol not found"
	return
    }

    # Extract numeric data
    set xdata {}
    set ydata {}
    set nums {}
    for {set r 1} {$r < [llength $lines]} {incr r} {
	set row [lindex $lines $r]
	if {[string trim $row] eq {}} continue
	set fields [split $row "\t"]
	set xv [lindex $fields $xi]
	if {![string is double -strict $xv]} continue
	if {$logx && $xv > 0} { set xv [expr {log10($xv)}] }
	lappend xdata $xv
	if {$yi >= 0} {
	    set yv [lindex $fields $yi]
	    if {![string is double -strict $yv]} { set yv 0 }
	    if {$logy && $yv > 0} { set yv [expr {log10($yv)}] }
	    lappend ydata $yv
	}
	if {$ni >= 0} { lappend nums [lindex $fields $ni] }
    }

    if {[llength $xdata] == 0} {
	set catpanel(status) "No numeric data in $xcol"
	return
    }

    # Create plot window
    incr catpanel(plot,counter)
    set n $catpanel(plot,counter)
    set w .catplot_$n
    toplevel $w
    wm geometry $w 560x440

    if {$ptype eq "Histogram"} {
	CatalogPanelPlotHistogram $w $xdata $xcol $nbins $logx
    } else {
	if {$yi < 0 || $ycol eq {}} {
	    set catpanel(status) "Select Y column for scatter"
	    destroy $w
	    return
	}
	CatalogPanelPlotScatter $w $xdata $ydata $xcol $ycol $nums $logx $logy
    }
}

proc CatalogPanelPlotScatter {w xdata ydata xcol ycol nums logx logy} {
    global catpanel

    set xlabel $xcol
    set ylabel $ycol
    if {$logx} { set xlabel "log10($xcol)" }
    if {$logy} { set ylabel "log10($ycol)" }
    wm title $w "$ylabel vs $xlabel"

    # Create BLT vectors
    set xvec catplot_xv_$catpanel(plot,counter)
    set yvec catplot_yv_$catpanel(plot,counter)
    blt::vector create $xvec $yvec
    $xvec set $xdata
    $yvec set $ydata

    # Create graph
    blt::graph $w.g -width 520 -height 400 \
	-plotpadx {60 10} -plotpady {10 40} \
	-title "$ylabel vs $xlabel"
    pack $w.g -fill both -expand 1

    $w.g element create data -xdata $xvec -ydata $yvec \
	-symbol circle -pixels 3 -color blue -linewidth 0 \
	-label ""

    $w.g axis configure x -title $xlabel
    $w.g axis configure y -title $ylabel

    # Store mapping for click → source
    set catpanel(plot,$w,nums) $nums
    set catpanel(plot,$w,graph) $w.g

    # Bind click for source selection
    bind $w.g <ButtonPress-1> [list CatalogPanelPlotClick $w %x %y]

    # Cleanup on close
    bind $w <Destroy> [list CatalogPanelPlotCleanup $w $xvec $yvec]
}

proc CatalogPanelPlotHistogram {w xdata xcol nbins logx} {
    global catpanel

    set xlabel $xcol
    if {$logx} { set xlabel "log10($xcol)" }
    wm title $w "Histogram of $xlabel"

    # Compute bins
    set xmin [lindex $xdata 0]
    set xmax $xmin
    foreach v $xdata {
	if {$v < $xmin} { set xmin $v }
	if {$v > $xmax} { set xmax $v }
    }
    set range [expr {$xmax - $xmin}]
    if {$range <= 0} { set range 1.0 }
    set binw [expr {$range / double($nbins)}]

    # Count per bin
    set counts {}
    set centers {}
    for {set i 0} {$i < $nbins} {incr i} {
	lappend counts 0
	lappend centers [expr {$xmin + ($i + 0.5) * $binw}]
    }
    foreach v $xdata {
	set bi [expr {int(($v - $xmin) / $binw)}]
	if {$bi >= $nbins} { set bi [expr {$nbins - 1}] }
	if {$bi < 0} { set bi 0 }
	lset counts $bi [expr {[lindex $counts $bi] + 1}]
    }

    # BLT barchart
    set xvec catplot_hx_$catpanel(plot,counter)
    set yvec catplot_hy_$catpanel(plot,counter)
    blt::vector create $xvec $yvec
    $xvec set $centers
    $yvec set $counts

    blt::barchart $w.g -width 520 -height 400 \
	-plotpadx {60 10} -plotpady {10 40} \
	-title "Histogram of $xlabel" -barwidth $binw
    pack $w.g -fill both -expand 1

    $w.g element create hist -xdata $xvec -ydata $yvec \
	-foreground steelblue -borderwidth 1 -label ""

    $w.g axis configure x -title $xlabel
    $w.g axis configure y -title "Count"

    bind $w <Destroy> [list CatalogPanelPlotCleanup $w $xvec $yvec]
}

proc CatalogPanelPlotClick {w sx sy} {
    global catpanel

    if {![info exists catpanel(plot,$w,graph)]} return
    set g $catpanel(plot,$w,graph)

    if {[catch {$g element closest $sx $sy info -halo 10}]} return
    if {![info exists info(index)]} return

    set idx $info(index)
    if {[info exists catpanel(plot,$w,nums)]} {
	set num [lindex $catpanel(plot,$w,nums) $idx]
	if {$num ne {}} {
	    CatalogPanelGotoSource $num
	}
    }
}

proc CatalogPanelPlotCleanup {w xvec yvec} {
    global catpanel
    catch {blt::vector destroy $xvec}
    catch {blt::vector destroy $yvec}
    catch {unset catpanel(plot,$w,nums)}
    catch {unset catpanel(plot,$w,graph)}
}

proc CatalogPanelAnalysisViewer {} {
    global catpanel

    set script [CatalogPanelGetScript ds9_analysis_gui.py]
    if {![file exists $script]} {
	set catpanel(status) "Script not found: ds9_analysis_gui.py"
	return
    }

    set args [list [OGFPython] $script]

    # Pass FITS file if available
    set fn [CatalogPanelGetFITS]
    if {$fn ne {} && [file exists $fn]} {
	lappend args --fits $fn
    }

    # Pass catalog if available
    if {[info exists catpanel(alldata)] && $catpanel(alldata) ne {}} {
	set tmpcat [CatalogPanelSaveTempCatalog "analysis_viewer"]
	if {$tmpcat ne {}} {
	    lappend args --catalog $tmpcat
	}
    }

    # Launch non-blocking
    set catpanel(status) "Launching Analysis Viewer..."
    update idletasks
    exec {*}$args &
    set catpanel(status) "Analysis Viewer launched"
}

