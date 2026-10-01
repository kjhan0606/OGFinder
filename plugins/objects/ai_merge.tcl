# objects/ai_merge.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/objects/plugin.json.

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
    if {[catch {set data [exec [OGFPython] $script $fn {*}$paramargs 2>$errfile]} err]} {
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

