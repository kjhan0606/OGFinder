# CCD reduction dialog: bias, dark, flat and snapshot.
# The AI button classifies by header type, file name, file size and exposure, then flat-fields.
# Manual assigns each type and runs one step (mean, median, sum, subtract, divide, normalize) or the recipe.

proc OGFCcdWork {} {
    return [file join [OGFSessWorkDir] ccdreduce]
}

proc OGFCcdOpen {} {
    OGFCcdDialog
}

proc OGFCcdEnsureTrace {} {
    if {[info exists ::ogfccd_traced]} return
    set ::ogfccd_traced 1
    trace add variable ::ogfccd_class write OGFCcdClassTrace
}

proc OGFCcdClassTrace {name idx op} {
    if {[info exists ::ogfccd_filling] || $idx eq {}} return
    set ::ogfccd_user($idx) 1
}

proc OGFCcdDialog {} {
    OGFCcdEnsureTrace
    set w .ogfccd
    if {[winfo exists $w]} {raise $w; focus $w; return}
    if {![info exists ::ogfccd_files]} {set ::ogfccd_files {}}
    if {[llength $::ogfccd_files] == 0} {
        foreach f [split [::ogf::params::get ccdreduce frames] |] {
            set f [string trim $f]
            if {$f ne {}} {lappend ::ogfccd_files $f}
        }
    }
    toplevel $w
    wm title $w "CCD reduction"
    wm resizable $w 1 1
    ttk::frame $w.top
    pack $w.top -fill x -padx 8 -pady 6
    ttk::button $w.top.add -text "Add FITS…" -command OGFCcdAdd
    ttk::button $w.top.ai -text "✦ AI reduce" -command OGFCcdAi
    ttk::button $w.top.man -text "Manual" -command OGFCcdManual
    pack $w.top.add $w.top.ai $w.top.man -side left -padx 3
    listbox $w.list -height 8 -width 78 -exportselection 0
    pack $w.list -fill both -expand 1 -padx 8
    ttk::label $w.recipe -wraplength 560 -justify left -text "Master bias (median or mean), dark scaled by exposure to a rate, flat bias/dark-subtracted and divided by its median, snapshot = (raw − bias − dark_rate×t) / flat."
    pack $w.recipe -fill x -padx 8 -pady 4
    ttk::frame $w.ed
    ttk::label $w.ed.h -wraplength 560 -justify left -text "Assign each frame, then apply the recipe or one step. Subtract and divide use the first two frames in this list: the first minus or over the second."
    pack $w.ed.h -anchor w -pady 2
    ttk::frame $w.ed.rows
    pack $w.ed.rows -fill x
    ttk::frame $w.ed.ops
    pack $w.ed.ops -fill x -pady 4
    ttk::button $w.ed.ops.sug -text "Suggest types" -command OGFCcdSuggest
    ttk::button $w.ed.ops.app -text "Apply recipe" -command OGFCcdApply
    ttk::label $w.ed.ops.cl -text "Combine"
    ttk::combobox $w.ed.ops.cb -values {median mean} -state readonly -width 8
    set comb [::ogf::params::get ccdreduce combine]
    if {$comb ni {median mean}} {set comb median}
    $w.ed.ops.cb set $comb
    pack $w.ed.ops.sug $w.ed.ops.app $w.ed.ops.cl $w.ed.ops.cb -side left -padx 3
    foreach op {mean median sum subtract divide normalize} {
        ttk::button $w.ed.ops.$op -text [string totitle $op] -command [list OGFCcdOp $op]
        pack $w.ed.ops.$op -side left -padx 2
    }
    OGFCcdPaint
}

proc OGFCcdPaint {} {
    set w .ogfccd
    if {![winfo exists $w]} return
    $w.list delete 0 end
    foreach f $::ogfccd_files {$w.list insert end $f}
    OGFCcdRows
}

proc OGFCcdAdd {} {
    set picked [tk_getOpenFile -title "Bias, dark, flat and snapshot FITS" -multiple 1 -filetypes {{FITS {.fits .fit .fts .gz}} {All *}}]
    if {$picked eq {}} return
    foreach f $picked {
        if {$f ni $::ogfccd_files} {lappend ::ogfccd_files $f}
    }
    ::ogf::params::put ccdreduce frames [join $::ogfccd_files |]
    OGFCcdPaint
}

proc OGFCcdManual {} {
    set w .ogfccd
    if {![winfo exists $w]} {OGFCcdDialog}
    if {![winfo ismapped $w.ed]} {pack $w.ed -fill x -padx 8 -pady 2}
    OGFCcdRows
}

proc OGFCcdRows {} {
    set w .ogfccd
    if {![winfo exists $w.ed.rows]} return
    foreach c [winfo children $w.ed.rows] {destroy $c}
    set ::ogfccd_filling 1
    set n 0
    foreach f $::ogfccd_files {
        if {![info exists ::ogfccd_class($f)]} {set ::ogfccd_class($f) snapshot}
        set row $w.ed.rows.r$n
        ttk::frame $row
        pack $row -fill x -pady 1
        ttk::label $row.n -text [file tail $f] -width 28 -anchor w
        ttk::combobox $row.c -textvariable ::ogfccd_class($f) -values {bias dark flat snapshot ignore} -state readonly -width 12
        set why ""
        if {[info exists ::ogfccd_why($f)]} {set why $::ogfccd_why($f)}
        ttk::label $row.w -text $why -anchor w
        pack $row.n $row.c -side left -padx 2
        pack $row.w -side left -fill x -expand 1
        incr n
    }
    unset -nocomplain ::ogfccd_filling
}

proc OGFCcdCombine {} {
    set c median
    if {[winfo exists .ogfccd.ed.ops.cb]} {set c [.ogfccd.ed.ops.cb get]}
    if {$c ni {median mean}} {set c median}
    ::ogf::params::put ccdreduce combine $c
    return $c
}

proc OGFCcdAssignment {} {
    set parts {}
    foreach f $::ogfccd_files {
        set c snapshot
        if {[info exists ::ogfccd_class($f)] && $::ogfccd_class($f) ne {}} {set c $::ogfccd_class($f)}
        lappend parts $f:$c
    }
    return [join $parts |]
}

proc OGFCcdLaunch {title kind tail} {
    if {![info exists ::ogfccd_files] || [llength $::ogfccd_files] < 1} {
        ::ogf::status "CCD reduction: add the bias, dark, flat and snapshot FITS first"
        return
    }
    set py [OGFPython]
    set drv [file join [::ogf::step::plugin_dir ccdreduce] run_ccd.py]
    set core [::ogf::params::get ccdreduce core-dir]
    set frames [join $::ogfccd_files |]
    ::ogf::params::put ccdreduce frames $frames
    file mkdir [OGFCcdWork]
    set argv [list $py $drv --core-dir $core {*}$tail --frames $frames]
    ::ogf::job::run $argv -step ccdreduce.dialog -class manual -title $title -plugin ccdreduce -done [list OGFCcdDone $kind]
}

proc OGFCcdAi {} {
    if {![winfo exists .ogfccd]} {OGFCcdDialog}
    set dir [OGFCcdWork]
    set comb [OGFCcdCombine]
    OGFCcdLaunch "AI CCD reduction" snapshot_flat.fits [list reduce --assignment auto --combine $comb --out [file join $dir snapshot_flat.fits] --summary [file join $dir summary.json]]
}

proc OGFCcdSuggest {} {
    set dir [OGFCcdWork]
    OGFCcdLaunch "Classify frames" classify [list inspect --summary [file join $dir summary.json]]
}

proc OGFCcdApply {} {
    set dir [OGFCcdWork]
    set comb [OGFCcdCombine]
    set asg [OGFCcdAssignment]
    ::ogf::params::put ccdreduce assignment $asg
    OGFCcdLaunch "CCD recipe" snapshot_flat.fits [list reduce --assignment $asg --combine $comb --out [file join $dir snapshot_flat.fits] --summary [file join $dir summary.json]]
}

proc OGFCcdOp {op} {
    if {![info exists ::ogfccd_files]} {set ::ogfccd_files {}}
    set all $::ogfccd_files
    if {$op in {subtract divide}} {
        if {[llength $all] < 2} {
            ::ogf::status "Subtract and divide use the first two frames, in list order: the first minus or over the second."
            return
        }
        set ::ogfccd_files [lrange $all 0 1]
    } elseif {$op eq "normalize"} {
        if {[llength $all] < 1} {::ogf::status "Normalize uses the first frame in the list."; return}
        set ::ogfccd_files [lrange $all 0 0]
    }
    set dir [OGFCcdWork]
    OGFCcdLaunch "CCD $op" manual.fits [list op --op $op --out [file join $dir manual.fits]]
    set ::ogfccd_files $all
    ::ogf::params::put ccdreduce frames [join $all |]
}

proc OGFCcdDone {kind ok output ms} {
    if {!$ok} return
    set dir [OGFCcdWork]
    set st [file join $dir status.txt]
    if {[file exists $st]} {
        set fh [open $st r]
        set txt [string trim [read $fh]]
        close $fh
        if {$txt ne {}} {::ogf::status $txt}
    }
    if {$kind eq "classify"} {
        OGFCcdApplyTsv [file join $dir classification.tsv]
        return
    }
    OGFCcdLoad [file join $dir $kind]
}

proc OGFCcdApplyTsv {path} {
    if {![file exists $path]} return
    set fh [open $path r]
    set body [read $fh]
    close $fh
    set ::ogfccd_filling 1
    foreach ln [lrange [split $body \n] 1 end] {
        if {$ln eq {}} continue
        set c [split $ln \t]
        set name [lindex $c 0]
        set label [lindex $c 2]
        if {$label ni {bias dark flat snapshot ignore}} {set label ignore}
        set why [lindex $c 7]
        foreach f $::ogfccd_files {
            if {[file tail $f] ne $name} continue
            if {[info exists ::ogfccd_user($f)]} continue
            set ::ogfccd_class($f) $label
            set ::ogfccd_why($f) $why
        }
    }
    unset -nocomplain ::ogfccd_filling
    OGFCcdRows
}

proc OGFCcdLoad {fits} {
    if {![file exists $fits]} return
    CreateFrame
    if {[catch {LoadFitsFile $fits {} {}} err]} {::ogf::log ERROR "ccdreduce: cannot load $fits: $err"}
}

proc OGFCcdAfter {} {
    set dir [file join [OGFSessWorkDir] ccdreduce]
    set st [file join $dir status.txt]
    if {[file exists $st]} {
        set fh [open $st r]
        set txt [string trim [read $fh]]
        close $fh
        if {$txt ne {}} {::ogf::status $txt}
    }
    OGFCcdLoad [file join $dir snapshot_flat.fits]
}
