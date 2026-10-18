# Guided recipes. Opening a tool is not proof that its scientific task succeeded.
namespace eval ::ogf::workflow {
    variable window {}
    variable catalogue {}
    variable choice {}
    variable recipe {}
    variable progress {}
    variable inputs {}
    variable outputs {}
    variable notes {}
    variable selected {}
    variable first_record 0
    variable busy 0
    variable message {}
    variable saved_records {}
    variable saved_parameters {}
    variable previous_parameters {}
    variable bindings {}
    variable refresh_pending 0
    variable restored_artifacts {}
}

proc ::ogf::workflow::script {} {
    return [file join [::ogf::step::plugin_dir workflow_designer] workflow_designer.py]
}

proc OGFWorkflowOpen {} {
    if {![info exists ::ogfui(root)] || ![winfo exists $::ogfui(root).info]} {
        ::ogf::status "Open the main analysis panel first."; return
    }
    set w $::ogfui(root).info.science_workflow
    set ::ogf::workflow::window $w
    if {[winfo exists $w]} {return}
    if {[catch {set raw [exec [OGFPython] [::ogf::workflow::script] catalog]} err]} {
        ::ogf::status "Workflow: $err"; return
    }
    set ::ogf::workflow::catalogue [dict get [::ogf::json::parse $raw] recipes]
    if {![llength [trace info variable ::ogfsess(steps)] ] || [lsearch -glob [trace info variable ::ogfsess(steps)] *::ogf::workflow::records_changed*] < 0} {
        trace add variable ::ogfsess(steps) write ::ogf::workflow::records_changed
    }
    ttk::frame $w
    set siblings [pack slaves $::ogfui(root).info]
    if {[llength $siblings]} {
        pack $w -side top -fill x -before [lindex $siblings 0]
    } else {
        pack $w -side top -fill x
    }
    ttk::frame $w.top
    ttk::label $w.top.label -text "분석 목적"
    set titles {}
    foreach r $::ogf::workflow::catalogue {lappend titles [dict get $r title]}
    ttk::combobox $w.top.pick -state readonly -width 25 -values $titles -textvariable ::ogf::workflow::choice
    pack $w.top.label $w.top.pick -side left -padx 5 -pady 5
    pack $w.top -fill x
    ttk::frame $w.controls
    ttk::button $w.controls.resume -text "저장 세션 열기" -command ::ogf::workflow::resume
    ttk::button $w.controls.save -text "현재 단계까지 저장" -command ::ogf::workflow::save
    pack $w.controls.resume $w.controls.save -side left -padx 5 -pady 3
    pack $w.controls -fill x
    bind $w.top.pick <<ComboboxSelected>> ::ogf::workflow::choose
    ttk::label $w.summary -wraplength 540 -justify left
    pack $w.summary -fill x -padx 12 -pady 5
    ttk::panedwindow $w.body -orient horizontal
    ttk::frame $w.body.left
    canvas $w.body.left.flow -width 235 -height 240 -background #f4f6fa -yscrollcommand [list $w.body.left.scroll set]
    ttk::scrollbar $w.body.left.scroll -orient vertical -command [list $w.body.left.flow yview]
    pack $w.body.left.scroll -side right -fill y
    pack $w.body.left.flow -fill both -expand 1
    ttk::frame $w.body.right
    text $w.body.right.detail -width 32 -height 6 -wrap word -state disabled
    ttk::label $w.body.right.nl -text "단계 메모 / 검토 근거"
    text $w.body.right.note -height 2 -width 32 -wrap word
    ttk::button $w.body.right.open -text "선택 단계 도구 열기" -command ::ogf::workflow::launch
    ttk::button $w.body.right.review -text "검토 완료 표시 → 다음" -command ::ogf::workflow::review
    ttk::button $w.body.right.attach -text "산출물 파일 추가…" -command ::ogf::workflow::attach
    ttk::button $w.body.right.settings -text "선택 단계 입력·매개변수" -command ::ogf::workflow::settings
    foreach child {detail nl note settings open review attach} {pack $w.body.right.$child -fill x -padx 8 -pady 5}
    $w.body add $w.body.left -weight 1
    $w.body add $w.body.right -weight 1
    pack $w.body -fill both -expand 1 -padx 10
    ttk::label $w.status -textvariable ::ogf::workflow::message -wraplength 540
    pack $w.status -fill x -padx 12 -pady 8
    if {$::ogf::workflow::recipe eq {}} {
        $w.top.pick current 0
        ::ogf::workflow::choose 0
    } else {
        ::ogf::workflow::draw
    }
}

proc ::ogf::workflow::close {} {
    variable busy
    if {$busy} {set ::ogf::workflow::message "저장 작업이 끝난 뒤 창을 닫을 수 있습니다."; return}
    ::ogf::workflow::keep_note
    destroy $::ogf::workflow::window
}

proc ::ogf::workflow::keep_note {} {
    variable selected; variable notes
    if {$selected ne {} && [winfo exists $::ogf::workflow::window.body.right.note]} {
        dict set notes $selected [string trim [$::ogf::workflow::window.body.right.note get 1.0 end]]
    }
}

proc ::ogf::workflow::choose {{ask 1}} {
    variable catalogue; variable choice; variable recipe; variable progress
    variable inputs; variable outputs; variable notes; variable selected; variable first_record
    variable saved_records; variable saved_parameters; variable busy
    if {$busy} {if {$recipe ne {}} {set choice [dict get $recipe title]}; return}
    if {$ask && $recipe ne {}} {
        set answer [tk_messageBox -parent $::ogf::workflow::window -type yesno -message "다른 연구를 새로 시작할까요? 현재 진행은 먼저 저장해야 다시 열 수 있습니다."]
        if {$answer ne "yes"} {set choice [dict get $recipe title]; return}
    }
    foreach r $catalogue {
        if {[dict get $r title] eq $choice} {set recipe $r; break}
    }
    set progress {}; set inputs {}; set outputs {}; set notes {}; set saved_records {}; set saved_parameters {}
    set ::ogf::workflow::bindings {}
    set ::ogf::workflow::previous_parameters {}
    set ::ogf::workflow::restored_artifacts {}
    set selected load
    set first_record [llength [::ogf::session::steps]]
    set ::ogf::workflow::message "연결된 도구를 열고 결과를 검토하세요. 검토 표시는 과학 검증 통과를 뜻하지 않습니다."
    draw
}

proc ::ogf::workflow::select {sid} {
    variable selected
    keep_note
    set selected $sid
    draw
}

proc ::ogf::workflow::step {} {
    variable recipe; variable selected
    foreach s [dict get $recipe steps] {if {[dict get $s id] eq $selected} {return $s}}
    error "Unknown workflow step"
}

proc ::ogf::workflow::draw {} {
    variable recipe; variable selected; variable progress; variable notes; variable choice
    if {![winfo exists $::ogf::workflow::window]} return
    set choice [dict get $recipe title]
    $::ogf::workflow::window.summary configure -text "입력: [dict get $recipe inputs]\n최종 산출물: [dict get $recipe outputs]"
    set c $::ogf::workflow::window.body.left.flow
    $c delete all
    set y 16
    foreach s [dict get $recipe steps] {
        set sid [dict get $s id]
        set st [::ogf::json::get $progress $sid pending]
        set fill [expr {$sid eq $selected ? "#dceaff" : ($st eq "reviewed" ? "#e1f2e7" : "white")}]
        set tag "step_$sid"
        $c create rectangle 10 $y 235 [expr {$y+48}] -fill $fill -outline #aeb9c9 -tags $tag
        set label "[dict get $s title]  ·  $st"
        if {[dict get $s optional]} {append label " (조건부)"}
        $c create text 20 [expr {$y+24}] -anchor w -width 205 -text $label -tags $tag
        $c bind $tag <Button-1> [list ::ogf::workflow::select $sid]
        incr y 48
        if {$sid ne "results"} {$c create line 120 $y 120 [expr {$y+18}] -arrow last -fill #75859b}
        incr y 24
    }
    $c configure -scrollregion [list 0 0 245 $y]
    set s [step]
    set a [dict get $s action]
    set availability [expr {[llength $a] ? "기존 도구/자료 입력 연결" : "검토 또는 추가 분석이 필요한 단계"}]
    set d $::ogf::workflow::window.body.right.detail
    $d configure -state normal
    $d delete 1.0 end
    $d insert end "[dict get $s title]\n\n[dict get $s detail]\n\n$availability\n\n현재 입력 [llength $::ogf::workflow::inputs]개 / 명시적으로 연결한 산출물 [llength $::ogf::workflow::outputs]개"
    if {[llength $a] == 2 && [info exists ::ogf::step::last([join $a .])]} {
        $d insert end "\n\n실행 보고서:\n[string range $::ogf::step::last([join $a .]) 0 4000]"
    }
    $d configure -state disabled
    $::ogf::workflow::window.body.right.note delete 1.0 end
    $::ogf::workflow::window.body.right.note insert end [::ogf::json::get $notes $selected]
    $::ogf::workflow::window.body.right.open configure -state [expr {[llength $a] ? "normal" : "disabled"}]
}

proc ::ogf::workflow::input_artifact {node artifact} {
    variable recipe; variable saved_records; variable first_record; variable restored_artifacts; variable bindings
    set action {}
    foreach item [dict get $recipe steps] {if {[dict get $item id] eq $node} {set action [dict get $item action]}}
    if {[llength $action] != 2} {error "선행 입력 단계가 없습니다."}
    set records {}
    foreach raw $saved_records {lappend records [::ogf::json::parse $raw]}
    foreach rec [lrange [::ogf::session::steps] $first_record end] {lappend records $rec}
    foreach rec [lreverse $records] {
        if {[dict get $rec step] ne [join $action .]} continue
        if {[dict exists $bindings $node] && [dict get $rec seq] ni [dict get $bindings $node]} continue
        # Use the latest matching execution; never fall back to an older success.
        if {[dict get $rec failed] || [dict get $rec undone] || [dict get $rec ms] eq {}} {error "선행 적합이 완료되지 않았습니다."}
        foreach source [dict get $rec outputs] {
            if {[file tail $source] ne $artifact} continue
            set path $source
            if {[dict exists $restored_artifacts $source]} {set path [dict get $restored_artifacts $source]}
            if {[file isfile $path]} {return $path}
        }
        error "선행 적합의 진단 입력 파일이 없습니다."
    }
    error "선행 적합을 먼저 실행하세요."
}

proc ::ogf::workflow::launch {} {
    variable progress; variable selected; variable inputs; variable recipe; variable busy
    if {$busy} return
    keep_note
    set a [dict get [step] action]
    if {[lindex $a 0] eq "checkpoint"} {save; return}
    set before [llength [::ogf::session::steps]]
    if {$a eq {}} return
    if {[catch {
        if {[lindex $a 0] eq "load"} {
            set files [tk_getOpenFile -parent $::ogf::workflow::window -multiple 1 -title "분석 자료 선택 (영상 / 스펙트럼 / 표)"]
            if {![llength $files]} return
            set inputs [lsort -unique $files]
            dict for {sid value} $progress {if {$value in {reviewed succeeded}} {dict set progress $sid needs_review}}
            if {[string match *_spectrum [dict get $recipe id]]} {
                if {[string equal -nocase [file extension [lindex $files 0]] .tsv]} {
                    if {[llength $files] != 1} {error "1D TSV 경로는 스펙트럼 하나씩 불러옵니다."}
                    if {[file size [lindex $files 0]] > 33554432} {error "1D TSV는 32 MiB 이하이어야 합니다."}
                    set fd [open [lindex $files 0] r]; fconfigure $fd -encoding utf-8
                    try {set data [read $fd]} finally {::close $fd}
                    ::ogf::cat::load_tsv $data [file tail [lindex $files 0]]
                } else {
                    ::ogf::params::put spectra spec-dir [file dirname [lindex $files 0]]
                    OGFParamDialog spectra
                }
            } else {
                foreach f $files {CreateFrame; LoadFitsFile $f {} {}}
            }
        } elseif {[lindex $a 0] eq "settings"} {
            OGFParamDialog [lindex $a 1]
        } else {
            set spec [step]
            if {[dict exists $spec personal_digest]} {
                set pm [::ogf::reg::get [lindex $a 0]]
                if {[dict get $pm personal_digest] ne [dict get $spec personal_digest]} {error "등록 버전이 변경되었습니다. 저장한 세션을 복원하거나 새 버전을 다시 편입하세요."}
            }
            if {[dict exists $spec desktop_input_param]} {
                set source [input_artifact [dict get $spec input_of] [dict get $spec input_artifact]]
                ::ogf::params::put [lindex $a 0] [dict get $spec desktop_input_param] $source
            }
            set registered [::ogf::reg::step {*}$a]
            if {[dict exists $registered variants]} {
                set menu $::ogf::workflow::window.variants
                catch {destroy $menu}
                menu $menu -tearoff 0
                foreach v [dict get $registered variants] {
                    $menu add command -label [dict get $v label] -command [list ::ogf::workflow::variant [lindex $a 0] [lindex $a 1] [dict get $v id]]
                }
                tk_popup $menu [winfo pointerx $::ogf::workflow::window] [winfo pointery $::ogf::workflow::window]
                return
            }
            if {![::ogf::step::run {*}$a]} {error "도구가 시작되지 않았습니다. 입력·설정과 상태 메시지를 확인하세요."}
        }
        dict set progress $selected opened
    } err]} {bind_new_records $selected $before; reconcile; set ::ogf::workflow::message $err; draw; return}
    set ::ogf::workflow::message "단계를 열었습니다. 완료된 결과·검증을 확인한 뒤 검토 완료를 표시하세요."
    bind_new_records $selected $before
    reconcile
    draw
}

proc ::ogf::workflow::variant {pid sid vid} {
    if {$::ogf::workflow::busy} return
    foreach v [dict get [::ogf::reg::step $pid $sid] variants] {
        if {[dict get $v id] eq $vid} {
            set before [llength [::ogf::session::steps]]
            if {[catch {uplevel #0 [dict get $v proc]} err]} {set ::ogf::workflow::message $err; return}
            dict set ::ogf::workflow::progress $::ogf::workflow::selected opened
            bind_new_records $::ogf::workflow::selected $before
            reconcile
            draw
            return
        }
    }
}

proc ::ogf::workflow::review {} {
    variable progress; variable selected; variable recipe; variable busy
    if {$busy || [::ogf::job::busy]} {set ::ogf::workflow::message "작업이 진행 중입니다. 결과가 완료된 뒤 검토하세요."; return}
    keep_note
    reconcile
    if {[::ogf::json::get $progress $selected] in {failed cancelled}} {
        set ::ogf::workflow::message "실패하거나 취소된 계산은 재실행한 뒤 검토하세요."; return
    }
    dict set progress $selected reviewed
    set steps [dict get $recipe steps]
    for {set i 0} {$i < [llength $steps]-1} {incr i} {
        if {[dict get [lindex $steps $i] id] eq $selected} {
            set selected [dict get [lindex $steps [expr {$i+1}]] id]; break
        }
    }
    draw
}

proc ::ogf::workflow::attach {} {
    variable outputs; variable busy
    if {$busy} return
    set files [tk_getOpenFile -parent $::ogf::workflow::window -multiple 1 -title "보존할 완료 산출물 파일 선택"]
    set outputs [lsort -unique [concat $outputs $files]]
    keep_note
    draw
}

proc ::ogf::workflow::bind_new_records {node before} {
    variable bindings
    set ids {}
    foreach rec [lrange [::ogf::session::steps] $before end] {
        if {[llength [dict get $rec argv]]} {lappend ids [dict get $rec seq]}
    }
    if {[llength $ids]} {dict set bindings $node $ids}
    if {[llength $ids]} {
        set downstream 0
        foreach s [dict get $::ogf::workflow::recipe steps] {
            set sid [dict get $s id]
            if {$downstream && [::ogf::json::get $::ogf::workflow::progress $sid] in {reviewed succeeded}} {
                dict set ::ogf::workflow::progress $sid needs_review
            }
            if {$sid eq $node} {set downstream 1}
        }
    }
}

proc ::ogf::workflow::execution_status {records active_seq} {
    if {![llength $records]} {return opened}
    set status succeeded
    foreach rec $records {
        if {[dict get $rec seq] eq $active_seq} {return running}
        if {[dict get $rec undone] || [string match -nocase *cancelled* [dict get $rec note]]} {return cancelled}
        if {[dict get $rec failed]} {return failed}
        if {[dict get $rec ms] eq {}} {set status opened}
    }
    return $status
}

proc ::ogf::workflow::reconcile {} {
    variable bindings; variable progress
    set active {}
    if {[::ogf::job::busy]} {set active $::ogf::job::cur(seq)}
    set records [::ogf::session::steps]
    dict for {node ids} $bindings {
        set matched {}
        foreach rec $records {if {[dict get $rec seq] in $ids} {lappend matched $rec}}
        if {![llength $matched]} continue
        set status [execution_status $matched $active]
        if {$status eq "succeeded" && [::ogf::json::get $progress $node] in {reviewed needs_review}} continue
        dict set progress $node $status
    }
}

proc ::ogf::workflow::records_changed {args} {
    variable refresh_pending
    if {$refresh_pending} return
    set refresh_pending 1
    after idle ::ogf::workflow::refresh_records
}

proc ::ogf::workflow::refresh_records {} {
    variable refresh_pending
    set refresh_pending 0
    if {$::ogf::workflow::recipe eq {}} return
    keep_note
    reconcile
    draw
}

proc ::ogf::workflow::save {} {
    variable busy; variable selected; variable recipe; variable progress; variable inputs; variable outputs
    variable notes; variable first_record; variable saved_records; variable saved_parameters
    if {$busy} return
    keep_note
    reconcile
    set dest [tk_chooseDirectory -parent $::ogf::workflow::window -title "연구 세션 저장 폴더"]
    if {$dest eq {}} return
    # The state/request contains only explicit domain parameters and recorded work.
    set draft [file join $dest .request-[pid]-[clock clicks]]
    file mkdir $draft
    set catalog [file join $draft current_catalog.tsv]
    set fd [open $catalog w]; fconfigure $fd -encoding utf-8
    puts -nonewline $fd [::ogf::cat::tsv]; ::close $fd
    set out [concat $outputs [list $catalog]]
    set records $saved_records
    set active [::ogf::job::busy]
    if {!$active} {OGFSessFinalize}
    foreach rec [lrange [::ogf::session::steps] $first_record end] {
        lappend records [OGFSessRecordJSON $rec]
        if {[dict get $rec ms] ne {} && ![dict get $rec failed] && ![dict get $rec undone]} {
            foreach f [dict get $rec outputs] {lappend out $f}
        }
    }
    set params $saved_parameters
    foreach s [dict get $recipe steps] {
        set a [dict get $s action]
        if {[llength $a] != 2} continue
        set pid [expr {[lindex $a 0] eq "settings" ? [lindex $a 1] : [lindex $a 0]}]
        set values {}
        foreach p [::ogf::json::get [::ogf::reg::get $pid] params] {
            set n [dict get $p name]
            dict set values $n [::ogf::params::get $pid $n]
        }
        dict set params ${pid}_json [OGFJDict $values]
    }
    if {$::ogf::workflow::previous_parameters ne {} && $params ne $::ogf::workflow::previous_parameters} {
        dict for {sid value} $progress {if {$value in {reviewed succeeded}} {dict set progress $sid needs_review}}
    }
    set ::ogf::workflow::previous_parameters $params
    set state [dict create personal_steps_json [personal_json] recipe [dict get $recipe id] selected $selected progress_json [OGFJDict $progress] \
        notes_json [OGFJDict $notes] parameters_json [OGFJDict $params] records_json "\[[join $records ,]\]" \
        catalog_source $catalog running_at_save $active input_policy external_references bindings_json [OGFJDict $::ogf::workflow::bindings]]
    set request [OGFJDict [dict create state_json [OGFJDict $state] outputs_list [lsort -unique $out] inputs_list $inputs]]
    set path [file join $draft request.json]
    set fd [open $path w]; fconfigure $fd -encoding utf-8; puts $fd $request; ::close $fd
    set busy 1
    set ::ogf::workflow::message "산출물 복사·체크섬 검증 중… 저장 완료 전까지 이 창을 유지합니다."
    set cmd [list [OGFPython] [script] save --request $path --destination $dest]
    if {[catch {set channel [open [concat | $cmd [list 2>@1]] r]} err]} {
        set busy 0; set ::ogf::workflow::message $err; return
    }
    fconfigure $channel -blocking 0 -encoding utf-8
    set ::ogf::workflow::buffer {}
    fileevent $channel readable [list ::ogf::workflow::saved $channel $draft]
}

proc ::ogf::workflow::saved {channel draft} {
    append ::ogf::workflow::buffer [read $channel]
    if {![eof $channel]} return
    fileevent $channel readable {}
    fconfigure $channel -blocking 1
    set rc [catch {::close $channel} err]
    set ::ogf::workflow::busy 0
    if {$rc || [catch {set obj [::ogf::json::parse $::ogf::workflow::buffer]; set manifest [dict get $obj manifest]}]} {
        set ::ogf::workflow::message "저장 실패: $::ogf::workflow::buffer (이전 스냅샷은 유지됩니다.)"
    } else {
        set ::ogf::workflow::message "저장 완료: $manifest"
        file delete -force $draft
    }
}

proc ::ogf::workflow::resume {} {
    variable busy
    if {$busy || [::ogf::job::busy]} return
    set manifest [tk_getOpenFile -parent $::ogf::workflow::window -title "저장한 session.json 선택" -filetypes {{{Research session} {.json}}}]
    if {$manifest eq {}} return
    set work [file join [OGFSessWorkDir] workflow-resume-[pid]-[clock clicks]]
    set cmd [list [OGFPython] [script] restore --manifest $manifest --work $work]
    if {[catch {set channel [open [concat | $cmd [list 2>@1]] r]} err]} {
        set ::ogf::workflow::message $err; return
    }
    set busy 1
    set ::ogf::workflow::buffer {}
    set ::ogf::workflow::message "스냅샷을 검증하고 재개용 작업 사본을 만드는 중…"
    fconfigure $channel -blocking 0 -encoding utf-8
    fileevent $channel readable [list ::ogf::workflow::resumed $channel]
}

proc ::ogf::workflow::resumed {channel} {
    append ::ogf::workflow::buffer [read $channel]
    if {![eof $channel]} return
    fileevent $channel readable {}
    fconfigure $channel -blocking 1
    set rc [catch {::close $channel} err]
    set ::ogf::workflow::busy 0
    if {$rc || [catch {set data [::ogf::json::parse $::ogf::workflow::buffer]}]} {
        set ::ogf::workflow::message "복원 실패: $::ogf::workflow::buffer"; return
    }
    if {[catch {::ogf::workflow::apply_restore $data} err]} {
        set ::ogf::workflow::message "복원 설정 적용 실패: $err"
    }
}

proc ::ogf::workflow::apply_restore {data} {
    variable busy; variable recipe; variable progress; variable selected; variable inputs; variable outputs; variable notes
    variable first_record; variable saved_records; variable saved_parameters
    set st [dict get $data state]
    if {[dict exists $st personal_steps] && [llength [dict get $st personal_steps]]} {OGFUIUserToolReload}
    set recipe [dict get $st recipe_definition]
    set progress [dict get $st progress]; set selected [dict get $st selected]; set notes [dict get $st notes]
    set ::ogf::workflow::bindings {}
    dict for {sid value} $progress {
        if {$value in {running queued}} {dict set progress $sid needs_review}
    }
    set mapping {}; set outputs {}
    foreach a [dict get $data artifacts] {
        dict set mapping [dict get $a source] [dict get $a resolved_path]
        lappend outputs [dict get $a resolved_path]
    }
    set ::ogf::workflow::restored_artifacts $mapping
    set inputs {}
    foreach r [dict get $data inputs] {lappend inputs [dict get $r path]}
    set saved_records {}
    # Preserve the original record dictionaries as JSON without executing their argv.
    foreach rec [dict get $st records] {lappend saved_records [OGFSessRecordJSON $rec]}
    set saved_parameters {}
    dict for {pid values} [dict get $st parameters] {
        dict for {k v} $values {
            if {[dict exists $mapping $v]} {set v [dict get $mapping $v]}
            catch {::ogf::params::put $pid $k $v}
            dict set values $k $v
        }
        dict set saved_parameters ${pid}_json [OGFJDict $values]
    }
    set ::ogf::workflow::previous_parameters $saved_parameters
    set first_record [llength [::ogf::session::steps]]
    if {[dict exists $mapping [dict get $st catalog_source]]} {
        set f [open [dict get $mapping [dict get $st catalog_source]] r]; fconfigure $f -encoding utf-8
        set tsv [read $f]; ::close $f
        if {[string trim $tsv] ne {}} {::ogf::cat::load_tsv $tsv "Restored science workflow"}
    }
    set changes [dict get $data input_changes]
    if {[llength $changes]} {
        dict for {sid value} $progress {if {$value in {reviewed succeeded}} {dict set progress $sid needs_review}}
        set ::ogf::workflow::message "산출물 복원 완료. 원본 입력 [llength $changes]개가 변경/누락되어 입력 연결과 후속 결과 재검토가 필요합니다."
    } else {
        set ::ogf::workflow::message "진행·설정·표·산출물을 복원했습니다. 원본 영상은 자료 불러오기에서 연결하세요. 실행 중이던 계산은 재시작해야 합니다."
    }
    draw
}

proc ::ogf::workflow::attach_personal {pid} {
    variable recipe; variable selected; variable progress; variable busy
    if {$busy || [::ogf::job::busy]} return
    if {$recipe eq {}} {set ::ogf::workflow::message "먼저 분석 목적을 선택하세요."; return}
    set m [::ogf::reg::get $pid]
    set sid personal_[clock clicks]
    set node [dict create id $sid title "$pid · 사용자 단계" detail [dict get $m description] action [list $pid run] optional 0 personal_digest [dict get $m personal_digest]]
    set steps [dict get $recipe steps]; set index 0
    foreach s $steps {incr index; if {[dict get $s id] eq $selected} break}
    dict set recipe steps [linsert $steps $index $node]
    set selected $sid;dict set progress $sid pending
    draw
}
proc ::ogf::workflow::personal_json {} {
    variable recipe
    set rows {};set previous {}
    foreach s [dict get $recipe steps] {
        if {[dict exists $s personal_digest]} {
            lappend rows [OGFJDict [dict create id [dict get $s id] after $previous tool [lindex [dict get $s action] 0] digest [dict get $s personal_digest]]]
        }
        set previous [dict get $s id]
    }
    return "\[[join $rows ,]\]"
}

proc ::ogf::workflow::settings {} {
    set a [dict get [step] action]
    if {[llength $a] != 2} return
    set pid [expr {[lindex $a 0] eq "settings" ? [lindex $a 1] : [lindex $a 0]}]
    OGFParamDialog $pid
}
