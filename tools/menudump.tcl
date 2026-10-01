# dump every command reachable from the new UI: Workflow/Tools menus and every plugin chip menu
proc dumpmenu {m path out} {
  set n [$m index end]
  if {$n eq "none"} return
  for {set i 0} {$i <= $n} {incr i} {
    set t [$m type $i]
    if {$t eq "separator" || $t eq "tearoff"} continue
    set lab [$m entrycget $i -label]
    set cmd {}; set var {}
    catch {set cmd [$m entrycget $i -command]}
    catch {set var [$m entrycget $i -variable]}
    puts $out [format "%s\t%s\t%s\t%s\t%s" $path $t $lab $cmd $var]
    if {$t eq "cascade"} {
      set sub [$m entrycget $i -menu]
      catch {eval [$sub cget -postcommand]}
      dumpmenu $sub "$path > $lab" $out
    }
  }
}
after 3500 {
  global ds9 ogfui
  set out [open $::env(MENU_OUT) w]
  set bar $ds9(catalog_frame).menubar
  foreach w {workflow tools} {dumpmenu $bar.$w.m [$bar.$w cget -text] $out}
  foreach t $::ogf::tabs {
    foreach id [::ogf::reg::on_tab $t] {
      if {[info exists ogfui(menu,$id)]} {dumpmenu $ogfui(menu,$id) "$t/$id" $out}
    }
  }
  close $out
  exit
}
