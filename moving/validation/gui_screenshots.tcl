proc snap {name} {
  update idletasks
  exec /workspace/ogf_venv/bin/python -c "from PIL import ImageGrab; ImageGrab.grab(xdisplay=':77').save('/workspace/fits/$name')"
}
proc step1 {} {
  global ogfmov ds9
  set ogfmov(workdir) /workspace/work/wd_real
  set ogfmov(files) [glob /home/box/.ds9/mast_cache/j8pu38{c7q,caq,ceq,ciq}_flc.fits]
  OGFMovShowMovers
  update idletasks
  after 2500 step2
}
proc step2 {} {
  global ds9
  snap moving_table.png
  set mb $ds9(catalog_frame).menubar.moving
  $mb.m post [winfo rootx $mb] [expr {[winfo rooty $mb]+[winfo height $mb]}]
  update
  after 1500 step3
}
proc step3 {} {
  global ds9 ogfmov
  snap moving_menu.png
  $ds9(catalog_frame).menubar.moving.m unpost
  set ogfmov(tracklet) 19
  OGFMovOrbitDialog
  set ogfmov(orbitlabel) tracklet19
  OGFMovOrbitDone 1 {}
  wm geometry .ogfmovorb +20+60
  update
  after 1500 step4
}
proc step4 {} { snap moving_orbit.png; exit }
proc bgerror {m} { set f [open /workspace/work/shots_dbg.txt a]; puts $f "ERR $m $::errorInfo"; close $f }
after 6000 step1
