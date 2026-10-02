# Background & noise model (plugins/noisemodel): outputs of noisemodel.py -> catalog keys noisemodel,*, plot / summary windows, frames.

proc OGFNoiseFile {name} {return [file join [OGFSessWorkDir] noisemodel $name]}

proc OGFNoiseAfter {} {
    foreach {k f} {bkg_file noisemodel_bkg.fits rms_file noisemodel_rms.fits sub_file noisemodel_sub.fits json_file noisemodel_summary.json plot_file noisemodel_noise_curve.png curve_file noisemodel_curve.tsv} {
	set p [OGFNoiseFile $f]
	if {[file exists $p]} {::ogf::cat::set noisemodel,$k $p}
    }
    if {[::ogf::params::get noisemodel show-frames]} {catch {OGFNoiseFrames}}
}

proc OGFNoiseFrames {} {
    global current scale
    set orig $current(frame)
    set n 0
    foreach nm {noisemodel_bkg.fits noisemodel_rms.fits} {
	set f [OGFNoiseFile $nm]
	if {![file exists $f]} continue
	CreateFrame
	if {[catch {LoadFitsFile $f {} {}} err]} {::ogf::log ERROR "noisemodel: cannot load $f: $err"; continue}
	set scale(mode) minmax
	ChangeScaleMode
	incr n
    }
    catch {GotoFrame $orig}
    ::ogf::status "Noise model: $n frame(s) opened (background, rms)"
    return $n
}

proc OGFNoiseSummary {} {
    set f [OGFNoiseFile noisemodel_curve.tsv]
    if {![file exists $f]} {::ogf::status "Noise model: nothing yet - run Measure Background and Noise first"; return}
    set fd [open $f r]; set txt [read $fd]; close $fd
    set j [OGFNoiseFile noisemodel_summary.json]
    set head {}
    if {[file exists $j]} {
	set fd [open $j r]; set js [read $fd]; close $fd
	foreach {key} {sigma_pix sky_median masked_fraction} {
	    if {[regexp "\"$key\": (\[-0-9.eE+\]+)" $js -> v]} {append head "$key = $v\n"}
	}
	if {[regexp {"noise_law": \{"alpha": ([-0-9.eE+]+), "beta": ([-0-9.eE+]+)} $js -> a b]} {append head "noise law: sigma_N = sigma_pix * $a * N^$b\n"}
	if {[regexp {"rho1": ([-0-9.eE+]+)} $js -> r1]} {append head "pixel correlation rho(1) = $r1\n"}
    }
    OGFTextWindow "Noise model summary" "$head\n$txt"
}

proc OGFNoisePlot {} {
    set png [OGFNoiseFile noisemodel_noise_curve.png]
    if {![file exists $png]} {::ogf::status "Noise model: no plot yet - run Measure Background and Noise first"; return}
    set w .ogfnoiseplot
    if {[winfo exists $w]} {destroy $w}
    toplevel $w
    wm title $w "Noise curve, pixel autocorrelation, sky flatness"
    catch {image delete ogfnoiseimg}
    image create photo ogfnoiseimg -file $png
    label $w.l -image ogfnoiseimg
    pack $w.l -fill both -expand 1
    ttk::button $w.close -text Close -command [list destroy $w]
    pack $w.close -pady 3
    return $w
}
