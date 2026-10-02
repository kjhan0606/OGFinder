# deconv/deconv.tcl -- moved from ds9/library/layout.tcl (procs unchanged; see docs/architecture.md section 6).
# Loaded through the "tcl" field of plugins/deconv/plugin.json.

proc CatalogPanelDeconvolve {algorithm} {
    global current

    if {![::ogf::cat::get psf,has_psf] || ![file exists [::ogf::cat::get psf,file]]} {
	::ogf::cat::set status "No PSF available — build or load PSF first"
	return
    }

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	::ogf::cat::set status "No FITS image loaded"
	return
    }

    set script [CatalogPanelPSFGetScript]
    if {![file exists $script]} {
	::ogf::cat::set status "ERROR: ds9_psf_deconv.py not found"
	return
    }

    set outfile [file join [file normalize ~] .ds9 deconv_result.fits]

    set paramargs {}
    lappend paramargs "--mode" "deconvolve"
    lappend paramargs "--psf" [::ogf::cat::get psf,file]
    lappend paramargs "--algorithm" $algorithm
    lappend paramargs "--output" $outfile

    # Algorithm-specific parameters
    switch $algorithm {
	rl - rl_accelerated {
	    lappend paramargs "--iterations" [::ogf::cat::get psf,param,rl-iterations]
	}
	rl_tv {
	    lappend paramargs "--iterations" [::ogf::cat::get psf,param,rl-iterations]
	    lappend paramargs "--tv-lambda" [::ogf::cat::get psf,param,tv-lambda]
	}
	wiener {
	    lappend paramargs "--wiener-nsr" [::ogf::cat::get psf,param,wiener-nsr]
	}
	tikhonov {
	    lappend paramargs "--tikhonov-lambda" [::ogf::cat::get psf,param,tikhonov-lambda]
	}
	clean {
	    lappend paramargs "--clean-gain" [::ogf::cat::get psf,param,clean-gain]
	    lappend paramargs "--clean-niter" [::ogf::cat::get psf,param,clean-niter]
	    lappend paramargs "--clean-threshold" [::ogf::cat::get psf,param,clean-threshold]
	}
	mem {
	    lappend paramargs "--mem-lambda" [::ogf::cat::get psf,param,mem-lambda]
	    lappend paramargs "--mem-niter" [::ogf::cat::get psf,param,mem-niter]
	}
    }

    ::ogf::cat::set status "Deconvolving ($algorithm) ..."
    update idletasks

    set errfile [file join [file normalize ~] .ds9 psf_stderr.txt]
    OGFSessLog deconv.run auto [list [OGFPython] $script $fn {*}$paramargs] -title {Deconvolution ($algorithm)} -requires [list psf]
    if {[catch {set data [exec [OGFPython] $script $fn {*}$paramargs 2>$errfile]} err]} {
	set stderr_msg ""
	catch {
	    set fd [open $errfile r]
	    set stderr_msg [read $fd]
	    close $fd
	}
	if {$stderr_msg ne ""} {
	    set stderr_lines [split [string trim $stderr_msg] \n]
	    set last_err [lindex $stderr_lines end]
	    ::ogf::cat::set status "Deconvolution error: $last_err"
	    puts "Deconvolution stderr:\n$stderr_msg"
	} else {
	    ::ogf::cat::set status "Deconvolution error: $err"
	}
	return
    }
    catch {file delete $errfile}

    # Load result in a new frame (preserve original)
    if {[file exists $outfile]} {
	CreateFrame
	if {[catch {LoadFitsFile $outfile {} {}} loaderr]} {
	    ::ogf::cat::set status "Deconvolution error: cannot load result: $loaderr"
	    return
	}
	# Apply zscale via DS9's standard scale API
	global scale
	set scale(mode) zscale
	ChangeScaleMode
	::ogf::cat::set status "Deconvolution complete ($algorithm) — result in new frame"
    } else {
	::ogf::cat::set status "Deconvolution complete but output file not found"
    }
}

proc CatalogPanelQuickDeconvolve {} {
    global current

    if {![::ogf::cat::has]} {
	::ogf::cat::set status "Extract sources first"
	return
    }

    set fn [CatalogPanelPSFGetFITS]
    if {$fn eq {} || ![file exists $fn]} {
	::ogf::cat::set status "No FITS image loaded"
	return
    }

    ::ogf::cat::set status "Quick Deconvolve: Step 1/3 — Finding stars ..."
    update idletasks

    # Step 1: Find stars
    CatalogPanelFindStars combined

    if {![::ogf::cat::exists psf,star_indices] || [::ogf::cat::get psf,star_indices] eq {}} {
	::ogf::cat::set status "Quick Deconvolve failed: no stars found"
	return
    }

    ::ogf::cat::set status "Quick Deconvolve: Step 2/3 — Building PSF ..."
    update idletasks

    # Step 2: Build PSF
    CatalogPanelBuildPSF median

    if {![::ogf::cat::get psf,has_psf]} {
	::ogf::cat::set status "Quick Deconvolve failed: PSF build failed"
	return
    }

    ::ogf::cat::set status "Quick Deconvolve: Step 3/3 — Richardson-Lucy deconvolution ..."
    update idletasks

    # Step 3: Deconvolve
    CatalogPanelDeconvolve rl
}


# after-hook of the headless deconvolution step: result into a new frame (as CatalogPanelDeconvolve)
proc OGFDeconvAfter {} {
    set outfile [file join [OGFSessWorkDir] deconv_result.fits]
    if {[file exists $outfile]} {
	CreateFrame
	if {[catch {LoadFitsFile $outfile {} {}} loaderr]} {
	    ::ogf::cat::set status "Deconvolution error: cannot load result: $loaderr"
	    return
	}
	global scale
	set scale(mode) zscale
	ChangeScaleMode
	::ogf::cat::set status "Deconvolution complete ([::ogf::params::get deconv algorithm]) - result in new frame"
    } else {
	::ogf::cat::set status "Deconvolution complete but output file not found"
    }
}
