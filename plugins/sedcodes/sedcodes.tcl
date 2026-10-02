# SED codes (plugins/sedcodes): adapters to EAZY / CIGALE / Bagpipes / Prospector; helper windows.

proc OGFSedcodesAfter {} {
    set f [file join [OGFSessWorkDir] sedcodes sedcodes_results.json]
    if {[file exists $f]} {::ogf::cat::set sedcodes,results_file $f}
}

proc OGFSedcodesScript {} {
    return [file join [::ogf::step::plugin_dir sedcodes] sedcodes.py]
}

proc OGFSedcodesPython {} {
    return [OGFPython]
}

proc OGFSedcodesCheck {} {
    set txt {}
    foreach code {eazy cigale bagpipes prospector} {
	set rc [catch {exec [OGFSedcodesPython] [OGFSedcodesScript] --check --code $code --python [::ogf::params::get sedcodes python]} out]
	append txt "== $code ==\n$out\n"
    }
    OGFTextWindow "SED codes: installed codes" $txt
}

proc OGFSedcodesProfiles {} {
    set f [file join [OGFSessWorkDir] sedcodes_ai_services.json]
    catch {file mkdir [file dirname $f]}
    if {[catch {exec [OGFSedcodesPython] [OGFSedcodesScript] --write-profiles $f} err]} {::ogf::status "SED codes: $err"; return}
    ::ogf::cat::set sedcodes,profile_file $f
    OGFTextWindow "SED codes: AI-services profiles" "Profiles written to\n$f\n\nUse it with Analysis > AI Services (registry: profile file) or\n  ds9_ai_bridge.py --mode run --services-file $f --service sedcodes_eazy --task photoz ...\nwith PYTHONPATH pointing at the OGFinder checkout (see docs/sed_codes.md)."
}
