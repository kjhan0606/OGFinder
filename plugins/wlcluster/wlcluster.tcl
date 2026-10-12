# Cluster weak lensing (plugins/wlcluster): profile and mass files from wlcluster.py.

proc OGFWlclusterFile {name} {return [file join [OGFSessWorkDir] wlcluster $name]}

proc OGFWlclusterAfter {} {
    foreach {k f} {profile_file wl_profile.tsv mass_file wl_mass.json mcal_file wl_mcal.json} {
	set p [OGFWlclusterFile $f]
	if {[file exists $p]} {::ogf::cat::set wlcluster,$k $p}
    }
}
