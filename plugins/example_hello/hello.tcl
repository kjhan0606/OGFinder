# Example plugin Tcl: only what the manifest cannot express.  Everything here goes through ::ogf::* services.
package provide DS9 1.0

proc OGFHelloAbout {} {
    set n [::ogf::cat::nrows]
    ::ogf::status "hello: catalog has $n rows; greeting is '[::ogf::params::get example_hello greeting]'"
}
