 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1041781.svg)](https://doi.org/10.5281/zenodo.1041781)

![alt-text](http://ds9.si.edu/doc/sun.gif "SAOImageDS9")
# SAOImageDS9

SAOImageDS9 is an astronomical imaging and data visualization application. DS9 supports FITS images and binary tables,  multiple frame buffers, region manipulation, and many scale algorithms and colormaps. It provides for easy communication with external analysis tasks and is highly configurable and extensible via XPA and SAMP.

DS9 is a stand-alone application. It requires no installation or support files. All versions and platforms support a consistent set of GUI and functional capabilities.

DS9 supports advanced features such as 2-D, 3-D and RGB frame buffers, mosaic images, tiling, blinking, geometric markers, colormap manipulation, scaling, arbitrary zoom, cropping, rotation, pan, and a variety of coordinate systems.

The GUI for DS9 is user configurable. GUI elements such as the coordinate display, panner, magnifier, horizontal and vertical graphs, button bar, and color bar can be configured via menus or the command line.

SAOImageDS9 is fully funded by the Chandra X-ray Science Center (CXC) and is licensed in part under the GNU General Public License, version 3. 

## Session recorder / Python pipeline export (OGFinder)

The catalog panel records the analysis steps of a GUI session. *Analysis > Save Session as Python Script...* writes a stand-alone script that re-runs the automatic steps on new FITS files (`--mode pipeline`, batch, `--jobs`, `--resume`) or replays the whole session on the original data (`--mode replay`). See [docs/session_python_script.md](docs/session_python_script.md).

*Analysis > AI Services* provides connection points (not models) for external AI / astronomy services: JSON service profiles (`ai_services.example.json`, `~/.ds9/ai_services.json`), task contracts, cutout generation, response caching, provenance, and a CLI (`ds9_ai_bridge.py`); results are added as catalogue columns and recorded by the session recorder. Keys are read from environment variables only and payloads leave the machine - see [docs/ai_services.md](docs/ai_services.md). Built-in `mock` service for offline tests; the example profiles are templates with placeholder URLs.

## License

Astrafex Desktop (based on the OGFinder tool set: the SAOImageDS9 8.7 based stand-alone with the OGFinder plugins; the product is named "Astrafex", the licensed work and its notices keep the names OGFinder / SAOImageDS9) is released under the **GNU General Public License, version 3** (tentative decision of the
owner); the full text is in [LICENSE](LICENSE).  Bundled third-party components keep their own licences (Tcl/Tk BSD-style, zlib/libpng/libtiff/IJG, LGPL for AST, PAL, wcslib-derived
code, libwcs and SEP, BSD for ERFA/MINPACK, ...); the list, the open items (e.g. the old OpenSSL) and the reasons are in [docs/licensing.md](docs/licensing.md) and
[docs/license_audit.md](docs/license_audit.md).  The shared compute core **Astrafex Core** (Python package `ds10core`) is a separate package of the Astrafex Web repository with its own licence; the
`ds10core` plugin runs it as a separate process and does not include or import it.
