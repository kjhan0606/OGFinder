# Building OGFinder on Windows and macOS - requirements and expected problems

**Status: NOT VERIFIED.**  OGFinder has been built and tested only on Linux/X11 (this box, Tcl/Tk 8.6 from the bundled sources,
Xvfb for GUI tests).  Nothing below was run on Windows or macOS.  It is derived from reading `BUILD.txt`, the `macos/` and `win/`
configure directories (inherited from SAOImageDS9 8.7) and the OGFinder sources.  Treat it as a checklist of what to look at, not
as a recipe known to work.

## What is inherited from SAOImageDS9 8.7 (upstream build)

`BUILD.txt` lists, unchanged from upstream:

| target | command | requirement stated upstream |
|---|---|---|
| Linux / X11 (also macOS X11, Cygwin) | `unix/configure && make` | autotools, X11, zlib, libxml2, libxslt, Xft, Tcl/Tk dev packages |
| macOS native | `macos/configure && make` -> `SAOImageDS9.app` | Xcode |
| Windows | `win/configure && make` | Cygwin with mingw compilers |

The OGFinder additions (C/C++ packages in `sep_src/`, `star_finder/`, `psf_phot/`, `fickle/`, ... and the Tcl/Python layers) are
wired into the top-level `Makefile`; only the Linux/X11 path of that Makefile has been exercised.

## Components that are Tcl only (should be portable in principle)

`ds9/library/*.tcl` (layout, plugins, table, tile, click chooser, session recorder) is plain Tcl/Tk 8.6.  Points where the code
is POSIX-specific and would need checking on other platforms:

* `ogf_core.tcl` stops background jobs with `exec kill $pid` (`OGFsess_exec kill`): Windows has no `kill` (use `taskkill /PID`).
* Python is started as `[OGFPython]` (`ogf_util.tcl`; env `OGFINDER_PYTHON`, else `python3`): on Windows the interpreter is
  usually `python`/`py -3`; set `OGFINDER_PYTHON` explicitly.
* Paths: the code builds paths with `/` and `file join`; plugin manifests use relative paths.  Not tested with drive letters or
  with spaces in the install path.
* Temp files: `/tmp` is hard-coded in some harness scripts (`scripts/*.sh`, `tools/callgraph.py`), not in the application.
* The shell scripts (`scripts/*.sh`, `scripts/run_all_checks.sh`) need bash; on Windows run them under Cygwin/MSYS2 or WSL.
* Tests need an X server (Xvfb) and xdotool: they do not run natively on macOS (XQuartz would be needed for the X11 build) or
  Windows.  The Python tests (`moving/tests`, `ai_bridge/tests`) are platform independent in principle.
* Multiprocessing: the Python backends use process pools in places (`morphometry`, `parallel`); on Windows and macOS the default
  start method is *spawn*, so every worker entry point must be importable and guarded by `if __name__ == "__main__"`.  Not audited.

## Python environment (all platforms)

Versions installed in the Linux venv used for all measurements: numpy 2.2.4, scipy 1.15.3, astropy 8.0.1, astroquery 0.4.11,
sep 1.4.1, matplotlib 3.11.2, pillow 11.1.0, scikit-learn 1.9.1, emcee 3.1.6, rebound 4.6.0, assist 1.2.3, jplephem 2.24,
reproject 0.21.0, requests 2.32.3, PyYAML 6.0.2, pytest 9.1.1.  Likely trouble spots elsewhere:

* `rebound` / `assist` (orbit fitting, `moving/orbit*.py`) are compiled extensions; wheels exist for the common platforms, but
  `assist` needs ephemeris data files downloaded separately.  Availability for Windows and Apple-silicon was not checked.
* `sep` (source extraction) has wheels for the main platforms; `sep_src/` in the repo is the bundled C code used by ds9 itself.
* Network features (MAST, SkyBoT, Horizons) need outbound HTTPS and write caches below `~/.ds9`.

## macOS specifics to check

* Native build (`macos/configure`) uses Cocoa Tk (`tkmacosx/`, `tk8.6/macosx`): look at menu accelerators (Command vs Control),
  the right-button / middle-button numbering (Tk on macOS reports button 2/3 swapped relative to X11), and that the
  `bind all <ButtonPress>` click-away mechanism in `ogf_pick.tcl` behaves with Cocoa event ordering.
* The OGFinder chip bar / panel layout invariants (catalog table y=181, h=769, info area 154 at `-geometry 1300x950`) were
  measured with X11 fonts; other platforms use different default font metrics, so the numbers will differ.
* Gatekeeper/notarisation for the `.app` bundle and for the Python interpreter it launches.

## Windows specifics to check

* Cygwin/mingw toolchain per `BUILD.txt`; the `win/` directory only carries upstream's configure files.
* Process control (`kill`), line endings in `*.tcl`/`*.sh` (use `git config core.autocrlf false`), `exec` of `.py` files, console windows
  popping up for each helper process, path length limits for the `~/.ds9` cache.

## How to verify when a Windows/macOS machine is available

1. Build ds9 with the platform procedure above; run `bin/ds9 /path/to/m51.fits`.
2. `OGFINDER_PYTHON=<python> python -m pytest moving/tests ai_bridge/tests`.
3. Open the Moving Objects and Catalog panels; run the menu-completeness check `tools/menudump.tcl` + `tools/compare_menus.py`.
4. Record a session and replay it (`scripts/verify_session_replay.py`) - needs bash.
