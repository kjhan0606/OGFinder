# Licensing of the OGFinder stand-alone (release note)

**Decision (tentative, by the owner, 2026-10-05): the OGFinder stand-alone is released under the GNU General Public License, version 3.**  The top-level
[`LICENSE`](../LICENSE) is the full GPL-3.0 text.  This is an engineering note, not legal advice; the audit behind it is [license_audit.md](license_audit.md)
(sections and file names below refer to it).  Have a lawyer review it before the first public release.

## What is under which licence

| part | licence | note |
|---|---|---|
| SAOImageDS9 8.7 core, tksao, fitsy, tclfitsy and the SAO-written packages | GPL-3.0-only (`copyright`: "version 3") | (audit 2.1, 3.1) |
| OGFinder-written code: `ds9/library/ogf_*.tcl`, modifications of ds9 files, `plugins/*/*.tcl`, `plugins/*/ds10.py`-style launchers, `ogfkit/`, `ai_bridge/`, `scripts/`, `tools/`, the other analysis packages | GPL-3.0 (tentative decision above; the files carry no per-file header yet - to be added) | (audit 3.5) |
| Tcl/Tk 8.6, tcllib, tklib, itcl, thread, tdbc, tkimg, tls, tkblt, Tktable, tkcon, tclxml, tksvg, tclzipfs, scidthemes | Tcl/Tk BSD-style `license.terms` and similar permissive terms; keep the notices | (audit 3.2) |
| zlib, libpng, libtiff, IJG libjpeg, SQLite (public domain), libtommath | permissive (zlib / PNG / BSD-like / IJG acknowledgement) | (audit 3.2) |
| awthemes (zlib-style), ttkthemes (mixed Tcl-style / GPL-2+ / GPL-3; GPL-3 is the one licence for all) | permissive / GPL-3 | (audit 3.2) |
| AST (LGPL-3.0+), PAL (LGPL-3.0+), wcslib-derived files and libwcs (LGPL-2+), ERFA/SOFA (BSD-3 + SOFA acknowledgement), cminpack/MINPACK (BSD-style with acknowledgement) | LGPL / BSD | linked into `bin/ds9`; source is shipped with the GPL-3 source (audit 2.3, 3.3) |
| Funtools (GPL-2+; `COPYING` says LGPL-2.1), XPA (MIT `LICENSE` / GPL-2+ `copyright`), tkhtml1 (LGPL-2+), tkmpeg and tkagif (GPL-2 / 2+) | GPL-2+/LGPL/MIT mix; all combinable under GPL-3 through the "or later" wording | the inconsistencies in the vendored licence files are listed in audit section 1 item 8 |
| SEP C library (`sep_src/`, LGPL-3.0+) used by `ds9_sextract` | LGPL-3.0+ | ship `sep_src` and `build_sextract.sh` (audit 3.3) |
| Python packages of the analysis tools (numpy, scipy, astropy, matplotlib, ...) | BSD / MIT / PSF / Apache-2.0 | installed in the user's environment, not part of this repository |
| `rebound`, `assist` (GPL-3.0+), `sep` Python wheel (LGPL-3.0+), PyMuPDF (AGPL-3.0 or commercial) | copyleft Python packages | imported by the `moving` tools or installed in the development venv; PyMuPDF is not imported by any OGFinder code (audit 1 item 6) |
| OpenSSL 1.0.2u (static in `bin/ds9`) | OpenSSL/SSLeay licence (advertising clause) | **open item**: end-of-life and GPL-incompatible without an exception; replace by OpenSSL 3.x or the system TLS library before a binary release (audit 1 item 5) |

Data sets and AI services have their own terms (audit sections 3.8 and 5); they are not covered by this note.

## ds10core is a separate package, run as a separate process

The shared compute core `ds10core` (regions, pixel table, calibration, map association, forced photometry, script runner, offline bundles) is **not part of this
repository**.  It is the Python package `server/ds10core` of the separate ds10-web repository, which keeps its **own, still undecided licence** (it is not GPL and
imports nothing from OGFinder or ds9).  OGFinder only ships the plugin `plugins/ds10core` (manifest, `ds10.py` launcher, a small Tcl hook): glue code under the
same GPL-3 as the rest of OGFinder that starts `python -m ds10core ...` in a **separate process** with files and argument vectors as the only interface.
Nothing of ds10core is imported, linked or copied into `bin/ds9`, and nothing of ds9, AST or SEP is imported by ds10core.  A user needs a checkout or install
of ds10core next to OGFinder (`DS10_CORE`, see [ds10core.md](ds10core.md)); the plugin reports "not found" otherwise.

## If you convey a binary

Give recipients the GPL-3 text, the complete corresponding source of `bin/ds9` (including the OGFinder changes and the build scripts) and the notices listed
in audit section 2.5; do not add restrictions.  The ds10-web web application is a different product and different repository; it is not covered by this note.
