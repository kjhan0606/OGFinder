# OGFinder license audit for commercialization

> **Naming note (2026-10):** the product family is now called Astrafex ("ds10-web" = Astrafex Web, the stand-alone OGFinder = Astrafex Desktop, package `ds10core` = Astrafex Core). This document is kept as written at the time; see Astrafex Web `docs/naming.md`.

**Status: engineering audit, NOT legal advice.** Prepared 2026-10-03 (Asia/Seoul) from the repository tree at `/workspace/OGFinder` (HEAD at start of audit: `dfd0593a5`; the tree had uncommitted edits by another worker, which were not touched), the Python environment `/workspace/ogf_venv` and a few public web pages. Nothing here is a legal opinion. Licence texts are interpreted by courts, not by this document; every "OK" below means "no blocker found by reading the files", not "cleared". **Have a qualified open-source/IP lawyer review before any commercial release.** Items marked **[UNVERIFIED]** were not confirmed from a primary source in this audit.

Goal of the audit: commercialization of OGFinder (customized SAOImageDS9 8.7 + OGFinder plugins + `ogfkit` Python + C code + `ai_bridge`) as (A) a **stand-alone binary** and (B) a **web SaaS** ("ds10-web", custom, not JS9).

---

## 0. Scope, method and limits

* Read: every `LICENSE`/`COPYING`/`copyright`/`license.terms` file in the tree (about 130 found with `find`), source-file headers of the vendored packages and of OGFinder-written code, `ds9/unix/Makefile` (what is linked), `ldd`/`strings` of `bin/ds9` and `bin/ds9_sextract`, `pip` metadata of `/workspace/ogf_venv` (132 distributions; the venv was created with `--system-site-packages`, `ogf_venv/pyvenv.cfg`) and of four side venvs (`eazy_venv`, `sed_venv`, `lens_venv`, `acs_venv`), `ai_bridge/agent_cli.py`, `docs/*.md`, git history.
* Web checks (public pages, fetched 2026-10-03): GALFIT FAQ/home page, SExtractor/PSFEx licence pages, ESA Gaia licence, MAST data-use policy, SDSS image-use policy, Rubin DP1/data-policy pages, Claude Code legal page, search summaries for Codex CLI / Gemini CLI / CDS / ZTF / DAOPHOT / ZOGY.
* **Not done:** no automated software-composition scan (ScanCode/FOSSA/Black Duck) was run; no per-file scan of the ~6300 tcllib / 2900 ttkthemes / 2300 tcl / 2700 tkimg files (only top-level and spot-checked headers); the Windows/macOS packaging recipes (`win/`, `macos/`, `docs/windows_macos_build.md`) were not audited; **`ds10-web` is not present in this tree or elsewhere on the box (searched `/workspace` and `find / -maxdepth 3 -iname '*ds10*'`), so the web code base itself is entirely [UNVERIFIED]**; only `docs/phase2_web_design.md` describes it.
* A "not found" in a search is evidence of absence only for the strings searched (listed in section 4).

---

## 1. Executive summary (top risks)

1. **GPL-3.0 core (blocker for closed-source binary distribution).** `LICENSE` is the GPL-3.0 text; `copyright` says SAO's code is under "version 3" of the GPL (no "or later" wording). OGFinder modifies ds9 core files (`ds9/library/ds9.tcl`, `frame.tcl`, `menu.tcl`, `layout.tcl`, `fits.tcl`, `mview.tcl`, git log of author Juhan Kim) and adds `ds9/library/ogf_*.tcl` that run inside ds9's Tcl interpreter. `bin/ds9` is therefore a GPL-3.0 combined work: if you **convey** it, you must give recipients the complete corresponding source under GPL-3.0 and may not add restrictions. SAO and others hold the copyright, so **dual-licensing / relicensing the core is not available to you**.
2. **SaaS is a different regime but not a free pass.** GPLv3 does not treat "mere interaction ... through a computer network" as conveying (`LICENSE` lines 99-101), so running ds9/OGF code on your own servers does not by itself trigger the source obligation. But (a) anything you **ship to the browser or the user's machine** (JS, WASM builds of SEP/`ds9_sextract`, Pyodide bundles, "local mode", companion app; see `docs/phase2_web_design.md` 13.10/13.13) *is* distribution; (b) any ds10-web code ported or translated from ds9/tksao C++/Tcl is a derivative work; (c) AGPL components (PyMuPDF is installed) would add a network clause.
3. **Data licences are the bigger SaaS risk.** Gaia is CC BY-NC 3.0 IGO per ESA; the Digitized Sky Survey and Guide Star Catalog are copyrighted with commercial use prohibited without permission (MAST); Rubin DP1 is only for data-rights holders and may not be served to others; VizieR catalogues carry per-catalogue terms; SDSS/MAST-mission data are mostly public domain. A commercial service that proxies, caches or bundles these needs per-source clearance.
4. **AI CLIs.** Claude Code may not be used by a third-party product with Free/Pro/Max subscription credentials on behalf of users, and hosting it requires Anthropic's Commercial Terms; Codex CLI and Gemini CLI are Apache-2.0 but their services have their own terms; `agy` and `grok` terms are [UNVERIFIED]. For SaaS use API keys under commercial terms, not CLI logins.
5. **OpenSSL 1.0.2u is statically linked into `bin/ds9`** (`strings bin/ds9`; `openssl/README`, `tls/Makefile`). It is end-of-life (Dec 2019, unpatched CVEs) and its OpenSSL/SSLeay licence (advertising clause) is generally considered GPL-incompatible without an exception. Replace with OpenSSL 3.x (Apache-2.0) or the platform TLS library.
6. **Copyleft Python dependencies:** `rebound` and `assist` (GPL-3.0+) are imported in-process by `moving/orbit.py`; `PyMuPDF` (AGPL-3.0 or commercial) is installed but no import was found in the repo; `sep` (LGPL-3.0+), `chardet`, `jwcrypto`, `websockify` (LGPL) installed. These constrain how `moving` can be packaged.
7. **GALFIT** is closed-source freeware: the repo correctly does **not** contain it (`docs/multifit.md` line 67 says it "is not in the repository and may not be redistributed"), but it is used as a black-box reference for tests; keep it out of any shipped artifact.
8. **Vendored third-party licence files are inconsistent** (funtools: LGPL-2.1 `COPYING` vs GPL-2+ `copyright`; XPA: MIT `LICENSE` vs GPL-2+ `copyright`; tkhtml1: LGPL-2+ `COPYRIGHT` vs GPL-2 `LICENSE`). Under a GPL-3 product this is workable, but it must be resolved before any attempt to replace the GPL core.
9. **No project-level notices exist:** no `THIRD_PARTY_NOTICES`, no licence header or `LICENSE` for OGFinder-written code (`ogfkit/`, `ai_bridge/`, `plugins/`, ...), no pinned Python requirements file. The ownership/licence of the OGF code (and of AI-assisted code) is undeclared.

Overall: **a closed-source stand-alone product built on this tree is not possible without replacing ds9 (and funtools/xpa/tclfitsy etc.); a GPL-compliant paid product, or a SaaS whose server side stays on your servers, is possible with the clean-ups in section 7.**

---

## 2. The determining questions

### 2.1 Which GPL, exactly?

* `LICENSE` (top level): "GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007".
* `copyright` (top level; ds9 source files say "see copyright notice in `copyright`"): "Copyright (C) 1999-2024 Smithsonian Astrophysical Observatory ... under the terms of the GNU General Public License as published by the Free Software Foundation, **version 3**." No "or any later version". Treat the ds9 core as **GPL-3.0-only**. (The embedded FSF postal address in this notice is outdated; irrelevant legally.)
* `README.md`: "licensed in part under the GNU General Public License, version 3." `ds9/doc/faq.html` (Copyright section): ds9 is ~20 open-source packages "usually GPL, LGPL, or BSD" plus SAO packages under GPL; "as long as you continue to adhere to the provisions of the licenses, you are free to distribute DS9 along with your software".
* `ds9/library/ds9.tcl` line 10: `set ds9(version) {8.7}`.
* Consequence: GPL-3.0-only cannot be combined into one work with GPL-2.0-**only** code. No GPL-2-only file was found (funtools, xpa, tkmpeg/ezmpeg say "version 2 or any later version"), so the mix is consistent, subject to the clarifications in 3.1.

### 2.2 Stand-alone binary distribution (GPLv3 implications)

`bin/ds9` is a single executable with a zip file system attached (`ds9/unix/ds9.zip`, `ds9/doc/faq.html`). It statically links `libast*.a`, `libfuntools.a`, `libxpa.a`, tksao, fitsy, tclfitsy, tkhtml1, tkmpeg, tkagif, tksvg, tkblt, tclxml, tkimg stubs, tls+OpenSSL (`ds9/unix/Makefile` lines 92-125; `strings bin/ds9`). If you convey it (sell, give away, put on a download page, ship in a device or installer) you must, among other things (GPLv3 sections 4-6, 10, 11):

* keep all copyright and licence notices; give a copy of the GPL-3.0 text;
* provide the **Complete Corresponding Source** (all ds9 modifications including every OGF change to `ds9/library/*.tcl`, `ogf_*.tcl`, the build scripts, and the Tcl/Python plugin glue that is part of the same work) with the binary or via a written offer; the source must be under GPL-3.0;
* **not impose further restrictions** (no EULA clause forbidding redistribution, reverse engineering, or competing use of the GPL part; no per-seat licence keys that bar redistribution of the GPL code; a "licence check" in code is possible only if recipients remain free to remove it);
* patent licence (section 11): you grant a licence to patents you hold in the contributed code;
* "Installation Information" (anti-tivoization) only if shipped inside a consumer "User Product" (section 6).

You **may** charge for copies, sell support/updates/hosting/training, ship a signed official build, and use your trademarks (not forced to license the name). You **cannot** keep the ds9-derived code proprietary or sublicense it commercially under different terms.

What can stay proprietary (only if kept legally and technically separate): programs that merely **communicate** with ds9 through argv/stdout/files/XPA and do not share its address space or intimate data structures. In this tree the Python tools are launched as subprocesses (`ds9/library/ogf_util.tcl` `OGFPython`, `ds9/library/ogf_core.tcl` step runner; contract in `docs/architecture.md`, `docs/plugins.md`), and `ogfkit/__init__.py` and `ai_bridge/__init__.py` state that they have no Tk and no global state. Those packages are *plausible* candidates for separate licensing. The Tcl side (`plugins/*/*.tcl`, `ds9/library/ogf_*.tcl`) runs in ds9's interpreter, calls ds9 procs and state, and should be treated as part of the GPL work. Whether the FSF "separate programs" line is satisfied is a legal judgement; this audit only notes the structure.

### 2.3 LGPL linking

| Library | In the binary | Obligation (LGPL-2.1 / LGPL-3.0) |
|---|---|---|
| AST 9.2.14 (`ast/`, LGPL-3.0+ headers) | static `libast.a` in `bin/ds9` | provide source of AST (+ your modifications to it), licence text, and a way to relink (object files or source of the combined work). Since the whole ds9 is GPL-3, shipping its source satisfies this in practice. |
| wcslib-derived files in `ast/wcslib/` (LGPL-2+ per `ast/wcslib/proj.h`); `funtools/wcs` libwcs (LGPL-2+ per `funtools/wcs/wcs.c`) | static | same, LGPL-2.x: allow modification + reverse engineering for debugging of the combined work. |
| SEP (`sep_src/`, LGPL-3.0+) | static in `bin/ds9_sextract` (built by `build_sextract.sh`); Python `sep` extension in venv | same as AST. `ds9_sextract` is a separate executable, but the LGPL still applies to the SEP part; ship `sep_src` sources and build script, or link `libsep` dynamically. |
| tkhtml1 (LGPL-2+ per `tkhtml1/COPYRIGHT`) | static | same. |
| TinyCC binary `compilers/tcc-0.9.25-win32-bin.zip` (LGPL, [UNVERIFIED]: the zip listing contains no licence file) | Windows package only (`ds9/library/util.tcl` line ~751 uses `tcc/tcc.exe` as the funtools filter compiler) | include licence text + source offer, or stop bundling. |
| Python LGPL packages (`sep`, `chardet`, `jwcrypto`, `websockify`) | only if you ship a Python environment | dynamic import is the normal LGPL-friendly pattern; keep notices; allow replacement of the library. |

LGPL becomes a real problem only if you want a **closed** main program. Because ds9 itself is GPL-3, LGPL adds no extra burden on the stand-alone; it matters if the GPL core is ever replaced (then AST/SEP must be dynamically linked or replaced).

### 2.4 SaaS (web) implications

| Topic | Finding |
|---|---|
| GPL-3.0 code run on your server, users interact remotely | Not "conveying" (`LICENSE` lines 99-101), so no source-offer duty under GPLv3 alone. Verify with counsel, in particular for any variant that sends users a program (e.g. a desktop client, plugin, WASM module). |
| AGPL | **PyMuPDF 1.25.4 = AGPL-3.0 or commercial** (metadata of `ogf_venv`). AGPL section 13 would require offering source of the *whole* network-facing work. No import of `fitz`/`pymupdf` was found in the repo (grep), so remove it from the server image or buy the commercial licence. |
| Code delivered to browsers | JS bundles, WASM (SEP, `ds9_sextract` Emscripten builds, Pyodide drivers, "local mode", "companion", `docs/phase2_web_design.md` 13.10) are **distribution**: SEP (LGPL-3), AST, numpy-style libs, and anything derived from ds9 carry their obligations. Source for LGPL/GPL parts must be offered; LGPL relinkability for WASM is awkward. |
| ds10-web derived from ds9 | `docs/phase2_web_design.md` line 490: "ds9 (GPL) code reused in the web viewer would make it GPL as well; clean-room the renderer or accept GPL." Keep a documented clean-room process (specs only: FITS standard, WCS papers, scale algorithms) if you want a non-GPL web viewer. |
| Third-party JS/CSS of ds10-web | [UNVERIFIED] (not in tree). Check for GPL/AGPL front-end libraries (e.g. sky viewers) and for fonts/icons. |
| ds9 `xpa` and `tkblt` etc. | never leave the server in a pure SaaS; only matter for desktop/companion. |
| Data served to customers | See section 5; for SaaS this is usually the harder constraint (NC and rights-holder-only datasets). |
| AI services | Server-side calls must use API keys under commercial terms (section 3.8). |

### 2.5 Attribution/notice obligations (minimum list)

Ship a `THIRD_PARTY_NOTICES` (and about-box/docs text) containing at least: GPL-3.0 text; LGPL-2.1/3.0 texts; Tcl/Tk BSD-style `license.terms` (retain "existing copyright notices ... included verbatim in any distributions"); IJG JPEG acknowledgement ("this software is based in part on the work of the Independent JPEG Group"); libpng, zlib, libtiff notices; SQLite (public domain, no notice needed); BLT/tkblt (`tkblt/LICENSE`) and Tktable (`tktable/license.txt`) texts; tclsignal (`tclsignal/LICENSE`, permission notice must appear in supporting documentation); ERFA/SOFA notice (`ast/erfa/LICENSE`: BSD 3-clause, "please acknowledge that you are using a library derived from SOFA"); MINPACK notice with the required end-user-documentation acknowledgement (`ast/cminpack/CopyrightMINPACK.txt`); awthemes zlib-style notice (`awthemes/LICENSE`); the OpenSSL/SSLeay advertising notices if OpenSSL 1.0.2 is kept (`openssl/LICENSE`); Python packages' BSD/MIT/Apache notices (numpy, scipy, astropy, matplotlib, torch, ... as shipped); Apache-2.0 NOTICE files where present. Data credits: section 5.

---

## 3. Inventory of third-party and bundled code

Legend: **Commercial OK?** Y = permissive, no blocker; Y* = OK with obligations/conditions; C = copyleft (OK only on a GPL-compatible distribution model); N = not for closed commercial use; ? = unclear/[UNVERIFIED]. **Risk**: L/M/H for the commercial plan in this document. Tracked-file counts are from `git ls-files`.

### 3.1 ds9 core and SAO-written packages (GPL-3.0)

| Component (path; tracked files) | Version | Licence and evidence | Comm. OK? | Binary distribution | SaaS | Risk | Mitigation / alternative |
|---|---|---|---|---|---|---|---|
| SAOImageDS9 core Tcl (`ds9/`, 1710 files; `ds9/library/*.tcl`) | 8.7 (`ds9/library/ds9.tcl:10`) | GPL-3.0-only wording: `LICENSE` (GPL-3 text), `copyright` ("version 3"), file headers "see copyright notice in copyright" | C | whole binary under GPL-3; full corresponding source; no extra restrictions | no duty if only run remotely; duty if code shipped to client | **H** | accept GPL (sell service/support/hosted/signed builds); or replace with clean-room viewer; do not dual-license |
| tksao (`tksao/`, 344 files: frame, colorbar, panner, magnifier, tkutil, widget) | 1.0 | SAO headers "copyright"; no own LICENSE; `tksao/tclconfig/license.terms` is only the TEA template | C (GPL-3) | as above | as above | **H** | as above. Renderer algorithms are the part a web port would copy; see 2.4 |
| fitsy (`fitsy/`, 112 files) and tclfitsy (`tclfitsy/`, 15) | 1.0 | SAO headers; `tclfitsy/LICENSE` = GPL-3 text | C | as above | as above | **H** | as above |
| vector, tkwin, tclsignal glue (`vector/`, `tkwin/`, `tclsignal/`) | 1.0 / 1.4.4.1 | `vector/*.C` SAO header; `tclsignal/LICENSE`: Schwartz Computer Consulting permission notice (MIT-like) modified by SAO | Y*/C | keep notice in documentation | n/a | L | keep notice |
| ds9 parsers (`ds9/parsers/*.tcl`) generated by taccle/fickle (`taccle/COPYING`, `fickle/COPYING` = GPL-2) | - | generated files say the output is not covered by the tool's GPL ("author's license"; e.g. `ds9/parsers/alignparser.tcl`) | Y (generated output) | none beyond ds9's | none | L | build tools are not shipped |
| OGFinder modifications of ds9 files (`ds9/library/ds9.tcl, frame.tcl, menu.tcl, layout.tcl, fits.tcl, mview.tcl`; git author Juhan Kim) and `ds9/library/ogf_*.tcl`, `ds9_*.py`, `ds9_sextract.c` | - | derivative of GPL-3 files; new files carry no licence header | C | must be released under GPL-3 with the binary | n/a unless shipped | **H** | add GPL-3 headers to the modified/new Tcl; put proprietary logic in separate-process Python tools (2.2) |
| Staging copy `ds9/unix/mntpt/` and `ds9/unix/ds9.zip` | - | generated copy of `ds9/library`, themes, tcllib subsets | - | same as source | - | L | exclude from the repository (build output) |

### 3.2 Tcl/Tk and Tcl extension ecosystem

| Component | Version | Licence and evidence | Comm. OK? | Binary distribution | SaaS | Risk | Mitigation |
|---|---|---|---|---|---|---|---|
| Tcl (`tcl8.6/`, 2325 files) | 8.6.15 (`tcl8.6/generic/tcl.h`) | Tcl/BSD-style `tcl8.6/license.terms` ("for any purpose ... notices retained ... included verbatim") | Y | include `license.terms` | none | L | - |
| Tk (`tk8.6/`, 933) | 8.6.15 | `tk8.6/license.terms` (same terms) | Y | include notice | none | L | - |
| tcllib (`tcllib/`, 6296 files; only base64, ftp, log, math, textutil, ... subset is in `ds9.zip`) | 1.20 per `tcllib/DESCRIPTION.txt` | `tcllib/license.terms` (Tcl terms); sub-module licences (`tcllib/modules/snit/license.txt`, `modules/amazon-s3/LICENSE.txt`, ...) not individually reviewed | Y* | notice for shipped modules | none | L | ship only needed modules; review `snit`/others if shipped |
| tklib (`tklib/`, 1127 files; tooltip) | unverified | `tklib/license.terms`; `modules/{wcb,tablelist,mentry,scrollutil}/COPYRIGHT.txt`, `ctext/LICENSE` | Y* | notice | none | L | ship only the tooltip module actually used |
| itcl 4.3.0, thread 2.8.10, tdbc 1.1.9, sqlite 3.45.3, libtommath | as named | `tcl8.6/pkgs/*/license.terms`; SQLite public domain (`tcl8.6/pkgs/sqlite3.45.3/compat/sqlite3/sqlite3.c`, version string 3.45.3); libtommath `tcl8.6/libtommath/LICENSE` (public domain/Unlicense) | Y | none/notice | none | L | - |
| tkimg (`tkimg/`, 2723) | 2.0.1 | `tkimg/license.terms` (Tcl-style) | Y | notice | none | L | - |
| zlib (`tkimg/compat/zlib`, `tcl8.6/compat/zlib`, `lib/zlibtcl1.3.1`) | 1.3.1 | zlib licence `tkimg/compat/zlib/LICENSE` | Y | none required (acknowledgement appreciated) | none | L | - |
| libpng | 1.6.44 (`tkimg/compat/libpng/png.h`) | PNG Reference Library License v2 (`tkimg/compat/libpng/LICENSE`) | Y | notice in docs | none | L | - |
| libtiff | 4.7.0 (`tiffvers.h`) | BSD-like (`tkimg/compat/libtiff/LICENSE.md`) | Y | notice | none | L | - |
| libjpeg (IJG) | "9" series, jpegtcl 9.6.0 | IJG licence (`tkimg/compat/libjpeg/README`) | Y | must acknowledge "based in part on the work of the IJG" | none | L | add acknowledgement |
| tls | 1.6.7 (`tls/configure`) | `tls/license.terms` (Tcl-style, Matt Newman) | Y | notice | none | L | - |
| **OpenSSL** (`openssl/`, 2287 files) linked statically into `bin/ds9` | **1.0.2u** (`openssl/README`; `strings bin/ds9`) | OpenSSL License + original SSLeay licence (`openssl/LICENSE`), advertising clause | **? / N for GPL combination** | licence text + advertising acknowledgement; **EOL**, GPL-incompatibility without exception | n/a (client side) | **H** | upgrade to OpenSSL 3.x (Apache-2.0, GPLv3-compatible) or system TLS; do not ship 1.0.2 |
| tkblt (`tkblt/`, 136) | 3.2 | `tkblt/LICENSE`: BLT/George Howlett permissive (MIT-style) + SAO modification | Y | notice | none | L | - |
| Tktable | 2.12 | `tktable/license.txt` (Tcl-style); `lib/Tktable2.12/license.txt` | Y | notice | none | L | - |
| tkcon | 2.7 | `tkcon/docs/license.terms` (Tcl-style, "bourbonware" note in `tkcon/tkcon.tcl`) | Y | notice | none | L | - |
| tclxml, tclxmlrpc | 3.2 / 1.0 | `tclxml/LICENSE` ("Explain": free for any purpose, notice must be kept). tclxmlrpc: only `tclxmlrpc/tclconfig/license.terms` (TEA template) found: **[UNVERIFIED]** | Y*/? | keep notice | none | L-M | confirm tclxmlrpc origin/licence |
| tksvg (`tksvg/`) incl. nanosvg | 0.7 | `tksvg/license.terms` (Tcl-style); `tksvg/generic/nanosvg.h` (Mikko Mononen, zlib-style permission), `nanosvgrast.h` (Shemanarev antigrain) | Y | notice | none | L | - |
| **tkhtml1** (`tkhtml1/`) | 1.0 | `tkhtml1/COPYRIGHT`: LGPL-2+ (D. R. Hipp); `tkhtml1/LICENSE`: GPL-2 text, inconsistent | C (LGPL/GPL) | LGPL obligations; fine under GPL-3 | none | M | decide which licence applies; consider removing if help-browser use is minor |
| **tkmpeg** + ezmpeg (`tkmpeg/`) | 1.0 | `tkmpeg/LICENSE`, `doc/license.txt`: GPL-2 text; `tkmpeg/ezmpeg.c` header: GPL **2 or later** | C | fine under GPL-3 (via "or later") | none | M | MPEG-1 encoder for movies; can be dropped for non-GPL builds |
| **tkagif** (`tkagif/`) | 1.0 | `tkagif/LICENSE`: GPL-2 text; SAO source header only | C/? | as above | none | M | verify upstream; drop for non-GPL builds |
| tclzipfs | 1.0 | `tclzipfs/license.terms` (Tcl-style) | Y | notice | none | L | - |
| awthemes | - | `awthemes/LICENSE`: zlib-style (Brad Lanam) for awthemes; per-theme note lists arc/breeze/breeze-dark as GPLv3 *originals* not used by awthemes | Y | notice | none | L | - |
| ttkthemes (`ttkthemes/`, 2956 files; subset in `ds9.zip`) | - | `ttkthemes/docs/licenses.rst`: mix of Tcl-style, GPLv2+, GPLv3; "only licence under which all themes are available together is GPLv3"; `radiance`, `clearlooks` themes are Tcl-style (`ttkthemes/ttkthemes/themes/*/LICENSE`) | C | fine under GPL-3; trim unused themes (only clearlooks, radiance used by `ds9/library/ds9.tcl` lines ~230-255) | none | L-M | ship only Tcl-licensed themes if GPL is later removed |
| scidthemes (`scidthemes/`) | - | `scidthemes/LICENSE`: Tcl-style (Uwe Klimmek) | Y | notice | none | L | - |

### 3.3 Astronomy C/C++ libraries

| Component | Version | Licence and evidence | Comm. OK? | Binary distribution | SaaS | Risk | Mitigation |
|---|---|---|---|---|---|---|---|
| AST (`ast/`, 625 files) static `libast*.a` | 9.2.14 (`ast/config.h`) | LGPL-3.0+ (`ast/ast.h` header; `ast/COPYING.LESSER`); `ast/COPYING` = GPL-3, `ast/COPYING.LIB` = LGPL-2 text also present | Y* (LGPL) | LGPL-3: source + relink | none server-side; WASM/client = distribution | M | keep dynamic or relinkable if GPL core is replaced; alternative WCS code: wcslib (LGPL-3), astropy `wcs` |
| ERFA/SOFA derivative (`ast/erfa/`) | - | BSD-3-clause style + SOFA acknowledgement (`ast/erfa/LICENSE`) | Y | notice + SOFA wording | none | L | - |
| PAL (`ast/pal/`) | - | LGPL-3.0+ headers (`ast/pal/palAddet.c`) | Y* | with AST | none | L | - |
| cminpack/MINPACK (`ast/cminpack/CopyrightMINPACK.txt`) | 1999 UChicago | BSD-style with end-user-documentation acknowledgement clause | Y* | add acknowledgement | none | L | - |
| wcslib-derived (`ast/wcslib/*`) | wcslib 2.9-based | LGPL-2+ (`ast/wcslib/proj.h`) | Y* | LGPL | none | L | - |
| libwcs (`funtools/wcs/`) | 2016 | LGPL-2+ (`funtools/wcs/wcs.c`; `funtools/wcs/COPYING` = LGPL-2.1) | Y* | LGPL | none | L | - |
| Funtools (`funtools/`, 807 files) static `libfuntools.a` | - | `funtools/copyright`: GPL-2+; `funtools/COPYING` = **LGPL-2.1** (inconsistent); bundled `funtools/util/zlib-1.2.3` (old zlib, `funtools/util/zlib-1.2.3/contrib/dotzlib/LICENSE_1_0.txt` Boost) and `funtools/gnu` (GPL-2+, `funtools/gnu/copyright`) | C | GPL-2+ is fine under GPL-3 | none | M | zlib 1.2.3 is outdated (CVE history); funtools is used for region/filter support: evaluate replacing |
| XPA (`xpa/`, 186) static `libxpa.a` | - | `xpa/LICENSE`: **MIT** (Smithsonian 2014-2016); `xpa/copyright`: GPL-2+ (1999-2013); source headers carry only "Copyright (c) 1999-2003 SAO" | ? (MIT likely newer, dual statements) | notice | none | L-M | confirm upstream licence in release; MIT would permit non-GPL use |
| Rice/HCOMPRESS/PLIO decompressors in `fitsy/` (`ricecomp.c`, `hdecompress.c`, `pliocomp.c`) | CFITSIO-derived | headers: R. White (STScI) code "made available for use in CFITSIO in July 1999"; PLIO by D. Tody (NRAO); **CFITSIO/HEASARC permissive licence text not included in the tree [UNVERIFIED]** | Y* | add the CFITSIO/HEASARC and STScI notices | none | L | add notices; CFITSIO itself is not bundled |
| HEALPix projection code (`fitsy/hpx.C`) | from WCSLIB 4.7 | header: "modified from the original authored by Dr. Mark Calabretta ... WCSLIB under GNU GPL version 3" | C | GPL-3 | none | L | covered by ds9 GPL-3 |
| **SEP** C library (`sep_src/`, 13 files) | [UNVERIFIED] upstream version (PyPI `sep` 1.4.1 in venv) | LGPL-3.0+ headers in `extract.c, lutz.c, background.c, ...` ("Copyright 1993-2011 Emmanuel Bertin ... Copyright 2014 SEP developers"; "was scan.c / extract.c in SExtractor"); `sep_src/overlap.h` BSD-3. Upstream also advertises MIT/BSD alternatives [UNVERIFIED] | Y* | LGPL-3 for `ds9_sextract`; provide sources + build script (`build_sextract.sh`) | server-side: none; WASM/client = distribution | M | request/confirm MIT/BSD option from upstream; or call `libsep` dynamically |
| `ds9_sextract.c` (`ds9/library/ds9_sextract.c`) | OGF | no licence header; links SEP (static objects) + system CFITSIO 4.6.2 dynamically (`ldd bin/ds9_sextract`; `build_sextract.sh`) | Y* | CFITSIO (HEASARC licence, permissive) is not bundled: if shipped, add its notice | none | L-M | add header; becomes LGPL-linked via SEP |

### 3.4 Other bundled binaries and assets

| Item | Licence/evidence | Comm. OK? | Risk | Note |
|---|---|---|---|---|
| `compilers/tcc-0.9.25-win32-bin.zip` (TinyCC 2009 binary) | LGPL (TinyCC), **no licence file in the zip listing** [UNVERIFIED] | Y* | M | used by Windows filter compiler (`ds9/library/util.tcl`); include LGPL text/source offer or drop |
| Colormaps (`ds9/cmaps/*.sao,*.lut`; `mpl_*`, `h5_*`, `cubehelix*`, `gist_*`, `inferno/magma`) | provenance per map not documented in tree; matplotlib maps are PSF-style/CC0-type [UNVERIFIED] | Y? | L | document source of each colormap |
| Icons (`ds9/icons/`, 217 files) and docs images | provenance not documented [UNVERIFIED] | ? | L-M | confirm they are SAO-authored |
| Fonts | none bundled in the tree | - | - | system fonts via fontconfig/freetype (dynamic) |
| Dynamic system libs of `bin/ds9` (`ldd`): X11, Xft, fontconfig, freetype (FTL/GPL-2 dual), libxml2 (MIT), libpng, zlib, libbz2, brotli, liblzma, expat | system packages | Y | L | when you ship a bundle (AppImage/DMG/installer), you also redistribute these; add their notices |
| `ast/sun210.pdf, sun211.pdf, sun211.htx_tar` (Starlink docs, 15+ MB) | Starlink/STFC documentation, terms [UNVERIFIED] | ? | L | exclude documentation from product repo/package unless licence confirmed |

### 3.5 OGFinder-written code (undeclared licence)

| Component | Evidence | Third-party content found | Comm. OK? | Risk | Action |
|---|---|---|---|---|---|
| `ogfkit/` (37 files), `ai_bridge/` (30), `plugins/` (224), `moving/`, `sed_fit/`, `sed_adapters/`, `psf_deconv/`, `lsbg/`, `photo_z/`, `galaxy_morph/`, `morphometry/`, `icl/`, `star_finder/`, `psf_phot/`, `crowded_phot/`, `sersic_fit/`, `bulge_disk/`, `parallel/`, `ai_merge/`, `scripts/`, `tools/` | git authors "Juhan Kim"/"kjhan0606"; **no copyright header and no LICENSE file in these trees** (grep of headers; only `sep_src/` and ds9-derived files carry notices); `ogfkit/__init__.py` and `ai_bridge/__init__.py` document no Tk/global state, `ai_bridge/__init__.py` "contains NO models" | algorithms from papers (docstrings cite Peng 2002/2010, ZOGY 2016, HelioLinC, Ciotti & Bertin 1999); no copied third-party code found (section 4) | Y (you own it, subject to employer/contributor and AI-assistance questions) | M | add LICENSE + SPDX headers; decide licence per package (GPL for Tcl glue; permissive/proprietary for pure-Python libs that stay separate); record authorship and any employer/AI-tool terms [UNVERIFIED] |
| Git history: 5262 commits by William Joye (SAO), 259 by Juhan Kim, few others (Diab Jerius, Ken Glotfelty, Doug Burke, liujunyan, oldherl) | `git log --format=%an` | - | - | - | contributors outside SAO to ds9 core hold their own copyright under GPL-3 |
| Trained models / data files | `ai_merge/data/checkpoints/mlp_best.pt`, `lsbg/data/checkpoints/svm_lsbg_best.joblib`, `lsbg/data/training/labeled.tsv` (709 rows, DES coordinates), `photo_z/data/checkpoints/photo_z_mdn_best.pt`, `photo_z/data/sdss_specphoto.h5` (1.1 MB, from `scripts/fetch_photoz_training.py`, SDSS DR18 SkyServer), `moving/data/realbogus_model.json`, `results/montes2022/*` | derived from public survey data; training-set provenance and licence not recorded in the files | Y* | M | record provenance; see section 5 |

### 3.6 Python dependencies (`/workspace/ogf_venv`, Python 3.13.5, system-site-packages enabled)

`pip list`: 128 entries (132 distributions seen by `importlib.metadata`, which also counts system packages); licence data below is the *package metadata* (License-Expression / License / Classifiers). The venv is **not** a product dependency manifest: there is no `requirements*.txt`/`pyproject.toml` in the repo, and Python is not bundled (Tcl finds it via `OGFINDER_PYTHON`, `ds9/library/ogf_util.tcl`). Which of these a shipped product needs is therefore undecided.

**Flagged packages**

| Package (version) | Licence (metadata) | Used by repo? | Comm. OK? | Binary distribution | SaaS | Risk | Mitigation |
|---|---|---|---|---|---|---|---|
| **rebound 4.6.0** | GPL-3.0+ (`GNU General Public License v3 or later`) | yes: `moving/orbit.py` line 67 `import rebound, assist`; `moving/util.py` | C | any bundle that includes `moving` + rebound becomes GPL-3 as a whole | server-side OK; no distribution | **M-H** | run orbit integration as a separate CLI process, keep optional, or replace (own integrator; Apache/BSD alternatives to be evaluated) |
| **assist 1.2.3** | GPL-3.0+ | yes (`moving/orbit.py` lines 32-83) | C | as above | as above | **M-H** | same; also needs JPL DE440 data files (public) |
| **PyMuPDF 1.25.4** | **AGPL-3.0 or Artifex commercial** | **no import found** (`grep fitz|pymupdf`); reason for installation unknown | N (closed) / commercial licence | AGPL covers the whole app | **AGPL section 13 network clause**: source of the whole service | **H if used** | uninstall from server images; use `pypdf`/Poppler or buy Artifex licence |
| **sep 1.4.1** | LGPL-3.0+ | yes (53 files import `sep`; many use `import sep_pjw as sep` with fallback) | Y* | dynamic import is OK; notices; allow replacement | none server-side | M | keep as separately replaceable module |
| chardet 5.2.0 | LGPL (v2+) | not imported by repo | Y* | dynamic | none | L | drop if unused (pulled by requests stack) |
| jwcrypto 1.5.6 | LGPL-3.0+ | not imported | Y* | dynamic | none | L | drop if unused |
| websockify 0.12.0 | LGPL-3.0 | not imported by repo (installed for a display/web test setup [UNVERIFIED]) | Y* | dynamic | none | L | if used as the web bridge, ship its source offer |
| certifi 2025.1.31 | MPL-2.0 | indirect | Y* | file-level copyleft: keep source of certifi files | none | L | - |
| tqdm 4.67.1 | MPL-2.0 AND MIT | indirect | Y* | as above | none | L | - |
| torch 2.14.1+cpu | Apache-2.0 AND BSD-2/3 AND BSL-1.0 AND MIT (+ LLVM exception) | yes (55 files: `ai_merge`, `galaxy_morph`, `photo_z`) | Y | ship NOTICE/third-party licences with wheels | none | L | no CUDA here; GPU builds add NVIDIA redistribution terms [UNVERIFIED] |
| playwright 1.63.0 | Apache-2.0 (downloads Chromium/Firefox/WebKit with their own licences) | not imported by repo scripts | Y | do not ship browsers unless needed | none | L | - |
| matplotlib 3.11.2 | PSF-style matplotlib licence (+ bundled DejaVu/STIX fonts) | yes | Y | include licence | none | L | - |
| Pillow 11.1.0 | MIT-CMU (HPND) | yes | Y | notice | none | L | - |
| pyerfa 2.0.1.5 / astropy-iers-data | BSD-3 (ERFA derived from SOFA) | indirect | Y | SOFA acknowledgement | none | L | - |

**Unflagged (permissive) packages**, from the same metadata: BSD-2/3: astropy 8.0.1, astropy-healpix, astroquery 0.4.11, astroscrappy, click, cloudpickle, contourpy, cycler, dask, dask-image, decorator, fsspec, h5py 3.16.0, ImageIO, Jinja2, joblib, kiwisolver, locket, lxml, lz4, MarkupSafe, mpmath, msgspec, networkx, numpy 2.2.4, olefile, openapi-schema-validator, partd, photutils 3.0.0, PIMS, pooch, pyerfa, Pygments, pyvo, python-dotenv, reproject 0.21.0, sbpy 0.6.0, scikit-learn 1.9.1, scipy 1.15.3, SecretStorage, slicerator, sympy, tifffile, toolz, threadpoolctl, webencodings, wrapt, zarr(MIT), idna, cssselect, gyp-next, html5lib; MIT: annotated-types, beautifulsoup4, Brotli, charset-normalizer, donfig, emcee 3.1.6, et_xmlfile, filelock, fonttools, fs, html5lib-modern, iniconfig, jaraco.*, jeepney, jplephem, jsonschema(+path/specifications), keyring, more-itertools, narwhals, numcodecs, openpyxl, packaging (Apache/BSD), platformdirs, pluggy, pydantic(+core/settings), PyAVM, pyee, pyparsing, pytest, PyYAML, redis, referencing, rfc3339-validator, rpds-py, setuptools, six, soupsieve, typing-inspection, urllib3, wheel; Apache-2.0: bcrypt, cryptography (Apache-2.0 OR BSD-3), google-crc32c, jsonschema-path, openapi-spec-validator, pathable, requests, ufoLib2, unicodedata2, zopfli, python-dateutil (BSD/Apache); PSF-2.0: typing_extensions, greenlet (MIT AND PSF-2.0). `numpy`/`scipy`/`torch`/`Pillow` wheels bundle further libraries (OpenBLAS, libgfortran with GCC runtime exception, etc.); collect their `*.libs`/`licenses` directories when you freeze a distribution. `pandas` is used only in a side venv.

**Optional science-code environments** (separate venvs; used through lazy imports or subprocess protocol, `docs/sed_codes.md`, `docs/lensmodel.md`):

| Package | Where | Licence (metadata) | Comm. OK? | Note |
|---|---|---|---|---|
| eazy-py 0.8.7, astro-sedpy 0.4.1 | `/workspace/eazy_venv`; `sed_adapters/native/eazy_native.py` | MIT | Y | templates/filters come from the separate `eazy-photoz` repository [licence UNVERIFIED] |
| astro-prospector 1.4.1, python-fsps 0.5.0, dynesty 2.1.5, corner | `/workspace/sed_venv`; `sed_adapters/native/prospector_fit.py` | MIT / BSD | Y | FSPS Fortran is under MIT per `/workspace/fsps_master/LICENSE`; stellar libraries/isochrones/spectra under `fsps_master/{SPECTRA,ISOCHRONES}` have their own terms [UNVERIFIED] |
| Bagpipes 1.3.6 | `sed_venv`; `sed_adapters/native/bagpipes_fit.py` | MIT upstream [metadata empty, UNVERIFIED] | Y? | optional **MultiNest** sampler: `pymultinest 2.12` is GPLv3 and the MultiNest Fortran code has non-commercial/academic terms [UNVERIFIED]; use Nautilus/dynesty instead |
| CIGALE (`pcigale` CLI) | `sed_adapters/cigale_adapter.py` (via `shutil.which`) | CeCILL-v2 [UNVERIFIED, not installed] | ? | called as a user-installed external program, not bundled |
| lenstronomy 1.14.2, numba, scikit-image | `lens_venv`; `ogfkit/lensmodel.py` (lazy import), `plugins/lensmodel/validation/` | BSD-3 / BSD | Y | - |
| acstools 3.8.2, photutils 3.0.0 | `/workspace/acs_venv`; `plugins/trails/validation/compare_acstools.py` (comparison only, "separate venv") | BSD-3 | Y | not a runtime dependency |
| pysersic, JAX, GalSim, WebbPSF, TinyTim | imported in `plugins/multifit/validation/pysersic_map.py`, `ai_merge/simulation/profiles.py`, `psf_deconv/psf/sim_psf.py` | pysersic/JAX/GalSim/WebbPSF: permissive (BSD/Apache) [UNVERIFIED, not installed]; **TinyTim: STScI software, redistribution terms [UNVERIFIED]** (`psf_deconv/psf/sim_psf.py` only checks for `tiny1/2/3` on PATH) | Y*/? | all optional and not bundled |

### 3.7 External programs invoked or optional (bundled vs called)

| Tool | How the repo uses it | Bundled? | Licence/terms (source) | Commercial implication |
|---|---|---|---|---|
| **GALFIT 3.0.5** (Peng) | tests/validation only: `GALFIT_BIN` env or `shutil.which('galfit')` (`plugins/multifit/tests/test_galfit.py`, `test_advanced.py`, `plugins/multifit/validation/galfit_*.py`); feedme import/export in `ogfkit/galfitio.py` | **No** in the repo (`docs/multifit.md:67`). A copy exists outside it: `/workspace/galfit_bin/` | Closed-source binary "provided to the public at large free of charge", personal modification allowed but modified versions may not be publicly redistributed; source withheld because of embedded copyrighted algorithms (GALFIT home page and FAQ, users.obs.carnegiescience.edu/peng/work/galfit/). Explicit commercial-use wording: [UNVERIFIED] (the "Disclaimers, User Agreement" section was truncated when fetched) | never ship the binary or its libs (`galfit_bin/lib`: ncurses); never make it a required dependency; marketing: say "reads/writes GALFIT-format feedme files", do not imply endorsement |
| **SExtractor 2.28** | not invoked: no `sex`/`source-extractor` call found. Built-in `ds9_sextract` is an own driver on SEP | no | GPL-3.0+ (astromatic.github.io/sextractor/License.html; Debian copyright file), CeCILL period before 2.13 | calling it as an external user-installed program is fine; bundling it makes the bundle GPL-3 |
| **SWarp** | not invoked (no match in repo) | no | licence not established in this audit [UNVERIFIED] (GPL-3 in recent releases, not confirmed here) | - |
| **PSFEx** | not invoked; `plugins/psfex` and `ogfkit/psfext.py` are an independent "PSFEx-like" implementation | no | GPL-3.0+ (psfex.readthedocs.io/en/latest/License.html) | do not copy code; calling it externally is fine |
| **DAOPHOT II / ALLSTAR** | not invoked; `ogfkit/daophot.py` is an independent "-style" implementation (FIND via photutils `DAOStarFinder`, BSD) | no | proprietary, Crown copyright; NRC agreement allows free use by scientists, **forbids distribution of original or modified copies without NRC written permission** (CADC `daospec/agreement.txt`, search summary); no commercial licence found | never bundle; do not port its source |
| acstools `findsat_mrt` | comparison script only | no | BSD-3 | fine |
| Prospector/EAZY/Bagpipes/CIGALE | optional adapters (3.6) | no | see 3.6 | fine if user-installed |
| **Claude Code** (`claude`) | `ai_bridge/agent_cli.py` backend `claude` (`-p --output-format json`), tests `ai_bridge/tests/test_claude_live.py` | no | proprietary; Consumer Terms (Free/Pro/Max) or Commercial Terms (Team/Enterprise/API). Anthropic's page: products that run Claude Code need Commercial Terms, unmodified binary, each end user authenticates with own key/plan; "does not permit third-party developers to offer Claude.ai login ... or to route requests through Free, Pro, or Max plan credentials on behalf of their users"; no pay/resell for end users (code.claude.com/docs/en/legal-and-compliance, 2026-10-03) | desktop: user logs in with own account is OK; **SaaS: do not drive a shared/subscription login from your servers**; use API key under Commercial Terms |
| **Codex CLI** (`codex`) | backend `codex` (`codex exec --sandbox read-only ...`) | no | CLI is Apache-2.0 (github.com/openai/codex, search summary); service use governed by ChatGPT/API terms | similar: per-user login or API key under OpenAI business terms [service terms UNVERIFIED] |
| **Gemini CLI / `agy` (Antigravity CLI)** | backend `agy`/`gemini` (`agy --print`, fallback `gemini -p`) | no | Gemini CLI Apache-2.0 and Google API/Code Assist terms (search summary); **`agy` licence/terms [UNVERIFIED]** | verify before offering as a feature |
| **Grok CLI** (`grok`) | backend `grok` (`--prompt-file ... --output-format json`), `ai_bridge/tests/test_grok_live.py` | no | **[UNVERIFIED]** (xAI terms for the CLI and API: docs.x.ai/build/cli/headless-scripting not reviewed) | verify before offering |
| Generic HTTP/TAP/local-command AI profiles | `ai_bridge/` (`ai_services.example.json`, placeholders `example.invalid`) | no | provider terms | secrets only via env var names (`ai_bridge/agent_cli.py` docstring); payloads leave the machine, which is a privacy/DPA matter for SaaS |

The CLIs are launched as argv lists in an empty temp dir with a minimal environment (`ai_bridge/agent_cli.py` docstring) and prompts never carry credentials; this is a good separation. Output of AI services is added as catalog columns: check each provider's terms on output ownership and on use of outputs to develop competing models [UNVERIFIED].

---

## 4. Provenance check: code copied or ported from restricted-licence sources?

Method: grep of own code trees (`ogfkit ai_bridge ai_merge plugins moving scripts tools psf_phot crowded_phot sersic_fit star_finder psf_deconv sed_fit sed_adapters photo_z galaxy_morph morphometry icl lsbg bulge_disk parallel regression ds9/library`) for "ported from / translated from / adapted from / copied from / derived from / clean room / reimplement / copyright / licence", for tool names (galfit, daophot, sextractor, psfex, swarp, zogy, statmorph, ppxf, find_orb, heliolinc, lenstool, GIM2D, Stetson, Bertin, Peng, Numerical Recipes, github URLs) and reading of the module docstrings and of `docs/{daophot,psfex,multifit}.md`. This is a keyword/heuristic review, **not a code-similarity analysis**; it cannot prove absence of copying and it makes no claim of certainty.

| Area | Evidence for independent work | Evidence of dependency/derivation | Assessment |
|---|---|---|---|
| **SExtractor** | `ds9/library/ds9_sextract.c` header: "SExtractor-like source extraction", "Uses: sep (C library), cfitsio". CLASS_STAR is a simple logistic function of a size ratio (`ds9_sextract.c` lines 1095-1103), not SExtractor's neural network; no `.nnw/.conv/.param/default.sex` handling found in repo; own `--conv-filter default|gauss5x5|mexhat|tophat` | `sep_src/*` is the SEP library, explicitly derived from SExtractor ("Note: was scan.c in SExtractor", "Copyright 1993-2011 Emmanuel Bertin"), distributed with LGPL-3.0+ headers by SEP's authors; output column names follow SExtractor catalogue names (names/format are interfaces) | Derivation is **declared and licensed** through SEP (LGPL-3+; SExtractor itself is GPL-3+ since 2.13, per Debian copyright file and upstream licence page). No evidence of additional SExtractor code in OGF files. Residual risk: confirm upstream SEP licence statement (LGPL-only vs MIT/BSD options) [UNVERIFIED] |
| **DAOPHOT** | `ogfkit/daophot.py` is a numpy/scipy "DAOPHOT/ALLSTAR/DAOMASTER-style" design; FIND uses photutils `DAOStarFinder` (BSD); `docs/daophot.md` records deviations: "`SHARP` is residual-based (not DAOPHOT's definition)"; no mention of Stetson or Fortran sources; no DAOPHOT-style fixed-format files | method names and workflow (PICKPSF, NSTAR, ALLSTAR, SUBSTAR, ADDSTAR, DAOMATCH/DAOMASTER) follow the published DAOPHOT II workflow; the "penny" PSF function is a DAOPHOT-lineage analytic function (published) | No evidence of copied DAOPHOT code. Original DAOPHOT II is proprietary and non-redistributable (see 3.7), so keep this independence documented (design notes, no source consulted) |
| **GALFIT** | `ogfkit/profiles.py`, `ogfkit/models.py`, `ogfkit/multifit.py` implement Sersic/Moffat/King/Nuker/Fourier/bending/truncation from Peng et al. 2002/2010 formulas (docstring citations: `profiles.py:11`, `profiles.py:55`); GALFIT binary is not in the repo | `ogfkit/galfitio.py` parses/writes GALFIT feedme and constraints; `docs/multifit.md` lines 54-130: conventions "measured against the GALFIT 3.0.5 binary (not from the paper where they differ)", e.g. truncation constant "binary 4.95, paper 4.98", pixel-parity tests against GALFIT-rendered images, GALFIT binary "downloaded from the author's page for these tests only" | Black-box behavioural matching, no source copying seen. Formulas and file formats are not copyrightable per se, but (i) GALFIT's own terms for using its binary in a commercial product comparison are unverified, (ii) numerical constants copied from observation of the binary are facts, still note them; (iii) the name GALFIT is the author's mark/identifier: use only descriptively |
| **PSFEx** | `ogfkit/psfext.py`: PCA rank reduction, wings, stacking, "PSFEx-like" labels; own JSON/FITS layout ("PSFEx-like layout" in `docs/daophot.md` line 27) | none beyond naming | No evidence of PSFEx (GPL-3+) code. Do not look at/paste its source if you intend a non-GPL result |
| **ZOGY** | `moving/zogy.py` states "Zackay, Ofek & Gal-Yam 2016, ApJ 830, 27 ... in numpy"; formulas follow paper notation and equation numbers | upstream ZOGY Python code (pmvreeswijk/ZOGY) is MIT-licensed (GitHub); no reference to it in the repo | Independent implementation from the paper; low risk. Patent status of ZOGY [UNVERIFIED] |
| HelioLinC / linking | `moving/tracklet.py`, `moving/nightlink.py` call themselves "simplified HelioLinC-style" and cite Holman 2018 / Heinze 2022 | - | algorithm from papers; no copied code seen |
| Isophote | `plugins/isophote/isophote.py` "IRAF `ellipse` style, engine photutils.isophote" | photutils (BSD) | fine |
| Sersic constants | `ogfkit/models.py:9`, `sersic_fit/utils.py:10`: Ciotti & Bertin 1999 series | math | fine |
| Other | grep hits for "copied from" (23 files) are comments about test decoys etc. (`moving/validation/nightlink_validate.py:98`), not provenance | - | none relevant |

Conclusion: **no evidence was found of source code copied from GALFIT, DAOPHOT, SExtractor or PSFEx in OGFinder-written code**; the one explicit derivative (SEP) is under its own LGPL notice. This is an absence-of-evidence statement from keyword review and reading of docstrings; a similarity scan (ScanCode, SCANOSS, or a code-clone tool against SExtractor/PSFEx/photutils/statmorph/astropy sources) is recommended before release.

---

## 5. Data, catalogues, services and test data

| Source (where used) | Terms found | Comm. OK? | Risk | Obligation / note |
|---|---|---|---|---|
| **Gaia** (`gea.esac.esa.int` TAP; `ogfkit/xmatch.py`, `plugins/xmatch`; `docs/testing.md` live Gaia DR3 test) | ESA: "Gaia data are distributed under the CC BY-NC 3.0 IGO licence"; credit "ESA, Gaia DPAC"; T&Cs give guidelines for commercial use (cosmos.esa.int/web/gaia-users/license, /web/esdc/terms-and-conditions) | **? NonCommercial** | **H (SaaS)** | do not bundle/cache/redistribute Gaia content in a paid product without ESA clearance; if users query Gaia themselves from their own client, the burden shifts but is not eliminated; carry the acknowledgement text |
| **Pan-STARRS1** (`ps1images.stsci.edu`, MAST) | PS1 archive at STScI; acknowledgement text and citations required (outerspace.stsci.edu/spaces/PANSTARRS/overview); MAST: most data public domain | Y* | L | acknowledgement |
| **MAST / HST / JWST / HLA** (`mast.stsci.edu`, `hla.stsci.edu`; `moving/mast.py`, tests) | Public domain mostly, acknowledgement expected; HLSP under CC BY 4.0; **DSS and Guide Star Catalogs are copyrighted: "Commercial, for-profit use ... is prohibited without written permission"** (archive.stsci.edu/publishing/data-use) | Y* / **N for DSS, GSC** | M-H | exclude DSS/GSC (and any DSS image survey via SkyView `skyview.gsfc.nasa.gov`) from commercial features or get permission; acknowledgements per mission |
| **SDSS** (`data.sdss.org`, `skyserver.sdss.org`; `scripts/fetch_photoz_training.py`, `photo_z/data/sdss_specphoto.h5`, `moving/validation/sdss_known_asteroids.py`, spectra tests) | sdss.org image-use policy: public data releases "considered in the public domain"; website images CC-BY; older credits pages mention ARC approval for commercial use [UNVERIFIED] | Y | L-M | keep credits; confirm with ARC if building a commercial product on SDSS-trained models |
| **Rubin/LSST DP1** (`moving/lsst.py`, `docs/moving_objects.md` section on DP1) | DP1 is proprietary: only Rubin data-rights holders (US/Chile scientists, named in-kind members); RSP accounts for rights holders only; derived products shareable only if they cannot recreate proprietary data (dp1.lsst.io/access, rubinobs.org data policy) | **N** for serving to customers | **H (SaaS)** | the code correctly never stores credentials and uses local files supplied by rights holders; a SaaS must **not** host or proxy DP1 images/catalogs for customers; commercial users generally have no data rights |
| **ZTF / ALeRCE / BTS** (`plugins/lightcurves/fetch_ztf.py`, `plugins/lightcurves/tests/data/ztf_sample.json`, `moving/tests/data/rb_heldout_sample.json`) | ZTF public releases require acknowledgement and citation (IRSA ZTF pages); ALeRCE API terms and ZTF commercial-use terms not found [UNVERIFIED]; BTS labels via VizieR J/ApJ/895/32 | ? | M | record origin of the sample JSON; avoid shipping raw light curves in the product |
| **VizieR / SIMBAD / CDS** (`tapvizier.cds.unistra.fr`, `simbad.cds.unistra.fr`; `ai_services.example.json`, `plugins/xmatch`) | SIMBAD: ODbL, cite/acknowledge; VizieR: free for science, **commercial rights depend on each catalogue's publisher/ReadMe** (cds.u-strasbg.fr/vizier-org/licences_vizier.html) | Y* / per catalogue | M | acknowledgement "SIMBAD database, operated at CDS" / VizieR DOI 10.26093/cds/vizier; restrict commercial features to catalogues whose ReadMe allows it |
| JPL Horizons/SBDB (`ssd-api.jpl.nasa.gov`), MPC (`data.minorplanetcenter.net`), IMCCE SkyBoT/SSP (`vo.imcce.fr`) | NASA/JPL data generally public; MPC and IMCCE terms and rate limits [UNVERIFIED] | ? | L-M | read each service's terms/ToS; cache politely; JPL DE440 kernels used by `assist`/`jplephem` are public JPL products |
| NED, IRSA (`ned.ipac.caltech.edu`, `irsa.ipac.caltech.edu`), DES (`desdr-server.ncsa.illinois.edu`), Legacy Survey, HSC/NAOJ, CADC, NVSS/FIRST, VO endpoints | not reviewed | ? | M | [UNVERIFIED]; DES/Legacy Survey/IRSA have acknowledgement requirements |
| **Trained models** (`ai_merge`, `lsbg`, `photo_z`, `moving/data/realbogus_model.json`) | derived from SDSS, DES Y3 LSBG catalogue (Tanoglidis+2021, `scripts/fetch_tanoglidis2021_data.py`), COSMOS HST (`docs/moving_objects.md` line 201) | Y* | M | keep a training-data provenance sheet; check whether a labelled catalogue (DES Y3 LSBG) permits commercial model training [UNVERIFIED] |
| **Test data** | Repo tracks only small samples: `plugins/lightcurves/tests/data/ztf_sample.json`, `moving/tests/data/rb_heldout_sample.json`, `lsbg/data/training/labeled.tsv`, `photo_z/data/sdss_specphoto.h5`, `results/montes2022/*` (derived from JWST SMACS0723 ERO, MAST public data). Large data are fetched at test time into `$OGF_DATA_CACHE` or `/workspace/fits` (HUDF12, M51, ACS `jc8m32j5q_flc.fits`, GALFIT example `gal.fits`/`psf.fits`, SDSS frames, HST J0946) | mostly public; **GALFIT example images belong to the GALFIT distribution** | L-M | do not ship test FITS in the product; do not redistribute the GALFIT example files; keep acknowledgement lines in docs |
| ds9's built-in image/catalogue servers (`ds9/library/cat.tcl`, DSS/SkyView/2MASS/VizieR menus, remote hosts such as `stdatu.stsci.edu`, `skyview.gsfc.nasa.gov`) | inherited from upstream ds9; each remote archive has its own terms (DSS: see MAST row) | ? | M | in the SaaS/commercial build, decide which built-in servers stay enabled; surface the acknowledgements |
| Rubin "First Look" HiPS, CDS HiPS | colour PNG art from press-release TIFFs, "NOT science data" (`moving/lsst.py` docstring) | ? | L | image credits apply; do not redistribute |

---

## 6. Prioritized risk list

| # | Risk | Severity | Applies to | Why |
|---|---|---|---|---|
| 1 | ds9 core GPL-3.0-only; OGF modifies core and adds in-process Tcl | **Critical** | stand-alone; WASM/companion | no closed-source/dual-licensed binary possible; source offer and no-further-restriction duties |
| 2 | Possible reuse/porting of ds9/tksao code into ds10-web (unverified; web code not in tree) | **Critical** if reused | SaaS client code | makes shipped JS/WASM GPL; needs clean-room discipline and evidence |
| 3 | Dataset terms: Gaia BY-NC, DSS/GSC commercial ban, Rubin DP1 rights-holders only, VizieR per-catalogue | **High** | SaaS (also paid desktop with bundled data) | contractual/copyright exposure independent of code licences |
| 4 | OpenSSL 1.0.2u EOL + GPL-incompatible licence | **High** | stand-alone (network TLS) | security and licence defect; easy to fix |
| 5 | AI CLI terms (Claude Code subscription use, hosted CLIs, agy/grok unknown) | **High** for SaaS, Medium desktop | SaaS | account suspension / terms breach; use API under commercial terms |
| 6 | GPL/AGPL Python deps: rebound/assist in-process, PyMuPDF installed | **Medium-High** | stand-alone bundle; SaaS (AGPL) | bundle becomes GPL; AGPL network clause |
| 7 | Undeclared licence/copyright of OGF code (no LICENSE, no headers; AI-assisted authorship; employer/contributor rights) | **Medium-High** | all | cannot grant licences you cannot show you own |
| 8 | Vendored licence inconsistencies (funtools, xpa, tkhtml1, tkmpeg/tkagif) | Medium | stand-alone | needs upstream clarification if GPL core is ever replaced; no notices file |
| 9 | SEP/AST LGPL obligations for client/WASM builds | Medium | local-mode web, companion | relinking, source offer |
| 10 | GALFIT binary use in tests; DAOPHOT/PSFEx/SExtractor naming and independence record | Medium | all (reputational/legal) | keep binary out, keep independence documented, use descriptive naming only |
| 11 | Training data/model provenance (SDSS, DES, COSMOS, ZTF) | Medium | SaaS, models in installer | undocumented rights |
| 12 | Bundled TinyCC LGPL binary, old zlib 1.2.3 in funtools, no CFITSIO/HEASARC notice | Low-Medium | Windows package | notices/updates |
| 13 | Notices and attribution package missing | Medium (easy) | all | compliance hygiene: no `THIRD_PARTY_NOTICES` |
| 14 | Icons/colormaps provenance undocumented | Low | all | verify |

---

## 7. Recommended next steps

**Decisions (week 0-2)**
1. Choose the product model: (a) **GPL-compliant commercial** (publish corresponding source of the desktop app; sell hosted service, support, signed installers, enterprise features as separate services); (b) **replace the GPL core** with a clean-room or permissively-licensed viewer (ds10-web as primary product) and keep desktop ds9-based builds as a GPL "community edition"; (c) negotiate a commercial licence/permission with SAO for non-GPL use of SAO-owned parts (tksao/fitsy etc.; third-party parts such as funtools/xpa/ast still need their own terms). Get legal advice on (a)-(c).
2. Write the architecture rule: **anything that ships to a user's computer or browser is either GPL-clean or separately licensed and separate-process**; server-only code may be proprietary (subject to AGPL/data/AI terms).

**Engineering (stand-alone)**
3. Replace OpenSSL 1.0.2u with OpenSSL 3.x or platform TLS (`tls/`, `openssl/`, `ds9/unix/Makefile`).
4. Add GPL-3 headers/`LICENSE` to OGF-modified ds9 files and to `ds9/library/ogf_*.tcl`, `plugins/*/*.tcl`; publish them with the binary; keep a build script that reproduces `bin/ds9` (Corresponding Source).
5. Keep `ogfkit/`, `ai_bridge/` and the plugin Python as **separate-process** tools with a documented argv/JSON contract (already the design); do not import ds9 code from them; ship them under a separate, explicit licence of your choice; keep optional GPL pieces (`rebound`, `assist`) outside that package (separate executable or user-installed plugin).
6. Remove PyMuPDF (and unused LGPL packages) from images; add a **pinned dependency manifest** (`requirements.lock`/`pyproject.toml`) with licence report (`pip-licenses`) in CI; fail the build on GPL/AGPL/NC/unknown.
7. Resolve upstream licences: funtools, XPA (MIT vs GPL-2+), tkhtml1, tkagif/tkmpeg, tclxmlrpc, tkcon; add the CFITSIO/HEASARC and STScI notices for the fitsy decompressors; replace `funtools/util/zlib-1.2.3`.
8. Drop optional GPL-only components (tkmpeg, tkagif, tkhtml1, GPL ttk themes, tcc) from any build intended to be relicensed later.
9. Provide `bin/ds9_sextract` as a separate executable with SEP either dynamically linked or with object files/sources and the build script (`build_sextract.sh`) shipped; confirm the SEP licence options with its maintainers.

**Web SaaS**
10. Keep ds10-web independent: clean-room process document (specs used, who wrote what, no ds9 source open during implementation); licence/SBOM scan of its JS dependencies (not in this tree).
11. Do not ship WASM/Pyodide bundles that include SEP/AST/`rebound` until the source-offer/relinking approach is settled (`docs/phase2_web_design.md` 13.13 already requires a licence review before M14).
12. Data policy layer: per-source allow-list with commercial flag; block or gate Gaia (NC), DSS/GSC, Rubin DP1, uncleared VizieR catalogues for paid tiers; user-supplied data stays the user's (terms of service: no licence claim on user uploads, DPA for AI providers); show required acknowledgements in the UI/reports.
13. AI features: API-key based providers under commercial terms; never route end-user subscription credentials; document data leaving the server; confirm `agy` and `grok` terms.

**Repository hygiene**
14. Add top-level `LICENSE` pointers, `THIRD_PARTY_NOTICES.md` (generated from section 3), SPDX headers for OGF code, `docs/provenance.md` (independent-implementation statements for DAOPHOT/PSFEx/GALFIT/ZOGY/SExtractor), data provenance sheet for models and test samples.
15. Run a full SCA scan (ScanCode toolkit or equivalent) on the tree and on the built artifacts; run a code-similarity check of `ogfkit/` and `moving/` against SExtractor, PSFEx, photutils, statmorph and ZOGY sources.
16. Re-run this audit before each release; extend it to `win/`, `macos/`, the installer, and ds10-web when it exists.

---

## 8. Open / unverified items (checklist)

* ds10-web source, front-end dependencies and deployment images (not available).
* Upstream SEP licence options and exact version of `sep_src/`.
* XPA, funtools, tkhtml1, tkagif, tclxmlrpc effective licences; CFITSIO/HEASARC notice text.
* TinyCC bundled zip licence text; Starlink PDF documentation terms; icon and colormap provenance.
* GALFIT explicit commercial-use and comparison clauses (page section not captured); DAOPHOT II terms are from a search summary of the CADC agreement.
* `agy` (Antigravity CLI) and Grok CLI terms; Codex service terms; AI output terms.
* SWarp licence; MultiNest terms; FSPS stellar-library/isochrone terms; CIGALE (CeCILL-v2) and TinyTim distribution terms.
* MPC, IMCCE, ALeRCE, ZTF commercial-use terms; DES/Legacy Survey/IRSA/NED terms; SDSS older "ARC approval for commercial use" statements.
* Whether PyMuPDF or websockify are used outside the repo (e.g. docs build, display bridge).
* Authorship/ownership: employer rights, contributors other than the committer, and licence terms of AI coding tools used to write OGF code.
* Windows/macOS bundles (extra libraries, installers, code signing) and patent exposure of bundled codecs (MPEG-1, GIF/LZW, others).

---

## 9. Disclaimer

This document is an informal technical inventory and risk assessment prepared by an automated assistant for Juhan Kim. **It is not legal advice, does not create an attorney-client relationship, and must not be relied on as a licence clearance.** Licence interpretations (copyleft scope, linking, "separate programs", network use, dataset terms, trademark/naming) are fact- and jurisdiction-dependent and need review by qualified counsel. Terms of third-party websites and services change; web facts were read on 2026-10-03 and may be outdated. Statements marked [UNVERIFIED] are unconfirmed. The absence of a finding does not imply the absence of an issue.
