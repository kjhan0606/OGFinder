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
| Funtools | `funtools/copyright` is GPL-2 or later (SAO, 1999-2007). `funtools/COPYING` is the LGPL-2.1 text. Both files stay. | see the notice section below |
| XPA | `xpa/LICENSE` is the MIT text (Smithsonian, 2014-2016). `xpa/copyright` is GPL-2 or later (1999-2013). Both files stay. | see the notice section below |
| tkhtml1 | `tkhtml1/COPYRIGHT` is the GNU Library General Public License, version 2 or later (D. Richard Hipp). `tkhtml1/LICENSE` is the GPL version 2 text. Both files stay. | help browser |
| tkmpeg and tkagif | GPL-2 / GPL-2 or later | combined in this GPL-3 build |
| SEP C library (`sep_src/`) | LGPL-3.0+ | sources stay in the tree. `build_sextract.sh` does not link them. The product extractor is `ogfmeas/sextract.py` |
| fitsy Rice / HCOMPRESS / PLIO | NASA CFITSIO permissive notice | `fitsy/NOTICE`. CFITSIO itself is not bundled |
| funtools `util/zlib-1.2.3` | zlib licence (Jean-loup Gailly and Mark Adler) | zlib 1.2.3, as shipped inside funtools. This note does not upgrade it |
| TinyCC `compilers/tcc-0.9.25-win32-bin.zip` | TinyCC is LGPL. The zip listing has no licence file (audit, unverified) | the Windows package rule does not unpack the zip |
| Python packages of the analysis tools (numpy, scipy, astropy, matplotlib, ...) | BSD / MIT / PSF / Apache-2.0 | installed in the user's environment, not part of this repository |
| `rebound`, `assist` (GPL-3.0+), Python `sep` (LGPL-3.0+), PyMuPDF (AGPL-3.0 or commercial) | copyleft Python packages | rebound and assist run only in the explicit `propagate_assist` child. Python sep is used only when `OGF_USE_SEP=1`. PyMuPDF is not imported |
| OpenSSL | Apache-2.0 for OpenSSL 3, which the build links from the system | the vendored `openssl/` tree is 1.0.2u and is not configured or linked. See below |

Data sets and AI services have their own terms (audit sections 3.8 and 5); they are not covered by this note.

## Notices for the 2026-10-09 build change

This section is an engineering record, not a conclusion that a lawyer has cleared the combination.

The product build stages the system OpenSSL 3 (`scripts/stage_openssl3.sh`) and links `libssl` and `libcrypto` from that stage (`$(prefix)/sys-openssl/lib`). `make` does not run the OpenSSL 1.0.2u configure in `openssl/`. TclTLS 1.6.7 is compiled against the OpenSSL 1.1/3 BIO and DH calls. A checkout still contains the `openssl/` 1.0.2u sources; they are not a build input.

`build_sextract.sh` exits without compiling `sep_src` or writing `bin/ds9_sextract`. Source extraction for the product is `python ogfmeas/sextract.py`. The LGPL sources remain in `sep_src/` and `ds9/library/ds9_sextract.c`.

`funtools/copyright` grants GPL version 2 or any later version. `funtools/COPYING` is the LGPL-2.1 document. `xpa/LICENSE` is the MIT grant dated 2014-2016, and `xpa/copyright` is the GPL-2-or-later grant dated 1999-2013. `tkhtml1/COPYRIGHT` names the Library General Public License version 2 or later, and `tkhtml1/LICENSE` is the GPL version 2 document. None of those files was rewritten. `funtools/util/zlib-1.2.3` remains zlib 1.2.3 under the zlib licence. `compilers/tcc-0.9.25-win32-bin.zip` stays in the source tree and is not unpacked into the Windows application. The Rice, HCOMPRESS and PLIO decompressors carry the NASA CFITSIO notice in `fitsy/NOTICE`.

## Product name and attribution

The product family is called **Astrafex**; the stand-alone program is called **Astrafex Desktop (based on the OGFinder tool set)**.  This is a product name only.  The
licensed work stays "OGFinder" here: the copyright and licence notices of SAOImageDS9 (Smithsonian Astrophysical Observatory), AST, SEP, tcllib, ... and the GPL-3.0
statement above name the works they apply to, are part of the source tree and of the GPL attribution/provenance, and are not edited by the rename.  The repository, the
directory names, the plugin infrastructure and the code that is called OGFinder in comments keep that name.  Nothing in "Astrafex Desktop" implies endorsement by the
SAOImageDS9 project; the ds9 name is used only to say what the program is derived from.  Trademark clearance for "Astrafex" (KR/EU) is pending; the rename is mechanical
and reversible (Astrafex Web `docs/naming.md`).

## ds10core is a separate package, run as a separate process

The shared compute core `ds10core` (regions, pixel table, calibration, map association, forced photometry, script runner, offline bundles) is **not part of this
repository**.  It is the Python package `server/ds10core` of the separate Astrafex Web repository, which keeps its **own, still undecided licence** (it is not GPL and
imports nothing from OGFinder or ds9).  OGFinder only ships the plugin `plugins/ds10core` (manifest, `ds10.py` launcher, a small Tcl hook): glue code under the
same GPL-3 as the rest of OGFinder that starts `python -m astrafex_core ...` in a **separate process** with files and argument vectors as the only interface.
Nothing of ds10core is imported, linked or copied into `bin/ds9`, and nothing of ds9, AST or SEP is imported by ds10core.  A user needs a checkout or install
of ds10core next to OGFinder (`ASTRAFEX_CORE` or `DS10_CORE`, see [ds10core.md](ds10core.md)); the plugin reports "not found" otherwise.

## If you convey a binary

Give recipients the GPL-3 text, the complete corresponding source of `bin/ds9` (including the OGFinder changes and the build scripts) and the notices listed
in audit section 2.5; do not add restrictions.  The Astrafex Web application is a different product and different repository; it is not covered by this note.
