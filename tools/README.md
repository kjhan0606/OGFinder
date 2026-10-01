# tools/ - generators and comparison scripts used for the architecture work

These are developer scripts, not part of the application.  They were written during the plugin restructuring (see
`docs/architecture.md`); paths are now relative to the repository (override with `OGF_ROOT=/path/to/OGFinder`).
Tcl-proc scanning needs `tclsh` on `PATH` (`tclprocs.tcl`).  None of them is run by the application.  `scripts/run_all_checks.sh` runs a syntax check of every script here (`py_compile`, `tclsh tclprocs.tcl` on one file) and the `compare_menus.py` data check.

| script | what it does | inputs -> outputs |
|---|---|---|
| `tclprocs.tcl` / `tclprocs.py` (library; imported by other tools, no CLI) | list top-level Tcl procs of a file with line ranges (comment block included) | `tclsh tclprocs.tcl FILE` -> `name start end` |
| `audit.py` | feature inventory from `layout.tcl` / `ogf_*.tcl` (proc-name groups, `catpanel(...)` readers/writers, cross-feature calls) | writes `tools/audit.json` |
| `callgraph.py` | proc reference graph over `layout.tcl`, `ogf_*.tcl`, `plugins/*/*.tcl` | writes `/tmp/refs.json`, `/tmp/owner.json` |
| `gen_arch.py` | regenerates the generated tables of `docs/architecture.md` from `audit.json`, `plugin.json` files and `data/menu_baseline.tsv` | **overwrites** `docs/architecture.md` - the hand-written sections of the current document are NOT in the script; do not run it on the current tree (it reproduces the 2026-09 version of the text) |
| `gen_manifests.py`, `gen_params.py` | one-off generators of `plugins/*/plugin.json` and the declarative `params` blocks; the JSON files are now the source of truth | writes `plugins/*/plugin.json` - do not re-run, it would overwrite hand edits |
| `migrate.py`, `rmprocs.py` | move / delete named procs between `layout.tcl` and `plugins/<id>/*.tcl` verbatim (used for the 7 migration stages) | in-place edit of Tcl files |
| `menudump.tcl` | dump every command reachable from the Workflow/Tools menus and plugin chip menus of a running ds9 (`MENU_OUT=file DISPLAY=:77 bin/ds9 img.fits -source tools/menudump.tcl`) | TSV: path, type, label, command, variable |
| `compare_menus.py` | compare two menu dumps (baseline vs new) and list entries that disappeared / are new / were replaced on purpose | `compare_menus.py data/menu_baseline.tsv data/menu_new3.tsv` |
| `steps_sig.py` | print the step signatures (`seq step class argv post`) of exported session scripts, to diff two GUI sessions | `steps_sig.py session1.py session2.py` |
| `data/menu_baseline.tsv` | the 154 menu entries of the original 11 top-level menus (measured on commit `d6d9dff82`) | |
| `data/menu_new3.tsv` | menu dump after the restructuring (measured) | |
| `feature_table.md`, `verification.md` | text blocks included by `gen_arch.py` | |

`steps_sig.py` and `compare_menus.py` are the useful ones to re-run after UI changes; the generators are kept for provenance.
