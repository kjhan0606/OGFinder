# Test infrastructure

`scripts/run_all_checks.sh [--quick] [--only a,b] [--list] [--display :77]` runs every harness and prints a table
(CHECK / PASS|FAIL|SKIP / seconds / last log line); exit status 1 if anything FAILs (SKIP is not a failure).  Logs: `$OGF_CHECK_OUT`
(default `/tmp/ogf_checks.<pid>/NAME.log`).  GUI checks start an Xvfb on the display if none runs there.

| check | what | needs | in `--quick` |
|---|---|---|---|
| `tools_syntax` | `py_compile` of `tools/*.py`, `scripts/*.py`; `tclprocs.tcl`; `compare_menus.py` on the stored dumps | python, tclsh | yes |
| `ai_bridge_tests` | `pytest ai_bridge/tests` | venv | yes |
| `moving_tests` | `pytest moving/tests` | venv | yes |
| `ai_gui` | `scripts/verify_ai_gui.py` (AI panels under Xvfb) | Xvfb | yes |
| `icl_export` | `scripts/verify_icl_export.tcl` smoke on m51 | Xvfb, `/workspace/fits/m51.fits` | yes |
| `click_chooser` | `scripts/verify_click_chooser.tcl` (click procs called directly, layout geometry) | Xvfb | yes |
| `click_xevent` | `scripts/verify_click_xevent.tcl` (real X events with xdotool; SKIP if xdotool absent) | Xvfb, xdotool | yes |
| `cat_api` | `scripts/verify_cat_api.tcl` (unit test of the `::ogf::cat` accessor: get/set/trace/registry) | Xvfb | yes |
| `cat_behavior` | `scripts/verify_cat_behavior.sh` (82 features: exec argv, `catpanel` keys, `.prf` files, session steps and status text, diffed against `scripts/golden/cat_behavior.golden`; `--update` only from deliberately reviewed code) | Xvfb | yes |
| `report_tests` | `plugins/report/tests` (pytest: HTML report contents on m51 and a HUDF F160W crop, cut-outs, escaping, optional PDF) | python deps, `bin/ds9_sextract` | yes |
| `geometry` | `scripts/verify_geometry.tcl` (table y/h and info-area height 181/769/154 in 21 panel states: 6 tabs, extract, search, all time-domain views, review filter, selection, tile/single, detached/reattached; main width 1300 -> 736 -> 1300) | Xvfb | yes |
| `mouse` | `scripts/verify_mouse.tcl` (real xdotool events: notebook tab, chip menu + entry, cascade submenu, table row click, header-click sort, wheel, Tools-menu detach/reattach; 25 checks) | Xvfb, xdotool | yes |
| `review_td` | `scripts/verify_review_td.tcl` (review on moving / transient / detection rows and the All view: columns, filters, tint, steps, save/load, re-run pruning, report from a time-domain view; 55 checks) | Xvfb | yes |
| `review_gui` | `scripts/verify_review_gui.tcl` (review columns, table filter and tint, session steps, save/load round trip, export job, layout invariants; 51 checks) | Xvfb | yes |
| `session_replay` | `scripts/verify_session_replay.sh`: record GUI sessions, export scripts, replay, compare (78 checks) | Xvfb, venv, ~10+ min | no |
| `link_bench` | `moving/validation/link_bench.py` on the injected sets `/workspace/work/inj1-6.json.pkl` (SKIP when missing); when the regenerated sets `/workspace/work/i1/sets/inj1..12.json.pkl` exist (incl. CR-heavy 9-12, held-out 7-8) all twelve run and per-group totals are printed | venv | no |
| `moving_options` | `scripts/verify_moving_options.tcl` (Xvfb): ZOGY parameters of the moving plugin -> CLI flags, default argv unchanged, recorded argv, layout invariants | X server | yes |
| `manifests` | `tools/validate_manifests.py` (static check of all plugin manifests and `cli` templates) + `pytest tools/tests` (13 mistake classes are detected) | venv | yes |
| `cli_templates` | `scripts/verify_cli_templates.tcl` (Xvfb, fake interpreter `scripts/fake_python_for_templates.sh`): the five converted steps vs. their legacy procs - argv, recorder record, resulting catalog, status | X server | yes |
| `moving_session` | self-contained by default (`scripts/verify_moving_session_auto.sh`): a real ds9 runs the Moving Objects GUI steps on the 4 cached BB89 exposures (`scripts/verify_moving_record.tcl`), saves the session script, and `verify_moving_session.sh` replays it and compares detections / movers / tracklets / identified / orbit / transients / light curves / export byte for byte; with `OGF_MOVING_SESSION`, `OGF_MOVING_REF`, `OGF_MOVING_FIELD` set it replays that session instead. SKIPs (77) without the cached exposures, an X server, or (network or `~/.ds9/moving_cache`); ~10 min | cache (+ network when uncached) | no (long) |

Test data location: `OGF_TEST_FITS` (default `/workspace/fits`).  `docs/windows_macos_build.md` lists what is *not* verified on
other platforms.  The generator/compare scripts are in `tools/` (see `tools/README.md`).
