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
| `review_gui` | `scripts/verify_review_gui.tcl` (review columns, table filter and tint, session steps, save/load round trip, export job, layout invariants; 51 checks) | Xvfb | yes |
| `session_replay` | `scripts/verify_session_replay.sh`: record GUI sessions, export scripts, replay, compare (78 checks) | Xvfb, venv, ~10+ min | no |
| `link_bench` | `moving/validation/link_bench.py` on the injected sets `/workspace/work/inj*.json.pkl` (SKIP when missing) | venv | no |
| `moving_session` | `scripts/verify_moving_session.sh` replay of a Moving Objects session; needs `OGF_MOVING_SESSION`, `OGF_MOVING_REF`, `OGF_MOVING_FIELD`; network + MAST cache | network | no |

Test data location: `OGF_TEST_FITS` (default `/workspace/fits`).  `docs/windows_macos_build.md` lists what is *not* verified on
other platforms.  The generator/compare scripts are in `tools/` (see `tools/README.md`).
