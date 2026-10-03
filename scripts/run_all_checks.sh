#!/bin/bash
# Run every test harness of OGFinder and print a pass/fail table.
#   scripts/run_all_checks.sh [--quick] [--only NAME[,NAME...]] [--list] [--display :77] [--keep-going|--fail-fast]
# --quick  skips the long and the network checks (session replay, link benchmark, Moving-Objects session replay).
# Needs: bin/ds9 (built), a python with the test deps ($OGFINDER_PYTHON, default /workspace/ogf_venv/bin/python3), and for the GUI
# checks an X server: an Xvfb on --display (default :77, 1400x1000x24) is started if none is running there.  xdotool is needed for
# the real-X-event check (SKIP if absent).  Logs go to $OGF_CHECK_OUT (default /tmp/ogf_checks.<pid>/NAME.log).
# Exit status: 0 if no check FAILed (SKIP does not fail), 1 otherwise.
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(dirname "$HERE")"
PY="${OGFINDER_PYTHON:-/workspace/ogf_venv/bin/python3}"; export OGFINDER_PYTHON="$PY"
DISP=:77; QUICK=0; ONLY=""; LIST=0; FAILFAST=0
while [ $# -gt 0 ]; do case "$1" in
  --quick) QUICK=1;; --only) ONLY="$2"; shift;; --list) LIST=1;; --display) DISP="$2"; shift;; --fail-fast) FAILFAST=1;; --keep-going) FAILFAST=0;;
  -h|--help) sed -n 2,10p "$0"; exit 0;; *) echo "unknown option $1" >&2; exit 2;; esac; shift; done
OUT="${OGF_CHECK_OUT:-/tmp/ogf_checks.$$}"; mkdir -p "$OUT"
cd "$ROOT" || exit 2
DS9=bin/ds9
FITS="${OGF_TEST_FITS:-/workspace/fits}"
declare -a NAMES KINDS RES SECS NOTE
need_x() { if ! DISPLAY=$DISP xdpyinfo >/dev/null 2>&1; then
    command -v Xvfb >/dev/null || return 1
    (setsid Xvfb $DISP -screen 0 1400x1000x24 >/dev/null 2>&1 &); sleep 2; DISPLAY=$DISP xdpyinfo >/dev/null 2>&1 || return 1; fi; return 0; }
gui() { rm -rf ~/ds9.auto ~/ds9.auto.dir; DISPLAY=$DISP timeout -s KILL "${T:-240}" "$@"; }

# ---- checks: each defines  chk_NAME  (print details on stdout, return 0 pass / 1 fail / 77 skip) and is registered in ALL / LONG
ALL="tools_syntax ai_bridge_tests moving_tests ai_gui agent_gui agent_real icl_export click_chooser click_xevent cat_api cat_behavior report_tests review_gui review_td geometry mouse"
LONG="session_replay link_bench moving_session"
chk_tools_syntax() { local bad=0 f
  for f in tools/*.py scripts/*.py; do "$PY" -m py_compile "$f" 2>&1 || bad=1; done
  command -v tclsh >/dev/null && { tclsh tools/tclprocs.tcl plugins/mask/mask.tcl >/dev/null 2>&1 || { echo "tclprocs.tcl failed"; bad=1; }; }
  "$PY" tools/compare_menus.py tools/data/menu_baseline.tsv tools/data/menu_new3.tsv >/dev/null 2>&1 || { echo "compare_menus.py failed"; bad=1; }
  "$PY" tools/steps_sig.py >/dev/null 2>&1; [ $bad = 0 ]; }
# the live-agent tests (grok / codex / claude / agy / TAP) talk to real services and are flaky under load: run_all_checks skips them unless OGF_LIVE=1
# (plain `pytest ai_bridge/tests` still runs them when the CLIs are installed and logged in); with OGF_LIVE=1 a failed run is repeated once for the failed tests only
chk_ai_bridge_tests() { local r
  if [ "${OGF_LIVE:-0}" = 1 ]; then "$PY" -m pytest -q ai_bridge/tests 2>&1 | tail -3; r=${PIPESTATUS[0]}
    [ $r = 0 ] || { echo "-- repeating the failed tests once"; "$PY" -m pytest -q --lf ai_bridge/tests 2>&1 | tail -3; r=${PIPESTATUS[0]}; }
  else OGF_AI_OFFLINE=1 "$PY" -m pytest -q ai_bridge/tests 2>&1 | tail -3; r=${PIPESTATUS[0]}; echo "(live-agent tests skipped; OGF_LIVE=1 includes them)"; fi
  [ $r = 0 ]; }
chk_moving_tests() { "$PY" -m pytest -q moving/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
chk_ai_gui() { need_x || { echo "no X server"; return 77; }
  "$PY" scripts/verify_ai_gui.py --python "$PY" --display $DISP --workdir "$OUT/ai_gui_work" 2>&1 | tail -6; [ ${PIPESTATUS[0]} = 0 ]; }
chk_agent_gui() { need_x || { echo "no X server"; return 77; }
  "$PY" scripts/verify_agent_gui.py --python "$PY" --display $DISP --workdir "$OUT/agent_gui_work" 2>&1 | tail -6; [ ${PIPESTATUS[0]} = 0 ]; }
# real agent-CLI binaries (if installed): flags accepted + not-logged-in error parsed; NO login, NO model answer, nothing is sent
chk_agent_real() { "$PY" scripts/verify_agent_cli_real.py 2>&1 | tail -12; return ${PIPESTATUS[0]}; }
chk_icl_export() { need_x || { echo "no X server"; return 77; }
  OGF_VERIFY_OUT="$OUT/icl" T=120 gui $DS9 $FITS/m51.fits -source scripts/verify_icl_export.tcl > "$OUT/icl_raw.txt" 2>&1
  local rc=$?; grep -E "VERIFY|rror" "$OUT/icl_raw.txt" | tail -8
  [ $rc = 0 ] && grep -q "VERIFY" "$OUT/icl_raw.txt" && ! grep -E "VERIFY.*(FAIL|fail)" "$OUT/icl_raw.txt" >/dev/null; }
tclsum() { # file -> prints counts, returns 0 when "SUMMARY failures=0" and no FAIL
  grep -c '^PASS' "$1" | sed 's/^/PASS lines: /'; grep -E '^FAIL|^ERROR|^BGERROR' "$1" | head -5; grep -E '^SUMMARY' "$1"
  grep -q '^SUMMARY failures=0' "$1" && ! grep -q '^FAIL' "$1"; }
chk_click_chooser() { need_x || { echo "no X server"; return 77; }
  OGF_CHK_OUT="$OUT/chooser.txt" gui $DS9 $FITS/m51.fits -geometry 1300x950 -source scripts/verify_click_chooser.tcl > "$OUT/chooser_raw.txt" 2>&1; tclsum "$OUT/chooser.txt"; }
chk_click_xevent() { need_x || { echo "no X server"; return 77; }; command -v xdotool >/dev/null || { echo "xdotool not installed"; return 77; }
  OGF_CHK_OUT="$OUT/xevent.txt" gui $DS9 $FITS/m51.fits -geometry 1300x950 -source scripts/verify_click_xevent.tcl > "$OUT/xevent_raw.txt" 2>&1; tclsum "$OUT/xevent.txt"; }
chk_session_replay() { T=1800; "$HERE/verify_session_replay.sh" --workdir "$OUT/replay_work" 2>&1 | tail -400 > "$OUT/replay_tail.txt"
  grep -E "checks,|ALL CHECKS|FAIL" "$OUT/replay_tail.txt" | tail -5; grep -q "ALL CHECKS PASSED\|0 failed" "$OUT/replay_tail.txt" && ! grep -q "^FAIL" "$OUT/replay_tail.txt"; }
chk_link_bench() { [ -f /workspace/work/inj1.json.pkl ] || { echo "injection sets /workspace/work/inj*.json.pkl missing"; return 77; }
  # regenerated + relabelled sets (validation/inject.py, incl. the CR-heavy inj9-12) when present, else the original inj1-6
  local D=/workspace/work/inj; [ -f /workspace/work/i1/sets/inj12.json.pkl ] && D=/workspace/work/i1/sets/inj
  local L="${D}1.json.pkl ${D}2.json.pkl ${D}3.json.pkl ${D}4.json.pkl ${D}5.json.pkl ${D}6.json.pkl"
  [ "$D" != /workspace/work/inj ] && L="$L ${D}7.json.pkl ${D}8.json.pkl ${D}9.json.pkl ${D}10.json.pkl ${D}11.json.pkl ${D}12.json.pkl"
  (cd moving/validation && "$PY" link_bench.py $L 2>&1 | tail -16); }
chk_moving_session() {
  # with OGF_MOVING_SESSION / _REF / _FIELD set: replay that recorded session against that reference (as before).  Otherwise the check is
  # self-contained: verify_moving_session_auto.sh records a GUI session on the cached BB89 exposures and replays it (77 = skipped when
  # the exposures, an X server or the network+cache are missing).
  if [ -n "$OGF_MOVING_SESSION" ]; then "$HERE/verify_moving_session.sh" "$OGF_MOVING_SESSION" "$OGF_MOVING_REF" "$OGF_MOVING_FIELD"
  else DISPLAY_OVERRIDE=$DISP "$HERE/verify_moving_session_auto.sh"; fi; }
chk_cat_api() { need_x || { echo "no X server"; return 77; }
  OGF_CAT_API_OUT="$OUT/cat_api.txt" gui $DS9 $FITS/m51.fits -geometry 1300x950 -source scripts/verify_cat_api.tcl > "$OUT/cat_api_raw.txt" 2>&1; tclsum "$OUT/cat_api.txt"; }
chk_cat_behavior() { need_x || { echo "no X server"; return 77; }; DISPLAY_OVERRIDE=$DISP "$HERE/verify_cat_behavior.sh" 2>&1 | tail -8; [ ${PIPESTATUS[0]} = 0 ]; }
chk_report_tests() { "$PY" -m pytest -q plugins/report/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
chk_review_gui() { need_x || { echo "no X server"; return 77; }
  OGF_REVIEW_OUT="$OUT/review.txt" OGF_REVIEW_DIR="$OUT/review_work" gui $DS9 $FITS/m51.fits -geometry 1300x950 -source scripts/verify_review_gui.tcl > "$OUT/review_raw.txt" 2>&1; tclsum "$OUT/review.txt"; }
chk_review_td() { need_x || { echo "no X server"; return 77; }
  OGF_RVTD_OUT="$OUT/review_td.txt" OGF_RVTD_DIR="$OUT/review_td_work" gui $DS9 $FITS/m51.fits -geometry 1300x950 -source scripts/verify_review_td.tcl > "$OUT/review_td_raw.txt" 2>&1; tclsum "$OUT/review_td.txt"; }
chk_geometry() { need_x || { echo "no X server"; return 77; }
  OGF_GEO_OUT="$OUT/geometry.txt" gui $DS9 $FITS/m51.fits -geometry 1300x950 -source scripts/verify_geometry.tcl > "$OUT/geometry_raw.txt" 2>&1; tclsum "$OUT/geometry.txt"; }
chk_mouse() { need_x || { echo "no X server"; return 77; }; command -v xdotool >/dev/null || { echo "xdotool missing"; return 77; }
  OGF_MOUSE_OUT="$OUT/mouse.txt" OGF_MOUSE_DIR="$OUT/mouse_work" gui $DS9 $FITS/m51.fits -geometry 1300x950 -source scripts/verify_mouse.tcl > "$OUT/mouse_raw.txt" 2>&1; tclsum "$OUT/mouse.txt"; }
# ---- extra checks registered by later work (appended below by the feature that adds them)
chk_moving_options() { need_x || { echo "no X server"; return 77; }
  OGF_MOVOPT_OUT="$OUT/movopt.txt" gui $DS9 $FITS/m51.fits -geometry 1300x950 -source scripts/verify_moving_options.tcl > "$OUT/movopt_raw.txt" 2>&1; tclsum "$OUT/movopt.txt"; }
EXTRA="$EXTRA moving_options"
chk_manifests() { "$PY" tools/validate_manifests.py 2>&1 | tail -5; [ ${PIPESTATUS[0]} = 0 ] || return 1
  "$PY" -m pytest -q tools/tests 2>&1 | tail -2; [ ${PIPESTATUS[0]} = 0 ]; }
chk_cli_templates() { need_x || { echo "no X server"; return 77; }
  local H; H=$(mktemp -d /tmp/ogf_clitpl_home.XXXXXX); rm -f "$OUT/clitpl_fake.log"
  HOME=$H OGF_CLITPL_OUT="$OUT/clitpl.txt" FAKE_LOG="$OUT/clitpl_fake.log" OGFINDER_PYTHON="$ROOT/scripts/fake_python_for_templates.sh" \
    gui $DS9 $FITS/m51.fits -geometry 1300x950 -source scripts/verify_cli_templates.tcl > "$OUT/clitpl_raw.txt" 2>&1
  tclsum "$OUT/clitpl.txt"; local rc=$?; rm -rf "$H"; return $rc; }
chk_cli_headless() { need_x || { echo "no X server"; return 77; }
  local H; H=$(mktemp -d /tmp/ogf_clihl_home.XXXXXX); rm -f "$OUT/clihl_fake.log"
  HOME=$H OGF_CLIHL_OUT="$OUT/clihl.txt" FAKE_LOG="$OUT/clihl_fake.log" OGFINDER_PYTHON="$ROOT/scripts/fake_python_for_templates.sh" \
    gui $DS9 $FITS/m51.fits -geometry 1300x950 -source scripts/verify_cli_headless.tcl > "$OUT/clihl_raw.txt" 2>&1
  tclsum "$OUT/clihl.txt"; local rc=$?; rm -rf "$H"; return $rc; }
EXTRA="$EXTRA manifests cli_templates cli_headless"
chk_newplugins() { need_x || { echo "no X server"; return 77; }
  DISPLAY_OVERRIDE=$DISP "$HERE/verify_newplugins.sh" "$OUT/newplugins_work" 2>&1 | tail -8; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA newplugins"
chk_isophote_tests() { "$PY" -m pytest -q plugins/isophote/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA isophote_tests"
chk_completeness_tests() { "$PY" -m pytest -q plugins/completeness/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA completeness_tests"
chk_daophot_tests() { "$PY" -m pytest -q plugins/daophot/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA daophot_tests"
chk_psfex_tests() { "$PY" -m pytest -q plugins/psfex/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA psfex_tests"
chk_multifit_tests() { "$PY" -m pytest -q plugins/multifit/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA multifit_tests"
chk_morphext_tests() { "$PY" -m pytest -q plugins/morph_ext/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA morphext_tests"
chk_noisemodel_tests() { "$PY" -m pytest -q plugins/noisemodel/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA noisemodel_tests"
chk_stacking_tests() { "$PY" -m pytest -q plugins/stacking/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA stacking_tests"
chk_photoz_tests() { "$PY" -m pytest -q plugins/photoz_sed/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA photoz_tests"
chk_sedcodes_tests() { "$PY" -m pytest -q plugins/sedcodes/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA sedcodes_tests"
chk_cluster_tests() { "$PY" -m pytest -q plugins/cluster/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA cluster_tests"
chk_spectra_tests() { "$PY" -m pytest -q plugins/spectra/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA spectra_tests"
chk_spectra_kin_gui() { need_x || { echo "no X server"; return 77; }
  DISPLAY_OVERRIDE=$DISP "$HERE/verify_spectra_kin.sh" "$OUT/spectra_kin_work" 2>&1 | tail -6; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA spectra_kin_gui"
chk_xmatch_tests() { "$PY" -m pytest -q plugins/xmatch/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA xmatch_tests"
chk_lightcurves_tests() { "$PY" -m pytest -q plugins/lightcurves/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA lightcurves_tests"
chk_batch_tests() { "$PY" -m pytest -q plugins/batch/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA batch_tests"
chk_repro_tests() { "$PY" -m pytest -q plugins/repro/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA repro_tests"
chk_lensmodel_tests() { "$PY" -m pytest -q plugins/lensmodel/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA lensmodel_tests"
chk_lensmodel_gui() { need_x || { echo "no X server"; return 77; }
  DISPLAY_OVERRIDE=$DISP "$HERE/verify_lensmodel.sh" "$OUT/lensmodel_work" 2>&1 | tail -6; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA lensmodel_gui"
chk_trails_tests() { "$PY" -m pytest -q plugins/trails/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA trails_tests"
chk_trails_gui() { need_x || { echo "no X server"; return 77; }
  DISPLAY_OVERRIDE=$DISP "$HERE/verify_trails.sh" "$OUT/trails_work" 2>&1 | tail -6; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA trails_gui"
# cross-night orbit-fit linking (moving/nightlink.py): synthetic injection + decoys, CLI mode, real MPC astrometry when the network is up (also part of moving_tests)
chk_nightlink_tests() { "$PY" -m pytest -q moving/tests/test_nightlink.py 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA nightlink_tests"
chk_regression_tests() { "$PY" -m pytest -q regression 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
EXTRA="$EXTRA regression_tests"
# per-tool regression set on small public data (HUDF F160W crop, SDSS Stripe 82 run 94, HST J0946+1006, SDSS photo-z set): data are fetched on demand
# into $OGF_DATA_CACHE (docs/testing.md); a case whose data are unreachable is SKIP, and the whole check is skipped (77) when nothing could run
chk_regression_data() { local o rc; o=$("$PY" regression/run_regression.py --report "$OUT/regression_report.json" 2>&1); rc=$?; echo "$o" | tail -14
  [ $rc != 0 ] && return 1; echo "$o" | grep -q 'regression: 0 pass, 0 fail' && return 77; return 0; }
LONG="$LONG regression_data"
#@EXTRA-CHECKS

run_one() { local n=$1 t0 t1 rc
  t0=$(date +%s); "chk_$n" > "$OUT/$n.log" 2>&1; rc=$?; t1=$(date +%s)
  NAMES+=("$n"); SECS+=($((t1-t0)))
  case $rc in 0) RES+=(PASS);; 77) RES+=(SKIP);; *) RES+=(FAIL);; esac
  NOTE+=("$(tail -1 "$OUT/$n.log" | cut -c1-70)"); }
LIST_ALL="$ALL $EXTRA"; [ $QUICK = 0 ] && LIST_ALL="$LIST_ALL $LONG"
[ $LIST = 1 ] && { echo "$ALL $EXTRA (long: $LONG)"; exit 0; }
for n in $LIST_ALL; do
  [ -n "$ONLY" ] && ! [[ ",$ONLY," == *",$n,"* ]] && continue
  printf '... %-18s' "$n" >&2; run_one "$n"; echo "${RES[-1]} (${SECS[-1]}s)" >&2
  [ $FAILFAST = 1 ] && [ "${RES[-1]}" = FAIL ] && break
done
if [ $QUICK = 1 ]; then for n in $LONG; do [ -z "$ONLY" ] && { NAMES+=("$n"); RES+=("SKIP"); SECS+=(0); NOTE+=("--quick"); }; done; fi
echo; printf '%-20s %-6s %6s  %s\n' CHECK RESULT SEC DETAIL; printf '%-20s %-6s %6s  %s\n' -------------------- ------ ------ ------
nf=0; for i in "${!NAMES[@]}"; do printf '%-20s %-6s %6s  %s\n' "${NAMES[$i]}" "${RES[$i]}" "${SECS[$i]}" "${NOTE[$i]}"; [ "${RES[$i]}" = FAIL ] && nf=$((nf+1)); done
echo; echo "logs: $OUT   failures: $nf"
[ $nf = 0 ]
