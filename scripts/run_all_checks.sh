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
ALL="tools_syntax ai_bridge_tests moving_tests ai_gui icl_export click_chooser click_xevent cat_api cat_behavior report_tests review_gui"
LONG="session_replay link_bench moving_session"
chk_tools_syntax() { local bad=0 f
  for f in tools/*.py scripts/*.py; do "$PY" -m py_compile "$f" 2>&1 || bad=1; done
  command -v tclsh >/dev/null && { tclsh tools/tclprocs.tcl plugins/mask/mask.tcl >/dev/null 2>&1 || { echo "tclprocs.tcl failed"; bad=1; }; }
  "$PY" tools/compare_menus.py tools/data/menu_baseline.tsv tools/data/menu_new3.tsv >/dev/null 2>&1 || { echo "compare_menus.py failed"; bad=1; }
  "$PY" tools/steps_sig.py >/dev/null 2>&1; [ $bad = 0 ]; }
chk_ai_bridge_tests() { "$PY" -m pytest -q ai_bridge/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
chk_moving_tests() { "$PY" -m pytest -q moving/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
chk_ai_gui() { need_x || { echo "no X server"; return 77; }
  "$PY" scripts/verify_ai_gui.py --python "$PY" --display $DISP --workdir "$OUT/ai_gui_work" 2>&1 | tail -6; [ ${PIPESTATUS[0]} = 0 ]; }
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
  (cd moving/validation && "$PY" link_bench.py /workspace/work/inj{1,2,3,4,5,6}.json.pkl 2>&1 | tail -12); }
chk_moving_session() { [ -n "$OGF_MOVING_SESSION" ] || { echo "set OGF_MOVING_SESSION=session.py REF=dir FIELD=... (needs network, MAST cache)"; return 77; }
  "$HERE/verify_moving_session.sh" "$OGF_MOVING_SESSION" "$OGF_MOVING_REF" "$OGF_MOVING_FIELD"; }
chk_cat_api() { need_x || { echo "no X server"; return 77; }
  OGF_CAT_API_OUT="$OUT/cat_api.txt" gui $DS9 $FITS/m51.fits -geometry 1300x950 -source scripts/verify_cat_api.tcl > "$OUT/cat_api_raw.txt" 2>&1; tclsum "$OUT/cat_api.txt"; }
chk_cat_behavior() { need_x || { echo "no X server"; return 77; }; DISPLAY_OVERRIDE=$DISP "$HERE/verify_cat_behavior.sh" 2>&1 | tail -8; [ ${PIPESTATUS[0]} = 0 ]; }
chk_report_tests() { "$PY" -m pytest -q plugins/report/tests 2>&1 | tail -3; [ ${PIPESTATUS[0]} = 0 ]; }
chk_review_gui() { need_x || { echo "no X server"; return 77; }
  OGF_REVIEW_OUT="$OUT/review.txt" OGF_REVIEW_DIR="$OUT/review_work" gui $DS9 $FITS/m51.fits -geometry 1300x950 -source scripts/verify_review_gui.tcl > "$OUT/review_raw.txt" 2>&1; tclsum "$OUT/review.txt"; }
# ---- extra checks registered by later work (appended below by the feature that adds them)
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
