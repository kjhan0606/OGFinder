#!/bin/bash
# Self-contained moving-object session check (run_all_checks: moving_session, long).
#   1. a real ds9 (Xvfb) runs the GUI steps on the cached BB89 exposures (scripts/verify_moving_record.tcl) in a scratch HOME and
#      saves the session script; its outputs under $HOME/.ds9/moving_work are the REFERENCE;
#   2. scripts/verify_moving_session.sh replays that script (--mode replay, empty HOME) and compares the deterministic outputs byte for byte.
# The reference is made by the CURRENT code on every run, so a change of the pipeline does not need a stored golden file; what is
# tested is "GUI run == replayed script run".  Takes ~10 minutes (two full pipeline runs) and needs network (SkyBoT, Horizons, MAST/Gaia).
# Exit 77 (skipped) when the 4 BB89 exposures are not in ~/.ds9/mast_cache, there is no X server / Xvfb, or the network is unreachable and
# there is no ~/.ds9/moving_cache.  The replay env also gets OGF_MOVING_CACHE/HOME of the scratch dir (the cache is a copy).
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(dirname "$HERE")"
DISP="${DISPLAY_OVERRIDE:-:77}"; REALHOME="${REAL_HOME:-$HOME}"
CACHE="$REALHOME/.ds9/mast_cache"; EPH="$REALHOME/.ds9/ephem"
FILES=""; for b in j8pu38c7q j8pu38caq j8pu38ceq j8pu38ciq; do f="$CACHE/${b}_flc.fits"; [ -f "$f" ] || { echo "SKIP: $f not cached (BB89 exposures, HST program 13758? see docs/moving_objects.md)"; exit 77; }; FILES="$FILES $f"; done
DISPLAY=$DISP xdpyinfo >/dev/null 2>&1 || { command -v Xvfb >/dev/null || { echo "SKIP: no X server"; exit 77; }
  # our Xvfb: started with fds >2 closed (no inherited lock fd) and killed on every exit (trap)
  ( for f in /proc/$BASHPID/fd/*; do f=${f##*/}; [ "$f" -gt 2 ] 2>/dev/null && eval "exec $f>&-"; done
    exec Xvfb $DISP -screen 0 1400x1000x24 </dev/null >/dev/null 2>&1 ) & XVFB_PID=$!
  trap 'kill $XVFB_PID 2>/dev/null; wait $XVFB_PID 2>/dev/null' EXIT; trap 'exit 130' INT; trap 'exit 143' TERM; trap 'exit 129' HUP
  sleep 2; DISPLAY=$DISP xdpyinfo >/dev/null 2>&1 || { echo "SKIP: no X server"; exit 77; }; }
NET=1; curl -s -m 15 -o /dev/null https://ssp.imcce.fr/ || NET=0   # SkyBoT/Horizons answers are cached in ~/.ds9/moving_cache (linked below)
[ $NET = 0 ] && [ ! -d "$REALHOME/.ds9/moving_cache" ] && { echo "SKIP: network unreachable and no moving_cache"; exit 77; }
PYBIN="${OGFINDER_PYTHON:-/workspace/ogf_venv/bin/python3}"; export OGFINDER_PYTHON="$PYBIN"
W=$(mktemp -d /tmp/ogf_movsess.XXXXXX); H="$W/home"; mkdir -p "$H/.ds9" "$W/rec"
ln -s "$CACHE" "$H/.ds9/mast_cache"; [ -d "$EPH" ] && ln -s "$EPH" "$H/.ds9/ephem"
[ -d "$REALHOME/.ds9/moving_cache" ] && cp -r "$REALHOME/.ds9/moving_cache" "$H/.ds9/moving_cache"   # a copy: the run may add entries; the real cache stays as it was
echo "network: $([ $NET = 1 ] && echo reachable || echo UNREACHABLE - relying on the cached SkyBoT/Horizons answers)"
cd "$ROOT"; rm -rf "$H/ds9.auto" "$H/ds9.auto.dir"
echo "recording the GUI session (work dir $W) ..."
DISPLAY=$DISP HOME="$H" MOV_OUT="$W/rec" timeout -s KILL 1500 bin/ds9 -geometry 1300x950 -source scripts/verify_moving_record.tcl > "$W/gui_stdout.txt" 2>&1
[ -f "$W/rec/ogfinder_session.py" ] || { echo "FAIL: the GUI run did not save a session script (see $W/gui_stdout.txt, $W/rec/gui.log)"; tail -5 "$W/rec/gui.log" 2>/dev/null; exit 1; }
echo "GUI run: $(tail -1 "$W/rec/gui.log")"
FIELD="BB89:main=$(echo $FILES | cut -d' ' -f1),img2=$(echo $FILES | cut -d' ' -f2),img3=$(echo $FILES | cut -d' ' -f3),img4=$(echo $FILES | cut -d' ' -f4)"
"$HERE/verify_moving_session.sh" "$W/rec/ogfinder_session.py" "$H/.ds9/moving_work" "$FIELD" "$W/replay_home"
rc=$?; echo "work dir: $W"; exit $rc
