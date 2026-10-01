#!/bin/bash
# Replay (and optionally compare) a recorded Moving Objects session script.
#   scripts/verify_moving_session.sh SESSION.py REF_OUTDIR  NAME:main=a.fits,img2=b.fits,...  [HOME_DIR]
# Runs `--mode replay` with an empty HOME (so ephemerides come from SESSION_META['ephem_dir']) and compares the
# deterministic outputs (detections/movers/tracklets/identified/orbit/transients/lightcurves/export) byte for byte
# with REF_OUTDIR (the outputs of the GUI run).  Needs network (SkyBoT/Horizons/MAST) and the cached exposures;
# exits 77 (skipped) if the exposures or the session script are absent.
set -u
S=${1:?session.py}; REF=${2:?reference outdir}; FIELD=${3:?FIELD=key=path,...}
H=${4:-/tmp/ogf_mov_home}
OUT=$(mktemp -d /tmp/ogf_mov_replay.XXXXXX)
[ -f "$S" ] || { echo "SKIP: no session script"; exit 77; }
for f in $(echo "${FIELD#*:}" | tr ',' '\n' | sed 's/^[^=]*=//'); do [ -f "$f" ] || { echo "SKIP: missing $f"; exit 77; }; done
rm -rf "$H"; mkdir -p "$H"
env -i HOME="$H" PATH="$PATH" python3 "$S" --mode replay --allow-network --python "${OGFINDER_PYTHON:-python3}" \
     --outdir "$OUT" --field "$FIELD" > "$OUT/run.log" 2>&1
rc=$?
echo "replay exit=$rc  $(grep -o 'steps ran=[0-9]* cached=[0-9]* skipped=[0-9]* failed=[0-9]*' "$OUT/run.log")  (log $OUT/run.log)"
bad=0; n=0
for f in detections.tsv movers.tsv tracklets.json identified.tsv orbit_tracklet0.json orbit_tracklet0.kv transients.tsv lightcurves.json elements_export.json tracklet0.obs80; do
  a=$(find "$REF" -name "$f" -not -path '*/pip*' | head -1); b=$(find "$OUT" -name "$f" | head -1)
  [ -n "$a" ] && [ -n "$b" ] || { echo "MISSING $f (ref=$a new=$b)"; bad=$((bad+1)); continue; }
  n=$((n+1)); cmp -s "$a" "$b" && echo "same  $f" || { echo "DIFF  $f"; bad=$((bad+1)); }
done
echo "compared $n files, $bad problems"; [ $rc -eq 0 ] && [ $bad -eq 0 ]
