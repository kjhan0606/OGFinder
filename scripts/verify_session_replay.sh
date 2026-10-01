#!/bin/bash
# Session-script replay verification (see scripts/verify_session_replay.py).
# Usage: scripts/verify_session_replay.sh [options]   (options are passed on, e.g. --skip-hudf)
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="${OGFINDER_PYTHON:-/workspace/ogf_venv/bin/python3}"
exec "$PY" "$HERE/verify_session_replay.py" --python "$PY" "$@"
