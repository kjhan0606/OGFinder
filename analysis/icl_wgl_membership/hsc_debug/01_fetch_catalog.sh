#!/bin/bash
# Download the public CAMIRA S20A Wide catalog with the bright-star mask.
# The file is the cluster list for the ICL stack debug. It has no mass column.
set -euo pipefail
ROOT=${ASTRAFEK_DATA:-/gpfs/kjhan/Astrafek/data}
URL=https://raw.githubusercontent.com/oguri/cluster_catalogs/main/hsc_s20a_camira/camira_s20a_wide_sm.tbl
SHA=d929985d7f6ee38c4435fc12099e9df1e5dbfd9ce9603cea319d5704129d22aa
DEST=$ROOT/catalog/camira_s20a_wide_sm.tbl
mkdir -p "$ROOT/catalog" "$ROOT/cutouts" "$ROOT/work" "$ROOT/logs" "$ROOT/scripts"
if [[ -f $DEST ]] && echo "$SHA  $DEST" | sha256sum -c --status; then
  echo "catalog present $DEST"
  exit 0
fi
tmp=$(mktemp)
curl -fsSL -o "$tmp" "$URL"
echo "$SHA  $tmp" | sha256sum -c --status
mv "$tmp" "$DEST"
echo "catalog written $DEST"
