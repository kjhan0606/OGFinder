#!/usr/bin/env python3
"""Audit the locally approved S23B CAMIRA cluster catalog's available fields.

This is a local-only catalog integrity and scope audit. It makes no network
requests and reads no HSC image pixels, shapes, or photo-z rows.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_BYTES = 827_339
EXPECTED_ROWS = 12_162
EXPECTED_SHA256 = "5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71"
EXPECTED_FIELDS = 7
EXPECTED_COLUMNS = (
    "ID", "Name", "RA [deg]", "Dec [deg]",
    "cluster photometric redshift z_cl", "richness N_mem",
    "spectroscopic redshift of BCG (-1.0 if not available)",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_catalog(path: Path) -> dict:
    digest = sha256(path)
    size = path.stat().st_size
    columns: list[str] = []
    nrows = 0
    malformed: list[dict[str, object]] = []
    with path.open("rt", encoding="ascii") as stream:
        for line_number, line in enumerate(stream, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                prefix = "# column "
                if stripped.startswith(prefix):
                    label = stripped[len(prefix):].split(":", 1)
                    if len(label) == 2:
                        columns.append(label[1].strip())
                continue
            fields = stripped.split()
            if len(fields) != EXPECTED_FIELDS:
                if len(malformed) < 20:
                    malformed.append({"line": line_number, "n_fields": len(fields)})
                continue
            nrows += 1

    integrity_ok = (
        size == EXPECTED_BYTES and digest == EXPECTED_SHA256
        and nrows == EXPECTED_ROWS and tuple(columns) == EXPECTED_COLUMNS
        and not malformed
    )
    return {
        "schema": "astrafek-s23b-camira-catalog-content-audit-v1",
        "status": "complete" if integrity_ok else "failed",
        "input": str(path.resolve()),
        "bytes": size,
        "sha256": digest,
        "data_rows": nrows,
        "column_count": len(columns),
        "columns_from_comments": columns,
        "expected_columns": list(EXPECTED_COLUMNS),
        "malformed_rows_sample": malformed,
        "contains_cluster_member_positions": False,
        "contains_object_level_shape_rows": False,
        "contains_individual_photoz_rows": False,
        "network_requests": False,
        "image_pixels_read": False,
        "science_measurement": False,
        "integrity_ok": integrity_ok,
        "audited_utc": datetime.now(timezone.utc).isoformat(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(args.input)
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    result = audit_catalog(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    if not result["integrity_ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
