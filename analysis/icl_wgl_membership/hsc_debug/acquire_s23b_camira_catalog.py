#!/usr/bin/env python3
"""Acquire only the approved, commit-pinned S23B CAMIRA Wide catalog."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


COMMIT = "d80734b10bd823d13b9cd5292dadc77c9c05f683"
FILENAME = "camira_s23b_wide_sm.tbl"
URL = (
    "https://raw.githubusercontent.com/oguri/cluster_catalogs/"
    f"{COMMIT}/hsc_s23b_camira/{FILENAME}"
)
EXPECTED_BYTES = 827_339
EXPECTED_GIT_BLOB_SHA1 = "066c2c8beca8f3a65186db7511fddca52d298908"
EXPECTED_SHA256 = "5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71"
MAX_BYTES = EXPECTED_BYTES


def git_blob_sha1(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def acquire(output: Path) -> dict[str, object]:
    if not platform.node().lower().startswith("lageunha"):
        raise RuntimeError("Run the approved catalog download on lageunha")
    if output.name != FILENAME:
        raise ValueError(f"output filename must be exactly {FILENAME}")
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing file: {output}")
    manifest_path = output.with_suffix(output.suffix + ".manifest.json")
    if manifest_path.exists():
        raise FileExistsError(f"refusing to overwrite manifest: {manifest_path}")
    output.parent.mkdir(parents=True, exist_ok=True)

    request = urllib.request.Request(
        URL,
        headers={"User-Agent": "Astrafek-S23B-catalog-audit/1.0"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        if response.status != 200:
            raise RuntimeError(f"catalog request returned HTTP {response.status}")
        chunks: list[bytes] = []
        received = 0
        while True:
            chunk = response.read(64 * 1024)
            if not chunk:
                break
            received += len(chunk)
            if received > MAX_BYTES:
                raise RuntimeError("response exceeded the approved byte limit")
            chunks.append(chunk)
    data = b"".join(chunks)
    if len(data) != EXPECTED_BYTES:
        raise RuntimeError(f"expected {EXPECTED_BYTES} bytes, received {len(data)}")
    blob_sha1 = git_blob_sha1(data)
    if blob_sha1 != EXPECTED_GIT_BLOB_SHA1:
        raise RuntimeError(f"Git blob SHA-1 mismatch: {blob_sha1}")
    sha256 = hashlib.sha256(data).hexdigest()
    if sha256 != EXPECTED_SHA256:
        raise RuntimeError(f"SHA-256 mismatch: {sha256}")

    # Stage on the same filesystem and publish only after all integrity gates pass.
    fd, temp_name = tempfile.mkstemp(prefix=f".{FILENAME}.", dir=output.parent)
    try:
        with os.fdopen(fd, "wb") as staged:
            staged.write(data)
            staged.flush()
            os.fsync(staged.fileno())
        if output.exists():
            raise FileExistsError(f"refusing to overwrite existing file: {output}")
        os.link(temp_name, output)
    finally:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass

    manifest: dict[str, object] = {
        "purpose": "catalog duplicate/coverage audit only",
        "filename": FILENAME,
        "source_url": URL,
        "repository_commit": COMMIT,
        "bytes": len(data),
        "git_blob_sha1": blob_sha1,
        "sha256": sha256,
        "downloaded_utc": datetime.now(timezone.utc).isoformat(),
        "images_shape_rows_photoz_included": False,
        "science_measurement": False,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = acquire(args.output)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
