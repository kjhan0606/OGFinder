#!/usr/bin/env python3
"""Fetch only the official HSC Y3 four-bin n(z) SACC calibration file.

The URL is pinned to the HSC release page's SACC link. Credentials come from
HSC_USER/HSC_PASSWORD or hidden interactive prompts and are never written.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import getpass
import hashlib
import json
import os
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import HTTPBasicAuthHandler, HTTPPasswordMgrWithDefaultRealm, Request, build_opener


OFFICIAL_URL = (
    "https://hsc-release.mtk.nao.ac.jp/archive/filetree/shape_catalog_y3/"
    "dalal23/hsc_y3_fourier_space_data_vector.sacc"
)
OFFICIAL_HOST = "hsc-release.mtk.nao.ac.jp"
MAX_BYTES = 64 * 1024 * 1024
FILENAME = "hsc_y3_fourier_space_data_vector.sacc"


def validate_url(url: str) -> None:
    parsed = urlparse(url)
    if (parsed.scheme != "https" or parsed.netloc != OFFICIAL_HOST
            or parsed.path != urlparse(OFFICIAL_URL).path or parsed.query or parsed.fragment):
        raise ValueError("refusing URL outside the pinned official HSC SACC resource")


def validate_payload(payload: bytes) -> None:
    if not payload:
        raise ValueError("downloaded SACC file is empty")
    if len(payload) > MAX_BYTES:
        raise ValueError("download exceeds the 64 MiB safety limit")
    if not payload.startswith(b"SIMPLE  ="):
        raise ValueError("download is not a FITS/SACC file; refusing HTML or error content")


def fetch(url: str, user: str, password: str) -> tuple[bytes, str]:
    validate_url(url)
    manager = HTTPPasswordMgrWithDefaultRealm()
    manager.add_password(None, f"https://{OFFICIAL_HOST}/", user, password)
    opener = build_opener(HTTPBasicAuthHandler(manager))
    request = Request(url, headers={"User-Agent": "Astrafek-HSC-calibration-fetch/1.0"})
    try:
        with opener.open(request, timeout=120) as response:
            final_url = response.geturl()
            validate_url(final_url)
            chunks, size = [], 0
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                size += len(chunk)
                if size > MAX_BYTES:
                    raise ValueError("download exceeds the 64 MiB safety limit")
                chunks.append(chunk)
            payload = b"".join(chunks)
    except HTTPError as error:
        raise RuntimeError(f"official HSC download returned HTTP {error.code}") from None
    except URLError as error:
        raise RuntimeError(f"official HSC download failed: {error.reason}") from None
    validate_payload(payload)
    return payload, final_url


def download(output_dir: Path, user: str, password: str, url: str = OFFICIAL_URL) -> dict:
    validate_url(url)
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing output directory: {output_dir}")
    if not user or not password:
        raise ValueError("HSC account and password are required")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    payload, final_url = fetch(url, user, password)
    digest = hashlib.sha256(payload).hexdigest()
    output_dir.mkdir(mode=0o700, parents=False, exist_ok=False)
    file_path = output_dir / FILENAME
    with file_path.open("xb") as stream:
        stream.write(payload)
    manifest = {
        "schema": "astrafek-hsc-nz-sacc-download-v1",
        "status": "downloaded",
        "requested_url": url,
        "final_url": final_url,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "file": FILENAME,
        "bytes": len(payload),
        "sha256": digest,
        "content_type": "FITS/SACC",
        "credential_values_saved": False,
        "science_measurement": False,
    }
    (output_dir / "download_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--user", default=os.environ.get("HSC_USER"))
    parser.add_argument("--url", default=OFFICIAL_URL)
    args = parser.parse_args()
    user = (args.user or input("HSC account name: ")).strip()
    password = os.environ.get("HSC_PASSWORD") or getpass.getpass("HSC password (hidden): ")
    try:
        manifest = download(args.output, user, password, args.url)
    except (OSError, ValueError, RuntimeError) as error:
        parser.error(str(error))
    print(json.dumps({"status": manifest["status"], "file": manifest["file"],
                      "bytes": manifest["bytes"], "sha256": manifest["sha256"],
                      "output": str(args.output.resolve())}, ensure_ascii=False))


if __name__ == "__main__":
    main()
