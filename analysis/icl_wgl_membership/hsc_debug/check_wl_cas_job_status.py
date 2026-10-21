#!/usr/bin/env python3
"""Read the status of an existing HSC CAS job; never submits or downloads."""
from __future__ import annotations

import argparse
import getpass
import json
import re
import sys
import urllib.error
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import Cas


def query_status(job_id: str, user: str = "kjhan0606") -> dict:
    if not re.fullmatch(r"\d+", str(job_id)):
        raise ValueError("CAS job ID must be a positive integer")
    password = getpass.getpass("HSC password (hidden): ")
    if not password:
        raise ValueError("Empty HSC password")
    client = Cas(user, password)
    try:
        result = client.json("status", id=str(job_id))
        error = result.get("error")
        return {"job_id": str(job_id), "status": str(result.get("status", "unknown")),
                "error": str(error)[:1000] if error else None,
                "query_submitted": False, "rows_downloaded": False,
                "credentials_saved": False}
    except urllib.error.HTTPError as error:
        # Do not print the response body or any credential-bearing details.
        return {"job_id": str(job_id), "status": "status_query_http_error",
                "http_status": int(error.code), "query_submitted": False,
                "rows_downloaded": False, "credentials_saved": False}
    finally:
        client.credential.clear()
        password = None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("job_id", help="existing HSC CAS catalog-job ID")
    parser.add_argument("--user", default="kjhan0606")
    args = parser.parse_args()
    print(json.dumps(query_status(args.job_id, args.user), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
