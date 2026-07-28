#!/usr/bin/env python3
"""DS9 bridge for preliminary optical orbit fitting with CODES.

The input CSV is intentionally explicit.  It must contain at least three
rows with ``obstime_utc`` or ``jd_tdb``, ``ra_deg`` and ``dec_deg``.  The
bridge writes the CODES fit products and a DS9 FK5 region file for review.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def _region_file(observations: Path, output: Path) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    with observations.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    with output.open("w", encoding="utf-8") as stream:
        stream.write("# Region file format: DS9 version 4.1\n")
        stream.write('global color=cyan width=2 font="helvetica 10 normal"\n')
        stream.write("fk5\n")
        for index, row in enumerate(rows, start=1):
            ra = row.get("ra_deg") or row.get("RA") or row.get("ra")
            dec = row.get("dec_deg") or row.get("DEC") or row.get("dec")
            if ra is None or dec is None:
                continue
            stream.write(f"point({float(ra)},{float(dec)}) # point=cross text={{obs {index}}}\n")
    return output


def run_codes(args: argparse.Namespace) -> dict:
    codes_root = Path(args.codes_root).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    command = [
        sys.executable,
        "-m",
        "neo_orbit_calculator.cli",
        "fit-observations",
        str(Path(args.observations).expanduser().resolve()),
        "--kernel-dir",
        str(Path(args.kernel_dir).expanduser().resolve()),
        "--output-dir",
        str(output_dir),
        "--ephemeris",
        args.ephemeris,
        "--backend",
        args.backend,
        "--openorb",
        "--observatory-code",
        args.observatory_code,
    ]
    if args.stop:
        command.extend(["--stop", args.stop])
    if args.major_bodies_only:
        command.append("--major-bodies-only")
    if args.jupiter_system:
        command.append("--jupiter-system")
    environment = dict(__import__("os").environ)
    environment["PYTHONPATH"] = str(codes_root) + ":" + environment.get("PYTHONPATH", "")
    completed = subprocess.run(
        command,
        cwd=codes_root,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
    )
    output_text = completed.stdout.strip()
    if not output_text:
        raise RuntimeError("CODES returned no JSON summary.")
    summary = json.loads(output_text)
    summary["regions"] = str(
        _region_file(
            Path(args.observations).expanduser().resolve(),
            output_dir / "observations.reg",
        )
    )
    (output_dir / "ds9_codes_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("observations", type=Path)
    parser.add_argument("--codes-root", type=Path, default=Path("/home/kjhan/BACKUP/CODES"))
    parser.add_argument("--kernel-dir", type=Path, default=Path("/home/kjhan/BACKUP/CODES/kernels"))
    parser.add_argument("--output-dir", type=Path, default=Path("~/.ds9/neo_codes"))
    parser.add_argument("--stop")
    parser.add_argument("--ephemeris", choices=("auto", "de442", "de441"), default="auto")
    parser.add_argument("--backend", choices=("fortran", "kepler-split", "scipy"), default="fortran")
    parser.add_argument("--major-bodies-only", action="store_true")
    parser.add_argument("--jupiter-system", action="store_true")
    parser.add_argument("--observatory-code", default="W84")
    args = parser.parse_args()
    print(json.dumps(run_codes(args), indent=2))


if __name__ == "__main__":
    main()
