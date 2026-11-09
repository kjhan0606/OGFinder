#!/usr/bin/env python3
"""Run an approved Gaia DR3 positional cross-match in bounded TAP batches."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import time
import urllib.parse
import urllib.request

import match_psf_gaia


GAIA_FIELDS = [
    "source_id", "ra", "dec", "ref_epoch", "parallax_over_error", "pmra", "pmra_error",
    "pmdec", "pmdec_error", "ruwe", "duplicated_source", "phot_g_mean_mag",
    "phot_variable_flag", "in_qso_candidates", "in_galaxy_candidates",
    "classprob_dsc_combmod_star",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_query(rows: list[dict]) -> str:
    """Query the union of approved 2-arcsec cones; associate IDs locally."""
    cones = [
        "1=CONTAINS(POINT('ICRS',g.ra,g.dec),"
        f"CIRCLE('ICRS',{float(row['ra_deg']):.12f},{float(row['dec_deg']):.12f},2.0/3600.0))"
        for row in rows
    ]
    return ("SELECT " + ", ".join(f"g.{name}" for name in GAIA_FIELDS)
            + " FROM gaiadr3.gaia_source AS g WHERE "
            + " OR ".join(f"({cone})" for cone in cones))


def angular_separation_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    """Stable great-circle separation for small Gaia/HSC match radii."""
    ra1, dec1, ra2, dec2 = map(math.radians, (ra1, dec1, ra2, dec2))
    dra = (ra2 - ra1 + math.pi) % (2 * math.pi) - math.pi
    ddec = dec2 - dec1
    hav = math.sin(ddec / 2) ** 2 + math.cos(dec1) * math.cos(dec2) * math.sin(dra / 2) ** 2
    return 2 * math.asin(math.sqrt(min(1.0, max(0.0, hav)))) * 206264.80624709636


def associate(targets: list[dict], sources: list[dict], offset: int) -> list[dict]:
    matches = []
    for local_index, target in enumerate(targets):
        for source in sources:
            separation = angular_separation_arcsec(
                float(target["ra_deg"]), float(target["dec_deg"]),
                float(source["ra"]), float(source["dec"]),
            )
            if separation <= 2.00001:
                matches.append({
                    "candidate_id": str(offset + local_index),
                    **source,
                    "separation_arcsec": f"{separation:.9f}",
                })
    return matches


def run(requests_path: Path, plan_path: Path, approved_plan_sha256: str,
        work: Path, batch_size: int, resume: bool = False) -> dict:
    if not platform.node().lower().startswith("lageunha"):
        raise RuntimeError("Run on lageunha CPU")
    if work.exists() and not resume:
        raise FileExistsError(f"refusing to overwrite Gaia output: {work}")
    if resume and not work.is_dir():
        raise FileNotFoundError(f"Gaia resume directory does not exist: {work}")
    if sha256(plan_path) != approved_plan_sha256:
        raise ValueError("plan SHA-256 does not match the user's approved plan")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    targets = json.loads(requests_path.read_text(encoding="utf-8"))["candidates"]
    expected_columns = ["candidate_id", *GAIA_FIELDS, "separation_arcsec"]
    if plan.get("schema") == "astrafek-s23b-matched-psf-expansion-gaia-plan-v1":
        planned_candidate_count = plan.get("n_candidate_positions")
        if plan.get("cas_query_requested") is not False:
            raise ValueError("expanded Gaia plan must explicitly exclude CAS requests")
    else:
        planned_candidate_count = plan.get("n_new_candidate_positions")
    if (plan.get("status") != "prepared_not_submitted"
            or planned_candidate_count != len(targets)
            or plan.get("candidate_list_sha256") != sha256(requests_path)
            or plan.get("match_radius_arcsec") != 2.0
            or plan.get("columns_returned") != expected_columns
            or plan.get("images_requested") is not False
            or plan.get("hsc_shape_rows_requested") is not False
            or plan.get("photoz_requested") is not False):
        raise ValueError("request rows or scope differ from the approved plan")
    planned_query_batches = plan.get("query_batches")
    if planned_query_batches is None:
        # Backward compatibility for the already approved target plan format.
        legacy_query_sha = hashlib.sha256((match_psf_gaia.QUERY + "\n").encode()).hexdigest()
        if plan.get("query_sha256") != legacy_query_sha:
            raise ValueError("legacy approved query template hash differs")
    else:
        expected_parameters = {
            "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "MAXREC": "20000",
        }
        if (batch_size != plan.get("batch_size_recommended")
                or plan.get("endpoint") != match_psf_gaia.ENDPOINT
                or plan.get("request_parameters") != expected_parameters
                or plan.get("maxrec_per_batch") != 20000
                or plan.get("gaia_columns_queried") != GAIA_FIELDS
                or len(planned_query_batches) != math.ceil(len(targets) / batch_size)
                or plan.get("expected_batches") != len(planned_query_batches)):
            raise ValueError("endpoint, request form, batch size, or Gaia column scope differs from the approved plan")
        manifest_bytes = json.dumps(
            planned_query_batches, sort_keys=True, separators=(",", ":")
        ).encode()
        if hashlib.sha256(manifest_bytes).hexdigest() != plan.get("query_batch_manifest_sha256"):
            raise ValueError("approved batch-query manifest hash differs")
        for batch_number, offset in enumerate(range(0, len(targets), batch_size)):
            batch_targets = targets[offset:offset + batch_size]
            entry = planned_query_batches[batch_number]
            expected_query = build_query(batch_targets) + "\n"
            query_path = Path(entry.get("query_file", ""))
            if (entry.get("batch") != batch_number
                    or entry.get("first_candidate_id") != offset
                    or entry.get("n_candidates") != len(batch_targets)
                    or hashlib.sha256(expected_query.encode()).hexdigest() != entry.get("query_sha256")
                    or not query_path.is_file()
                    or query_path.read_text(encoding="utf-8") != expected_query
                    or sha256(query_path) != entry.get("query_sha256")):
                raise ValueError(f"generated Gaia batch {batch_number} differs from the approved SQL")
    completed_batches = 0
    all_matches = []
    batch_count = (len(targets) + batch_size - 1) // batch_size
    if resume:
        if planned_query_batches is None:
            raise ValueError("resume is supported only for pinned batch-query plans")
        manifest_path = work / "status.json"
        candidates_snapshot = work / "candidates.json"
        if not manifest_path.is_file() or not candidates_snapshot.is_file():
            raise FileNotFoundError("resume output lacks its status or candidate snapshot")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if json.loads(candidates_snapshot.read_text(encoding="utf-8")) != targets:
            raise ValueError("resume candidate snapshot differs from the approved request list")
        if (manifest.get("status") != "running"
                or manifest.get("plan_sha256") != sha256(plan_path)
                or manifest.get("request_sha256") != sha256(requests_path)
                or manifest.get("approved_plan_query_sha256") != plan.get(
                    "query_batch_manifest_sha256", plan.get("query_sha256"))
                or manifest.get("n_candidates") != len(targets)
                or manifest.get("batch_size") != batch_size
                or manifest.get("match_radius_arcsec") != 2.0
                or manifest.get("columns") != expected_columns
                or manifest.get("images_requested") is not False
                or manifest.get("hsc_shape_rows_requested") is not False
                or manifest.get("photoz_requested") is not False):
            raise ValueError("resume status does not match the approved Gaia plan")
        completed_batches = int(manifest.get("completed_batches", -1))
        if not 0 <= completed_batches <= batch_count:
            raise ValueError("resume checkpoint has an invalid completed-batch count")
        if len(manifest.get("batches", [])) != completed_batches:
            raise ValueError("resume checkpoint batch list is incomplete")
        for batch_number in range(completed_batches):
            offset = batch_number * batch_size
            batch_targets = targets[offset:offset + batch_size]
            batch_dir = work / f"batch_{batch_number:03d}"
            batch = json.loads((batch_dir / "status.json").read_text(encoding="utf-8"))
            query_path = batch_dir / "query.adql"
            response_path = batch_dir / "gaia_cone_union.csv"
            matched_path = batch_dir / "matched_pairs.csv"
            planned_entry = planned_query_batches[batch_number]
            expected_query = build_query(batch_targets) + "\n"
            if (batch.get("status") != "complete"
                    or batch.get("batch") != batch_number
                    or batch.get("first_candidate_id") != offset
                    or batch.get("n_candidates") != len(batch_targets)
                    or batch.get("query_sha256") != planned_entry.get("query_sha256")
                    or not query_path.is_file()
                    or query_path.read_text(encoding="utf-8") != expected_query
                    or sha256(query_path) != planned_entry.get("query_sha256")
                    or not response_path.is_file()
                    or sha256(response_path) != batch.get("response_sha256")
                    or not matched_path.is_file()):
                raise ValueError(f"completed resume batch {batch_number} failed receipt checks")
            with response_path.open(encoding="utf-8-sig", newline="") as stream:
                reader = csv.DictReader(stream)
                if not set(GAIA_FIELDS).issubset(reader.fieldnames or []):
                    raise ValueError(f"unexpected Gaia CSV schema in completed batch {batch_number}")
                sources = list(reader)
            if (len(sources) != batch.get("n_gaia_union_rows")
                    or len(sources) >= int(plan.get("maxrec_per_batch", 20000))):
                raise ValueError(f"completed resume batch {batch_number} row count is invalid")
            source_ids = [row["source_id"] for row in sources]
            if len(source_ids) != len(set(source_ids)):
                raise ValueError(f"duplicate Gaia source in completed batch {batch_number}")
            expected_matches = associate(batch_targets, sources, offset)
            with matched_path.open(encoding="utf-8", newline="") as stream:
                stored_matches = list(csv.DictReader(stream))
            if stored_matches != expected_matches or len(stored_matches) != batch.get("n_match_rows"):
                raise ValueError(f"completed resume batch {batch_number} match table is inconsistent")
            all_matches.extend(expected_matches)
        if len(all_matches) != manifest.get("n_match_rows"):
            raise ValueError("resume checkpoint total matched-row count is inconsistent")
    else:
        work.mkdir(parents=True)
        (work / "candidates.json").write_text(
            json.dumps(targets, indent=2) + "\n", encoding="utf-8")
        manifest = dict(
            status="running", plan_sha256=sha256(plan_path), request_sha256=sha256(requests_path),
            approved_plan_query_sha256=plan.get(
                "query_batch_manifest_sha256", plan.get("query_sha256")
            ),
            execution_method="synchronous Gaia TAP cone-union queries; candidate IDs and exact separations associated locally",
            n_candidates=len(targets), batch_size=batch_size, match_radius_arcsec=2.0,
            columns=expected_columns, images_requested=False, hsc_shape_rows_requested=False,
            photoz_requested=False, science_measurement=False, psf_approved=False,
            completed_batches=0, n_match_rows=0, batches=[])
        write_json(work / "status.json", manifest)
    for batch_number, offset in enumerate(range(0, len(targets), batch_size)):
        if batch_number < completed_batches:
            continue
        batch_targets = targets[offset:offset + batch_size]
        batch_dir = work / f"batch_{batch_number:03d}"
        query = build_query(batch_targets)
        query_path = batch_dir / "query.adql"
        expected_query = query + "\n"
        if batch_dir.exists():
            if (not resume
                    or set(path.name for path in batch_dir.iterdir()) != {"query.adql"}
                    or query_path.read_text(encoding="utf-8") != expected_query):
                raise ValueError(f"incomplete resume batch {batch_number} has unexpected files")
        else:
            batch_dir.mkdir()
        if not query_path.exists():
            query_path.write_text(expected_query, encoding="utf-8")
        maxrec = int(plan.get("maxrec_per_batch", 20000))
        form = urllib.parse.urlencode({
            "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
            "MAXREC": str(maxrec), "QUERY": query,
        }).encode()
        req = urllib.request.Request(match_psf_gaia.ENDPOINT, form)
        started = time.monotonic()
        try:
            with urllib.request.urlopen(req, timeout=180) as response:
                raw = response.read(20 * 1024 * 1024 + 1)
        except Exception as error:
            batch = dict(batch=batch_number, first_candidate_id=offset,
                         n_candidates=len(batch_targets), query_sha256=sha256(query_path),
                         status="failed", error_type=type(error).__name__, error=str(error))
            write_json(batch_dir / "status.json", batch)
            manifest.update(status="failed", failed_batch=batch_number)
            write_json(work / "status.json", manifest)
            raise
        if len(raw) > 20 * 1024 * 1024:
            raise ValueError("Gaia response exceeds the 20 MiB safety limit")
        (batch_dir / "gaia_cone_union.csv").write_bytes(raw)
        reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
        if not set(GAIA_FIELDS).issubset(reader.fieldnames or []):
            raise ValueError(f"unexpected Gaia CSV schema in batch {batch_number}")
        sources = list(reader)
        if len(sources) >= maxrec:
            raise ValueError("possible MAXREC truncation in Gaia response")
        source_ids = [row["source_id"] for row in sources]
        if len(source_ids) != len(set(source_ids)):
            raise ValueError(f"duplicate Gaia sources in cone-union batch {batch_number}")
        matches = associate(batch_targets, sources, offset)
        with (batch_dir / "matched_pairs.csv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=expected_columns)
            writer.writeheader()
            writer.writerows(matches)
        batch = dict(
            batch=batch_number, first_candidate_id=offset, n_candidates=len(batch_targets),
            n_gaia_union_rows=len(sources), n_match_rows=len(matches),
            query=str(query_path.resolve()), query_sha256=sha256(query_path),
            response_sha256=hashlib.sha256(raw).hexdigest(),
            elapsed_seconds=round(time.monotonic() - started, 3), status="complete")
        write_json(batch_dir / "status.json", batch)
        all_matches.extend(matches)
        manifest["completed_batches"] = batch_number + 1
        manifest["n_match_rows"] = len(all_matches)
        manifest["batches"].append(batch)
        write_json(work / "status.json", manifest)
        print(f"Gaia batch {batch_number + 1}/{batch_count}: {len(batch_targets)} candidates, "
              f"{len(sources)} Gaia sources, {len(matches)} matched pairs", flush=True)
        time.sleep(1)
    output = work / "gaia_matches.csv"
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=expected_columns)
        writer.writeheader()
        writer.writerows(all_matches)
    manifest.update(status="complete", response_sha256=sha256(output),
                    n_match_rows=len(all_matches),
                    completed_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    write_json(work / "status.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--requests", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--approved-plan-sha256", required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=250)
    parser.add_argument("--resume", action="store_true",
                        help="resume a running, hash-matched plan from its completed batches")
    args = parser.parse_args()
    if not 1 <= args.batch_size <= 500:
        parser.error("batch size must be between 1 and 500")
    result = run(args.requests, args.plan, args.approved_plan_sha256,
                 args.work, args.batch_size, resume=args.resume)
    print(json.dumps({key: result[key] for key in (
        "status", "n_candidates", "batch_size", "completed_batches", "n_match_rows",
        "response_sha256", "science_measurement", "psf_approved")}, indent=2))


if __name__ == "__main__":
    main()
