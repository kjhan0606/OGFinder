#!/usr/bin/env python3
"""Parse and query official HSC PDR3 Wide tract/patch sky polygons.

The input files are the four coordinate tables linked from the HSC PDR3
Available Data page. This module never downloads files. It reads only the
explicit files in ``--directory`` and assigns coordinates by spherical
polygon containment in a local gnomonic projection. Overlapping polygons are
reported and resolved deterministically to the nearest containing center;
cross-field overlaps remain unassigned and are flagged.

Official sources:
https://hsc-release.mtk.nao.ac.jp/doc/index.php/available-data__pdr3/
https://hsc-release.mtk.nao.ac.jp/rsrc/pdr3/tract_patches/info/tracts_patches_W-AEGIS.txt
https://hsc-release.mtk.nao.ac.jp/rsrc/pdr3/tract_patches/info/tracts_patches_W-hectomap.txt
https://hsc-release.mtk.nao.ac.jp/rsrc/pdr3/tract_patches/info/tracts_patches_W-autumn.txt
https://hsc-release.mtk.nao.ac.jp/rsrc/pdr3/tract_patches/info/tracts_patches_W-spring.txt
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Iterable

import numpy as np
from scipy.spatial import cKDTree


WIDE_FILES = {
    "tracts_patches_W-AEGIS.txt": "AEGIS",
    "tracts_patches_W-hectomap.txt": "HECTOMAP",
    "tracts_patches_W-autumn.txt": "AUTUMN",
    "tracts_patches_W-spring.txt": "SPRING",
}
WIDE_FILE_BYTES = {
    "tracts_patches_W-AEGIS.txt": 156_887,
    "tracts_patches_W-hectomap.txt": 4_762_372,
    "tracts_patches_W-autumn.txt": 26_107_950,
    "tracts_patches_W-spring.txt": 18_198_617,
}
TRACT_MARKER = re.compile(r"^\*+\s*Tract\s*\*+$", re.IGNORECASE)
PATCH_MARKER = re.compile(r"^\*+\s*Patches\s*\*+$", re.IGNORECASE)
TRACT_ID = re.compile(r"\bTract:\s*(\d+)")
PATCH_ID = re.compile(r"\bPatch:\s*(\d+)\s*,\s*(\d+)")
SKY_POINT = re.compile(
    r"\b(Center|Corner(\d+))\s*\(RA,\s*Dec\):\s*"
    r"\(\s*([-+\d.eE]+)\s*,\s*([-+\d.eE]+)\s*\)",
)


@dataclass(frozen=True)
class SkyPolygon:
    field: str
    tract: int
    patch: tuple[int, int] | None
    center_ra_deg: float
    center_dec_deg: float
    corners_ra_dec: tuple[tuple[float, float], ...]


@dataclass
class PointAssignments:
    field: np.ndarray
    tract: np.ndarray
    patch: np.ndarray
    containing_polygon_count: np.ndarray
    cross_field_ambiguous: np.ndarray

    @property
    def unmatched_count(self) -> int:
        return int(np.count_nonzero(self.containing_polygon_count == 0))

    @property
    def overlap_count(self) -> int:
        return int(np.count_nonzero(self.containing_polygon_count > 1))

    @property
    def cross_field_ambiguous_count(self) -> int:
        return int(np.count_nonzero(self.cross_field_ambiguous))


def unit_vectors(ra_deg, dec_deg) -> np.ndarray:
    ra = np.deg2rad(np.asarray(ra_deg, dtype=float))
    dec = np.deg2rad(np.asarray(dec_deg, dtype=float))
    if (ra.ndim != 1 or dec.shape != ra.shape or not np.isfinite(ra).all()
            or not np.isfinite(dec).all() or np.any((dec < -np.pi / 2) | (dec > np.pi / 2))):
        raise ValueError("RA/Dec must be finite one-dimensional arrays with valid declinations")
    return np.column_stack((np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)))


def _point_xy(ra_deg: float, dec_deg: float, center_ra_deg: float,
              center_dec_deg: float) -> tuple[float, float] | None:
    ra, dec, ra0, dec0 = np.deg2rad([ra_deg, dec_deg, center_ra_deg, center_dec_deg])
    dra = math.atan2(math.sin(ra - ra0), math.cos(ra - ra0))
    sin_dec, cos_dec = math.sin(dec), math.cos(dec)
    sin_dec0, cos_dec0 = math.sin(dec0), math.cos(dec0)
    denominator = sin_dec0 * sin_dec + cos_dec0 * cos_dec * math.cos(dra)
    if denominator <= 0:
        return None
    return (
        cos_dec * math.sin(dra) / denominator,
        (cos_dec0 * sin_dec - sin_dec0 * cos_dec * math.cos(dra)) / denominator,
    )


def _point_in_xy_polygon(point: tuple[float, float], vertices: np.ndarray,
                         tolerance: float = 2e-9) -> bool:
    x, y = point
    inside = False
    previous = vertices[-1]
    for current in vertices:
        x0, y0 = float(previous[0]), float(previous[1])
        x1, y1 = float(current[0]), float(current[1])
        dx, dy = x1 - x0, y1 - y0
        cross = (x - x0) * dy - (y - y0) * dx
        dot = (x - x0) * dx + (y - y0) * dy
        length_sq = dx * dx + dy * dy
        if abs(cross) <= tolerance * max(math.sqrt(length_sq), 1.0) and (
                -tolerance <= dot <= length_sq + tolerance):
            return True
        if (y0 > y) != (y1 > y):
            crossing_x = x0 + (y - y0) * dx / dy
            if x < crossing_x:
                inside = not inside
        previous = current
    return inside


def polygon_contains(polygon: SkyPolygon, ra_deg: float, dec_deg: float) -> bool:
    point = _point_xy(ra_deg, dec_deg, polygon.center_ra_deg, polygon.center_dec_deg)
    if point is None:
        return False
    vertices = [
        _point_xy(ra, dec, polygon.center_ra_deg, polygon.center_dec_deg)
        for ra, dec in polygon.corners_ra_dec
    ]
    if any(vertex is None for vertex in vertices):
        raise ValueError(f"polygon for tract {polygon.tract} spans a gnomonic horizon")
    return _point_in_xy_polygon(point, np.asarray(vertices, dtype=float))


def _normalize_corners(corners: dict[int, tuple[float, float]], label: str) -> tuple[tuple[float, float], ...]:
    if sorted(corners) != list(range(len(corners))):
        raise ValueError(f"non-consecutive corner numbering in {label}")
    values = [corners[index] for index in sorted(corners)]
    if len(values) == 5 and np.allclose(values[0], values[-1], rtol=0, atol=1e-9):
        values.pop()
    if len(values) != 4:
        raise ValueError(f"expected four polygon corners in {label}; found {len(values)}")
    return tuple(values)


def _validate_sky_point(point: tuple[float, float], label: str) -> tuple[float, float]:
    ra, dec = map(float, point)
    if not (math.isfinite(ra) and 0.0 <= ra < 360.0
            and math.isfinite(dec) and -90.0 <= dec <= 90.0):
        raise ValueError(f"invalid sky coordinate in {label}: {point}")
    return ra, dec


def parse_geometry_file(path: Path) -> tuple[list[SkyPolygon], list[SkyPolygon]]:
    path = Path(path)
    try:
        field = WIDE_FILES[path.name]
    except KeyError as exc:
        raise ValueError(f"unrecognized HSC PDR3 Wide geometry filename: {path.name}") from exc

    section: str | None = None
    tract_records: dict[int, dict] = {}
    patch_records: dict[tuple[int, tuple[int, int]], dict] = {}
    current_tract: int | None = None
    current_patch: tuple[int, int] | None = None
    saw_tract_section = saw_patch_section = False

    with path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            text = line.strip()
            if not text:
                continue
            if TRACT_MARKER.match(text):
                section = "tract"
                saw_tract_section = True
                current_tract = None
                current_patch = None
                continue
            if PATCH_MARKER.match(text):
                section = "patch"
                saw_patch_section = True
                current_patch = None
                continue

            tract_match = TRACT_ID.search(text)
            patch_match = PATCH_ID.search(text) if section == "patch" else None
            if tract_match:
                next_tract = int(tract_match.group(1))
                if next_tract != current_tract:
                    current_patch = None
                current_tract = next_tract
            if section == "tract" and current_tract is not None:
                record = tract_records.setdefault(current_tract, {"center": None, "corners": {}})
            elif section == "patch" and current_tract is not None and patch_match:
                current_patch = (int(patch_match.group(1)), int(patch_match.group(2)))
                key = (current_tract, current_patch)
                record = patch_records.setdefault(key, {"center": None, "corners": {}})
            elif section == "patch" and current_patch is not None and current_tract is not None:
                record = patch_records[(current_tract, current_patch)]
            else:
                continue

            point_match = SKY_POINT.search(text)
            if not point_match:
                continue
            kind, corner_number, ra_text, dec_text = point_match.groups()
            point = _validate_sky_point((float(ra_text), float(dec_text)), f"{path.name}:{line_number}")
            if kind == "Center":
                if record["center"] is not None:
                    raise ValueError(f"duplicate center in {path.name}:{line_number}")
                record["center"] = point
            else:
                index = int(corner_number)
                if index in record["corners"]:
                    raise ValueError(f"duplicate corner in {path.name}:{line_number}")
                record["corners"][index] = point

    if not (saw_tract_section and saw_patch_section and tract_records and patch_records):
        raise ValueError(f"incomplete tract/patch geometry file: {path}")

    tracts = []
    for tract_id, record in sorted(tract_records.items()):
        if record["center"] is None:
            raise ValueError(f"missing center for {field} tract {tract_id}")
        polygon = SkyPolygon(
            field, tract_id, None, *record["center"],
            _normalize_corners(record["corners"], f"{field} tract {tract_id}"),
        )
        if not polygon_contains(polygon, *record["center"]):
            raise ValueError(f"center lies outside polygon for {field} tract {tract_id}")
        tracts.append(polygon)

    tract_ids = set(tract_records)
    patches = []
    for (tract_id, patch_xy), record in sorted(patch_records.items()):
        if tract_id not in tract_ids:
            raise ValueError(f"patch refers to absent {field} tract {tract_id}")
        if record["center"] is None:
            raise ValueError(f"missing center for {field} tract {tract_id} patch {patch_xy}")
        polygon = SkyPolygon(
            field, tract_id, patch_xy, *record["center"],
            _normalize_corners(record["corners"], f"{field} tract {tract_id} patch {patch_xy}"),
        )
        if not polygon_contains(polygon, *record["center"]):
            raise ValueError(f"center lies outside polygon for {field} tract {tract_id} patch {patch_xy}")
        patches.append(polygon)

    return tracts, patches


def load_wide_geometry(directory: Path) -> tuple[list[SkyPolygon], list[SkyPolygon], dict]:
    directory = Path(directory)
    tracts: list[SkyPolygon] = []
    patches: list[SkyPolygon] = []
    files = {}
    for filename, field in WIDE_FILES.items():
        path = directory / filename
        if not path.is_file():
            raise FileNotFoundError(f"required official HSC PDR3 Wide geometry file is missing: {path}")
        observed_bytes = path.stat().st_size
        if observed_bytes != WIDE_FILE_BYTES[filename]:
            raise ValueError(
                f"official geometry file size differs from the approved scope for {filename}: "
                f"{observed_bytes} != {WIDE_FILE_BYTES[filename]}"
            )
        tract_rows, patch_rows = parse_geometry_file(path)
        if any(row.field != field for row in tract_rows + patch_rows):
            raise ValueError(f"field label differs from filename for {path}")
        tracts.extend(tract_rows)
        patches.extend(patch_rows)
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1 << 20), b""):
                digest.update(chunk)
        files[filename] = {
            "bytes": observed_bytes, "sha256": digest.hexdigest(),
        }

    tract_keys = [(row.field, row.tract) for row in tracts]
    patch_keys = [(row.field, row.tract, row.patch) for row in patches]
    if len(tract_keys) != len(set(tract_keys)) or len(patch_keys) != len(set(patch_keys)):
        raise ValueError("duplicate tract or patch identifiers across official field files")
    summary = {
        "source": "official HSC PDR3 Wide tract/patch coordinate files",
        "n_fields": len(WIDE_FILES),
        "expected_total_bytes": sum(WIDE_FILE_BYTES.values()),
        "n_tract_polygons": len(tracts),
        "n_patch_polygons": len(patches),
        "files": files,
    }
    return tracts, patches, summary


class PolygonIndex:
    """Spatial index over HSC polygons with deterministic overlap handling.

    ``field`` mode uses tract polygons only to identify a field and leaves the
    tract and patch outputs unset. ``tract`` and ``patch`` modes report the
    nearest containing polygon within an unambiguous field.
    """

    def __init__(self, polygons: Iterable[SkyPolygon]):
        self.polygons = list(polygons)
        if not self.polygons:
            raise ValueError("cannot index an empty polygon set")
        centers = unit_vectors(
            [row.center_ra_deg for row in self.polygons],
            [row.center_dec_deg for row in self.polygons],
        )
        self.tree = cKDTree(centers)
        self.center_vectors = centers
        max_extent = 0.0
        for polygon, vector in zip(self.polygons, centers):
            corners = unit_vectors(
                [point[0] for point in polygon.corners_ra_dec],
                [point[1] for point in polygon.corners_ra_dec],
            )
            dots = np.clip(corners @ vector, -1.0, 1.0)
            max_extent = max(max_extent, float(np.max(np.arccos(dots))))
        self.search_radius_rad = max_extent + math.radians(1.0 / 3600.0)
        self.search_chord = 2.0 * math.sin(self.search_radius_rad / 2.0)

    def assign(self, ra_deg, dec_deg, *, level: str = "tract",
               chunk_size: int = 50_000) -> PointAssignments:
        if level not in {"field", "tract", "patch"}:
            raise ValueError("level must be field, tract, or patch")
        ra = np.asarray(ra_deg, dtype=float)
        dec = np.asarray(dec_deg, dtype=float)
        points = unit_vectors(ra, dec)
        n = len(ra)
        fields = np.full(n, None, dtype=object)
        tracts = np.full(n, -1, dtype=np.int64)
        patches = np.full((n, 2), -1, dtype=np.int32)
        match_counts = np.zeros(n, dtype=np.int16)
        cross_field = np.zeros(n, dtype=bool)

        for start in range(0, n, chunk_size):
            stop = min(start + chunk_size, n)
            candidate_groups = self.tree.query_ball_point(points[start:stop], self.search_chord)
            for offset, candidate_indices in enumerate(candidate_groups):
                row_index = start + offset
                if not candidate_indices:
                    continue
                point_ra, point_dec = float(ra[row_index]), float(dec[row_index])
                matches = [
                    index for index in candidate_indices
                    if polygon_contains(self.polygons[index], point_ra, point_dec)
                ]
                match_counts[row_index] = min(len(matches), np.iinfo(np.int16).max)
                if not matches:
                    continue
                match_fields = {self.polygons[index].field for index in matches}
                if len(match_fields) > 1:
                    cross_field[row_index] = True
                    continue
                if level == "field":
                    fields[row_index] = next(iter(match_fields))
                    continue
                point_vector = points[row_index]
                chosen = min(matches, key=lambda index: (
                    math.acos(float(np.clip(
                        np.dot(point_vector, self.center_vectors[index]), -1.0, 1.0,
                    ))),
                    self.polygons[index].tract,
                    self.polygons[index].patch or (-1, -1),
                ))
                polygon = self.polygons[chosen]
                fields[row_index] = polygon.field
                tracts[row_index] = polygon.tract
                if polygon.patch is not None:
                    patches[row_index] = polygon.patch
        return PointAssignments(fields, tracts, patches, match_counts, cross_field)


def load_index(directory: Path, *, level: str = "tract") -> tuple[PolygonIndex, dict]:
    if level not in {"field", "tract", "patch"}:
        raise ValueError("level must be field, tract, or patch")
    tracts, patches, summary = load_wide_geometry(directory)
    polygons = patches if level == "patch" else tracts
    summary["indexed_level"] = level
    return PolygonIndex(polygons), summary


def _write_assignments(path: Path, identifiers: list[str], assignments: PointAssignments) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite coordinate assignment output: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("id", "field", "tract", "patch_x", "patch_y", "containing_polygons", "cross_field_ambiguous"))
        for index, identifier in enumerate(identifiers):
            writer.writerow((
                identifier,
                assignments.field[index] or "",
                int(assignments.tract[index]),
                int(assignments.patch[index, 0]),
                int(assignments.patch[index, 1]),
                int(assignments.containing_polygon_count[index]),
                bool(assignments.cross_field_ambiguous[index]),
            ))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True,
                        help="directory containing the four explicit HSC PDR3 Wide coordinate files")
    parser.add_argument("--level", choices=("field", "tract", "patch"), default="tract")
    parser.add_argument("--points-tsv", type=Path,
                        help="optional TSV with ra_deg and dec_deg columns, and optional id")
    parser.add_argument("--output", type=Path,
                        help="new TSV path for assignments; refuses to overwrite")
    args = parser.parse_args()
    index, summary = load_index(args.directory, level=args.level)
    if args.points_tsv:
        if args.output is None:
            parser.error("--output is required with --points-tsv")
        with args.points_tsv.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            required = {"ra_deg", "dec_deg"}
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise ValueError("point TSV must include ra_deg and dec_deg columns")
            rows = list(reader)
        identifiers = [row.get("id", str(position)) for position, row in enumerate(rows)]
        assigned = index.assign(
            [float(row["ra_deg"]) for row in rows],
            [float(row["dec_deg"]) for row in rows], level=args.level,
        )
        _write_assignments(args.output, identifiers, assigned)
        summary["points"] = {
            "count": len(rows),
            "unmatched": assigned.unmatched_count,
            "overlapping_polygon_matches": assigned.overlap_count,
            "cross_field_ambiguous": assigned.cross_field_ambiguous_count,
            "assignment_output": str(args.output.resolve()),
        }
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
