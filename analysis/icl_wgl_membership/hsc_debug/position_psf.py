#!/usr/bin/env python3
"""Load and apply the approved position-dependent HSC model PSF cores."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.signal import fftconvolve


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _coord_key(ra: float, dec: float) -> tuple[float, float]:
    return round(float(ra), 9), round(float(dec), 9)


class PositionPSFBank:
    """Resolve planned cluster/member roles to audited PSF FITS files.

    ``plan_dir`` and ``download_dir`` must point to the immutable approved
    744-position request and its audited response. The already-approved
    cluster-10 center test supplies the one center omitted from that batch.
    """

    def __init__(self, plan_dir: Path, download_dir: Path,
                 test_center_psf: Path):
        plan_dir = Path(plan_dir)
        download_dir = Path(download_dir)
        plan = json.loads((plan_dir / "psf_picker_plan.json").read_text())
        mapping_path = plan_dir / "position_mapping.json"
        coordinate_path = plan_dir / "psf_picker_coordinates.txt"
        mapping = json.loads(mapping_path.read_text())
        audit_path = download_dir / "psf_picker_audit.json"
        audit = json.loads(audit_path.read_text())
        if plan.get("schema") != "astrafek-hsc-psf-picker-injection-plan-v1":
            raise ValueError("unexpected PSF plan schema")
        if plan.get("user_approved_requests") != 744 or plan.get("filter") != "HSC-I":
            raise ValueError("PSF plan does not match the approved HSC-I request")
        if sha256(mapping_path) != plan.get("position_mapping_sha256"):
            raise ValueError("PSF mapping checksum differs from approved plan")
        if sha256(coordinate_path) != plan.get("coordinate_list_sha256"):
            raise ValueError("PSF coordinate-list checksum differs from approved plan")
        if audit.get("status") != "verified" or audit.get("n_returned") != 744:
            raise ValueError("PSF batch has not passed the 744-file audit")
        if audit.get("position_mapping_sha256") != plan["position_mapping_sha256"]:
            raise ValueError("PSF audit was produced from a different position mapping")

        self._path_by_coord: dict[tuple[float, float], tuple[Path, str]] = {}
        for record in audit["files"]:
            key = _coord_key(record["ra_deg_filename"], record["dec_deg_filename"])
            if key in self._path_by_coord:
                raise ValueError("duplicate coordinate in audited PSF manifest")
            self._path_by_coord[key] = (download_dir / record["file"], record["sha256"])
        if len(self._path_by_coord) != 744:
            raise ValueError("expected exactly 744 audited PSF coordinates")

        self._role_coord: dict[tuple[str, str, str, int | None], tuple[float, float]] = {}
        for position in mapping:
            key = _coord_key(position["ra_deg"], position["dec_deg"])
            if key not in self._path_by_coord:
                raise ValueError("planned role has no returned PSF file")
            for role in position["roles"]:
                role_key = (str(role["cluster_id"]), str(role["frame_filter"]),
                            str(role["role"]),
                            None if role.get("target_index") is None
                            else int(role["target_index"]))
                if role_key in self._role_coord:
                    raise ValueError(f"duplicate PSF role mapping: {role_key}")
                self._role_coord[role_key] = key
        self.supported_frames = {
            (cluster, filt) for cluster, filt, role, _ in self._role_coord
            if role == "cluster_center"
        }
        # The already-approved center test is intentionally excluded from
        # the remaining-744 coordinate mapping.
        self.supported_frames.add(("10", "i2"))
        for cluster, filt, role, _ in self._role_coord:
            if role == "satellite_injection_center" and (cluster, filt) not in self.supported_frames:
                raise ValueError(f"satellite PSF roles lack a center PSF: {cluster}/{filt}")

        # The approved one-position test is the cluster-10/i2 center omitted
        # from the remaining-744 request plan.
        test_center_psf = Path(test_center_psf)
        if not test_center_psf.is_file():
            raise FileNotFoundError(test_center_psf)
        if sha256(test_center_psf) != "1276d46113d389116ded6f63bae4c9865bb42d24244b22247a84d2c6929cb202":
            raise ValueError("cluster-10 test PSF checksum differs from its audit record")
        test_key = ("10", "i2", "cluster_center", None)
        if test_key in self._role_coord:
            raise ValueError("test PSF role unexpectedly appears in remainder mapping")
        self._test_center_psf = test_center_psf
        self._cache: dict[Path, np.ndarray] = {}
        self.plan_path = plan_dir / "psf_picker_plan.json"
        self.mapping_path = mapping_path
        self.audit_path = audit_path
        self._test_center_sha = sha256(test_center_psf)

    def _load(self, path: Path, expected_sha: str | None = None) -> np.ndarray:
        if path not in self._cache:
            if expected_sha is not None and sha256(path) != expected_sha:
                raise ValueError(f"PSF file checksum mismatch: {path.name}")
            with fits.open(path, memmap=False, checksum=True) as hdus:
                if len(hdus) != 1 or hdus[0].data is None:
                    raise ValueError(f"unexpected PSF FITS structure: {path.name}")
                psf = np.asarray(hdus[0].data, dtype=np.float64)
            if psf.ndim != 2 or any(size % 2 != 1 for size in psf.shape):
                raise ValueError(f"PSF must have odd two-dimensional shape: {path.name}")
            if not np.isfinite(psf).all() or abs(float(psf.sum()) - 1.0) > 1e-3:
                raise ValueError(f"PSF is non-finite or not unit-normalized: {path.name}")
            positive = np.clip(psf, 0, None)
            yy, xx = np.indices(psf.shape)
            cx = float((positive * xx).sum() / positive.sum())
            cy = float((positive * yy).sum() / positive.sum())
            if np.hypot(cx - (psf.shape[1] - 1) / 2,
                        cy - (psf.shape[0] - 1) / 2) > 0.5:
                raise ValueError(f"PSF centroid is not centered: {path.name}")
            self._cache[path] = psf.astype(np.float32)
        return self._cache[path]

    def frame(self, cluster_id: str, frame_filter: str,
              n_satellites: int) -> tuple[np.ndarray, list[np.ndarray], list[str]]:
        """Return center and ordered per-satellite PSFs plus their file names."""
        ident, filt = str(cluster_id), str(frame_filter).lower()
        center_key = (ident, filt, "cluster_center", None)
        if center_key == ("10", "i2", "cluster_center", None):
            center_path, center_hash = self._test_center_psf, None
        else:
            try:
                coordinate = self._role_coord[center_key]
                center_path, center_hash = self._path_by_coord[coordinate]
            except KeyError as exc:
                raise ValueError(f"no approved center PSF for cluster {ident}/{filt}") from exc
        center = self._load(center_path, center_hash)
        satellites = []
        names = [center_path.name]
        for index in range(n_satellites):
            role_key = (ident, filt, "satellite_injection_center", index)
            try:
                coordinate = self._role_coord[role_key]
                path, digest = self._path_by_coord[coordinate]
            except KeyError as exc:
                raise ValueError(f"no approved PSF for {role_key}") from exc
            satellites.append(self._load(path, digest))
            names.append(path.name)
        return center, satellites, names

    def file_hashes(self, names: list[str]) -> dict[str, str]:
        """Return hashes for the named PSFs used by one frame."""
        wanted = set(names)
        result = {}
        for path, digest in self._path_by_coord.values():
            if path.name in wanted:
                result[str(path)] = digest
        if self._test_center_psf.name in wanted:
            result[str(self._test_center_psf)] = self._test_center_sha
        if len(result) != len(wanted):
            raise ValueError("could not resolve every used PSF filename to a checksum")
        return result


class SupportedPairsPositionPSFBank(PositionPSFBank):
    """Combine audited PSFs for the 224-pair injection-support sample.

    The base 744-position injection batch supplies PSFs for the original
    57-frame subset. A separate center batch supplies target centers, and an
    approved supported-pairs batch supplies the remaining satellite positions.
    The constructor verifies each immutable plan, its coordinate mapping, the
    audit manifest, and every file checksum before exposing a frame.
    """

    def __init__(self, injection_plan_dir: Path, injection_download_dir: Path,
                 test_center_psf: Path, center_plan_dir: Path,
                 center_download_dir: Path, satellite_plan_dir: Path,
                 satellite_download_dir: Path,
                 approved_satellite_plan_sha256: str):
        super().__init__(injection_plan_dir, injection_download_dir, test_center_psf)
        center_plan_dir, center_download_dir = Path(center_plan_dir), Path(center_download_dir)
        satellite_plan_dir, satellite_download_dir = Path(satellite_plan_dir), Path(satellite_download_dir)
        satellite_plan_path = satellite_plan_dir / "psf_request_plan.json"
        actual_approval_hash = sha256(satellite_plan_path)
        if actual_approval_hash != str(approved_satellite_plan_sha256).lower():
            raise ValueError(
                "supported-satellite PSF plan hash differs from the explicitly approved plan"
            )
        self.approved_satellite_plan_sha256 = actual_approval_hash
        self.provenance_paths = [self.plan_path, self.mapping_path, self.audit_path]
        self._center_coord_by_cluster: dict[str, tuple[float, float]] = {}
        self._load_center_batch(center_plan_dir, center_download_dir)
        self._load_satellite_batch(satellite_plan_dir, satellite_download_dir)

        satellite_indices: dict[tuple[str, str], set[int]] = {}
        for cluster, filt, role, target_index in self._role_coord:
            if role == "satellite_injection_center":
                if target_index is None:
                    raise ValueError("satellite PSF role is missing its target index")
                satellite_indices.setdefault((cluster, filt), set()).add(target_index)
        incomplete = {
            frame: sorted(indices)
            for frame, indices in satellite_indices.items()
            if indices != set(range(max(indices) + 1))
        }
        if incomplete:
            first = next(iter(sorted(incomplete.items())))
            raise ValueError(f"combined satellite PSF indices are incomplete for {first[0]}: {first[1]}")
        self.supported_frames = set(satellite_indices)
        self.request_approval_verified = True
        if len(self.supported_frames) != 224:
            raise ValueError(f"expected 224 PSF-complete ICL support frames, found {len(self.supported_frames)}")
        missing_centers = [frame for frame in self.supported_frames
                           if (frame[0], frame[1], "cluster_center", None) not in self._role_coord]
        if missing_centers:
            raise ValueError(f"supported injection frames lack center PSFs: {missing_centers[:10]}")

    @staticmethod
    def _read_mapping(plan_dir: Path, plan_name: str, audit_path: Path,
                      expected_schema: str, count_field: str) -> tuple[dict, list, dict]:
        plan_path = plan_dir / plan_name
        mapping_path = plan_dir / "position_mapping.json"
        coordinate_path = plan_dir / "psf_picker_coordinates.txt"
        if not all(path.is_file() for path in (plan_path, mapping_path, coordinate_path, audit_path)):
            raise FileNotFoundError(f"PSF plan, mapping, coordinate list, or audit is missing under {plan_dir}")
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
        mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if plan.get("schema") != expected_schema or plan.get("status") != "prepared_not_submitted":
            raise ValueError(f"unexpected or modified PSF request plan: {plan_path}")
        if plan.get("request_submitted") is True or plan.get("network_used") is True:
            raise ValueError(f"plan provenance unexpectedly records a submitted request: {plan_path}")
        request = plan.get("request", plan)
        if (request.get("rerun") != "pdr3_wide" or request.get("type") != "coadd"
                or request.get("filter") != "HSC-I"
                or request.get("centered", plan.get("centered")) is not True):
            raise ValueError(f"PSF plan is outside the centered PDR3 Wide HSC-I coadd scope: {plan_path}")
        if (plan.get("images_requested") is not False
                or plan.get("shape_rows_requested") is not False
                or plan.get("photoz_requested") is not False
                or plan.get("credentials_saved") is not False):
            raise ValueError(f"PSF plan includes out-of-scope products: {plan_path}")
        expected = int(plan.get(count_field, -1))
        if expected <= 0 or len(mapping) != expected:
            raise ValueError(f"PSF plan/mapping count mismatch: {plan_path}")
        if sha256(mapping_path) != plan.get("position_mapping_sha256"):
            raise ValueError(f"PSF position mapping checksum mismatch: {mapping_path}")
        if sha256(coordinate_path) != plan.get("coordinate_list_sha256"):
            raise ValueError(f"PSF coordinate list checksum mismatch: {coordinate_path}")
        if (audit.get("status") != "verified"
                or audit.get("n_returned") != expected
                or audit.get("n_matched_coordinates") != expected
                or audit.get("position_mapping_sha256") != sha256(mapping_path)):
            raise ValueError(f"PSF response audit does not match plan: {audit_path}")
        mapped = {_coord_key(row["ra_deg"], row["dec_deg"]): row for row in mapping}
        if len(mapped) != expected:
            raise ValueError(f"duplicate PSF coordinates in mapping: {mapping_path}")
        lines = coordinate_path.read_text(encoding="ascii").splitlines()
        if not lines or lines[0].strip() != "#? ra dec":
            raise ValueError(f"unexpected PSF coordinate list header: {coordinate_path}")
        listed = []
        for line in lines[1:]:
            fields = line.split()
            if len(fields) != 2:
                raise ValueError(f"malformed PSF coordinate row: {coordinate_path}")
            listed.append(_coord_key(fields[0], fields[1]))
        if len(listed) != expected or len(set(listed)) != expected or set(listed) != set(mapped):
            raise ValueError(f"PSF coordinate list differs from mapping: {coordinate_path}")
        return {"plan": plan_path, "mapping": mapping_path, "coordinates": coordinate_path}, mapping, audit

    def _merge_batch(self, plan_dir: Path, download_dir: Path, plan_name: str,
                     audit_name: str, expected_schema: str, count_field: str) -> list[dict]:
        audit_path = download_dir / audit_name
        paths, mapping, audit = self._read_mapping(
            plan_dir, plan_name, audit_path, expected_schema, count_field)
        plan_path = paths["plan"]
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
        if audit.get("plan_sha256") != sha256(plan_path):
            raise ValueError(f"PSF response audit plan checksum mismatch: {audit_path}")
        if (audit.get("credentials_saved") is not False
                or audit.get("images_returned") is not False
                or audit.get("shape_rows_returned") is not False
                or audit.get("photoz_rows_returned") is not False
                or audit.get("science_measurement") is not False):
            raise ValueError(f"unexpected returned data or gate state in {audit_path}")
        file_rows = audit.get("files")
        if not isinstance(file_rows, list) or len(file_rows) != len(mapping):
            raise ValueError(f"PSF file audit count mismatch: {audit_path}")
        expected_names = {str(row["file"]) for row in file_rows}
        actual_names = {path.name for path in download_dir.glob("*.fits")}
        if actual_names != expected_names:
            raise ValueError(f"download directory FITS set differs from its audit: {download_dir}")
        response_by_coord = {}
        for row in file_rows:
            coord = _coord_key(row["ra_deg_filename"], row["dec_deg_filename"])
            if coord in response_by_coord:
                raise ValueError(f"duplicate PSF coordinate in response audit: {audit_path}")
            path = download_dir / row["file"]
            if not path.is_file() or sha256(path) != row.get("sha256"):
                raise ValueError(f"missing or changed PSF file: {path}")
            response_by_coord[coord] = (path, row["sha256"])
        mapped_coords = {_coord_key(row["ra_deg"], row["dec_deg"])
                         for row in mapping}
        if set(response_by_coord) != mapped_coords:
            raise ValueError(f"returned PSF coordinates differ from plan mapping: {audit_path}")

        for coord, (path, digest) in response_by_coord.items():
            if coord in self._path_by_coord:
                raise ValueError(f"supplemental PSF batch duplicates an existing coordinate: {coord}")
            self._path_by_coord[coord] = (path, digest)
        self.provenance_paths.extend([paths["plan"], paths["mapping"], paths["coordinates"], audit_path])
        return mapping

    def _load_center_batch(self, plan_dir: Path, download_dir: Path) -> None:
        mapping = self._merge_batch(
            plan_dir, download_dir, "psf_picker_plan.json",
            "psf_picker_audit_reproducible.json",
            "astrafek-hsc-psf-picker-target-control-center-plan-v1",
            "n_requests_remaining")
        for row in mapping:
            coord = _coord_key(row["ra_deg"], row["dec_deg"])
            for role in row.get("roles", []):
                if role.get("role") != "target_cluster_center":
                    continue
                cluster = str(role["cluster_id"])
                previous = self._center_coord_by_cluster.setdefault(cluster, coord)
                if previous != coord:
                    raise ValueError(f"multiple target-center PSFs for cluster {cluster}")

    def _load_satellite_batch(self, plan_dir: Path, download_dir: Path) -> None:
        mapping = self._merge_batch(
            plan_dir, download_dir, "psf_request_plan.json", "psf_picker_audit.json",
            "astrafek-icl-supported-pairs-psf-request-plan-v1", "n_requests_remaining")
        satellite_roles = 0
        for row in mapping:
            coord = _coord_key(row["ra_deg"], row["dec_deg"])
            for role in row.get("roles", []):
                if role.get("role") != "satellite_injection_center":
                    continue
                key = (str(role["cluster_id"]), str(role["frame_filter"]).lower(),
                       "satellite_injection_center", int(role["target_index"]))
                previous = self._role_coord.setdefault(key, coord)
                if previous != coord:
                    raise ValueError(f"conflicting satellite PSF role mapping: {key}")
                satellite_roles += 1
        if satellite_roles != 2005:
            raise ValueError(f"expected 2,005 new satellite PSF roles, found {satellite_roles}")

        for cluster, filt, role, _ in list(self._role_coord):
            if role != "satellite_injection_center":
                continue
            center_key = (cluster, filt, "cluster_center", None)
            if center_key in self._role_coord:
                continue
            try:
                coord = self._center_coord_by_cluster[cluster]
            except KeyError as exc:
                raise ValueError(f"no center PSF for supported cluster {cluster}") from exc
            self._role_coord[center_key] = coord

    def frame(self, cluster_id: str, frame_filter: str,
              n_satellites: int) -> tuple[np.ndarray, list[np.ndarray], list[str]]:
        if (str(cluster_id), str(frame_filter).lower()) not in self.supported_frames:
            raise ValueError(f"frame is outside the verified 224-pair PSF support: {cluster_id}/{frame_filter}")
        return super().frame(cluster_id, frame_filter, n_satellites)


class S23BCamiraMemberProfilePSFBank:
    """Resolve the separately approved 427-position S23B member PSF batch."""

    def __init__(self, plan_dir: Path, download_dir: Path,
                 approved_plan_sha256: str):
        plan_dir, download_dir = Path(plan_dir), Path(download_dir)
        self.plan_path = plan_dir / "psf_request_plan.json"
        self.mapping_path = plan_dir / "position_mapping.json"
        self.coordinate_path = plan_dir / "psf_picker_coordinates.txt"
        self.audit_path = download_dir / "psf_picker_audit.json"
        plan = json.loads(self.plan_path.read_text(encoding="utf-8"))
        mapping = json.loads(self.mapping_path.read_text(encoding="utf-8"))
        audit = json.loads(self.audit_path.read_text(encoding="utf-8"))
        actual_plan_sha256 = sha256(self.plan_path)
        if actual_plan_sha256 != str(approved_plan_sha256).lower():
            raise ValueError("member-profile PSF plan differs from the explicitly approved plan")
        if (plan.get("schema") != "astrafek-s23b-camira-member-profile-psf-plan-v1"
                or plan.get("status") != "prepared_not_submitted"
                or plan.get("n_unique_member_positions") != 427
                or plan.get("requires_explicit_user_approval") is not True
                or plan.get("science_measurement") is not False):
            raise ValueError("unexpected or modified S23B member-profile PSF plan")
        request = plan.get("request", {})
        if (request.get("rerun") != "pdr3_wide" or request.get("type") != "coadd"
                or request.get("filter") != "HSC-I" or request.get("centered") is not True
                or request.get("product") != "position-specific coadd model PSF FITS core"
                or plan.get("images_requested") is not False
                or plan.get("coadd_images_requested") is not False
                or plan.get("shape_rows_requested") is not False
                or plan.get("photoz_requested") is not False):
            raise ValueError("member-profile PSF plan falls outside its approved data scope")
        if (sha256(self.mapping_path) != plan.get("position_mapping_sha256")
                or sha256(self.coordinate_path) != plan.get("coordinate_list_sha256")):
            raise ValueError("member-profile PSF mapping or coordinate checksum mismatch")
        if (audit.get("status") != "verified"
                or audit.get("plan_sha256") != actual_plan_sha256
                or audit.get("approved_plan_sha256") != actual_plan_sha256
                or audit.get("n_returned") != 427
                or audit.get("n_matched_coordinates") != 427
                or audit.get("position_mapping_sha256") != plan["position_mapping_sha256"]
                or audit.get("total_bytes", plan["max_total_bytes_for_review"] + 1)
                   > plan["max_total_bytes_for_review"]
                or audit.get("images_returned") is not False
                or audit.get("shape_rows_returned") is not False
                or audit.get("photoz_rows_returned") is not False
                or audit.get("science_measurement") is not False):
            raise ValueError("member-profile PSF response has not passed its scoped audit")

        self._path_by_coord: dict[tuple[float, float], tuple[Path, str]] = {}
        for record in audit.get("files", []):
            coord = _coord_key(record["ra_deg_filename"], record["dec_deg_filename"])
            if coord in self._path_by_coord:
                raise ValueError("duplicate position in member-profile PSF audit")
            path = download_dir / record["file"]
            if not path.is_file() or path.stat().st_size != int(record["bytes"]):
                raise ValueError(f"missing or size-mismatched member PSF: {path}")
            self._path_by_coord[coord] = (path, str(record["sha256"]))
        if len(self._path_by_coord) != 427:
            raise ValueError("expected exactly 427 audited member-profile PSF files")

        self._member_coord: dict[tuple[str, str, int, str], tuple[float, float]] = {}
        if len(mapping) != 427:
            raise ValueError("expected exactly 427 unique member-profile PSF requests")
        for position in mapping:
            coord = _coord_key(position["ra_deg"], position["dec_deg"])
            if coord not in self._path_by_coord:
                raise ValueError("planned member coordinate has no audited PSF file")
            for role in position.get("roles", []):
                key = (str(role["cluster_id"]), str(role["frame_filter"]).lower(),
                       int(role["member_index"]), str(role["source_row"]))
                previous = self._member_coord.setdefault(key, coord)
                if previous != coord:
                    raise ValueError(f"conflicting PSF coordinate for member {key}")
        if len(self._member_coord) != 427:
            raise ValueError("expected exactly 427 selected CAMIRA member PSF roles")

        self.approved_plan_sha256 = actual_plan_sha256
        self.provenance_paths = [self.plan_path, self.mapping_path,
                                 self.coordinate_path, self.audit_path]
        self._cache: dict[Path, np.ndarray] = {}

    def member(self, cluster_id: str, frame_filter: str, member_index: int,
               source_row: str) -> tuple[np.ndarray, str, str]:
        """Return the verified local PSF, filename, and recorded file digest."""
        key = (str(cluster_id), str(frame_filter).lower(), int(member_index), str(source_row))
        try:
            coord = self._member_coord[key]
            path, expected_sha = self._path_by_coord[coord]
        except KeyError as exc:
            raise ValueError(f"no approved member-profile PSF for {key}") from exc
        if not path.is_file() or sha256(path) != expected_sha:
            raise ValueError(f"member-profile PSF checksum mismatch: {path.name}")
        if path not in self._cache:
            with fits.open(path, memmap=False, checksum=True) as hdus:
                if len(hdus) != 1 or hdus[0].data is None:
                    raise ValueError(f"unexpected PSF FITS structure: {path.name}")
                psf = np.asarray(hdus[0].data, dtype=np.float64)
            if (psf.ndim != 2 or any(size % 2 != 1 for size in psf.shape)
                    or not np.isfinite(psf).all()
                    or abs(float(psf.sum()) - 1.0) > 1e-3):
                raise ValueError(f"invalid signed unit-normalized PSF core: {path.name}")
            positive = np.clip(psf, 0, None)
            yy, xx = np.indices(psf.shape)
            cx = float((positive * xx).sum() / positive.sum())
            cy = float((positive * yy).sum() / positive.sum())
            if np.hypot(cx - (psf.shape[1] - 1) / 2,
                        cy - (psf.shape[0] - 1) / 2) > 0.5:
                raise ValueError(f"member PSF is not centered: {path.name}")
            self._cache[path] = psf.astype(np.float32)
        return self._cache[path], path.name, expected_sha

    def manifest(self) -> dict:
        """Return immutable request, audit, and file provenance for a run."""
        return {
            "plan": str(self.plan_path.resolve()),
            "plan_sha256": self.approved_plan_sha256,
            "mapping": str(self.mapping_path.resolve()),
            "mapping_sha256": sha256(self.mapping_path),
            "coordinates": str(self.coordinate_path.resolve()),
            "coordinate_list_sha256": sha256(self.coordinate_path),
            "audit": str(self.audit_path.resolve()),
            "audit_sha256": sha256(self.audit_path),
            "files_sha256": {
                str(path.resolve()): digest
                for path, digest in sorted(self._path_by_coord.values(), key=lambda item: str(item[0]))
            },
            "n_files": len(self._path_by_coord),
            "science_measurement": False,
        }


def convolve_surface_brightness(template: np.ndarray, psf: np.ndarray) -> np.ndarray:
    """Convolve SB with a signed, unit-normalized centered model PSF core."""
    image = np.asarray(template, dtype=np.float32)
    kernel = np.asarray(psf, dtype=np.float32)
    if image.ndim != 2 or kernel.ndim != 2 or any(size % 2 != 1 for size in kernel.shape):
        raise ValueError("template and odd-sized 2-D PSF are required")
    if not np.isfinite(image).all() or not np.isfinite(kernel).all():
        raise ValueError("template and PSF must be finite")
    if abs(float(kernel.sum()) - 1.0) > 1e-3:
        raise ValueError("PSF must retain its signed, unit-normalized values")
    return fftconvolve(image, kernel, mode="same").astype(np.float32)


def build_core_wing_hybrid_psf(local_core: np.ndarray, public_wing: np.ndarray,
                               splice_radius_arcsec: float,
                               blend_width_arcsec: float,
                               wing_amplitude: float,
                               pixel_scale_arcsec: float = 0.168
                               ) -> tuple[np.ndarray, dict]:
    """Join a position-specific PSF core to a public extended wing for sensitivity tests.

    The local signed core is retained inside the transition annulus; a smooth
    radial taper blends it into the scaled public model outside.  The returned
    kernel is unit-sum normalized.  This is an exploratory model constructor,
    not an observational PSF calibration or a correction to science images.
    """
    core = np.asarray(local_core, dtype=np.float64)
    wing = np.asarray(public_wing, dtype=np.float64)
    if core.ndim != 2 or wing.ndim != 2 or any(n % 2 != 1 for n in (*core.shape, *wing.shape)):
        raise ValueError("local core and public wing must be odd-sized 2-D arrays")
    if wing.shape[0] != wing.shape[1]:
        raise ValueError("public wing must be square")
    if core.shape[0] > wing.shape[0] or core.shape[1] > wing.shape[1]:
        raise ValueError("public wing must be at least as large as the local core")
    if not np.isfinite(core).all() or not np.isfinite(wing).all():
        raise ValueError("PSF arrays must contain only finite values")
    if not np.isfinite([splice_radius_arcsec, blend_width_arcsec,
                        wing_amplitude, pixel_scale_arcsec]).all():
        raise ValueError("hybrid parameters must be finite")
    if splice_radius_arcsec <= 0 or blend_width_arcsec <= 0 or wing_amplitude <= 0 or pixel_scale_arcsec <= 0:
        raise ValueError("splice radius, blend width, wing amplitude, and pixel scale must be positive")
    if not np.isclose(core.sum(), 1.0, rtol=0.0, atol=1e-3):
        raise ValueError("local core must retain its signed unit normalization")
    if not np.isclose(wing.sum(), 1.0, rtol=0.0, atol=2e-5):
        raise ValueError("public wing must retain its signed unit normalization")

    n = wing.shape[0]
    c = n // 2
    half_y, half_x = core.shape[0] // 2, core.shape[1] // 2
    core_padded = np.zeros_like(wing)
    core_padded[c-half_y:c+half_y+1, c-half_x:c+half_x+1] = core
    yy, xx = np.ogrid[-c:c+1, -c:c+1]
    radius_arcsec = np.hypot(xx, yy) * pixel_scale_arcsec
    inner = splice_radius_arcsec - blend_width_arcsec / 2.0
    outer = splice_radius_arcsec + blend_width_arcsec / 2.0
    if inner < 0 or outer > min(half_y, half_x) * pixel_scale_arcsec:
        raise ValueError("splice transition must fit inside the local-core radius")
    t = np.clip((radius_arcsec - inner) / blend_width_arcsec, 0.0, 1.0)
    smooth = t * t * (3.0 - 2.0 * t)
    core_weight = 1.0 - smooth
    hybrid = core_weight * core_padded + (1.0 - core_weight) * wing_amplitude * wing
    raw_sum = float(hybrid.sum(dtype=np.float64))
    if not np.isfinite(raw_sum) or abs(raw_sum) < 1e-12:
        raise ValueError("hybrid PSF has an invalid signed sum")
    hybrid /= raw_sum
    metadata = {
        "splice_radius_arcsec": float(splice_radius_arcsec),
        "blend_width_arcsec": float(blend_width_arcsec),
        "wing_amplitude_pre_normalization": float(wing_amplitude),
        "pixel_scale_arcsec": float(pixel_scale_arcsec),
        "local_core_shape": list(core.shape),
        "public_wing_shape": list(wing.shape),
        "local_core_signed_sum": float(core.sum()),
        "public_wing_signed_sum": float(wing.sum()),
        "hybrid_signed_sum_before_normalization": raw_sum,
        "hybrid_signed_sum": float(hybrid.sum(dtype=np.float64)),
        "negative_pixels_preserved": bool(np.any(hybrid < 0)),
        "status": "exploratory_core_wing_hybrid_sensitivity_only",
        "psf_approved": False,
        "science_measurement": False,
    }
    return hybrid.astype(np.float32), metadata


def satellite_halos_with_psfs(shape: tuple[int, int], centers_xy: np.ndarray,
                              kpc_per_pixel: float, amplitude_at_50: float,
                              scale_kpc: float, psfs: list[np.ndarray],
                              truncation_scale_lengths: float = 15.0) -> np.ndarray:
    """Build and locally convolve each satellite halo with its own PSF core."""
    centers = np.asarray(centers_xy, dtype=float).reshape((-1, 2))
    if len(centers) != len(psfs):
        raise ValueError("one local PSF is required for each satellite center")
    if kpc_per_pixel <= 0 or scale_kpc <= 0 or amplitude_at_50 < 0:
        raise ValueError("scales must be positive and amplitude non-negative")
    if truncation_scale_lengths < 5:
        raise ValueError("halo truncation must be at least five scale lengths")
    output = np.zeros(shape, dtype=np.float32)
    radius = int(np.ceil(truncation_scale_lengths * scale_kpc / kpc_per_pixel))
    axis = np.arange(-radius, radius + 1, dtype=np.float32)
    yy, xx = np.meshgrid(axis, axis, indexing="ij")
    for (x, y), psf in zip(centers, psfs):
        cx, cy = int(round(float(x))), int(round(float(y)))
        frac_x, frac_y = float(x) - cx, float(y) - cy
        shifted_radius = np.hypot(xx - frac_x, yy - frac_y) * kpc_per_pixel
        shifted_template = amplitude_at_50 * np.exp((50.0 - shifted_radius) / scale_kpc)
        convolved = convolve_surface_brightness(shifted_template, psf)
        x0, x1 = max(0, cx - radius), min(shape[1], cx + radius + 1)
        y0, y1 = max(0, cy - radius), min(shape[0], cy + radius + 1)
        if x0 >= x1 or y0 >= y1:
            continue
        sx0, sx1 = x0 - (cx - radius), x1 - (cx - radius)
        sy0, sy1 = y0 - (cy - radius), y1 - (cy - radius)
        output[y0:y1, x0:x1] += convolved[sy0:sy1, sx0:sx1]
    return output
