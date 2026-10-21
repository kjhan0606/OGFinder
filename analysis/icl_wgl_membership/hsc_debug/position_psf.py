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
