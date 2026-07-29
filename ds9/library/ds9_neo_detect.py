#!/usr/bin/env python3
"""Detect and link moving-source candidates in a FITS frame sequence."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy.ndimage import maximum_filter


@dataclass
class Detection:
    frame: int
    path: Path
    x: float
    y: float
    ra_deg: float
    dec_deg: float
    jd_tdb: float
    snr: float


def _detections(path: Path, frame: int, threshold_sigma: float) -> list[Detection]:
    with fits.open(path) as hdul:
        data = np.asarray(hdul[0].data, dtype=float)
        header = hdul[0].header
    data = np.squeeze(data)
    finite = np.isfinite(data)
    if data.ndim != 2 or not finite.any():
        return []
    background = float(np.nanmedian(data[finite]))
    mad = float(np.nanmedian(np.abs(data[finite] - background)))
    sigma = max(1.4826 * mad, float(np.nanstd(data[finite])) * 0.05, 1.0e-12)
    peaks = maximum_filter(data, size=5, mode="nearest") == data
    peaks &= data > background + threshold_sigma * sigma
    yy, xx = np.nonzero(peaks)
    wcs = WCS(header)
    date_obs = header.get("DATE-OBS")
    if date_obs is None:
        raise ValueError(f"{path} has no DATE-OBS header.")
    epoch = Time(str(date_obs), scale="utc").tdb.jd
    output = []
    for x, y in zip(xx, yy, strict=True):
        ra, dec = wcs.pixel_to_world_values(float(x), float(y))
        output.append(Detection(frame, path, float(x), float(y), float(ra), float(dec), epoch, float((data[y, x] - background) / sigma)))
    return output


def _link(detections: list[Detection], max_step: float) -> list[list[Detection]]:
    by_frame: dict[int, list[Detection]] = {}
    for item in detections:
        by_frame.setdefault(item.frame, []).append(item)
    tracks: list[list[Detection]] = []
    for frame in sorted(by_frame):
        for detection in by_frame[frame]:
            best = None
            best_distance = float("inf")
            for track in tracks:
                if track[-1].frame != frame - 1:
                    continue
                distance = float(np.hypot(detection.x - track[-1].x, detection.y - track[-1].y))
                if distance < best_distance and distance <= max_step:
                    best = track
                    best_distance = distance
            if best is None:
                tracks.append([detection])
            else:
                best.append(detection)
    return tracks


def _tracklet_row(track_id: int, track: list[Detection]) -> dict[str, str | int]:
    ordered = sorted(track, key=lambda item: item.jd_tdb)
    jd = np.asarray([item.jd_tdb for item in ordered], dtype=float)
    ra_rad = np.unwrap(np.deg2rad([item.ra_deg for item in ordered]))
    dec_deg = np.asarray([item.dec_deg for item in ordered], dtype=float)
    jd_mid = float(np.median(jd))
    ra_mid_rad = float(np.median(ra_rad))
    dec_mid_deg = float(np.median(dec_deg))
    hours = (jd - jd_mid) * 24.0
    east_arcsec = (ra_rad - ra_mid_rad) * np.cos(np.deg2rad(dec_mid_deg)) * 206264.806247
    north_arcsec = (dec_deg - dec_mid_deg) * 3600.0
    design = np.column_stack((hours, np.ones_like(hours)))
    east_rate, east_offset = np.linalg.lstsq(design, east_arcsec, rcond=None)[0]
    north_rate, north_offset = np.linalg.lstsq(design, north_arcsec, rcond=None)[0]
    east_residual = east_arcsec - (east_rate * hours + east_offset)
    north_residual = north_arcsec - (north_rate * hours + north_offset)
    rms_arcsec = float(np.sqrt(np.mean(east_residual**2 + north_residual**2)))
    rate_arcsec_hour = float(np.hypot(east_rate, north_rate))
    position_angle_deg = float(np.degrees(np.arctan2(east_rate, north_rate)) % 360.0)
    ra_mid_deg = float(np.degrees(ra_mid_rad) % 360.0)
    utc_mid = Time(jd_mid, format="jd", scale="tdb").utc.isot
    return {
        "ID": f"NEO-{track_id:02d}",
        "NDET": len(ordered),
        "RA_DEG": f"{ra_mid_deg:.6f}",
        "DEC_DEG": f"{dec_mid_deg:.6f}",
        "UTC_MID": utc_mid,
        "RATE_AS_H": f"{rate_arcsec_hour:.2f}",
        "PA_DEG": f"{position_angle_deg:.1f}",
        "SNR_MED": f"{np.median([item.snr for item in ordered]):.1f}",
        "RMS_AS": f"{rms_arcsec:.2f}",
        "STATUS": "REVIEW",
    }


def _write_tracklet_catalog(tracks: list[list[Detection]], path: Path) -> None:
    fieldnames = [
        "ID",
        "NDET",
        "RA_DEG",
        "DEC_DEG",
        "UTC_MID",
        "RATE_AS_H",
        "PA_DEG",
        "SNR_MED",
        "RMS_AS",
        "STATUS",
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for track_id, track in enumerate(tracks, start=1):
            writer.writerow(_tracklet_row(track_id, track))


def run(paths: list[Path], output_dir: Path, threshold_sigma: float, max_step: float, min_track_length: int) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    all_detections = []
    for frame, path in enumerate(paths, start=1):
        all_detections.extend(_detections(path, frame, threshold_sigma))
    tracks = [track for track in _link(all_detections, max_step) if len(track) >= min_track_length]
    csv_path = output_dir / "neo_candidates.csv"
    tracklet_path = output_dir / "neo_tracklets.tsv"
    region_path = output_dir / "neo_candidates.reg"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["track_id", "frame", "obstime_utc", "ra_deg", "dec_deg", "sigma_arcsec", "observer", "snr"])
        writer.writeheader()
        for track_id, track in enumerate(tracks, start=1):
            for item in track:
                writer.writerow({"track_id": track_id, "frame": item.frame, "obstime_utc": Time(item.jd_tdb, format="jd", scale="tdb").utc.isot + "Z", "ra_deg": item.ra_deg, "dec_deg": item.dec_deg, "sigma_arcsec": 0.5, "observer": "EARTH", "snr": item.snr})
    _write_tracklet_catalog(tracks, tracklet_path)
    with region_path.open("w", encoding="utf-8") as stream:
        stream.write("# Region file format: DS9 version 4.1\n")
        stream.write('global color=cyan width=3 font="helvetica 10 bold"\n')
        stream.write("fk5\n")
        for track_id, track in enumerate(tracks, start=1):
            for item in track:
                stream.write(f"point({item.ra_deg},{item.dec_deg}) # point=circle text={{NEO {track_id} f{item.frame}}}\n")
    frame_regions = []
    for frame in range(1, len(paths) + 1):
        frame_region_path = output_dir / f"neo_candidates_f{frame}.reg"
        with frame_region_path.open("w", encoding="utf-8") as stream:
            stream.write("# Region file format: DS9 version 4.1\n")
            stream.write('global color=cyan width=3 font="helvetica 10 bold"\n')
            stream.write("fk5\n")
            for track_id, track in enumerate(tracks, start=1):
                for item in track:
                    if item.frame == frame:
                        stream.write(f"point({item.ra_deg},{item.dec_deg}) # point=circle text={{NEO {track_id}}}\n")
        frame_regions.append(str(frame_region_path.resolve()))
    summary = {"frames": len(paths), "detections": len(all_detections), "tracks": len(tracks), "csv": str(csv_path.resolve()), "tracklets": str(tracklet_path.resolve()), "regions": str(region_path.resolve()), "frame_regions": frame_regions, "status": "candidate tracks ready for CODES/OpenOrb"}
    (output_dir / "neo_detection_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("frames", nargs="+", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("~/.ds9/neo_codes/detection").expanduser())
    parser.add_argument("--threshold-sigma", type=float, default=4.0)
    parser.add_argument("--max-step-pixel", type=float, default=100.0)
    parser.add_argument("--min-track-length", type=int, default=3)
    args = parser.parse_args()
    print(json.dumps(run([path.expanduser().resolve() for path in args.frames], args.output_dir.expanduser().resolve(), args.threshold_sigma, args.max_step_pixel, args.min_track_length), indent=2))


if __name__ == "__main__":
    main()
