#!/usr/bin/env python3
"""In-tree detection and forced photometry for the desktop catalog tools.

Run ``python ogfmeas/sextract.py IMAGE`` or ``python -m ogfmeas.sextract IMAGE``.
X_IMAGE and Y_IMAGE are 1-based. THETA_IMAGE is degrees. FLUX_AUTO is the
elliptical sum at 2.5 times the Kron radius. MAG_ISOCOR is the magnitude of
the isophotal sum. CLASS_STAR is 0: this program does not apply a star/galaxy
classifier. It does not wrap the C extractor and it does not link SEP.

``catalog_table`` keeps the short column list used by the library test.
The command line writes the desktop column list, ``--info``, and forced
photometry (``--forced-catalog`` with ``--measure-image``).
"""
import argparse
import math
import pathlib
import sys

if __package__ in (None, ""):
    _root = str(pathlib.Path(__file__).resolve().parents[1])
    if _root not in sys.path:
        sys.path.insert(0, _root)

import numpy as np

from ogfmeas import Background, extract, flux_radius, kron_radius, sum_circle, sum_ellipse

_MAGERR = 2.5 / math.log(10.0)
_DETECTION_HEADER = (
    "NUMBER\tX_IMAGE\tY_IMAGE\tALPHA_J2000\tDELTA_J2000\t"
    "MAG_AUTO\tMAG_ISOCOR\tMAG_APER\tFLUX_AUTO\t"
    "FLUXERR_AUTO\tFLUXERR_APER\tMAGERR_AUTO\tMAGERR_APER\t"
    "FLUX_APER_2\tFLUX_APER_3\tFLUX_APER_5\t"
    "FLUXERR_APER_2\tFLUXERR_APER_3\tFLUXERR_APER_5\t"
    "MAG_APER_2\tMAG_APER_3\tMAG_APER_5\t"
    "FLUX_RADIUS\tA_IMAGE\tB_IMAGE\tTHETA_IMAGE\tELLIPTICITY\t"
    "KRON_RADIUS\tFWHM_IMAGE\tISO_RADIUS\tNPIX_ISO\t"
    "MU_MAX\tMU_THRESHOLD\tCLASS_STAR\tFLAGS"
)
_OUTSIDE = 4096


def load_fits(path):
    """Return ``(image, header)`` for the first 2-D image HDU."""
    from astropy.io import fits
    with fits.open(path) as hdul:
        for hdu in hdul:
            if hdu.data is None:
                continue
            data = np.squeeze(np.asarray(hdu.data))
            if data.ndim == 3:
                data = data[0]
            if data.ndim == 2:
                return np.asarray(data, dtype=np.float64), hdu.header
    raise SystemExit("no 2-D image in %s" % path)


def _load(path):
    data, _header = load_fits(path)
    return data


def _header_float(header, name):
    if header is None:
        return None
    try:
        return float(header[name])
    except Exception:
        return None


def _header_text(header, *names):
    if header is None:
        return ""
    for name in names:
        try:
            value = header[name]
        except Exception:
            continue
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def _cd_matrix(header):
    cd11 = _header_float(header, "CD1_1")
    if cd11 is None:
        cdelt1 = _header_float(header, "CDELT1") or 0.0
        cdelt2 = _header_float(header, "CDELT2") or 0.0
        pc11 = _header_float(header, "PC1_1")
        pc12 = _header_float(header, "PC1_2")
        pc21 = _header_float(header, "PC2_1")
        pc22 = _header_float(header, "PC2_2")
        if pc11 is None:
            pc11 = 1.0
        if pc22 is None:
            pc22 = 1.0
        pc12 = 0.0 if pc12 is None else pc12
        pc21 = 0.0 if pc21 is None else pc21
        return pc11 * cdelt1, pc12 * cdelt1, pc21 * cdelt2, pc22 * cdelt2
    return (
        cd11,
        _header_float(header, "CD1_2") or 0.0,
        _header_float(header, "CD2_1") or 0.0,
        _header_float(header, "CD2_2") or 0.0,
    )


def _pixel_scale_arcsec(header):
    """Plate scale in arcsec/pixel from the CD determinant, or None."""
    cd11, cd12, cd21, cd22 = _cd_matrix(header)
    det = cd11 * cd22 - cd12 * cd21
    if det == 0 or not math.isfinite(det):
        return None
    return math.sqrt(abs(det)) * 3600.0


def _ab_zeropoint(header):
    """AB zeropoint from PHOTFLAM/PHOTPLAM, PIXAR_SR, or PHOTFNU. ``(zp, pivot)``."""
    photflam = _header_float(header, "PHOTFLAM") or 0.0
    photplam = _header_float(header, "PHOTPLAM") or 0.0
    if photflam > 0 and photplam > 0:
        zp = -2.5 * math.log10(photflam) - 5.0 * math.log10(photplam) - 2.408
        return zp, photplam
    pixar = _header_float(header, "PIXAR_SR") or 0.0
    if pixar > 0:
        return -6.10 - 2.5 * math.log10(pixar), None
    photfnu = _header_float(header, "PHOTFNU") or 0.0
    if photfnu > 0:
        return -2.5 * math.log10(photfnu) + 8.90, None
    return None, None


def info_text(path):
    """KEY=VALUE block read by the band registry."""
    data, header = load_fits(path)
    ny, nx = data.shape
    lines = [
        "FILE=%s" % path,
        "NAXIS1=%d" % nx,
        "NAXIS2=%d" % ny,
        "FILTER=%s" % _header_text(header, "FILTER", "FILTER1"),
    ]
    pupil = _header_text(header, "PUPIL")
    if pupil:
        lines.append("PUPIL=%s" % pupil)
    lines.append("INSTRUME=%s" % _header_text(header, "INSTRUME"))
    lines.append("BUNIT=%s" % _header_text(header, "BUNIT"))
    zp, pivot = _ab_zeropoint(header)
    if zp is not None:
        lines.append("ZP_AB=%.4f" % zp)
    if pivot is not None:
        lines.append("PIVOT=%.1f" % pivot)
    exptime = _header_float(header, "EXPTIME")
    if exptime is not None and exptime > 0:
        lines.append("EXPTIME=%.1f" % exptime)
    crpix1 = _header_float(header, "CRPIX1")
    crpix2 = _header_float(header, "CRPIX2")
    crval1 = _header_float(header, "CRVAL1")
    crval2 = _header_float(header, "CRVAL2")
    scale = _pixel_scale_arcsec(header)
    valid = all(v is not None for v in (crpix1, crpix2, crval1, crval2)) and scale is not None and scale > 0
    if valid:
        cd11, cd12, cd21, cd22 = _cd_matrix(header)
        lines.append("PIXSCALE=%.5f" % scale)
        lines.append("CRPIX1=%.6f" % crpix1)
        lines.append("CRPIX2=%.6f" % crpix2)
        lines.append("CRVAL1=%.9f" % crval1)
        lines.append("CRVAL2=%.9f" % crval2)
        lines.append("CD1_1=%.12e" % cd11)
        lines.append("CD1_2=%.12e" % cd12)
        lines.append("CD2_1=%.12e" % cd21)
        lines.append("CD2_2=%.12e" % cd22)
    lines.append("WCSVALID=%d" % (1 if valid else 0))
    return "\n".join(lines) + "\n"


def convolution_kernel(name):
    """Detection kernel. ``None`` keeps the extractor's built-in 3x3."""
    name = (name or "default").strip().lower()
    if name in ("", "default"):
        return None
    if name == "gauss5x5":
        axis = np.arange(5.0) - 2.0
        yy, xx = np.meshgrid(axis, axis, indexing="ij")
        kernel = np.exp(-0.5 * (xx * xx + yy * yy))
        return kernel / kernel.sum()
    if name == "mexhat":
        axis = np.arange(5.0) - 2.0
        yy, xx = np.meshgrid(axis, axis, indexing="ij")
        radius2 = xx * xx + yy * yy
        kernel = (2.0 - radius2) * np.exp(-0.5 * radius2)
        return kernel - kernel.mean()
    if name == "tophat":
        return np.ones((5, 5), dtype=np.float64) / 25.0
    return None


def _subtract(data, mask, bw, fw):
    image = np.asarray(data, dtype=np.float64)
    finite = np.isfinite(image)
    bad = ~finite
    if mask is not None:
        bad = bad | np.asarray(mask, dtype=bool)
    work = np.where(bad, 0.0, image)
    mesh = max(int(bw), 1)
    filt = max(int(fw), 1)
    bkg = Background(work, mask=bad, bw=mesh, bh=mesh, fw=filt, fh=filt)
    return work - bkg.back(), bkg.rms(), bad


def _mag(flux, zp):
    flux = float(flux)
    if flux > 0.0 and math.isfinite(flux):
        return zp - 2.5 * math.log10(flux)
    return 99.0


def _magerr(flux, err):
    flux = float(flux)
    err = float(err)
    if flux > 0.0 and err > 0.0 and math.isfinite(flux) and math.isfinite(err):
        return _MAGERR * err / flux
    return 99.0


def _mu(value, zp, pixel_scale):
    scale = float(pixel_scale) if pixel_scale and pixel_scale > 0 else 1.0
    per_pixel = float(value) / (scale * scale)
    return _mag(per_pixel, zp)


def _sky(header, x_one, y_one):
    if header is None or len(x_one) == 0:
        return np.zeros(len(x_one)), np.zeros(len(y_one))
    if _header_float(header, "CRVAL1") is None or _pixel_scale_arcsec(header) is None:
        return np.zeros(len(x_one)), np.zeros(len(y_one))
    try:
        from astropy.wcs import WCS
        world = WCS(header)
        if not world.has_celestial:
            return np.zeros(len(x_one)), np.zeros(len(y_one))
        ra, dec = world.wcs_pix2world(np.asarray(x_one, dtype=np.float64), np.asarray(y_one, dtype=np.float64), 1)
        return np.asarray(ra, dtype=np.float64), np.asarray(dec, dtype=np.float64)
    except Exception:
        return np.zeros(len(x_one)), np.zeros(len(y_one))


def detect_sources(data, thresh=1.5, minarea=5, bw=64, fw=3, deblend_nthresh=32,
                   deblend_cont=0.005, zp=25.0, gain=0.0, pixel_scale=1.0,
                   phot_aperture=5.0, phot_aperture_2=4.0, phot_aperture_3=6.0,
                   phot_aperture_5=10.0, conv="default", header=None, mask=None,
                   seeing_fwhm=3.0):
    """Detect on a background-subtracted image and measure the desktop columns.

    ``seeing_fwhm`` is accepted from the extract panel. CLASS_STAR stays 0.
    Returned ``x`` and ``y`` are 0-based.
    """
    del seeing_fwhm
    sub, rms, bad = _subtract(data, mask, bw, fw)
    kernel = convolution_kernel(conv)
    options = {}
    if kernel is not None:
        options["filter_kernel"] = kernel
    gain_arg = float(gain) if gain and float(gain) > 0 else None
    objects = extract(
        sub, thresh, err=rms, minarea=minarea, mask=bad,
        deblend_nthresh=deblend_nthresh, deblend_cont=deblend_cont, **options)
    n = len(objects)
    empty = {
        "x": np.zeros(0), "y": np.zeros(0), "a": np.zeros(0), "b": np.zeros(0),
        "mag_auto": np.zeros(0),
    }
    if n == 0:
        return empty
    x = np.asarray(objects["x"], dtype=np.float64)
    y = np.asarray(objects["y"], dtype=np.float64)
    a = np.maximum(np.asarray(objects["a"], dtype=np.float64), 0.3)
    b = np.clip(np.asarray(objects["b"], dtype=np.float64), 0.3, a)
    theta = np.asarray(objects["theta"], dtype=np.float64)
    kron, _kron_flag = kron_radius(sub, x, y, a, b, theta, 6.0, mask=bad)
    kron = np.where(np.isfinite(kron) & (kron > 0), kron, 2.5)
    flux_auto, err_auto, flag_auto = sum_ellipse(
        sub, x, y, a, b, theta, 2.5 * kron, err=rms, gain=gain_arg, mask=bad, subpix=5)
    half, _half_flag = flux_radius(sub, x, y, np.maximum(6.0 * a, 4.0), 0.5, mask=bad, subpix=1)
    diameters = (phot_aperture, phot_aperture_2, phot_aperture_3, phot_aperture_5)
    circles = []
    for diameter in diameters:
        radius = np.full(n, 0.5 * float(diameter))
        circles.append(sum_circle(sub, x, y, radius, err=rms, gain=gain_arg, mask=bad, subpix=5))
    flux_iso = np.asarray(objects["flux"], dtype=np.float64)
    peak = np.asarray(objects["peak"], dtype=np.float64)
    npix = np.asarray(objects["npix"], dtype=np.float64)
    flags = np.asarray(objects["flag"], dtype=np.int32) | flag_auto | circles[0][2]
    order = np.argsort(np.where(flux_auto > 0, zp - 2.5 * np.log10(np.maximum(flux_auto, 1e-30)), 99.0), kind="mergesort")
    alpha, delta = _sky(header, x[order] + 1.0, y[order] + 1.0)
    scale = float(pixel_scale) if pixel_scale and float(pixel_scale) > 0 else 1.0
    thresh_level = float(thresh) * np.asarray(rms[np.clip(np.rint(y).astype(int), 0, rms.shape[0] - 1),
                                                   np.clip(np.rint(x).astype(int), 0, rms.shape[1] - 1)])
    rows = {
        "x": x[order],
        "y": y[order],
        "a": a[order],
        "b": b[order],
        "alpha": alpha,
        "delta": delta,
        "mag_auto": np.array([_mag(flux_auto[i], zp) for i in order]),
        "mag_isocor": np.array([_mag(flux_iso[i], zp) for i in order]),
        "mag_aper": np.array([_mag(circles[0][0][i], zp) for i in order]),
        "flux_auto": flux_auto[order],
        "fluxerr_auto": err_auto[order],
        "fluxerr_aper": circles[0][1][order],
        "magerr_auto": np.array([_magerr(flux_auto[i], err_auto[i]) for i in order]),
        "magerr_aper": np.array([_magerr(circles[0][0][i], circles[0][1][i]) for i in order]),
        "kron": kron[order],
        "half": half[order],
        "theta": theta[order],
        "npix": npix[order],
        "mu_max": np.array([_mu(peak[i], zp, scale) for i in order]),
        "mu_threshold": np.array([_mu(thresh_level[i], zp, scale) for i in order]),
        "flags": flags[order],
    }
    for slot, (flux, err, _flag) in zip((2, 3, 5), circles[1:]):
        rows["flux_aper_%d" % slot] = flux[order]
        rows["fluxerr_aper_%d" % slot] = err[order]
        rows["mag_aper_%d" % slot] = np.array([_mag(flux[i], zp) for i in order])
    return rows


def _fmt(value, spec):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "nan"
    if not math.isfinite(number):
        return "nan"
    return spec % number


def format_detection(rows):
    lines = [_DETECTION_HEADER]
    n = len(rows.get("x", ()))
    for i in range(n):
        a = float(rows["a"][i])
        b = float(rows["b"][i])
        ellipticity = 0.0 if not (a > 0 and math.isfinite(a) and math.isfinite(b)) else max(0.0, 1.0 - b / a)
        fwhm = 2.0 * math.sqrt(2.0 * math.log(2.0)) * math.sqrt(max(a * b, 0.0))
        iso = math.sqrt(max(float(rows["npix"][i]), 0.0) / math.pi)
        lines.append("\t".join([
            str(i + 1),
            _fmt(rows["x"][i] + 1.0, "%.4f"),
            _fmt(rows["y"][i] + 1.0, "%.4f"),
            _fmt(rows["alpha"][i], "%.6f"),
            _fmt(rows["delta"][i], "%.6f"),
            _fmt(rows["mag_auto"][i], "%.4f"),
            _fmt(rows["mag_isocor"][i], "%.4f"),
            _fmt(rows["mag_aper"][i], "%.4f"),
            _fmt(rows["flux_auto"][i], "%.6g"),
            _fmt(rows["fluxerr_auto"][i], "%.6g"),
            _fmt(rows["fluxerr_aper"][i], "%.6g"),
            _fmt(rows["magerr_auto"][i], "%.4f"),
            _fmt(rows["magerr_aper"][i], "%.4f"),
            _fmt(rows["flux_aper_2"][i], "%.6g"),
            _fmt(rows["flux_aper_3"][i], "%.6g"),
            _fmt(rows["flux_aper_5"][i], "%.6g"),
            _fmt(rows["fluxerr_aper_2"][i], "%.6g"),
            _fmt(rows["fluxerr_aper_3"][i], "%.6g"),
            _fmt(rows["fluxerr_aper_5"][i], "%.6g"),
            _fmt(rows["mag_aper_2"][i], "%.4f"),
            _fmt(rows["mag_aper_3"][i], "%.4f"),
            _fmt(rows["mag_aper_5"][i], "%.4f"),
            _fmt(rows["half"][i], "%.4f"),
            _fmt(a, "%.4f"),
            _fmt(b, "%.4f"),
            _fmt(math.degrees(float(rows["theta"][i])), "%.4f"),
            _fmt(ellipticity, "%.4f"),
            _fmt(rows["kron"][i], "%.4f"),
            _fmt(fwhm, "%.4f"),
            _fmt(iso, "%.4f"),
            str(int(rows["npix"][i])),
            _fmt(rows["mu_max"][i], "%.4f"),
            _fmt(rows["mu_threshold"][i], "%.4f"),
            "0.000",
            str(int(rows["flags"][i])),
        ]))
    return "\n".join(lines) + "\n"


def detection_text(data, header=None, **options):
    return format_detection(detect_sources(data, header=header, **options))


def catalog_table(data, thresh=1.5, minarea=5, bw=64, deblend_nthresh=32, deblend_cont=0.005):
    """Short catalog. Column 2 is X_IMAGE. Positions are 1-based."""
    image = np.asarray(data, dtype=np.float64)
    finite = np.isfinite(image)
    work = np.where(finite, image, 0.0)
    bkg = Background(work, mask=~finite, bw=bw, bh=bw)
    sub = work - bkg.back()
    rms = bkg.rms()
    objects = extract(
        sub, thresh, err=rms, minarea=minarea, mask=~finite,
        deblend_nthresh=deblend_nthresh, deblend_cont=deblend_cont)
    n = len(objects)
    lines = ["NUMBER\tX_IMAGE\tY_IMAGE\tA_IMAGE\tB_IMAGE\tTHETA_IMAGE\tFLUX_ISO\tFLUX_AUTO\tKRON_RADIUS\tFLUX_RADIUS\tFLAGS"]
    if n == 0:
        return "\n".join(lines) + "\n"
    a = np.maximum(objects["a"], 0.3)
    b = np.clip(objects["b"], 0.3, a)
    theta = objects["theta"]
    kron, _ = kron_radius(sub, objects["x"], objects["y"], a, b, theta, 6.0, mask=~finite)
    kron = np.where(np.isfinite(kron) & (kron > 0), kron, 2.5)
    flux_auto, _, auto_flag = sum_ellipse(
        sub, objects["x"], objects["y"], a, b, theta, 2.5 * kron, err=rms, mask=~finite, subpix=5)
    half, _ = flux_radius(sub, objects["x"], objects["y"], np.maximum(6.0 * a, 4.0), 0.5, mask=~finite)
    for i, obj in enumerate(objects):
        lines.append("\t".join([
            str(i + 1),
            "%.4f" % (obj["x"] + 1.0),
            "%.4f" % (obj["y"] + 1.0),
            "%.4f" % obj["a"],
            "%.4f" % obj["b"],
            "%.4f" % (np.degrees(obj["theta"])),
            "%.6g" % obj["flux"],
            "%.6g" % flux_auto[i],
            "%.4f" % kron[i],
            "%.4f" % half[i],
            str(int(obj["flag"]) | int(auto_flag[i])),
        ]))
    return "\n".join(lines) + "\n"


def _read_catalog(path):
    rows = []
    header = None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = [item.strip() for item in fields]
            else:
                rows.append(fields)
    if not header:
        raise ValueError("empty catalog")
    index = {name: i for i, name in enumerate(header)}
    for name in ("NUMBER", "X_IMAGE", "Y_IMAGE"):
        if name not in index:
            raise ValueError("catalog needs NUMBER, X_IMAGE, Y_IMAGE columns")
    return index, rows


def _cell(row, index, name, default):
    slot = index.get(name)
    if slot is None or slot >= len(row) or row[slot].strip() == "":
        return default
    try:
        return float(row[slot])
    except ValueError:
        return default


def _wcs_of(header):
    if header is None or _header_float(header, "CRVAL1") is None or _pixel_scale_arcsec(header) is None:
        return None
    try:
        from astropy.wcs import WCS
        world = WCS(header)
        if world.has_celestial:
            return world
    except Exception:
        return None
    return None


def _map_pixel(w_from, w_to, x_one, y_one):
    """1-based detection pixel -> 0-based measurement pixel."""
    sky = w_from.pixel_to_world(x_one - 1.0, y_one - 1.0)
    mapped = w_to.world_to_pixel(sky)
    return float(np.asarray(mapped[0]).reshape(-1)[0]), float(np.asarray(mapped[1]).reshape(-1)[0])


def _transformed_ellipse(a, b, theta, jacobian):
    j00, j01, j10, j11 = jacobian
    cosine, sine = math.cos(theta), math.sin(theta)
    s11 = a * a * cosine * cosine + b * b * sine * sine
    s22 = a * a * sine * sine + b * b * cosine * cosine
    s12 = (a * a - b * b) * cosine * sine
    t11 = j00 * (j00 * s11 + j01 * s12) + j01 * (j00 * s12 + j01 * s22)
    t12 = j00 * (j10 * s11 + j11 * s12) + j01 * (j10 * s12 + j11 * s22)
    t22 = j10 * (j10 * s11 + j11 * s12) + j11 * (j10 * s12 + j11 * s22)
    trace = 0.5 * (t11 + t22)
    root = math.sqrt(max(0.0, 0.25 * (t11 - t22) * (t11 - t22) + t12 * t12))
    major = math.sqrt(max(trace + root, 1e-8))
    minor = math.sqrt(max(trace - root, 1e-8))
    if minor > major:
        major, minor = minor, major
    position = 0.5 * math.atan2(2.0 * t12, t11 - t22)
    scale = math.sqrt(max(abs(j00 * j11 - j01 * j10), 1e-12))
    return major, minor, position, scale


def forced_text(detection_path, catalog_path, measure_path, band="B", zp=25.0, snr_min=1.0,
                phot_aperture=5.0, bw=64, fw=3, gain=0.0):
    """Measure catalog positions on another image. Magnitudes below the S/N floor are 99."""
    det_image, det_header = load_fits(detection_path)
    meas_image, meas_header = load_fits(measure_path)
    if gain is None or float(gain) <= 0:
        header_gain = _header_float(meas_header, "GAIN") or _header_float(meas_header, "EGAIN")
        gain = header_gain if header_gain and header_gain > 0 else 0.0
    sub, rms, bad = _subtract(meas_image, None, bw, fw)
    gain_arg = float(gain) if float(gain) > 0 else None
    index, rows = _read_catalog(catalog_path)
    w_det = _wcs_of(det_header)
    w_meas = _wcs_of(meas_header)
    same = det_image.shape == meas_image.shape
    if not same and (w_det is None or w_meas is None):
        raise ValueError("images differ in size and lack WCS")
    ny, nx = sub.shape
    prefix = str(band)
    header = (
        "NUMBER\tX_%s\tY_%s\tSCALE_%s\tFLUX_AUTO_%s\tFLUXERR_AUTO_%s\tMAG_AUTO_%s\tMAGERR_AUTO_%s\t"
        "FLUX_APER_%s\tFLUXERR_APER_%s\tMAG_APER_%s\tMAGERR_APER_%s\tFLAGS_%s"
        % (prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix, prefix)
    )
    lines = [header]
    for row in rows:
        number = int(_cell(row, index, "NUMBER", 0))
        x_one = _cell(row, index, "X_IMAGE", float("nan"))
        y_one = _cell(row, index, "Y_IMAGE", float("nan"))
        if not math.isfinite(x_one) or not math.isfinite(y_one):
            continue
        semi_a = _cell(row, index, "A_IMAGE", 0.0)
        semi_b = _cell(row, index, "B_IMAGE", 0.0)
        theta = math.radians(_cell(row, index, "THETA_IMAGE", 0.0))
        kron = _cell(row, index, "KRON_RADIUS", 2.5)
        if not math.isfinite(kron) or kron <= 0:
            kron = 2.5
        flags = 0
        if same or w_det is None or w_meas is None:
            xm, ym = x_one - 1.0, y_one - 1.0
            ap, bp, tp, scale = max(semi_a, 0.3), max(semi_b, 0.3), theta, 1.0
        else:
            try:
                xm, ym = _map_pixel(w_det, w_meas, x_one, y_one)
                x1, y1 = _map_pixel(w_det, w_meas, x_one + 1.0, y_one)
                x2, y2 = _map_pixel(w_det, w_meas, x_one, y_one + 1.0)
                jacobian = (x1 - xm, x2 - xm, y1 - ym, y2 - ym)
                if semi_a > 0 and semi_b > 0:
                    ap, bp, tp, scale = _transformed_ellipse(semi_a, semi_b, theta, jacobian)
                else:
                    ap, bp, tp, scale = 2.0, 2.0, 0.0, math.sqrt(max(abs(jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2]), 1.0))
            except Exception:
                xm = ym = float("nan")
                ap = bp = tp = scale = 1.0
                flags |= _OUTSIDE
        inside = math.isfinite(xm) and math.isfinite(ym) and 0.0 <= xm < nx and 0.0 <= ym < ny
        if not inside:
            flags |= _OUTSIDE
            flux_auto = err_auto = flux_ap = err_ap = 0.0
        else:
            if semi_a <= 0 or semi_b <= 0:
                ap, bp, tp = 2.0, 2.0, 0.0
                flags |= 0x100
            flux_auto, err_auto, flag_auto = sum_ellipse(
                sub, xm, ym, max(ap, 0.3), max(min(bp, ap), 0.3), tp, 2.5 * kron,
                err=rms, gain=gain_arg, mask=bad, subpix=5)
            flux_ap, err_ap, flag_ap = sum_circle(
                sub, xm, ym, 0.5 * float(phot_aperture) * scale,
                err=rms, gain=gain_arg, mask=bad, subpix=5)
            flux_auto, err_auto = float(np.asarray(flux_auto).reshape(-1)[0]), float(np.asarray(err_auto).reshape(-1)[0])
            flux_ap, err_ap = float(np.asarray(flux_ap).reshape(-1)[0]), float(np.asarray(err_ap).reshape(-1)[0])
            flags |= int(np.asarray(flag_auto).reshape(-1)[0]) | int(np.asarray(flag_ap).reshape(-1)[0])
        detected_auto = inside and flux_auto > 0 and (err_auto <= 0 or flux_auto >= float(snr_min) * err_auto)
        detected_ap = inside and flux_ap > 0 and (err_ap <= 0 or flux_ap >= float(snr_min) * err_ap)
        mag_auto = _mag(flux_auto, zp) if detected_auto else 99.0
        mer_auto = _magerr(flux_auto, err_auto) if detected_auto else 99.0
        mag_ap = _mag(flux_ap, zp) if detected_ap else 99.0
        mer_ap = _magerr(flux_ap, err_ap) if detected_ap else 99.0
        lines.append("\t".join([
            str(number),
            _fmt(xm + 1.0 if math.isfinite(xm) else float("nan"), "%.4f"),
            _fmt(ym + 1.0 if math.isfinite(ym) else float("nan"), "%.4f"),
            _fmt(scale, "%.4f"),
            _fmt(flux_auto, "%.6g"),
            _fmt(err_auto, "%.6g"),
            _fmt(mag_auto, "%.4f"),
            _fmt(mer_auto, "%.4f"),
            _fmt(flux_ap, "%.6g"),
            _fmt(err_ap, "%.6g"),
            _fmt(mag_ap, "%.4f"),
            _fmt(mer_ap, "%.4f"),
            str(int(flags)),
        ]))
    return "\n".join(lines) + "\n"


def _parser():
    parser = argparse.ArgumentParser(description="Detect sources and write a TSV catalog")
    parser.add_argument("image", nargs="?")
    parser.add_argument("--info", metavar="FITS")
    parser.add_argument("--out", default="-")
    parser.add_argument("--thresh", type=float, default=None)
    parser.add_argument("--detect-thresh", type=float, default=None)
    parser.add_argument("--minarea", type=int, default=None)
    parser.add_argument("--detect-minarea", type=int, default=None)
    parser.add_argument("--bw", type=int, default=None)
    parser.add_argument("--back-size", type=int, default=None)
    parser.add_argument("--back-filtersize", type=int, default=3)
    parser.add_argument("--deblend-nthresh", type=int, default=32)
    parser.add_argument("--deblend-cont", type=float, default=None)
    parser.add_argument("--deblend-mincont", type=float, default=None)
    parser.add_argument("--phot-aperture", type=float, default=5.0)
    parser.add_argument("--phot-aperture-2", type=float, default=4.0)
    parser.add_argument("--phot-aperture-3", type=float, default=6.0)
    parser.add_argument("--phot-aperture-5", type=float, default=10.0)
    parser.add_argument("--mag-zeropoint", type=float, default=25.0)
    parser.add_argument("--gain", type=float, default=0.0)
    parser.add_argument("--pixel-scale", type=float, default=1.0)
    parser.add_argument("--seeing-fwhm", type=float, default=3.0)
    parser.add_argument("--conv-filter", default="default")
    parser.add_argument("--forced-catalog")
    parser.add_argument("--measure-image")
    parser.add_argument("--band", default="B")
    parser.add_argument("--snr-min", type=float, default=1.0)
    return parser


def _first(primary, secondary, default):
    if primary is not None:
        return primary
    if secondary is not None:
        return secondary
    return default


def _resolve_scale_and_gain(header, pixel_scale, gain):
    scale = float(pixel_scale)
    if scale == 1.0:
        measured = _pixel_scale_arcsec(header)
        if measured is not None and measured > 0:
            scale = measured
    if gain is None or float(gain) <= 0:
        header_gain = _header_float(header, "GAIN")
        if header_gain is None or header_gain <= 0:
            header_gain = _header_float(header, "EGAIN")
        gain = header_gain if header_gain and header_gain > 0 else 0.0
    return scale, float(gain)


def main(argv=None):
    args = _parser().parse_args(argv)
    try:
        if args.info:
            text = info_text(args.info)
        elif not args.image:
            sys.stderr.write("usage: sextract.py IMAGE [options] | sextract.py --info IMAGE\n")
            return 1
        elif args.forced_catalog or args.measure_image:
            if not args.forced_catalog or not args.measure_image:
                sys.stderr.write("ERROR: --forced-catalog and --measure-image go together\n")
                return 1
            text = forced_text(
                args.image, args.forced_catalog, args.measure_image, band=args.band,
                zp=args.mag_zeropoint, snr_min=args.snr_min, phot_aperture=args.phot_aperture,
                bw=_first(args.bw, args.back_size, 64), fw=args.back_filtersize, gain=args.gain)
        else:
            _image, header = load_fits(args.image)
            scale, gain = _resolve_scale_and_gain(header, args.pixel_scale, args.gain)
            text = detection_text(
                _image, header=header,
                thresh=_first(args.thresh, args.detect_thresh, 1.5),
                minarea=_first(args.minarea, args.detect_minarea, 5),
                bw=_first(args.bw, args.back_size, 64),
                fw=args.back_filtersize,
                deblend_nthresh=args.deblend_nthresh,
                deblend_cont=_first(args.deblend_cont, args.deblend_mincont, 0.005),
                zp=args.mag_zeropoint, gain=gain, pixel_scale=scale,
                phot_aperture=args.phot_aperture, phot_aperture_2=args.phot_aperture_2,
                phot_aperture_3=args.phot_aperture_3, phot_aperture_5=args.phot_aperture_5,
                conv=args.conv_filter, seeing_fwhm=args.seeing_fwhm)
    except (OSError, ValueError) as exc:
        sys.stderr.write("ERROR: %s\n" % exc)
        return 1
    if args.out == "-":
        sys.stdout.write(text)
    else:
        with open(args.out, "w", encoding="utf-8") as handle:
            handle.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
