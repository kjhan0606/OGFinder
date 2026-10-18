"""Names for JWST and HST simulation PSFs.

The in-tree model is a circular aperture with a central obscuration, drawn
at one wavelength and smoothed by a Gaussian jitter. Its header method is
``ogfmeas``. WebbPSF is imported only when a JWST simulation PSF is requested
and the user-installed package imports. TinyTim is a separate program:
tiny1, tiny2, and tiny3 are used only when they are already on PATH. This
module does not ship those programs and does not label the circular pupil
with their names.
"""

import os
import re
import shutil
import sys
import tempfile

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.special import j1

# arcsec per radian
_ARCSEC_PER_RAD = 180.0 * 3600.0 / np.pi
_FILTER_RE = re.compile(r"^F(\d+)([A-Z]+)?$")


# ---- Availability checks ----

def check_webbpsf_available():
    """True when the user-installed webbpsf package imports."""
    try:
        import webbpsf  # noqa: F401
        return True
    except Exception:
        return False


def check_tinytim_available():
    """True when the user-installed TinyTim programs are on PATH."""
    return (shutil.which('tiny1') is not None and
            shutil.which('tiny2') is not None and
            shutil.which('tiny3') is not None)


# ---- FITS header detection ----

def detect_telescope_instrument(header):
    """Detect telescope, instrument, and filter from FITS header keywords.

    Parameters
    ----------
    header : astropy.io.fits.Header
        FITS header.

    Returns
    -------
    info : dict
        Keys: telescope, instrument, filter, detector, pixel_scale.
        Values are uppercase strings or 'unknown'.
    """
    info = {
        'telescope': 'unknown',
        'instrument': 'unknown',
        'filter': 'unknown',
        'detector': 'unknown',
        'pixel_scale': None,
    }

    # Telescope
    telescop = str(header.get('TELESCOP', '')).strip().upper()
    if 'JWST' in telescop:
        info['telescope'] = 'JWST'
    elif 'HST' in telescop or 'HUBBLE' in telescop:
        info['telescope'] = 'HST'

    # Instrument
    instrume = str(header.get('INSTRUME', '')).strip().upper()
    if instrume:
        info['instrument'] = instrume

    # Filter — try FILTER, FILTER1, FILTER2
    filt = ''
    for key in ('FILTER', 'FILTER1', 'FILTER2'):
        val = str(header.get(key, '')).strip().upper()
        if val and val != 'CLEAR' and val != 'NONE' and not val.startswith('N/A'):
            filt = val
            break
    if filt:
        info['filter'] = filt

    # Detector
    detector = str(header.get('DETECTOR', '')).strip().upper()
    if detector:
        info['detector'] = detector

    # Pixel scale from CD matrix or CDELT
    try:
        if 'CD1_1' in header:
            info['pixel_scale'] = abs(float(header['CD1_1'])) * 3600.0  # deg -> arcsec
        elif 'CDELT1' in header:
            info['pixel_scale'] = abs(float(header['CDELT1'])) * 3600.0
    except (ValueError, TypeError):
        pass

    return info


# ---- WebbPSF generation ----

WEBBPSF_INSTRUMENTS = {
    'NIRCAM': [
        'F070W', 'F090W', 'F115W', 'F140M', 'F150W', 'F162M', 'F164N',
        'F150W2', 'F182M', 'F187N', 'F200W', 'F210M', 'F212N',
        'F250M', 'F277W', 'F300M', 'F322W2', 'F323N', 'F335M',
        'F356W', 'F360M', 'F405N', 'F410M', 'F430M', 'F444W',
        'F460M', 'F466N', 'F470N', 'F480M',
    ],
    'MIRI': [
        'F560W', 'F770W', 'F1000W', 'F1065C', 'F1130W', 'F1140C',
        'F1280W', 'F1500W', 'F1550C', 'F1800W', 'F2100W', 'F2300C',
        'F2550W',
    ],
    'NIRISS': [
        'F090W', 'F115W', 'F140M', 'F150W', 'F158M', 'F200W',
        'F277W', 'F356W', 'F380M', 'F430M', 'F444W', 'F480M',
    ],
    'NIRSPEC': [
        'F070LP', 'F100LP', 'F170LP', 'F290LP', 'CLEAR',
    ],
    'FGS': [
        'FGS',
    ],
}


def generate_webbpsf(instrument, filter_name, psf_size=201,
                     oversample=1, jitter_sigma=0.007,
                     focus_offset=0.0, output=None):
    """Generate a JWST PSF with the user-installed WebbPSF package.

    A missing package raises and does not write a FITS file. There is no
    in-tree instrument PSF to substitute.
    """
    if not check_webbpsf_available():
        raise RuntimeError(
            "webbpsf was not found. the external package is not called. "
            "An instrument PSF is not produced here."
        )
    import webbpsf

    inst_name = instrument.upper()
    inst_cls = getattr(webbpsf, {
        'NIRCAM': 'NIRCam',
        'MIRI': 'MIRI',
        'NIRISS': 'NIRISS',
        'NIRSPEC': 'NIRSpec',
        'FGS': 'FGS',
    }.get(inst_name, inst_name), None)

    if inst_cls is None:
        raise ValueError(f"Unknown WebbPSF instrument: {instrument}")

    inst = inst_cls()
    inst.filter = filter_name.upper()

    if jitter_sigma > 0:
        inst.options['jitter'] = 'gaussian'
        inst.options['jitter_sigma'] = jitter_sigma

    if focus_offset != 0.0:
        inst.options['defocus_waves'] = focus_offset

    hdulist = inst.calc_psf(oversample=oversample, fov_pixels=psf_size)

    # Use the detector-sampled extension (ext 0).
    psf = hdulist[0].data.astype(np.float64)

    # Ensure correct size
    if psf.shape[0] > psf_size or psf.shape[1] > psf_size:
        cy, cx = psf.shape[0] // 2, psf.shape[1] // 2
        half = psf_size // 2
        psf = psf[cy - half:cy - half + psf_size,
                  cx - half:cx - half + psf_size]

    # Normalize
    total = psf.sum()
    if total > 0:
        psf /= total

    if output:
        from astropy.io import fits
        hdu = fits.PrimaryHDU(psf.astype(np.float32))
        hdu.header['PSFMETH'] = 'webbpsf'
        hdu.header['TELESCOP'] = 'JWST'
        hdu.header['INSTRUME'] = inst_name
        hdu.header['FILTER'] = filter_name.upper()
        hdu.header['PSFSIZE'] = psf_size
        hdu.header['OVERSAMP'] = oversample
        hdu.header['JITTER'] = jitter_sigma
        hdu.header['FOCUSOFF'] = focus_offset
        hdu.writeto(output, overwrite=True)
        print(f"Saved WebbPSF: {output}", file=sys.stderr)

    return psf


# ---- TinyTim generation ----

TINYTIM_INSTRUMENTS = {
    'ACS_WFC': {
        'chip': 15,
        'filters': [
            'F435W', 'F475W', 'F502N', 'F550M', 'F555W', 'F606W',
            'F625W', 'F658N', 'F660N', 'F775W', 'F814W', 'F850LP',
        ],
    },
    'ACS_HRC': {
        'chip': 16,
        'filters': [
            'F220W', 'F250W', 'F330W', 'F344N', 'F435W', 'F475W',
            'F502N', 'F550M', 'F555W', 'F606W', 'F625W', 'F658N',
            'F775W', 'F814W', 'F850LP',
        ],
    },
    'WFC3_UVIS': {
        'chip': 22,
        'filters': [
            'F200LP', 'F218W', 'F225W', 'F275W', 'F280N', 'F336W',
            'F350LP', 'F390M', 'F390W', 'F395N', 'F410M', 'F438W',
            'F467M', 'F475W', 'F475X', 'F487N', 'F502N', 'F547M',
            'F555W', 'F600LP', 'F606W', 'F621M', 'F625W', 'F631N',
            'F645N', 'F656N', 'F657N', 'F658N', 'F665N', 'F673N',
            'F680N', 'F689M', 'F763M', 'F775W', 'F814W', 'F845M',
            'F850LP', 'F953N',
        ],
    },
    'WFC3_IR': {
        'chip': 25,
        'filters': [
            'F098M', 'F105W', 'F110W', 'F125W', 'F126N', 'F127M',
            'F128N', 'F130N', 'F132N', 'F139M', 'F140W', 'F153M',
            'F160W', 'F164N', 'F167N',
        ],
    },
    'WFPC2': {
        'chip': 6,
        'filters': [
            'F218W', 'F255W', 'F300W', 'F336W', 'F380W', 'F439W',
            'F450W', 'F555W', 'F606W', 'F622W', 'F675W', 'F702W',
            'F791W', 'F814W', 'F850LP',
        ],
    },
}


def _tinytim_filter_number(instrument, filter_name):
    """Map filter name to TinyTim filter number (1-based index)."""
    inst = instrument.upper().replace('/', '_')
    if inst not in TINYTIM_INSTRUMENTS:
        raise ValueError(f"Unknown TinyTim instrument: {instrument}")

    filters = TINYTIM_INSTRUMENTS[inst]['filters']
    fname = filter_name.upper()
    if fname in filters:
        return filters.index(fname) + 1
    raise ValueError(
        f"Filter {filter_name} not available for {instrument}. "
        f"Available: {', '.join(filters)}")


def generate_tinytim(instrument, filter_name, psf_size=201,
                     oversample=1, focus_offset=0.0, output=None):
    """Generate an HST PSF by running the user-installed TinyTim programs.

    tiny1 writes the parameter file and tiny2 computes the PSF. The programs
    are not bundled. If they are not on PATH, this raises and does not write
    ``output``.
    """
    if not check_tinytim_available():
        raise RuntimeError(
            "TinyTim executables (tiny1, tiny2, tiny3) were not found on PATH. "
            "An instrument PSF is not produced here."
        )

    import subprocess
    from astropy.io import fits

    inst = instrument.upper().replace('/', '_')
    if inst not in TINYTIM_INSTRUMENTS:
        raise ValueError(f"Unknown TinyTim instrument: {instrument}")

    chip = TINYTIM_INSTRUMENTS[inst]['chip']
    filt_num = _tinytim_filter_number(instrument, filter_name)

    tmpdir = tempfile.mkdtemp(prefix='tinytim_')
    try:
        paramfile = os.path.join(tmpdir, 'tinytim')
        tiny1_input = (
            f"{paramfile}\n{chip}\n{filt_num}\n{psf_size}\n"
            f"{focus_offset}\n{paramfile}.tt3\n"
        )
        result = subprocess.run(
            ['tiny1', paramfile],
            input=tiny1_input, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            raise RuntimeError(f"tiny1 failed: {result.stderr}")

        result = subprocess.run(
            ['tiny2', paramfile + '.tt3'],
            capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            raise RuntimeError(f"tiny2 failed: {result.stderr}")

        psf_file = paramfile + '00.fits'
        if not os.path.exists(psf_file):
            for name in os.listdir(tmpdir):
                if name.endswith('.fits') and 'tinytim' in name:
                    psf_file = os.path.join(tmpdir, name)
                    break

        if not os.path.exists(psf_file):
            raise RuntimeError(f"TinyTim output not found in {tmpdir}")

        with fits.open(psf_file) as hdul:
            psf = hdul[0].data.astype(np.float64)

        if psf.shape[0] != psf_size or psf.shape[1] != psf_size:
            out = np.zeros((psf_size, psf_size), dtype=np.float64)
            sy, sx = psf.shape
            if sy >= psf_size and sx >= psf_size:
                cy, cx = sy // 2, sx // 2
                half = psf_size // 2
                psf = psf[cy - half:cy - half + psf_size,
                          cx - half:cx - half + psf_size]
            else:
                oy = (psf_size - sy) // 2
                ox = (psf_size - sx) // 2
                out[oy:oy + sy, ox:ox + sx] = psf
                psf = out

        total = psf.sum()
        if total > 0:
            psf /= total

        if output:
            hdu = fits.PrimaryHDU(psf.astype(np.float32))
            hdu.header['PSFMETH'] = 'tinytim'
            hdu.header['TELESCOP'] = 'HST'
            hdu.header['INSTRUME'] = inst
            hdu.header['FILTER'] = filter_name.upper()
            hdu.header['PSFSIZE'] = psf_size
            hdu.header['OVERSAMP'] = oversample
            hdu.header['FOCUSOFF'] = focus_offset
            hdu.writeto(output, overwrite=True)
            print(f"Saved TinyTim PSF: {output}", file=sys.stderr)

        return psf
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ---- Circular diffraction pupil ----
#
# JWST entries use a filled annulus: primary diameter 6.5 m and secondary
# diameter 0.74 m (linear obscuration 0.74/6.5). The real JWST pupil is
# segmented and is not what this draws. HST entries use a 2.4 m aperture
# and linear obscuration 0.33. Filter digits are microns times 100 for the
# JWST names (F150W -> 1.50 um) and nanometres for the HST names
# (F814W -> 814 nm).

_JWST_PUPIL = {"diameter_m": 6.5, "obscuration": 0.114, "wave": "um100", "telescope": "JWST"}
_HST_PUPIL = {"diameter_m": 2.4, "obscuration": 0.33, "wave": "nm", "telescope": "HST"}

PUPILS = {
    "NIRCAM": dict(_JWST_PUPIL, pixel_scale="nircam"),
    "MIRI": dict(_JWST_PUPIL, pixel_scale=0.11),
    "NIRISS": dict(_JWST_PUPIL, pixel_scale=0.065),
    "FGS": dict(_JWST_PUPIL, pixel_scale=0.069, wavelength_um=0.65),
    "ACS_WFC": dict(_HST_PUPIL, pixel_scale=0.05),
    "ACS_HRC": dict(_HST_PUPIL, pixel_scale=0.025),
    "WFC3_UVIS": dict(_HST_PUPIL, pixel_scale=0.04),
    "WFC3_IR": dict(_HST_PUPIL, pixel_scale=0.13),
}


def _jinc(x):
    x = np.asarray(x, dtype=np.float64)
    out = np.ones(x.shape, dtype=np.float64)
    mask = np.abs(x) >= 1e-8
    out[mask] = 2.0 * j1(x[mask]) / x[mask]
    return out


def _annulus_intensity(rho, obscuration):
    """Intensity of a circular pupil with linear obscuration. The centre is 1."""
    eps = float(obscuration)
    amp = (_jinc(rho) - eps * eps * _jinc(eps * rho)) / (1.0 - eps * eps)
    return amp * amp


def filter_wavelength_um(filter_name, scheme):
    """Approximate wavelength in microns from a filter name such as F814W."""
    match = _FILTER_RE.match(str(filter_name).strip().upper())
    if match is None:
        raise ValueError("wavelength was not set")
    number = int(match.group(1))
    kind = "cutoff" if (match.group(2) or "").upper() == "LP" else "pivot"
    if scheme == "um100":
        return number / 100.0, kind
    if scheme == "nm":
        return number / 1000.0, kind
    raise ValueError("wavelength was not set")


def resolve_diffraction(instrument=None, filter_name=None, diameter_m=None,
                        obscuration=None, wavelength_um=None, pixel_scale=None):
    """Fill missing pupil numbers from a known instrument. Explicit numbers win."""
    inst = None
    if instrument and str(instrument).lower() not in ("auto", "unknown", ""):
        inst = normalize_instrument_name(instrument)
    spec = PUPILS.get(inst, {}) if inst else {}
    if diameter_m is None:
        diameter_m = spec.get("diameter_m")
    if obscuration is None:
        obscuration = spec.get("obscuration", 0.0)
    wavekind = "given"
    filt = None if filter_name is None else str(filter_name).strip()
    filt_missing = filt is None or filt.lower() in ("auto", "unknown", "")
    if wavelength_um is None:
        if "wavelength_um" in spec and (filt_missing or filt.upper() == "FGS"):
            wavelength_um = spec["wavelength_um"]
            wavekind = "pivot"
        elif not filt_missing and spec.get("wave"):
            wavelength_um, wavekind = filter_wavelength_um(filt, spec["wave"])
        else:
            raise ValueError("wavelength was not set")
    if pixel_scale is None:
        scale = spec.get("pixel_scale")
        if scale == "nircam":
            pixel_scale = 0.031 if float(wavelength_um) < 2.4 else 0.063
        else:
            pixel_scale = scale
    if diameter_m is None:
        raise ValueError("diameter was not set")
    if pixel_scale is None:
        raise ValueError("pixel scale was not set")
    diameter_m = float(diameter_m)
    obscuration = float(obscuration)
    wavelength_um = float(wavelength_um)
    pixel_scale = float(pixel_scale)
    if diameter_m <= 0 or wavelength_um <= 0 or pixel_scale <= 0:
        raise ValueError("diameter, wavelength, and pixel scale must be positive")
    if not (0.0 <= obscuration < 1.0):
        raise ValueError("obscuration must be in [0, 1)")
    return {
        "instrument": inst,
        "filter": None if filt_missing else filt.upper(),
        "telescope": spec.get("telescope"),
        "diameter_m": diameter_m,
        "obscuration": obscuration,
        "wavelength_um": wavelength_um,
        "wavekind": wavekind,
        "pixel_scale": pixel_scale,
    }


def generate_diffraction(psf_size=201, oversample=1, jitter_sigma=0.0,
                         diameter_m=None, obscuration=None, wavelength_um=None,
                         pixel_scale=None, instrument=None, filter_name=None,
                         output=None):
    """Circular-aperture PSF. The result is summed to 1 and labeled ogfmeas.

    ``jitter_sigma`` is a Gaussian sigma in arcsec. ``oversample`` draws on
    a finer grid and bins back to ``psf_size``. A missing diameter,
    wavelength, or pixel scale raises before any file is written.
    """
    spec = resolve_diffraction(
        instrument, filter_name, diameter_m, obscuration, wavelength_um, pixel_scale)
    psf_size = int(psf_size)
    oversample = int(oversample)
    if psf_size < 3 or oversample < 1:
        raise ValueError("psf_size must be at least 3 and oversample at least 1")
    fine_n = psf_size * oversample
    fine_scale = spec["pixel_scale"] / oversample
    axis = np.arange(fine_n, dtype=np.float64) - (fine_n - 1) / 2.0
    yy, xx = np.meshgrid(axis, axis, indexing="ij")
    theta_rad = np.hypot(xx, yy) * fine_scale / _ARCSEC_PER_RAD
    rho = np.pi * spec["diameter_m"] / (spec["wavelength_um"] * 1e-6) * theta_rad
    image = _annulus_intensity(rho, spec["obscuration"])
    if float(jitter_sigma) > 0:
        sigma = float(jitter_sigma) / fine_scale
        if sigma >= 0.05:
            image = gaussian_filter(image, sigma=sigma, mode="constant", cval=0.0)
    if oversample > 1:
        image = image.reshape(psf_size, oversample, psf_size, oversample).sum(axis=(1, 3))
    total = float(image.sum())
    if total <= 0 or not np.isfinite(total):
        raise RuntimeError("diffraction PSF is zero")
    image = image / total
    if output:
        from astropy.io import fits
        hdu = fits.PrimaryHDU(image.astype(np.float32))
        hdu.header["PSFMETH"] = "ogfmeas"
        hdu.header["PUPIL"] = "circular"
        if spec["telescope"]:
            hdu.header["TELESCOP"] = spec["telescope"]
        if spec["instrument"]:
            hdu.header["INSTRUME"] = spec["instrument"]
        if spec["filter"]:
            hdu.header["FILTER"] = spec["filter"]
        hdu.header["DIAM_M"] = spec["diameter_m"]
        hdu.header["OBSCFRAC"] = spec["obscuration"]
        hdu.header["WAVE_UM"] = spec["wavelength_um"]
        hdu.header["WAVEKIND"] = spec["wavekind"]
        hdu.header["PIXSCALE"] = spec["pixel_scale"]
        hdu.header["JITTER"] = float(jitter_sigma)
        hdu.header["PSFSIZE"] = psf_size
        hdu.header["OVERSAMP"] = oversample
        hdu.writeto(output, overwrite=True)
        print(f"Saved diffraction PSF: {output}", file=sys.stderr)
    return image


# ---- Instrument name normalization ----

def normalize_instrument_name(instrument, telescope='auto'):
    """Normalize instrument name for lookup.

    Maps FITS header values like 'ACS', 'ACS/WFC' to lookup keys like 'ACS_WFC'.
    """
    inst = instrument.upper().strip().replace('/', '_').replace(' ', '_')

    # Common aliases
    aliases = {
        'ACS': 'ACS_WFC',
        'WFC3': 'WFC3_UVIS',
        'WFPC': 'WFPC2',
    }
    return aliases.get(inst, inst)
