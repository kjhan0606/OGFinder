"""Simulation PSF generation via WebbPSF (JWST) and TinyTim (HST)."""

import os
import sys
import shutil
import tempfile
import numpy as np


# ---- Availability checks ----

def check_webbpsf_available():
    """Check if webbpsf package is importable."""
    try:
        import webbpsf  # noqa: F401
        return True
    except ImportError:
        return False


def check_tinytim_available():
    """Check if TinyTim executables (tiny1, tiny2, tiny3) are on PATH."""
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
    """Generate JWST PSF using WebbPSF.

    Parameters
    ----------
    instrument : str
        JWST instrument name (NIRCAM, MIRI, NIRISS, NIRSPEC, FGS).
    filter_name : str
        Filter name (e.g., F150W).
    psf_size : int
        Output PSF size in pixels.
    oversample : int
        Oversampling factor.
    jitter_sigma : float
        Jitter in arcsec (Gaussian sigma).
    focus_offset : float
        Focus offset in waves of optical path difference.
    output : str, optional
        Output FITS path.

    Returns
    -------
    psf : 2D ndarray
        PSF normalized to unit sum.
    """
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

    fov = psf_size * oversample
    hdulist = inst.calc_psf(oversample=oversample, fov_pixels=psf_size)

    # Use the detector-sampled extension (ext 0 = oversampled, ext 1 = detector)
    if oversample == 1:
        psf = hdulist[0].data.astype(np.float64)
    else:
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
    """Generate HST PSF using TinyTim.

    Parameters
    ----------
    instrument : str
        HST instrument (ACS_WFC, WFC3_UVIS, WFC3_IR, WFPC2, etc.).
    filter_name : str
        Filter name (e.g., F814W).
    psf_size : int
        Output PSF size in pixels.
    oversample : int
        Oversampling factor.
    focus_offset : float
        Focus offset in microns of secondary mirror despace.
    output : str, optional
        Output FITS path.

    Returns
    -------
    psf : 2D ndarray
        PSF normalized to unit sum.
    """
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

        # tiny1: create parameter file
        tiny1_input = f"{paramfile}\n{chip}\n{filt_num}\n{psf_size}\n{focus_offset}\n{paramfile}.tt3\n"
        result = subprocess.run(
            ['tiny1', paramfile],
            input=tiny1_input, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            raise RuntimeError(f"tiny1 failed: {result.stderr}")

        # tiny2: compute PSF
        result = subprocess.run(
            ['tiny2', paramfile + '.tt3'],
            capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            raise RuntimeError(f"tiny2 failed: {result.stderr}")

        # Find output FITS
        psf_file = paramfile + '00.fits'
        if not os.path.exists(psf_file):
            # Try alternative naming
            for f in os.listdir(tmpdir):
                if f.endswith('.fits') and 'tinytim' in f:
                    psf_file = os.path.join(tmpdir, f)
                    break

        if not os.path.exists(psf_file):
            raise RuntimeError(f"TinyTim output not found in {tmpdir}")

        with fits.open(psf_file) as hdul:
            psf = hdul[0].data.astype(np.float64)

        # Ensure correct size
        if psf.shape[0] != psf_size or psf.shape[1] != psf_size:
            out = np.zeros((psf_size, psf_size), dtype=np.float64)
            sy, sx = psf.shape
            # Center crop or pad
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

        # Normalize
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
