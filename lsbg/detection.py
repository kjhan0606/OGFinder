"""Low-threshold source detection for LSBG candidates.

Uses SEP with low detection threshold and large convolution kernels
to enhance diffuse, extended sources that standard SExtractor misses.
"""

import sys
import numpy as np
try:
    import sep
except ImportError:
    import sep_pjw as sep


def build_convolution_kernel(name='gauss5x5'):
    """Build a convolution kernel for source detection.

    Larger kernels enhance diffuse, low-surface-brightness features.

    Parameters
    ----------
    name : str
        Kernel name: 'gauss3x3', 'gauss5x5', 'gauss7x7', 'gauss9x9',
        'tophat5', 'tophat7', 'mexhat'.

    Returns
    -------
    kernel : 2D array or None
        Convolution kernel (normalized). None if name is 'none'.
    """
    if name == 'none' or name is None:
        return None

    if name == 'gauss3x3':
        k = np.array([[0.0625, 0.125, 0.0625],
                       [0.125, 0.25, 0.125],
                       [0.0625, 0.125, 0.0625]])
    elif name == 'gauss5x5':
        y, x = np.mgrid[-2:3, -2:3]
        k = np.exp(-(x * x + y * y) / (2 * 1.5 ** 2))
    elif name == 'gauss7x7':
        y, x = np.mgrid[-3:4, -3:4]
        k = np.exp(-(x * x + y * y) / (2 * 2.0 ** 2))
    elif name == 'gauss9x9':
        y, x = np.mgrid[-4:5, -4:5]
        k = np.exp(-(x * x + y * y) / (2 * 2.5 ** 2))
    elif name == 'tophat5':
        y, x = np.mgrid[-2:3, -2:3]
        k = ((x * x + y * y) <= 4).astype(float)
    elif name == 'tophat7':
        y, x = np.mgrid[-3:4, -3:4]
        k = ((x * x + y * y) <= 9).astype(float)
    elif name == 'mexhat':
        y, x = np.mgrid[-4:5, -4:5]
        r2 = x * x + y * y
        sigma = 2.0
        k = (1.0 - r2 / sigma ** 2) * np.exp(-r2 / (2 * sigma ** 2))
    else:
        # Default to gauss5x5
        y, x = np.mgrid[-2:3, -2:3]
        k = np.exp(-(x * x + y * y) / (2 * 1.5 ** 2))

    k = k / k.sum()
    return k.astype(np.float64)


def detect_lsbg_candidates(data, bkg_rms=None, detect_thresh=0.8,
                            detect_minarea=50, filter_kernel='gauss5x5',
                            deblend_nthresh=32, deblend_mincont=0.005):
    """Detect LSBG candidates with low threshold on cleaned image.

    Parameters
    ----------
    data : 2D array
        Background-subtracted (cleaned) image.
    bkg_rms : 2D array or None
        RMS map. If None, computed from data.
    detect_thresh : float
        Detection threshold in units of sigma.
    detect_minarea : int
        Minimum area in pixels for a detection.
    filter_kernel : str
        Convolution kernel name for detection.
    deblend_nthresh : int
        Number of deblending thresholds.
    deblend_mincont : float
        Minimum deblending contrast.

    Returns
    -------
    objects : structured array
        SEP detection catalog.
    segmap : 2D int array
        Segmentation map.
    rms : 2D array
        RMS map used for detection.
    """
    data_c = np.ascontiguousarray(data, dtype=np.float64)

    # Compute RMS if not provided
    if bkg_rms is None:
        try:
            bkg = sep.Background(data_c)
            rms = bkg.rms()
        except Exception:
            rms = np.full_like(data_c,
                               np.std(data_c[data_c != 0]) if
                               np.any(data_c != 0) else 1.0)
    else:
        rms = np.ascontiguousarray(bkg_rms, dtype=np.float64)

    # Build convolution kernel
    kernel = build_convolution_kernel(filter_kernel)

    # Run SEP extraction
    kwargs = dict(
        thresh=detect_thresh,
        err=rms,
        minarea=detect_minarea,
        deblend_nthresh=deblend_nthresh,
        deblend_cont=deblend_mincont,
        segmentation_map=True,
    )
    if kernel is not None:
        kwargs['filter_kernel'] = kernel

    try:
        objects, segmap = sep.extract(data_c, **kwargs)
    except Exception as e:
        print(f"LSBG detection failed: {e}", file=sys.stderr)
        objects = np.array([])
        segmap = np.zeros(data.shape, dtype=np.int32)

    print(f"LSBG detect: {len(objects)} candidates at {detect_thresh}sigma, "
          f"minarea={detect_minarea}", file=sys.stderr)

    return objects, segmap, rms
