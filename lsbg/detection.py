"""Low-threshold source detection for LSBG candidates.

Uses SEP with low detection threshold and large convolution kernels
to enhance diffuse, extended sources that standard SExtractor misses.

Multi-scale detection (Greco+2018, Prole+2018): rebin at 1×, 2×, 4×
scales, detect at each scale, merge catalogs.
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
    elif name == 'gauss15x15':
        y, x = np.mgrid[-7:8, -7:8]
        k = np.exp(-(x * x + y * y) / (2 * 4.0 ** 2))
    elif name == 'gauss21x21':
        y, x = np.mgrid[-10:11, -10:11]
        k = np.exp(-(x * x + y * y) / (2 * 5.0 ** 2))
    elif name == 'gauss33x33':
        y, x = np.mgrid[-16:17, -16:17]
        k = np.exp(-(x * x + y * y) / (2 * 8.0 ** 2))
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
    # Ensure contiguous float64 (avoid copy if already correct)
    if data.dtype == np.float64 and data.flags['C_CONTIGUOUS']:
        data_c = data
    else:
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
        if bkg_rms.dtype == np.float64 and bkg_rms.flags['C_CONTIGUOUS']:
            rms = bkg_rms
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

    # Scale pixstack to image size: at low thresholds (~0.8σ), up to
    # ~21% of noise pixels exceed threshold, plus real signal pixels.
    # Use 50% of total pixels as buffer, with 10M minimum.
    npixels = data_c.shape[0] * data_c.shape[1]
    pixstack = max(10_000_000, int(npixels * 0.5))
    sep.set_extract_pixstack(pixstack)
    sep.set_sub_object_limit(4096)

    try:
        objects, segmap = sep.extract(data_c, **kwargs)
    except Exception as e:
        print(f"LSBG detection failed: {e}", file=sys.stderr)
        # Retry with higher threshold if pixstack/deblend overflow
        if 'pixel buffer full' in str(e) or 'deblending overflow' in str(e):
            print("  retrying with detect_thresh * 1.5...", file=sys.stderr)
            kwargs['thresh'] = detect_thresh * 1.5
            try:
                objects, segmap = sep.extract(data_c, **kwargs)
            except Exception as e2:
                print(f"  retry also failed: {e2}", file=sys.stderr)
                objects = np.array([])
                segmap = np.zeros(data.shape, dtype=np.int32)
        else:
            objects = np.array([])
            segmap = np.zeros(data.shape, dtype=np.int32)

    print(f"LSBG detect: {len(objects)} candidates at {detect_thresh}sigma, "
          f"minarea={detect_minarea}", file=sys.stderr)

    return objects, segmap, rms


def rebin_image(data, factor):
    """Rebin image by block-averaging factor×factor pixels.

    Parameters
    ----------
    data : 2D array
        Input image.
    factor : int
        Rebinning factor (2, 4, etc.).

    Returns
    -------
    rebinned : 2D array
        Rebinned image.
    """
    ny, nx = data.shape
    # Trim to multiple of factor
    ny_trim = (ny // factor) * factor
    nx_trim = (nx // factor) * factor
    trimmed = data[:ny_trim, :nx_trim]
    # Block average
    rebinned = trimmed.reshape(ny_trim // factor, factor,
                               nx_trim // factor, factor).mean(axis=(1, 3))
    return rebinned


def detect_multiscale(data, bkg_rms=None, detect_thresh=0.8,
                      detect_minarea=50, filter_kernel='gauss5x5',
                      deblend_nthresh=32, deblend_mincont=0.005,
                      scale_factors='1,2,4', match_radius=5.0):
    """Multi-scale LSBG detection.

    Detects sources at multiple rebinning scales to capture both
    compact and very diffuse LSBGs. Merges catalogs, removing duplicates.

    Parameters
    ----------
    data : 2D array
        Background-subtracted (cleaned) image.
    bkg_rms : 2D array or None
        RMS map at native resolution.
    detect_thresh : float
        Detection threshold in sigma.
    detect_minarea : int
        Minimum area at native scale.
    filter_kernel : str
        Convolution kernel name.
    deblend_nthresh, deblend_mincont : int, float
        Deblending parameters.
    scale_factors : str
        Comma-separated rebinning factors (e.g., '1,2,4').
    match_radius : float
        Cross-match radius in native pixels for duplicate removal.

    Returns
    -------
    objects : structured array
        Merged detection catalog (coordinates in native pixels).
    segmap : 2D int array
        Segmentation map from native-scale detection.
    rms : 2D array
        RMS map from native-scale detection.
    """
    if isinstance(scale_factors, str):
        factors = [int(f) for f in scale_factors.split(',')]
    else:
        factors = [int(f) for f in scale_factors]
    if 1 not in factors:
        factors = [1] + factors

    all_detections = []  # list of (x, y, a, b, theta, flux, flag, scale)
    native_segmap = None
    native_rms = None

    for factor in sorted(factors):
        if factor == 1:
            # Native-scale detection
            objects, segmap, rms = detect_lsbg_candidates(
                data, bkg_rms=bkg_rms,
                detect_thresh=detect_thresh,
                detect_minarea=detect_minarea,
                filter_kernel=filter_kernel,
                deblend_nthresh=deblend_nthresh,
                deblend_mincont=deblend_mincont,
            )
            native_segmap = segmap
            native_rms = rms

            for obj in objects:
                all_detections.append({
                    'x': float(obj['x']),
                    'y': float(obj['y']),
                    'a': float(obj['a']),
                    'b': float(obj['b']),
                    'theta': float(obj['theta']),
                    'flux': float(obj['flux']),
                    'flag': int(obj['flag']),
                    'scale': 1,
                })
        else:
            # Rebinned detection
            rebinned = rebin_image(data, factor)
            # Scale minimum area down by factor^2
            scaled_minarea = max(5, detect_minarea // (factor * factor))
            # Use slightly lower threshold for rebinned (noise averages down)
            scaled_thresh = detect_thresh * 0.9

            rms_rebin = None
            if bkg_rms is not None:
                # Rebin RMS: RMS scales as 1/sqrt(N) for block average
                rms_rebin = rebin_image(bkg_rms, factor) / np.sqrt(factor * factor)

            objects_r, _, _ = detect_lsbg_candidates(
                rebinned, bkg_rms=rms_rebin,
                detect_thresh=scaled_thresh,
                detect_minarea=scaled_minarea,
                filter_kernel=filter_kernel,
                deblend_nthresh=deblend_nthresh,
                deblend_mincont=deblend_mincont,
            )

            # Scale coordinates back to native resolution
            for obj in objects_r:
                all_detections.append({
                    'x': float(obj['x']) * factor + (factor - 1) * 0.5,
                    'y': float(obj['y']) * factor + (factor - 1) * 0.5,
                    'a': float(obj['a']) * factor,
                    'b': float(obj['b']) * factor,
                    'theta': float(obj['theta']),
                    'flux': float(obj['flux']) * factor * factor,
                    'flag': int(obj['flag']),
                    'scale': factor,
                })

    # Ensure valid fallbacks for segmap and rms
    if native_segmap is None:
        native_segmap = np.zeros(data.shape, dtype=np.int32)
    if native_rms is None:
        native_rms = np.ones_like(data, dtype=np.float64)

    if not all_detections:
        return np.array([]), native_segmap, native_rms

    # Merge: remove duplicates (keep native-scale when overlapping)
    merged = _merge_detections(all_detections, match_radius)

    print(f"Multi-scale: {len(all_detections)} total → "
          f"{len(merged)} merged", file=sys.stderr)

    # Convert to structured array matching SEP format
    objects_merged = _dicts_to_sep_array(merged)

    return objects_merged, native_segmap, native_rms


def _merge_detections(detections, match_radius):
    """Merge multi-scale detections, removing duplicates.

    Priority: native scale > 2× > 4× (prefer higher resolution).
    Uses KDTree for O(N log N) when N > 500, brute-force otherwise.
    """
    # Sort by scale (native first)
    sorted_dets = sorted(detections, key=lambda d: d['scale'])

    n = len(sorted_dets)
    if n <= 500:
        # Brute-force for small catalogs
        merged = []
        for det in sorted_dets:
            is_dup = False
            for existing in merged:
                dx = det['x'] - existing['x']
                dy = det['y'] - existing['y']
                dist = np.sqrt(dx * dx + dy * dy)
                effective_radius = match_radius * max(det['scale'],
                                                       existing['scale'])
                if dist < effective_radius:
                    is_dup = True
                    break
            if not is_dup:
                merged.append(det)
        return merged

    # KDTree approach for large catalogs
    from scipy.spatial import cKDTree

    # Maximum possible match radius (for tree query)
    max_scale = max(d['scale'] for d in sorted_dets)
    max_match = match_radius * max_scale

    merged = []
    merged_xy = []  # separate list for building tree incrementally

    for det in sorted_dets:
        if merged_xy:
            # Query existing merged detections
            tree = cKDTree(np.array(merged_xy))
            dists, idxs = tree.query([det['x'], det['y']], k=1)
            if dists < max_match:
                # Check exact effective radius
                existing = merged[idxs]
                effective_radius = match_radius * max(det['scale'],
                                                       existing['scale'])
                if dists < effective_radius:
                    continue
        merged.append(det)
        merged_xy.append([det['x'], det['y']])

    return merged


def _dicts_to_sep_array(detections):
    """Convert list of dicts to numpy structured array matching SEP format."""
    n = len(detections)
    dtype = np.dtype([
        ('x', np.float64), ('y', np.float64),
        ('a', np.float64), ('b', np.float64),
        ('theta', np.float64), ('flux', np.float64),
        ('flag', np.int32),
    ])
    arr = np.zeros(n, dtype=dtype)
    for i, d in enumerate(detections):
        arr[i]['x'] = d['x']
        arr[i]['y'] = d['y']
        arr[i]['a'] = d['a']
        arr[i]['b'] = d['b']
        arr[i]['theta'] = d['theta']
        arr[i]['flux'] = d['flux']
        arr[i]['flag'] = d['flag']
    return arr
