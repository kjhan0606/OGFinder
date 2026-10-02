"""GALFIT (Peng et al. 2002/2010) feedme / constraints import and export for ``ogfkit.multifit``.

    cfg = parse_feedme(text_or_path, exptime=1.0)     # -> config dict (the JSON config of plugins/multifit, plus psf/mask/input file names)
    text = write_feedme(cfg, ...)                      # config dict (or a fitted config) -> feedme text;  write_constraints(cfg) -> constraints text

Mapping (GALFIT -> multifit):  sersic -> sersic;  devauc -> dev;  expdisk -> exp (R_e = 1.67835 R_s);  psf -> psf;  gaussian -> sersic n = 0.5 with R_e = FWHM/2 (exact);
sky -> shared sky (value at the centre of the fitting region (xmin+xmax)/2, (ymin+ymax)/2, converted to the array centre of multifit, plus dsky/dx, dsky/dy).
Position angle: GALFIT measures PA from +y towards -x (up = 0, left = 90), multifit counter-clockwise from +x:  pa_multifit = pa_galfit + 90.
Coordinates: 1-based image pixels in both the feedme and the config dict.  Magnitudes: m_multifit = m_galfit - 2.5 log10(exptime) at the same zero point
(GALFIT: F = exptime 10^(-0.4 (m - zp)) ).  Fix flags: 0 = fixed, 1 = free (2 is treated as fixed).
Constraints file: `N param lo hi` (relative to the start value) and `N param lo to hi` (absolute) for x y mag re rs n q pa; `N_M param offset` ties N and M (accepted only if the start values
are equal, i.e. a common centre / PA / ...; other offsets and `ratio` are reported in cfg['warnings'] and ignored).  Object numbers count every feedme object including the sky.
Not supported (ValueError, or skipped with a warning when strict=False): moffat, nuker, ferrer, king, edgedisk, powsersic, fourier / bending / truncation modes, 'diffusion kernel'
in D), non-unit PSF fine sampling is passed through as cfg['psf_sampling'] (the fitter does not oversample)."""
import math
import os
import re

RS_TO_RE = 1.678347                     # b_1 : R_e = b_1 R_s for the exponential disc
SUPPORTED = ('sersic', 'devauc', 'expdisk', 'psf', 'gaussian', 'sky')
_PNAMES = {'sersic': {1: 'xy', 3: 'mag', 4: 're', 5: 'n', 9: 'q', 10: 'pa'}, 'devauc': {1: 'xy', 3: 'mag', 4: 're', 9: 'q', 10: 'pa'},
           'expdisk': {1: 'xy', 3: 'mag', 4: 'rs', 9: 'q', 10: 'pa'}, 'gaussian': {1: 'xy', 3: 'mag', 4: 'fwhm', 9: 'q', 10: 'pa'},
           'psf': {1: 'xy', 3: 'mag'}, 'sky': {1: 'sky', 2: 'dsky_dx', 3: 'dsky_dy'}}
HEADER_KEYS = 'ABCDEFGHIJKLMNOPQRSTUVWXY'


class FeedmeError(ValueError):
    pass


def _tokens(rest):
    return rest.split('#', 1)[0].split()


def parse_feedme_raw(text):
    """-> (header {letter: [tokens]}, objects [{type, params {n: [tokens]}, Z}])"""
    header, objs, cur = {}, [], None
    for ln, line in enumerate(text.splitlines(), 1):
        s = line.strip()
        if not s or s.startswith('#') or s.startswith('==='):
            continue
        m = re.match(r'^([A-Za-z]|[0-9]{1,2})\)\s*(.*)$', s)
        if not m:
            continue
        key, rest = m.group(1), m.group(2)
        toks = _tokens(rest)
        if key == '0':
            cur = dict(type=toks[0].lower() if toks else '', params={}, line=ln)
            objs.append(cur)
        elif key.isdigit():
            if cur is None:
                raise FeedmeError('line %d: parameter %s) before any object' % (ln, key))
            cur['params'][int(key)] = toks
        elif key == 'Z' and cur is not None:
            cur['Z'] = toks[0] if toks else '0'
        else:
            if cur is not None and key not in HEADER_KEYS:
                continue
            header[key.upper()] = toks
    return header, objs


def _f(tok, what):
    try:
        return float(tok)
    except (ValueError, TypeError):
        raise FeedmeError('cannot read %s from %r' % (what, tok))


def _val_flag(toks, what):
    if not toks:
        raise FeedmeError('missing value for %s' % what)
    return _f(toks[0], what), (int(float(toks[1])) if len(toks) > 1 else 1)


def parse_feedme(src, exptime=1.0, strict=True, base_dir=None):
    """Feedme text or file path -> config dict (see module doc)."""
    if '\n' not in src and os.path.isfile(src):
        base_dir = base_dir or os.path.dirname(os.path.abspath(src))
        src = open(src).read()
    header, objs = parse_feedme_raw(src)
    warn = []
    base_dir = base_dir or '.'

    def path(k):
        t = header.get(k, [])
        if not t or t[0].lower() in ('none', 'null', ''):
            return ''
        return t[0] if os.path.isabs(t[0]) else os.path.join(base_dir, t[0])
    cfg = dict(one_based=True, galfit=True, input=header.get('A', [''])[0] if header.get('A') else '', output=header.get('B', [''])[0] if header.get('B') else '',
               sigma=path('C'), psf=path('D'), mask=path('F'), constraints=path('G'), zp=_f(header['J'][0], 'zero point J)') if header.get('J') else 25.0, exptime=float(exptime))
    cfg['input'] = path('A') or cfg['input']
    if header.get('D') and len(header['D']) > 1 and header['D'][1].lower() != 'none':
        warn.append('diffusion kernel in D) is ignored')
    cfg['psf_sampling'] = int(float(header['E'][0])) if header.get('E') else 1
    if header.get('H'):
        t = [int(float(v)) for v in header['H'][:4]]
        if len(t) < 4:
            raise FeedmeError('H) needs xmin xmax ymin ymax')
        cfg['region'] = t
        cfg['bbox'] = [t[0] - 1, t[1], t[2] - 1, t[3]]
    if header.get('I'):
        cfg['conv_box'] = [int(float(v)) for v in header['I'][:2]]
    if header.get('K'):
        cfg['plate_scale'] = [float(v) for v in header['K'][:2]]
    comps, sky, gobj = [], None, []
    lm = 2.5 * math.log10(max(cfg['exptime'], 1e-30))
    for k, o in enumerate(objs, 1):
        t = o['type']
        if t not in SUPPORTED:
            if strict:
                raise FeedmeError('object %d: type %r is not supported by multifit (supported: %s)' % (k, t, ', '.join(SUPPORTED)))
            warn.append('object %d (%s) skipped' % (k, t))
            continue
        P = o['params']
        if t == 'sky':
            s, fs = _val_flag(P.get(1), 'sky')
            gx, fgx = _val_flag(P.get(2, ['0', '0']), 'dsky/dx')
            gy, fgy = _val_flag(P.get(3, ['0', '0']), 'dsky/dy')
            sky = dict(value=s, grad=[gx, gy], free=[fs == 1, fgx == 1, fgy == 1], obj=k)
            gobj.append(('sky', k))
            continue
        pn = _PNAMES[t]
        if 1 not in P or len(P[1]) < 2:
            raise FeedmeError('object %d (%s): missing position' % (k, t))
        c = dict(x=_f(P[1][0], 'x'), y=_f(P[1][1], 'y'), fixed=[], bounds={}, galfit_object=k, galfit_type=t)
        fl = [int(float(v)) for v in P[1][2:4]] + [1, 1]
        if fl[0] != 1:
            c['fixed'].append('x')
        if fl[1] != 1:
            c['fixed'].append('y')
        if 3 not in P:
            raise FeedmeError('object %d (%s): missing magnitude' % (k, t))
        v, f = _val_flag(P[3], 'magnitude')
        c['mag'] = v - lm
        if f != 1:
            c['fixed'].append('flux')
        if t == 'psf':
            c['kind'] = 'psf'
        else:
            c['kind'] = {'sersic': 'sersic', 'devauc': 'dev', 'expdisk': 'exp', 'gaussian': 'sersic'}[t]
            v, f = _val_flag(P.get(4), 'R_e / R_s / FWHM')
            c['re'] = v * (RS_TO_RE if t == 'expdisk' else 0.5 if t == 'gaussian' else 1.0)
            if f != 1:
                c['fixed'].append('re')
            if t == 'sersic':
                v, f = _val_flag(P.get(5), 'n')
                c['n'] = v
                if f != 1:
                    c['fixed'].append('n')
            elif t == 'gaussian':
                c['n'] = 0.5
                c['fixed'].append('n')
            v, f = _val_flag(P.get(9, ['1', '1']), 'axis ratio')
            c['q'] = v
            if f != 1:
                c['fixed'].append('q')
            v, f = _val_flag(P.get(10, ['0', '1']), 'PA')
            c['pa'] = v + 90.0
            if f != 1:
                c['fixed'].append('pa')
        comps.append(c)
        gobj.append(('comp', k))
    cfg['components'] = comps
    cfg['n_objects'] = len(objs)
    # sky mode / value at the multifit array centre
    if sky is None:
        cfg['sky'], cfg['sky_value'], cfg['sky_grad'] = 'fixed', 0.0, [0.0, 0.0]
    else:
        fs, fgx, fgy = sky['free']
        cfg['sky_grad'] = list(sky['grad'])
        cfg['sky_value_galfit'] = sky['value']
        if fgx or fgy:
            cfg['sky'] = 'plane'
            if not (fgx and fgy):
                warn.append('only one of dsky/dx, dsky/dy is free: both are fitted')
            if not fs:
                warn.append('sky fixed but gradient free: the sky level is fitted too')
        else:
            cfg['sky'] = 'const' if fs else 'fixed'
        cfg['sky_value'] = sky['value']
        if 'region' in cfg:                       # value at the GALFIT region centre -> value at multifit's array centre (nx/2, ny/2 in 0-based = nx/2+1 in 1-based)
            x0, x1, y0, y1 = cfg['region']
            nx, ny = x1 - x0 + 1, y1 - y0 + 1
            xc_g, yc_g = (x0 + x1) / 2.0, (y0 + y1) / 2.0
            xc_m, yc_m = x0 - 1 + nx / 2.0 + 1.0, y0 - 1 + ny / 2.0 + 1.0
            cfg['sky_value'] = sky['value'] + sky['grad'][0] * (xc_m - xc_g) + sky['grad'][1] * (yc_m - yc_g)
    if cfg.get('constraints'):
        if os.path.isfile(cfg['constraints']):
            _apply_constraints(cfg, open(cfg['constraints']).read(), warn)
        else:
            warn.append('constraints file %s not found' % cfg['constraints'])
    cfg['warnings'] = warn
    return cfg


# ------------------------------------------------------------------------------------------------------------------ constraints
_CNAME = {'x': 'x', 'y': 'y', 'mag': 'flux', 're': 're', 'rs': 're', 'n': 'n', 'q': 'q', 'pa': 'pa'}


def _bound_to_multifit(c, name, lo, hi):
    """GALFIT-unit absolute (lo, hi) -> multifit bound entry (name, lo, hi) for component dict c (config units)."""
    if name == 'mag':
        m_lo, m_hi = lo - 0.0, hi - 0.0
        return 'flux', 'mag', (m_lo, m_hi)
    if name == 'rs':
        return 're', None, (lo * RS_TO_RE, hi * RS_TO_RE)
    if name == 'pa':
        return 'pa', None, (lo + 90.0, hi + 90.0)
    return _CNAME[name], None, (lo, hi)


def _apply_constraints(cfg, text, warn):
    comps = {c['galfit_object']: c for c in cfg['components']}
    ties = cfg.setdefault('tie', [])
    index = {c['galfit_object']: i for i, c in enumerate(cfg['components'])}
    lm = 2.5 * math.log10(max(cfg.get('exptime', 1.0), 1e-30))
    for line in text.splitlines():
        s = line.split('#', 1)[0].strip()
        if not s:
            continue
        t = s.replace(',', ' ').split()
        if len(t) < 3:
            continue
        ids = t[0].split('_')
        try:
            ids = [int(v) for v in ids]
        except ValueError:
            warn.append('constraint line not understood: %s' % s)
            continue
        pname = t[1].lower()
        if pname in ('x', 'y', 'mag', 're', 'rs', 'n', 'q', 'pa') and any(w in t[2:] for w in ('offset', 'ratio')):
            if 'ratio' in t[2:] or len(ids) < 2:
                warn.append('constraint "%s" (ratio / single object) is not supported' % s)
                continue
            if pname in ('mag',):
                warn.append('constraint "%s": magnitude tie is not supported' % s)
                continue
            m = ids[0]
            for o in ids[1:]:
                if m not in comps or o not in comps:
                    warn.append('constraint "%s" refers to a skipped object' % s)
                    continue
                key = _CNAME[pname]
                v1 = comps[m].get(key if key != 'flux' else 'mag')
                v2 = comps[o].get(key if key != 'flux' else 'mag')
                if v1 is None or v2 is None or abs(v1 - v2) > 1e-6 * max(1.0, abs(v1)):
                    warn.append('constraint "%s": start values differ (%s vs %s), constant offsets are not supported' % (s, v1, v2))
                    continue
                ties.append(['%d.%s' % (index[o], key), '%d.%s' % (index[m], key)])
            continue
        if len(ids) != 1 or ids[0] not in comps or pname not in _CNAME:
            if len(ids) == 1 and ids[0] not in comps:
                warn.append('constraint "%s" refers to a sky / skipped object' % s)
            else:
                warn.append('constraint not understood: %s' % s)
            continue
        c = comps[ids[0]]
        try:
            if 'to' in t[2:]:
                j = t.index('to')
                lo, hi = float(t[2]), float(t[j + 1])
                rel = False
            else:
                lo, hi = float(t[2]), float(t[3])
                rel = True
        except (ValueError, IndexError):
            warn.append('constraint not understood: %s' % s)
            continue
        if rel:                                     # relative to the start value, in GALFIT units
            g0 = {'x': c['x'], 'y': c['y'], 'mag': c['mag'] + lm, 're': c.get('re', 0.0) / (RS_TO_RE if c['galfit_type'] == 'expdisk' else 0.5 if c['galfit_type'] == 'gaussian' else 1.0),
                  'rs': c.get('re', 0.0) / RS_TO_RE, 'n': c.get('n', 0.0), 'q': c.get('q', 1.0), 'pa': c.get('pa', 90.0) - 90.0}[pname]
            lo, hi = g0 + lo, g0 + hi
        key, kind, (blo, bhi) = _bound_to_multifit(c, pname, lo, hi)
        if kind == 'mag':                            # mag bound -> flux bound (counts at the config zero point)
            zp = cfg['zp']
            blo, bhi = 10 ** (-0.4 * (bhi - lm - zp)), 10 ** (-0.4 * (blo - lm - zp))
        c['bounds'][key] = [blo, bhi]


# ------------------------------------------------------------------------------------------------------------------ export
def _fmt(v, p=4):
    return ('%.' + str(p) + 'f') % v


def write_feedme(cfg, image='image.fits', output='imgblock.fits', sigma='none', psf='none', mask='none', constraints='none', region=None, shape=None, conv_box=None, plate_scale=(1.0, 1.0),
                 zp=None, exptime=None, mode=0, comment=''):
    """Config dict (as returned by parse_feedme / used by plugins/multifit --config) -> GALFIT feedme text.
    Component values are written as the start values (use the fitted config to export a result).  Object order: components, then the sky."""
    zp = cfg.get('zp', 25.0) if zp is None else zp
    exptime = cfg.get('exptime', 1.0) if exptime is None else exptime
    lm = 2.5 * math.log10(max(exptime, 1e-30))
    if region is None:
        if cfg.get('region'):
            region = cfg['region']
        elif cfg.get('bbox'):
            b = cfg['bbox']; region = [b[0] + 1, b[1], b[2] + 1, b[3]]
        elif shape is not None:
            region = [1, shape[1], 1, shape[0]]
        else:
            raise ValueError('region / bbox / shape needed')
    conv_box = conv_box or cfg.get('conv_box') or [min(100, region[1] - region[0] + 1), min(100, region[3] - region[2] + 1)]
    out = ['# GALFIT feedme written by OGFinder ogfkit.galfitio' + (' - ' + comment if comment else ''), '#' + '=' * 78,
           '# IMAGE and GALFIT CONTROL PARAMETERS',
           'A) %-18s # Input data image (FITS file)' % image, 'B) %-18s # Output data image block' % output, 'C) %-18s # Sigma image name' % sigma,
           'D) %-18s # Input PSF image and (optional) diffusion kernel' % psf, 'E) %-18d # PSF fine sampling factor relative to data' % int(cfg.get('psf_sampling', 1)),
           'F) %-18s # Bad pixel mask' % mask, 'G) %-18s # File with parameter constraints' % constraints,
           'H) %d %d %d %d   # Image region to fit (xmin xmax ymin ymax)' % tuple(region), 'I) %d %d   # Size of the convolution box (x y)' % tuple(conv_box),
           'J) %-18s # Magnitude photometric zeropoint' % _fmt(zp, 4), 'K) %s %s   # Plate scale (dx dy) [arcsec per pixel]' % (plate_scale[0], plate_scale[1]),
           'O) regular             # Display type (regular, curses, both)', 'P) %d                   # 0=optimize, 1=model, 2=imgblock, 3=subcomps' % mode, '']
    fl = lambda c, n: 0 if n in c.get('fixed', ()) else 1
    for c in cfg['components']:
        k = c['kind']
        out.append(' 0) %-18s # object type' % {'sersic': 'sersic', 'exp': 'expdisk', 'dev': 'devauc', 'psf': 'psf'}[k])
        out.append(' 1) %s %s %d %d   # position x, y' % (_fmt(c['x']), _fmt(c['y']), fl(c, 'x'), fl(c, 'y')))
        mag = c['mag'] if 'mag' in c else zp - 2.5 * math.log10(c['flux'])
        out.append(' 3) %s %d   # integrated magnitude' % (_fmt(mag + lm), fl(c, 'flux')))
        if k != 'psf':
            out.append(' 4) %s %d   # %s' % (_fmt(c['re'] / (RS_TO_RE if k == 'exp' else 1.0)), fl(c, 're'), 'R_s (disc scale length) [pix]' if k == 'exp' else 'R_e (half-light radius) [pix]'))
            if k == 'sersic':
                out.append(' 5) %s %d   # Sersic index n' % (_fmt(c['n']), fl(c, 'n')))
            out.append(' 9) %s %d   # axis ratio (b/a)' % (_fmt(c['q']), fl(c, 'q')))
            out.append('10) %s %d   # position angle (PA) [deg: Up=0, Left=90]' % (_fmt(((c['pa'] - 90.0 + 90.0) % 180.0) - 90.0), fl(c, 'pa')))
        out.append(' Z) 0                      # output option (0 = resid., 1 = do not subtract)')
        out.append('')
    mode_sky = cfg.get('sky', 'const')
    s = cfg.get('sky_value', 0.0) or 0.0
    g = cfg.get('sky_grad') or [0.0, 0.0]
    if shape is None and region:
        nx, ny = region[1] - region[0] + 1, region[3] - region[2] + 1
    else:
        ny, nx = shape
    xc_g, yc_g = (region[0] + region[1]) / 2.0, (region[2] + region[3]) / 2.0
    xc_m, yc_m = region[0] - 1 + nx / 2.0 + 1.0, region[2] - 1 + ny / 2.0 + 1.0
    s_g = s - g[0] * (xc_m - xc_g) - g[1] * (yc_m - yc_g)          # multifit value at its array centre -> GALFIT value at the region centre
    out += [' 0) sky                  # object type', ' 1) %s %d   # sky background at the centre of the fitting region [ADU]' % (_fmt(s_g, 5), 0 if mode_sky == 'fixed' else 1),
            ' 2) %s %d   # dsky/dx' % (_fmt(g[0], 6), 1 if mode_sky == 'plane' else 0), ' 3) %s %d   # dsky/dy' % (_fmt(g[1], 6), 1 if mode_sky == 'plane' else 0),
            ' Z) 0                      # output option', '', '#' + '=' * 78]
    return '\n'.join(out) + '\n'


def write_constraints(cfg):
    """Bounds and ties of a config dict -> GALFIT constraints file text (objects numbered in the order of write_feedme, the sky last)."""
    lines = ['# Component  parameter  constraint        written by ogfkit.galfitio']
    zp = cfg.get('zp', 25.0)
    lm = 2.5 * math.log10(max(cfg.get('exptime', 1.0), 1e-30))
    for i, c in enumerate(cfg['components'], 1):
        for k, (lo, hi) in (c.get('bounds') or {}).items():
            if k == 'flux':
                lines.append('%d  mag  %s to %s' % (i, _fmt(-2.5 * math.log10(hi) + zp + lm), _fmt(-2.5 * math.log10(max(lo, 1e-30)) + zp + lm)))
            elif k == 're':
                f = RS_TO_RE if c['kind'] == 'exp' else 1.0
                lines.append('%d  %s  %s to %s' % (i, 'rs' if c['kind'] == 'exp' else 're', _fmt(lo / f), _fmt(hi / f)))
            elif k == 'pa':
                lines.append('%d  pa  %s to %s' % (i, _fmt(lo - 90.0), _fmt(hi - 90.0)))
            elif k in ('x', 'y', 'n', 'q'):
                lines.append('%d  %s  %s to %s' % (i, k, _fmt(lo), _fmt(hi)))
    for a, b in cfg.get('tie') or ():
        (ia, na), (ib, nb) = a.split('.'), b.split('.')
        if na == nb and na in ('x', 'y', 'n', 'q', 're', 'pa'):
            lines.append('%d_%d  %s  offset' % (int(ib) + 1, int(ia) + 1, 'rs' if na == 're' and cfg['components'][int(ia)]['kind'] == 'exp' else na))
        else:
            lines.append('# tie %s = %s has no GALFIT equivalent' % (a, b))
    return '\n'.join(lines) + '\n'
