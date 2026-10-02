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
Advanced objects and modifiers (GALFIT 3.0.x; see docs/multifit.md for the measured agreement with the binary): moffat, ferrer, king, nuker, edgedisk (surface-brightness parameters mu [mag/arcsec^2]
are converted to counts per pixel with the plate scale K) and the flux normalisations sersic1 / sersic2 / sersic3 (also expdisk1.., devauc1..: `norm` center / re / break); C0 (diskiness / boxiness), F<m> (Fourier
modes), B<m> (bending modes), R0..R10 (coordinate rotation, `power` and `log`), truncation objects (`T0) radial`, T1 T4 T5 T9 T10, referenced by Ti) / To) of the object they modify).  Constraints:
bounds, `offset` and `ratio` ties.  Not supported (ValueError, or skipped with a warning when strict=False): powsersic, nuker/other objects not listed, truncation types other than `radial`, 'diffusion kernel'
in D); non-unit PSF fine sampling is passed through as cfg['psf_sampling'] (the fitter does not oversample)."""
import math
import os
import re

RS_TO_RE = 1.678347                     # b_1 : R_e = b_1 R_s for the exponential disc
SUPPORTED = ('sersic', 'devauc', 'expdisk', 'psf', 'gaussian', 'sky', 'moffat', 'ferrer', 'king', 'nuker', 'edgedisk', 'trunc')
_PNAMES = {'sersic': {1: 'xy', 3: 'mag', 4: 're', 5: 'n', 9: 'q', 10: 'pa'}, 'devauc': {1: 'xy', 3: 'mag', 4: 're', 9: 'q', 10: 'pa'},
           'expdisk': {1: 'xy', 3: 'mag', 4: 'rs', 9: 'q', 10: 'pa'}, 'gaussian': {1: 'xy', 3: 'mag', 4: 'fwhm', 9: 'q', 10: 'pa'},
           'psf': {1: 'xy', 3: 'mag'}, 'sky': {1: 'sky', 2: 'dsky_dx', 3: 'dsky_dy'},
           'moffat': {1: 'xy', 3: 'mag', 4: 'fwhm', 5: 'beta', 9: 'q', 10: 'pa'}, 'ferrer': {1: 'xy', 3: 'mu', 4: 'rout', 5: 'alpha', 6: 'beta', 9: 'q', 10: 'pa'},
           'king': {1: 'xy', 3: 'mu', 4: 'rc', 5: 'rt', 6: 'alpha', 9: 'q', 10: 'pa'}, 'nuker': {1: 'xy', 3: 'mu', 4: 'rb', 5: 'alpha', 6: 'beta', 7: 'gamma', 9: 'q', 10: 'pa'},
           'edgedisk': {1: 'xy', 3: 'mu', 4: 'hs', 5: 'rs', 10: 'pa'}}
# GALFIT type -> (multifit kind, key of the surface-brightness parameter or None)
_KIND = {'sersic': 'sersic', 'devauc': 'dev', 'expdisk': 'exp', 'gaussian': 'sersic', 'psf': 'psf', 'moffat': 'moffat', 'ferrer': 'ferrer', 'king': 'king', 'nuker': 'nuker', 'edgedisk': 'edgedisk'}
_SBKEY = {'ferrer': 'i0', 'king': 'i0', 'nuker': 'ib', 'edgedisk': 'i0'}
_NORM_SUFFIX = {'1': 'center', '2': 're', '3': 'break'}
_SIZE_KEY = {'sersic': 're', 'devauc': 're', 'expdisk': 're'}
HEADER_KEYS = 'ABCDEFGHIJKLMNOPQRSTUVWXY'


class FeedmeError(ValueError):
    pass


def _tokens(rest):
    return rest.split('#', 1)[0].split()


def parse_feedme_raw(text):
    """-> (header {letter: [tokens]}, objects [{type, params {n: [tokens]}, mods {'C0'|'F1'|'B2'|'R3'|'TI'|'TO'|'T4'...: [tokens]}, Z}]).  A `T0)` line starts a truncation object (type 'trunc')."""
    header, objs, cur = {}, [], None
    for ln, line in enumerate(text.splitlines(), 1):
        s = line.strip()
        if not s or s.startswith('#') or s.startswith('==='):
            continue
        m = re.match(r'^([A-Za-z][0-9]{0,2}|[Tt][IiOo]|[0-9]{1,2})\)\s*(.*)$', s)
        if not m:
            continue
        key, rest = m.group(1), m.group(2)
        toks = _tokens(rest)
        if key == '0':
            cur = dict(type=toks[0].lower() if toks else '', params={}, mods={}, line=ln)
            objs.append(cur)
        elif key.upper() == 'T0':
            cur = dict(type='trunc', ttype=(toks[0].lower() if toks else 'radial'), params={}, mods={}, line=ln)
            objs.append(cur)
        elif key.isdigit():
            if cur is None:
                raise FeedmeError('line %d: parameter %s) before any object' % (ln, key))
            cur['params'][int(key)] = toks
        elif key == 'Z' and cur is not None:
            cur['Z'] = toks[0] if toks else '0'
        elif len(key) >= 2 and key[0].upper() in 'CFBRT' and (key[1:].isdigit() or key.upper() in ('TI', 'TO')):
            if cur is None:
                raise FeedmeError('line %d: %s) before any object' % (ln, key))
            cur['mods'][key.upper()] = toks
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
    ps = cfg.get('plate_scale') or [1.0, 1.0]
    pixarea = float(ps[0]) * float(ps[1]) * cfg['exptime']
    zp_ = cfg['zp']
    comp_of_obj, pairs = {}, []
    for k, o in enumerate(objs, 1):
        t = o['type']
        norm = None
        m_ = re.match(r'^(sersic|expdisk|devauc)([123])$', t)
        if m_:
            t, norm = m_.group(1), _NORM_SUFFIX[m_.group(2)]
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
        if t == 'trunc':
            if o.get('ttype', 'radial') != 'radial':
                msg = 'object %d: truncation type %r is not supported (only "radial")' % (k, o.get('ttype'))
                if strict:
                    raise FeedmeError(msg)
                warn.append(msg + ', skipped')
                continue
            M_ = o['mods']
            c = dict(kind='trunc', fixed=[], bounds={}, galfit_object=k, galfit_type='trunc')
            if 'T1' in M_ and len(M_['T1']) >= 2:
                c['x'], c['y'] = _f(M_['T1'][0], 'x'), _f(M_['T1'][1], 'y')
                fl = [int(float(v)) for v in M_['T1'][2:4]] + [1, 1]
                c['fixed'] += [n_ for n_, f_ in (('x', fl[0]), ('y', fl[1])) if f_ != 1]
            for key, nm, dflt in (('T4', 'rbreak', None), ('T5', 'dsoft', None), ('T9', 'q', None), ('T10', 'pa', None)):
                if key in M_:
                    v, f = _val_flag(M_[key], key)
                    c[nm] = v + 90.0 if nm == 'pa' else v
                    if f != 1:
                        c['fixed'].append(nm)
            if 'rbreak' not in c or 'dsoft' not in c:
                raise FeedmeError('object %d: truncation needs T4 (break radius) and T5 (softening length)' % k)
            comps.append(c)
            comp_of_obj[k] = len(comps) - 1
            gobj.append(('comp', k))
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
        v, f = _val_flag(P[3], 'magnitude / surface brightness')
        sbkey = _SBKEY.get(t) or ('i0' if norm else None)
        if sbkey:                                               # surface brightness [mag/arcsec^2] -> counts per pixel (exptime and plate scale included)
            c['mu'] = v
            c[sbkey] = pixarea * 10 ** (-0.4 * (v - zp_))
            if f != 1:
                c['fixed'].append(sbkey)
        else:
            c['mag'] = v - lm
            if f != 1:
                c['fixed'].append('flux')
        if norm:
            c['norm'] = norm
        c['kind'] = _KIND[t]
        if t != 'psf':
            if t in ('sersic', 'devauc', 'expdisk', 'gaussian'):
                v, f = _val_flag(P.get(4), 'R_e / R_s / FWHM')
                c['re'] = v * (RS_TO_RE if t == 'expdisk' else 0.5 if t == 'gaussian' else 1.0)
                if f != 1:
                    c['fixed'].append('re')
            else:
                for num, nm in pn.items():
                    if num in (1, 3, 9, 10):
                        continue
                    if num not in P:
                        raise FeedmeError('object %d (%s): missing parameter %d) %s' % (k, t, num, nm))
                    v, f = _val_flag(P[num], nm)
                    c[nm] = v
                    if f != 1:
                        c['fixed'].append(nm)
            if t == 'sersic':
                v, f = _val_flag(P.get(5), 'n')
                c['n'] = v
                if f != 1:
                    c['fixed'].append('n')
            elif t == 'gaussian':
                c['n'] = 0.5
                c['fixed'].append('n')
            if t != 'edgedisk':
                v, f = _val_flag(P.get(9, ['1', '1']), 'axis ratio')
                c['q'] = v
                if f != 1:
                    c['fixed'].append('q')
            v, f = _val_flag(P.get(10, ['0', '1']), 'PA')
            c['pa'] = v + 90.0
            if f != 1:
                c['fixed'].append('pa')
            _read_modifiers(c, o['mods'], k, warn)
        elif o['mods']:
            warn.append('object %d (psf): shape modifiers are ignored' % k)
        comps.append(c)
        comp_of_obj[k] = len(comps) - 1
        pairs.append((c, o))
        gobj.append(('comp', k))
    for c, o in pairs:
        for key, nm in (('TI', 'trunc_in'), ('TO', 'trunc_out')):
            if key in o['mods']:
                try:
                    c[nm] = [comp_of_obj[int(float(v))] for v in o['mods'][key]]
                except (KeyError, ValueError):
                    raise FeedmeError('object %d: %s) refers to an object that is not a (supported) truncation component: %s' % (c['galfit_object'], key.capitalize(), ' '.join(o['mods'][key])))
        if c.get('norm') == 'break' and not (c.get('trunc_in') or c.get('trunc_out')):
            raise FeedmeError('object %d: normalisation at the break radius (suffix 3) needs a truncation' % c['galfit_object'])
        if (c.get('trunc_in') or c.get('trunc_out')) and 'norm' not in c:
            raise FeedmeError('object %d: GALFIT writes the flux of a truncated profile as a surface brightness: use sersic1/2/3 (centre / R_e / break radius)' % c['galfit_object'])
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


def _read_modifiers(c, M_, k, warn):
    """C0 / F<m> / B<m> / R0..R10 lines of one object -> flat keys of the component dict (see ogfkit.profiles)."""
    if 'C0' in M_:
        v, f = _val_flag(M_['C0'], 'C0')
        c['c0'] = v
        if f != 1:
            c['fixed'].append('c0')
    for key in sorted((q for q in M_ if re.match(r'^F[0-9]+$', q)), key=lambda q: int(q[1:])):
        m = int(key[1:])
        t = M_[key]
        if len(t) < 2:
            raise FeedmeError('object %d: %s) needs amplitude and phase' % (k, key))
        c['f%da' % m], c['f%dp' % m] = _f(t[0], key), _f(t[1], key)
        fl = [int(float(v)) for v in t[2:4]] + [1, 1]
        if fl[0] != 1:
            c['fixed'].append('f%da' % m)
        if fl[1] != 1:
            c['fixed'].append('f%dp' % m)
    for key in sorted((q for q in M_ if re.match(r'^B[0-9]+$', q)), key=lambda q: int(q[1:])):
        m = int(key[1:])
        v, f = _val_flag(M_[key], key)
        c['b%d' % m] = v
        if f != 1:
            c['fixed'].append('b%d' % m)
    if 'R0' in M_:
        fn = (M_['R0'][0] if M_['R0'] else 'power').lower()
        if fn in ('power', 'powerlaw'):
            fn = 'power'
        elif fn in ('log', 'logarithmic'):
            fn = 'log'
        else:
            if fn != 'none':
                raise FeedmeError('object %d: rotation function %r not supported (power, log)' % (k, fn))
            fn = ''
        if fn:
            c['rot_func'] = fn
            for key, nm in (('R1', 'rot_in'), ('R2', 'rot_out'), ('R3', 'rot_theta'), ('R4', 'rot_alpha' if fn == 'power' else 'rot_ws'), ('R9', 'rot_incl'), ('R10', 'rot_pa')):
                if key not in M_:
                    if key in ('R9', 'R10'):
                        c[nm] = 0.0
                        c['fixed'].append(nm)
                        continue
                    raise FeedmeError('object %d: rotation needs %s)' % (k, key))
                v, f = _val_flag(M_[key], key)
                c[nm] = v
                if f != 1:
                    c['fixed'].append(nm)


# ------------------------------------------------------------------------------------------------------------------ constraints
_CNAME = {'x': 'x', 'y': 'y', 'mag': 'flux', 're': 're', 'rs': 're', 'n': 'n', 'q': 'q', 'pa': 'pa', 'fwhm': 'fwhm', 'c0': 'c0', 'mu': 'i0', 'rb': 'rb', 'rout': 'rout', 'rc': 'rc', 'rt': 'rt',
          'alpha': 'alpha', 'beta': 'beta', 'gamma': 'gamma', 'hs': 'hs', 'rbreak': 'rbreak', 'dsoft': 'dsoft', 'rsoft': 'dsoft'}


def _cname(c, pname):
    """GALFIT constraint parameter name -> (multifit key, kind of conversion) for component c; None if unknown."""
    t = c['galfit_type']
    p = pname.lower()
    if p in ('mag', 'mu'):
        if t in _SBKEY or c.get('norm'):
            return (_SBKEY.get(t, 'i0'), 'mu')
        return ('flux', 'mag')
    if p == 'rs' and t == 'edgedisk':
        return ('rs', None)
    if p == 're' and t in ('gaussian',):
        return ('re', 'fwhm')
    if p == 'fwhm' and t == 'gaussian':
        return ('re', 'fwhm')
    if p == 'n' and t == 'moffat':
        return ('beta', None)
    if p in ('rs',) and t == 'expdisk':
        return ('re', 'rs')
    if p == 'pa':
        return ('pa', 'pa')
    if p in ('f%s' % q for q in range(1, 40)) or re.match(r'^f[0-9]+$', p):
        return None
    if re.match(r'^b[0-9]+$', p):
        return (p, None)
    if p in _CNAME:
        return (_CNAME[p], None)
    return None


def _bound_to_multifit(c, name, lo, hi):
    """GALFIT-unit absolute (lo, hi) -> multifit bound entry (name, kind, (lo, hi)) for component dict c (config units)."""
    key, conv = _cname(c, name)
    if conv == 'rs':
        return key, None, (lo * RS_TO_RE, hi * RS_TO_RE)
    if conv == 'fwhm':
        return key, None, (lo * 0.5, hi * 0.5)
    if conv == 'pa':
        return key, None, (lo + 90.0, hi + 90.0)
    return key, conv, (lo, hi)


def _start_in_galfit_units(c, name, lm, key, conv):
    v = c.get(key)
    if conv == 'mag':
        return c['mag'] + lm
    if conv == 'mu':
        return c['mu']
    if conv == 'rs':
        return v / RS_TO_RE
    if conv == 'fwhm':
        return v * 2.0
    if conv == 'pa':
        return v - 90.0
    return v


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
        if any(w in t[2:] for w in ('offset', 'ratio')):
            mode = 'ratio' if 'ratio' in t[2:] else 'offset'
            if len(ids) < 2:
                warn.append('constraint "%s": %s needs at least two objects' % (s, mode))
                continue
            m = ids[0]
            for o in ids[1:]:
                if m not in comps or o not in comps:
                    warn.append('constraint "%s" refers to a skipped object' % s)
                    continue
                cm, co = _cname(comps[m], pname), _cname(comps[o], pname)
                if cm is None or co is None or cm[0] != co[0] and pname not in ('mag', 'mu'):
                    warn.append('constraint "%s": parameter %r is not supported' % (s, pname))
                    continue
                key = cm[0]
                if cm[1] in ('mag', 'mu') or key in ('flux', 'i0', 'ib'):        # magnitude offset = flux ratio (surface brightness likewise)
                    if cm[1] == 'mag' or co[1] == 'mag':
                        ties.append(['%d.flux' % index[o], '%d.flux' % index[m], 'ratio'])
                    else:
                        ties.append(['%d.%s' % (index[o], key), '%d.%s' % (index[m], key), 'ratio'])
                else:
                    if mode == 'offset' and abs(comps[o].get(key, 0.0) - comps[m].get(key, 0.0)) < 1e-12:
                        ties.append(['%d.%s' % (index[o], key), '%d.%s' % (index[m], key)])       # plain equality
                    else:
                        ties.append(['%d.%s' % (index[o], key), '%d.%s' % (index[m], key), mode])
            continue
        if len(ids) != 1 or ids[0] not in comps or _cname(comps[ids[0]], pname) is None:
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
        key, conv = _cname(c, pname)
        if rel:                                     # relative to the start value, in GALFIT units
            g0 = _start_in_galfit_units(c, pname, lm, key, conv)
            lo, hi = g0 + lo, g0 + hi
        key, kind, (blo, bhi) = _bound_to_multifit(c, pname, lo, hi)
        if kind == 'mag':                            # mag bound -> flux bound (counts at the config zero point)
            zp = cfg['zp']
            blo, bhi = 10 ** (-0.4 * (bhi - lm - zp)), 10 ** (-0.4 * (blo - lm - zp))
        elif kind == 'mu':
            ps = cfg.get('plate_scale') or [1.0, 1.0]
            pixarea = float(ps[0]) * float(ps[1]) * cfg.get('exptime', 1.0)
            blo, bhi = pixarea * 10 ** (-0.4 * (bhi - cfg['zp'])), pixarea * 10 ** (-0.4 * (blo - cfg['zp']))
        c['bounds'][key] = [blo, bhi]


# ------------------------------------------------------------------------------------------------------------------ export
def _fmt(v, p=4):
    return ('%.' + str(p) + 'f') % v


def _write_modifiers(out, c, fl, ic):
    if 'c0' in c:
        out.append('C0) %s %d   # diskyness / boxyness' % (_fmt(c['c0']), fl(c, 'c0')))
    for key in sorted((q for q in c if re.match(r'^f[0-9]+a$', q)), key=lambda q: int(q[1:-1])):
        m = int(key[1:-1])
        out.append('F%d) %s %s %d %d   # Fourier mode %d: amplitude, phase [deg]' % (m, _fmt(c[key]), _fmt(c.get('f%dp' % m, 0.0)), fl(c, key), fl(c, 'f%dp' % m), m))
    for key in sorted((q for q in c if re.match(r'^b[0-9]+$', q)), key=lambda q: int(q[1:])):
        out.append('B%s) %s %d   # bending mode %s' % (key[1:], _fmt(c[key]), fl(c, key), key[1:]))
    if c.get('rot_func'):
        out.append('R0) %-12s # coordinate rotation function' % c['rot_func'])
        for key, nm in (('R1', 'rot_in'), ('R2', 'rot_out'), ('R3', 'rot_theta'), ('R4', 'rot_alpha' if c['rot_func'] == 'power' else 'rot_ws'), ('R9', 'rot_incl'), ('R10', 'rot_pa')):
            out.append('%s) %s %d   # %s' % (key, _fmt(c.get(nm, 0.0)), fl(c, nm) if nm in c else 0, nm))
    for key, nm in (('Ti', 'trunc_in'), ('To', 'trunc_out')):
        if c.get(nm):
            out.append('%s) %s   # truncation object(s)' % (key, ' '.join(str(i + 1) for i in c[nm])))


def write_feedme(cfg, image='image.fits', output='imgblock.fits', sigma='none', psf='none', mask='none', constraints='none', region=None, shape=None, conv_box=None, plate_scale=None,
                 zp=None, exptime=None, mode=0, comment=''):
    """Config dict (as returned by parse_feedme / used by plugins/multifit --config) -> GALFIT feedme text.
    Component values are written as the start values (use the fitted config to export a result).  Object order: components, then the sky."""
    plate_scale = plate_scale or cfg.get('plate_scale') or (1.0, 1.0)
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
    pixarea = float(plate_scale[0]) * float(plate_scale[1]) * exptime
    gtype = {v: k for k, v in _KIND.items() if k != 'gaussian'}
    for ic, c in enumerate(cfg['components']):
        k = c['kind']
        if k == 'trunc':
            out.append('T0) radial                # truncation (object %d)' % (ic + 1))
            if 'x' in c:
                out.append('T1) %s %s %d %d   # centre x, y' % (_fmt(c['x']), _fmt(c['y']), fl(c, 'x'), fl(c, 'y')))
            out.append('T4) %s %d   # break radius [pix]' % (_fmt(c['rbreak']), fl(c, 'rbreak')))
            out.append('T5) %s %d   # softening length [pix]' % (_fmt(c['dsoft']), fl(c, 'dsoft')))
            if 'q' in c:
                out.append('T9) %s %d   # axis ratio' % (_fmt(c['q']), fl(c, 'q')))
            if 'pa' in c:
                out.append('T10) %s %d   # position angle' % (_fmt(((c['pa'] - 90.0 + 90.0) % 180.0) - 90.0), fl(c, 'pa')))
            out.append(' Z) 0                      # output option')
            out.append('')
            continue
        if k not in gtype:
            raise ValueError('component %d (%s) has no GALFIT equivalent' % (ic, k))
        t = gtype[k]
        if c.get('norm') and t in ('sersic', 'expdisk', 'devauc'):
            t += {v: q for q, v in _NORM_SUFFIX.items()}[c['norm']]
        out.append(' 0) %-18s # object type' % t)
        out.append(' 1) %s %s %d %d   # position x, y' % (_fmt(c['x']), _fmt(c['y']), fl(c, 'x'), fl(c, 'y')))
        sbk = _SBKEY.get(k) or ('i0' if c.get('norm') else None)
        if sbk:
            mu = zp - 2.5 * math.log10(max(c[sbk], 1e-300) / pixarea)
            out.append(' 3) %s %d   # surface brightness [mag/arcsec^2]' % (_fmt(mu), fl(c, sbk)))
        else:
            mag = c['mag'] if 'mag' in c else zp - 2.5 * math.log10(c['flux'])
            out.append(' 3) %s %d   # integrated magnitude' % (_fmt(mag + lm), fl(c, 'flux')))
        if k != 'psf':
            tt = 'sersic' if t.startswith('sersic') else 'expdisk' if t.startswith('expdisk') else 'devauc' if t.startswith('devauc') else t
            for num, nm in _PNAMES[tt].items():
                if num in (1, 3):
                    continue
                if nm == 'q':
                    out.append(' 9) %s %d   # axis ratio (b/a)' % (_fmt(c['q']), fl(c, 'q')))
                elif nm == 'pa':
                    out.append('10) %s %d   # position angle (PA) [deg: Up=0, Left=90]' % (_fmt(((c['pa'] - 90.0 + 90.0) % 180.0) - 90.0), fl(c, 'pa')))
                elif nm in ('re', 'rs'):
                    out.append(' 4) %s %d   # %s' % (_fmt(c['re'] / (RS_TO_RE if k == 'exp' else 1.0)), fl(c, 're'), 'R_s (disc scale length) [pix]' if k == 'exp' else 'R_e (half-light radius) [pix]'))
                elif tt == 'edgedisk' and nm == 'rs':
                    out.append(' 5) %s %d   # R_s' % (_fmt(c['rs']), fl(c, 'rs')))
                else:
                    out.append('%2d) %s %d   # %s' % (num, _fmt(c[nm]), fl(c, nm), nm))
            _write_modifiers(out, c, fl, ic)
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


def _gname(c, key):
    """multifit parameter key -> (GALFIT constraint name, conversion)."""
    if key in ('flux', 'i0', 'ib') and (c['kind'] in ('sersic', 'exp', 'dev', 'psf', 'moffat') and key == 'flux' or key in ('i0', 'ib')):
        return ('mag' if key == 'flux' else 'mu'), key
    if key == 're':
        return ('rs' if c['kind'] == 'exp' else 're'), key
    if key == 'beta' and c['kind'] == 'moffat':
        return 'n', key
    return key, key


def write_constraints(cfg):
    """Bounds and ties of a config dict -> GALFIT constraints file text (objects numbered in the order of write_feedme, the sky last)."""
    lines = ['# Component  parameter  constraint        written by ogfkit.galfitio']
    zp = cfg.get('zp', 25.0)
    lm = 2.5 * math.log10(max(cfg.get('exptime', 1.0), 1e-30))
    ps = cfg.get('plate_scale') or (1.0, 1.0)
    pixarea = float(ps[0]) * float(ps[1]) * cfg.get('exptime', 1.0)
    for i, c in enumerate(cfg['components'], 1):
        for k, (lo, hi) in (c.get('bounds') or {}).items():
            nm, _ = _gname(c, k)
            if k == 'flux':
                lines.append('%d  mag  %s to %s' % (i, _fmt(-2.5 * math.log10(hi) + zp + lm), _fmt(-2.5 * math.log10(max(lo, 1e-30)) + zp + lm)))
            elif k in ('i0', 'ib'):
                lines.append('%d  mag  %s to %s' % (i, _fmt(zp - 2.5 * math.log10(hi / pixarea)), _fmt(zp - 2.5 * math.log10(max(lo, 1e-300) / pixarea))))
            elif k == 're':
                f = RS_TO_RE if c['kind'] == 'exp' else 1.0
                lines.append('%d  %s  %s to %s' % (i, nm, _fmt(lo / f), _fmt(hi / f)))
            elif k == 'pa':
                lines.append('%d  pa  %s to %s' % (i, _fmt(lo - 90.0), _fmt(hi - 90.0)))
            else:
                lines.append('%d  %s  %s to %s' % (i, nm.upper() if re.match(r'^[fb][0-9]', nm) else nm, _fmt(lo), _fmt(hi)))
    for t in cfg.get('tie') or ():
        a, b = t[0], t[1]
        mode = t[2] if len(t) > 2 else 'offset'
        (ia, na), (ib, nb) = a.split('.'), b.split('.')
        ca, cb = cfg['components'][int(ia)], cfg['components'][int(ib)]
        if na == nb:
            nm, _ = _gname(ca, na)
            if na in ('flux', 'i0', 'ib') and mode == 'ratio':
                lines.append('# flux ratio %s = ratio * %s: GALFIT ties magnitudes by offset (mag offset = -2.5 log10 ratio), the value follows the start values' % (a, b))
                mode = 'offset'
            lines.append('%d_%d  %s  %s' % (int(ib) + 1, int(ia) + 1, nm, mode))
        else:
            lines.append('# tie %s = %s has no GALFIT equivalent' % (a, b))
    return '\n'.join(lines) + '\n'
