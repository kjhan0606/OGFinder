"""Inject synthetic movers into the real HST ACS BB89 chips (j8pu38{c7,a,e,i}q), detect, and pickle (dets, truth, offsets).

usage: inject.py SEED NOBJ OUT.json [--nodq] [--cr-extra N] [--detect-opts JSON] [--save-chips CHIPS.pkl] [--chips CHIPS.pkl] [--field STEM,STEM,..]
  --field       use another set of exposures of the mast cache instead of BB89 (j8pu38*): the injection centre is the middle of the [SCI,1] chip of
                the first exposure (BB89 keeps its original centre).  Used to build the TRAINING sets of the real/bogus classifier on fields other
                than the evaluation field.

  --nodq        drop the archive cosmic-ray flags (DQ bit 4096) of all chips: the "no DQ available" case; real CRs are then unflagged.
  --cr-extra N  add N synthetic, unflagged cosmic-ray tracks per chip inside the analysed region (CR-heavy sets).
  --detect-opts JSON   keyword arguments for pipeline.detect_in_region (e.g. '{"cr_reject": false}'); default = pipeline defaults.
  --save-chips / --chips  cache the injected chips (pickle) so that detection variants can be re-run without re-injecting.

With no option the output is bit-identical to the original script (same RNG order).  Injected objects: Gaussian PSF (FWHM = max(estimated, 0.085")),
trailed over the exposure, Poisson noise of the source; rate 1.5-40 arcsec/h, mag 21-26.5, parallax Delta 1.5-4.5 AU."""
import sys, time, pickle, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
from scipy import ndimage as ndi
from moving import imaging as I, align as A, pipeline as P, util
from moving.align import tangent_inverse

C0 = os.path.expanduser('~/.ds9/mast_cache/')
NAMES = ['j8pu38c7q_flc', 'j8pu38caq_flc', 'j8pu38ceq_flc', 'j8pu38ciq_flc']
CENTER = (150.1375, 2.3610)


def gauss_psf(fwhm, size=41):
    r = size // 2; yy, xx = np.mgrid[-r:r + 1, -r:r + 1]; s = fwhm / 2.355
    g = np.exp(-(xx ** 2 + yy ** 2) / (2 * s * s)); return g / g.sum()


def add_cosmic_rays(c, n, rng, half=700, center=CENTER):
    """Add `n` synthetic cosmic-ray tracks (unflagged) to chip `c` around the field centre (e-/s).  A track is a random walk of 1-12
    pixels with a total charge drawn log-normally (median 700 e-, 0.35 dex), concentrated in 1-2 pixels per step, NOT convolved with the PSF
    (CRs are sharper than stars), Poisson noise on top."""
    x, y = c.wcs.all_world2pix([center[0]], [center[1]], 0); x = float(x[0]); y = float(y[0])
    ny, nx = c.shape
    add = np.zeros(c.shape, np.float32)
    for _ in range(n):
        xs = rng.uniform(max(2, x - half), min(nx - 3, x + half)); ys = rng.uniform(max(2, y - half), min(ny - 3, y + half))
        L = 1 + int(rng.exponential(3.0)); L = min(L, 12); ang = rng.uniform(0, 2 * np.pi)
        tot = 10 ** rng.normal(np.log10(700.0), 0.35)
        per = tot / L
        for k in range(L):
            xi = int(round(xs + k * np.cos(ang))); yi = int(round(ys + k * np.sin(ang)))
            if 1 <= xi < nx - 1 and 1 <= yi < ny - 1:
                v = per * rng.uniform(0.5, 1.5)
                add[yi, xi] += v
                if rng.random() < 0.3:                      # charge diffusion into one neighbour
                    add[yi + int(rng.integers(-1, 2)), xi + int(rng.integers(-1, 2))] += 0.15 * v
    e = add * c.texp
    e = e + rng.normal(0, np.sqrt(np.clip(e, 0, None)))
    c.data = c.data + (e / c.texp).astype(np.float32)


def make(seed, nobj, nodq=False, cr_extra=0, names=None):
    rng = np.random.default_rng(seed)
    chips = []
    for n in (names or NAMES):
        chips += I.load_chips(C0 + n + ".fits")
    A.align_field(chips)
    psf = I.estimate_psf(chips, snr_min=6.0)
    chips = [c for c in chips if c.name.endswith('[SCI,1]')]
    offs, tm = P.hst_offsets(chips)
    center = CENTER
    if names:
        cx, cy = chips[0].wcs.all_pix2world([chips[0].shape[1] / 2.0], [chips[0].shape[0] / 2.0], 0)
        center = (float(cx[0]), float(cy[0]))
    c0 = chips[0]
    ra0, de0 = center
    a0, d0 = np.radians(ra0), np.radians(de0)
    ex_ = np.array([-np.sin(a0), np.cos(a0), 0]); ey_ = np.array([-np.sin(d0) * np.cos(a0), -np.sin(d0) * np.sin(a0), np.cos(d0)])
    fw = max(psf[1] if psf[0] is not None else 1.7, 0.085 / chips[0].pixscale)
    truth = []
    for k in range(nobj):
        mag = rng.uniform(21.0, 26.5)
        rate = np.exp(rng.uniform(np.log(1.5), np.log(40)))      # arcsec/h geocentric
        pa = rng.uniform(0, 360); delta = rng.uniform(1.5, 4.5)
        xi0 = rng.uniform(-18, 18); eta0 = rng.uniform(-18, 18)   # arcsec offsets from center, at t_ref
        truth.append(dict(mag=mag, rate=rate, pa=pa, delta=delta, xi0=xi0, eta0=eta0))
    tref = np.mean(tm)
    for ci, c in enumerate(chips):
        img = np.zeros(c.shape, np.float32)
        for k, tk in enumerate(truth):
            vx = tk['rate'] * np.sin(np.radians(tk['pa'])) * 24.0; vy = tk['rate'] * np.cos(np.radians(tk['pa'])) * 24.0
            o = offs[ci] * 206264.806
            nsub = 40
            tt = np.linspace(-0.5, 0.5, nsub) * c.texp / 86400.0
            xs = []; ys = []
            for dt in tt:
                tday = tm[ci] - tref + dt
                ox = o @ ex_; oy = o @ ey_
                xs.append(tk['xi0'] + vx * tday - ox / tk['delta']); ys.append(tk['eta0'] + vy * tday - oy / tk['delta'])
            ra_s, de_s = tangent_inverse(np.array(xs) / 3600.0 * 0 + np.array(xs), np.array(ys), ra0, de0)
            px, py = c.wcs.all_world2pix(ra_s, de_s, 0)
            flux = 10 ** (-0.4 * (tk['mag'] - c.zp_ab))   # e-/s total
            g = gauss_psf(fw)
            tr = np.zeros(c.shape, np.float32)
            for xp, yp in zip(px, py):
                xi_, yi_ = int(round(xp)), int(round(yp)); r = 20
                if xi_ < r or yi_ < r or xi_ >= c.shape[1] - r or yi_ >= c.shape[0] - r: continue
                sh = ndi.shift(g, (yp - yi_, xp - xi_), order=1)
                tr[yi_ - r:yi_ + r + 1, xi_ - r:xi_ + r + 1] += sh * (flux / nsub)
            n_e = tr * c.texp
            noise = rng.normal(0, np.sqrt(np.clip(n_e, 0, None))) / c.texp
            img += tr + noise.astype(np.float32)
            if ci == 0:
                tk['x0'] = float(px[nsub // 2]); tk['y0'] = float(py[nsub // 2])
            tk.setdefault('pos', {})[ci] = (float(np.mean(px)), float(np.mean(py)), float(np.mean(ra_s)), float(np.mean(de_s)))
            # extra truth (not used by the older benchmark): trail end points and orientation in pixels
            tk.setdefault('trail', {})[ci] = (float(px[0]), float(py[0]), float(px[-1]), float(py[-1]))
        c.data = c.data + img
        c.inj = img                                   # injected-only image (e-/s), kept for the exact real/bogus labels (inj_labels)
    if nodq or cr_extra:
        rng2 = np.random.default_rng(seed + 1000003)             # separate stream: the injected movers are identical with/without CRs
        for c in chips:
            if cr_extra:
                add_cosmic_rays(c, cr_extra, rng2, center=center)
            if nodq:
                c.cr = np.zeros(c.shape, bool)
    return chips, truth, offs, psf, center


def inj_labels(dets, chips):
    """Exact labels from the injected-only image `chip.inj` (e-/s): `inj_box` = injected flux inside the detection's own footprint (a box rotated
    along the detection's major axis: half-length a_pix + 3 px (trail channel: a_pix is half the trail length), half-width 3 px) and
    `real` = that flux is a significant part of the detected flux (inj_box >= 0.3 |flux| and >= 2 flux_err; trail channel: >= 30 % of the injected
    flux within reach).  A cosmic ray that only lies within
    1 arcsec of an injected object is NOT real; an injected object found with a CR-contaminated centroid still is."""
    byname = {}
    for d in dets:
        d['inj_box'] = 0.0; d['real'] = False
    for c in chips:
        if not hasattr(c, 'inj'):
            continue
        for d in dets:
            if not d.get('chip', '').startswith(c.name):
                continue
            x = float(d['x_chip']); y = float(d['y_chip'])
            th = d.get('theta_pix', 0.0); th = float(th) if np.isfinite(th) else 0.0
            a = d.get('a_pix', 3.0); a = float(a) if np.isfinite(a) else 3.0
            hl = a + 3.0 if d.get('channel') == 'trail' else 3.0
            R = int(np.ceil(hl)) + 3
            xi, yi = int(round(x)), int(round(y))
            y0, y1, x0, x1 = max(0, yi - R), min(c.shape[0], yi + R + 1), max(0, xi - R), min(c.shape[1], xi + R + 1)
            if y1 <= y0 or x1 <= x0:
                continue
            yy, xx = np.mgrid[y0:y1, x0:x1]
            u = (xx - x) * np.cos(th) + (yy - y) * np.sin(th); v = -(xx - x) * np.sin(th) + (yy - y) * np.cos(th)
            sub = c.inj[y0:y1, x0:x1]
            box = float(sub[(np.abs(u) <= hl) & (np.abs(v) <= 3.0)].sum())
            fl = abs(float(d.get('flux_e_s', 0.0) or 0.0)); fe = float(d.get('flux_err', 0.0) or 0.0)
            d['inj_box'] = box
            ok = box > 0 and box >= 2.0 * fe and box >= 0.3 * fl
            if not ok and d.get('channel') == 'trail' and box >= 2.0 * fe:
                # the trail-channel flux is a sum of matched-filter (point-source) fluxes over the component, ~PSF-area times the true flux, so
                # it is not comparable with the injected total: require the footprint to hold >= 30 % of the injected flux within its reach
                ok = box >= 0.3 * float(sub[(np.abs(u) <= hl + 10.0) & (np.abs(v) <= 8.0)].sum())
            d['real'] = bool(d['sign'] > 0 and ok)
    return dets


def detect_and_save(chips, truth, offs, psf, center, out, detect_opts=None):
    t0 = time.time()
    res = P.detect_in_region(chips, center=center, half_pix=620, psf=psf, **(detect_opts or {}))
    dets, infos = res
    inj_labels(dets, chips)
    pickle.dump((dets, truth, offs), open(out + '.pkl', 'wb'))
    print(len(dets), 'dets', round(time.time() - t0, 1), 's detect')
    return dets


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = {}
    def pop(flag, val=False):
        if flag in a:
            i = a.index(flag)
            if val:
                opt[flag] = a[i + 1]; del a[i:i + 2]
            else:
                opt[flag] = True; del a[i]
    for f in ("--nodq",): pop(f)
    for f in ("--cr-extra", "--detect-opts", "--save-chips", "--chips", "--field"): pop(f, True)
    seed = int(a[0]); nobj = int(a[1]); out = a[2]
    if "--chips" in opt:
        chips, truth, offs, psf, center = pickle.load(open(opt["--chips"], "rb"))
    else:
        chips, truth, offs, psf, center = make(seed, nobj, nodq=bool(opt.get("--nodq")), cr_extra=int(opt.get("--cr-extra", 0)),
                                                  names=opt["--field"].split(",") if "--field" in opt else None)
        if "--save-chips" in opt:
            pickle.dump((chips, truth, offs, psf, center), open(opt["--save-chips"], "wb"), protocol=4)
    dets = detect_and_save(chips, truth, offs, psf, center, out, json.loads(opt.get("--detect-opts", "{}")))
    trs = P.link_detections(dets, offs, snr_min=8.0, tol_arcsec=1.0, min_exposures=3, max_per_exposure=600, max_tracklets=400)
    print(len(trs), 'tracklets')
    rec = []
    for k, tk in enumerate(truth):
        pos = tk['pos']
        det_hit = []
        for ci in range(len(pos)):
            px, py, ra, de = pos[ci]
            best = None
            for d in dets:
                if d['ex'] != ci or d['sign'] < 0: continue
                s = np.hypot((d['ra'] - ra) * np.cos(np.radians(de)), d['dec'] - de) * 3600
                if s < 1.0 and (best is None or s < best[0]): best = (s, d)
            det_hit.append(best)
        nd = sum(b is not None for b in det_hit)
        tr_hit = None
        for t in trs:
            m = 0
            for e, r, d in zip(t['ex'], t['ra'], t['dec']):
                px, py, ra, de = pos[e]
                if np.hypot((r - ra) * np.cos(np.radians(de)), d - de) * 3600 < 1.0: m += 1
            if m >= 3: tr_hit = t; break
        rec.append(dict(mag=tk['mag'], rate=tk['rate'], pa=tk['pa'], delta=tk['delta'], n_det=nd,
                        det_astrom=[round(b[0], 3) if b else None for b in det_hit],
                        det_cls=[b[1]['cls'] if b else None for b in det_hit],
                        det_snr=[round(b[1]['snr'], 1) if b else None for b in det_hit],
                        linked=tr_hit is not None,
                        rate_fit=tr_hit['rate_ash'] if tr_hit else None, pa_fit=tr_hit['pa_deg'] if tr_hit else None))
    json.dump(dict(rec=rec, n_tracklets=len(trs), n_false_tracklets=len(trs) - sum(r['linked'] for r in rec), n_dets=len(dets)), open(out, 'w'), indent=1)
    for r in rec: print(r)
