"""Native-grid diffuse-light extraction and empirical stack covariance diagnostics."""
import numpy as np

EDGES = np.array([50., 100., 150., 200., 300., 400., 500.])
MODELS = [(h, q, a) for h, q in ((100., 1.), (250., .65)) for a in (1., 3., 10.)]
FACTORS = (2.5, 3.5, 5.)


def geometry(shape, center, kpp):
    y, x = np.ogrid[:shape[0], :shape[1]]
    dx = (x-center[0])*kpp; dy = (y-center[1])*kpp
    radius = np.hypot(dx, dy)
    bins = np.searchsorted(EDGES, radius, side='right')-1
    bins[(radius < EDGES[0]) | (radius >= EDGES[-1])] = -1
    outer = (radius >= 600) & (radius < 900)
    tiles = []
    for y in range(0, shape[0], 64):
        for x in range(0, shape[1], 64):
            sl = np.s_[y:y+64, x:x+64]; area = int(outer[sl].sum())
            if area >= 128:
                tiles.append((y, x, sl, outer[sl], area))
    return dx, dy, radius, bins, tiles


def diffuse_model(dx, dy, scale, axis_ratio, amplitude):
    angle = np.deg2rad(37.)
    major = dx*np.cos(angle)+dy*np.sin(angle)
    minor = -dx*np.sin(angle)+dy*np.cos(angle)
    # amplitude is observed nJy/arcsec2 at elliptical radius 100 kpc.
    return amplitude*np.exp((100.-np.hypot(major, minor/axis_ratio))/scale)


def measure(image, variance, excluded, geom, center, kpp):
    dx, dy, radius, bins, tile_geometry = geom
    valid = ~excluded & np.isfinite(image) & np.isfinite(variance) & (variance > 0)
    use = valid & (bins >= 0); labels = bins[use]
    count = np.bincount(labels, minlength=6)
    sums = np.bincount(labels, weights=image[use], minlength=6)
    vsums = np.bincount(labels, weights=variance[use], minlength=6)
    yy, xx = np.nonzero(use)
    means = np.divide(sums, count, out=np.full(6, np.nan), where=count > 0)
    design = np.ones((6, 3))
    for j, values in [(1, (xx-center[0])*kpp/1000), (2, (yy-center[1])*kpp/1000)]:
        design[:, j] = np.divide(np.bincount(labels, weights=values, minlength=6), count,
                                  out=np.zeros(6), where=count > 0)
    area = np.pi*np.diff(EDGES**2)/kpp**2
    coverage = np.minimum(1., count/area)
    diag = np.divide(vsums, count**2, out=np.full(6, np.nan), where=count > 0)
    tiles = []
    for y, x, sl, region_mask, region in tile_geometry:
        keep = valid[sl] & region_mask; n = int(keep.sum())
        if n < 128 or 2*n < region:
            continue
        iy, ix = np.nonzero(keep)
        tiles.append([1., (x+ix.mean()-center[0])*kpp/1000,
                      (y+iy.mean()-center[1])*kpp/1000,
                      image[sl][keep].astype(float).mean(), variance[sl][keep].astype(float).sum()/n**2])
    result = dict(raw=means, corrected=np.full(6, np.nan), coverage=coverage, count=count,
                  nominal_raw=np.diag(diag), nominal_corrected=np.full((6, 6), np.nan),
                  plane=np.full(3, np.nan), n_tiles=len(tiles), design=design)
    if len(tiles) >= 10:
        tiles = np.array(tiles); matrix = tiles[:, :3]; sigma = np.sqrt(tiles[:, 4])
        weighted = matrix/sigma[:, None]
        coef, _, rank, _ = np.linalg.lstsq(weighted, tiles[:, 3]/sigma, rcond=None)
        if rank == 3:
            covariance = np.linalg.inv(weighted.T@weighted)
            result.update(corrected=means-design@coef,
                          nominal_corrected=np.diag(diag)+design@covariance@design.T, plane=coef)
    return result


def spatial_labels(ra, dec, width=5.):
    return np.array([f'{int(r//width)}_{int((d+90)//width)}' for r, d in zip(ra, dec)])


def covariance_diagnostics(values, labels, nominal, seed=20261015, replicas=2000):
    """Complete-case covariance; block resampling keeps every cell's fields together."""
    values = np.asarray(values, float); labels = np.asarray(labels)
    n, p = values.shape
    if n < 2 or not np.isfinite(values).all():
        return dict(status='insufficient_complete_fields', n=n, p=p), {}
    sample = np.atleast_2d(np.cov(values, rowvar=False, ddof=1)); independent = sample/n
    rng = np.random.default_rng(seed)
    counts = rng.multinomial(n, np.full(n, 1/n), size=replicas)
    iid_means = counts@values/n
    iid_cov = np.atleast_2d(np.cov(iid_means, rowvar=False, ddof=1))
    groups = np.unique(labels); group_n = np.array([(labels == g).sum() for g in groups])
    group_sum = np.array([values[labels == g].sum(axis=0) for g in groups])
    block_cov = np.full((p, p), np.nan); jk_cov = block_cov.copy(); folds = []
    if len(groups) >= 3:
        draw = rng.multinomial(len(groups), np.full(len(groups), 1/len(groups)), size=replicas)
        block_mean = (draw@group_sum)/(draw@group_n)[:, None]
        block_cov = np.atleast_2d(np.cov(block_mean, rowvar=False, ddof=1))
        jack = np.array([(values.sum(axis=0)-s)/(n-m) for s, m in zip(group_sum, group_n)])
        # Unequal delete-cell counts make this an explicit sensitivity diagnostic.
        centered = jack-jack.mean(axis=0)
        jk_cov = (len(groups)-1)/len(groups)*(centered.T@centered)
        # Fixed sky-cell folds; each held-out cell is excluded from training completely.
        for fold in range(5):
            test = np.isin(labels, groups[fold::5]); train = ~test
            nt, nh = int(train.sum()), int(test.sum())
            if nh < 2 or nt < max(10, p+2):
                continue
            train_cov = np.atleast_2d(np.cov(values[train], rowvar=False, ddof=1))
            pred = train_cov*(1/nt+1/nh)
            delta = values[test].mean(axis=0)-values[train].mean(axis=0)
            z = delta/np.sqrt(np.diag(pred))
            folds.append(dict(fold=fold, n_train=nt, n_test=nh, delta=delta.tolist(),
                              predicted_sigma=np.sqrt(np.diag(pred)).tolist(), z=z.tolist()))
    eig = np.linalg.eigvalsh(sample)
    ratio = np.divide(np.diag(iid_cov), np.diag(independent), out=np.full(p, np.nan), where=np.diag(independent)>0)
    z = np.array([f['z'] for f in folds])
    mean_nominal = np.mean(nominal, axis=0)/n
    result = dict(status='pilot_estimate', n=n, p=p, n_sky_cells=len(groups),
                  rank=int(np.linalg.matrix_rank(sample)), minimum_eigenvalue=float(eig.min()),
                  condition_number=float(np.linalg.cond(sample)),
                  iid_bootstrap_diagonal_ratio=ratio.tolist(),
                  empirical_to_nominal_variance_ratio=(np.diag(independent)/np.diag(mean_nominal)).tolist(),
                  folds=folds, cv_fraction_within_1p96=float(np.mean(abs(z) <= 1.96)) if z.size else None,
                  covariance_sample_size_gate=bool(n >= 5*p and len(groups) >= 2*p+1),
                  sky_correction_approved=False, science_measurement=False)
    arrays = dict(sample_covariance=sample, independent_stack_covariance=independent,
                  iid_bootstrap_covariance=iid_cov, spatial_bootstrap_covariance=block_cov,
                  delete_cell_covariance=jk_cov, nominal_independent_stack_covariance=mean_nominal,
                  mean=values.mean(axis=0), n_group_fields=group_n)
    return result, arrays
