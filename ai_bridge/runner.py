"""The run engine: catalog -> records -> (cutouts) -> batches -> cache / send with retry -> output TSV + provenance."""
import hashlib
import json
import os
import platform
import sys
import time

from . import CONTRACT_VERSION, __version__, adapters, cache as cachemod, contracts, profile as profmod, records
from .errors import AuthError, BridgeError, ProfileError, RequestFailed, TemplateError

_sleep = time.sleep          # replaced in tests
_clock = time.monotonic


def now_iso():
    return time.strftime('%Y-%m-%dT%H:%M:%S%z')


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


class RateLimiter:
    def __init__(self, per_s):
        self.gap = (1.0 / per_s) if per_s and per_s > 0 else 0.0
        self.last = None

    def wait(self):
        if self.gap <= 0:
            return
        t = _clock()
        if self.last is not None and t - self.last < self.gap:
            _sleep(self.gap - (t - self.last))
        self.last = _clock()


class Stats:
    def __init__(self):
        self.requests_sent = 0
        self.cache_hits = 0
        self.retries = 0
        self.request_failures = 0
        self.bytes_sent = 0
        self.dry_run_requests = 0


def effective_profile(prof, task):
    p = profmod.with_defaults(prof)
    if p.get('task') in (None, 'any'):
        p['task'] = task
    return p


def send_with_retry(adapter, prepared, prof, stats, limiter, log):
    attempts = int(prof['retries']) + 1
    delay = float(prof['backoff_s'])
    last = None
    for a in range(attempts):
        limiter.wait()
        stats.requests_sent += 1
        stats.bytes_sent += prepared.size_bytes
        try:
            return adapter.send(prepared, float(prof['timeout_s']))
        except AuthError:
            stats.request_failures += 1
            raise
        except RequestFailed as e:
            last = e
            stats.request_failures += 1
            if not e.retryable or a == attempts - 1:
                raise
            wait = min(float(prof['backoff_max_s']), delay)
            if e.retry_after is not None:
                wait = min(float(prof['backoff_max_s']), max(wait, e.retry_after))
            stats.retries += 1
            log('  retry %d/%d in %.2fs after: %s' % (a + 1, attempts - 1, wait, e))
            _sleep(wait)
            delay *= float(prof['backoff_factor'])
    raise last


def run(profile, task, catalog, images, output=None, opts=None, log=None):
    """Execute a run.  Returns dict(summary, exit_code).  opts keys (all optional):
       numbers, id_col, mag_columns, size_pix, size_arcsec, normalize, cutout_format, cutout_dir, bands,
       pixel_scale, psf_fwhm, params(dict), dry_run, no_cache, cache_dir, resume, allow_network, provenance,
       max_objects, mock_fail_ids, limit_preview"""
    opts = dict(opts or {})
    log = log or (lambda m: print(m, file=sys.stderr, flush=True))
    t_start = time.time()
    if task not in contracts.TASKS:
        raise BridgeError('unknown task %r (known: %s)' % (task, ', '.join(contracts.TASK_NAMES)))
    prof = effective_profile(profile, task)
    errs, warns = profmod.validate(prof, task)
    if errs:
        raise ProfileError('profile %r is invalid:\n  ' % prof.get('name') + '\n  '.join(errs))
    for w in warns:
        log('warning: ' + w)
    network = prof['transport'] in profmod.NETWORK_TRANSPORTS
    agent = prof['transport'] == 'agent_cli'
    dry = bool(opts.get('dry_run'))
    if network and not dry and not opts.get('allow_network'):
        raise BridgeError('service %r sends data over the network; re-run with --allow-network '
                          '(or use --mode dry-run to inspect the request)' % prof['name'])
    if agent and not dry and not opts.get('allow_agent_cli'):
        raise BridgeError('service %r hands catalog data to the agent CLI %r, which forwards it to its provider\'s cloud '
                          'model; re-run with --allow-agent-cli (and --send-images only if images may leave the machine), '
                          'or use --mode dry-run to see the exact prompt and argv' % (prof['name'], prof.get('backend')))
    id_col = opts.get('id_col') or 'NUMBER'
    header, rows = records.read_catalog(catalog)
    recs = records.build_records(header, rows, id_col=id_col, mag_columns=opts.get('mag_columns'),
                                 pixel_scale=opts.get('pixel_scale'), psf_fwhm=opts.get('psf_fwhm'),
                                 numbers=opts.get('numbers'))
    if opts.get('max_objects'):
        recs = recs[:int(opts['max_objects'])]
    if not recs:
        raise BridgeError('no objects selected')
    ids = [r['id'] for r in recs]

    # ---- images / cutouts
    cu = dict(prof.get('cutouts') or {})
    imgs = list(images or [])
    bands = list(opts.get('bands') or [])
    maker_holder = {}
    if imgs:
        if not bands:
            bands = ['band%d' % (i + 1) for i in range(len(imgs))] if len(imgs) > 1 else ['main']
        if len(bands) != len(imgs):
            raise BridgeError('--bands has %d names for %d images' % (len(bands), len(imgs)))
        ps0 = opts.get('pixel_scale')
        if ps0 is None:
            try:
                from . import cutouts as cut_mod
                ps0 = cut_mod.header_pixel_scale(imgs[0])
            except Exception:            # numpy/astropy missing or unreadable header: only needed for cutouts
                ps0 = None
        for r in recs:
            if r['pixel_scale_arcsec'] is None:
                r['pixel_scale_arcsec'] = ps0
        sp, sa = opts.get('size_pix'), opts.get('size_arcsec')
        if sp is None and sa is None:
            sp, sa = cu.get('size_pix'), cu.get('size_arcsec')
        cfg = dict(size_pix=sp, size_arcsec=sa, norm=opts.get('normalize') or cu.get('normalize', 'asinh'),
                   fmt=opts.get('cutout_format') or cu.get('format', 'png'), pixel_scale=opts.get('pixel_scale'))
        cu['format'] = cfg['fmt']
        prof['cutouts'] = dict(cu, format=cfg['fmt'])

        def get_maker():          # lazy: images are only opened when a cutout is really needed
            if 'm' not in maker_holder:
                from . import cutouts as cut_mod
                cdir = opts.get('cutout_dir')
                if not cdir:
                    import tempfile
                    cdir = tempfile.mkdtemp(prefix='ogf_ai_cutouts_')
                    opts['_tmp_cutout_dir'] = cdir
                ims = cut_mod.ImageSet(dict(zip(bands, imgs)))
                maker_holder['im'] = ims
                maker_holder['m'] = cut_mod.CutoutMaker(ims, cdir, **cfg)
            return maker_holder['m']
    else:
        get_maker = None
    params = dict(prof.get('params') or {})
    params.update(opts.get('params') or {})
    ctx = {'params': params, 'bands': bands, 'mock_fail_ids': set(opts.get('mock_fail_ids') or []),
           'send_images': bool(opts.get('send_images'))}
    ad = adapters.make_adapter(prof, ctx)
    phash = profmod.profile_hash(prof)

    needs = []
    for fmt, band, kind in (ad.cutout_needs() or []):
        bl = bands if band == '*' else [band or (bands[0] if bands else None)]
        for b in bl:
            if (fmt, b) not in needs:
                needs.append((fmt, b))
    cutcache, cutfail = {}, {}

    def cut(rec, fmt, band):
        k = (rec['id'], fmt, band)
        if k in cutcache:
            return cutcache[k]
        if get_maker is None:
            raise ValueError('this service needs image cutouts but no --image was given')
        if band not in bands:
            raise ValueError('band %r not among the images (%s)' % (band, ', '.join(bands)))
        info = get_maker().make(rec, band, fmt)
        if info is None:
            raise ValueError('object %s is outside image %s (or has no position)' % (rec['id'], band))
        cutcache[k] = info
        return info

    ok_recs = []
    table = {}
    res_extra = {'tmp_cutout_dir': opts.get('_tmp_cutout_dir')}
    prefix = contracts.prov_prefix(task, prof)
    pc = {k: 'AI_%s_%s' % (prefix, k) for k in contracts.PROVENANCE}

    def fail(oid, msg):
        table[oid] = {pc['SERVICE']: prof['name'], pc['ERROR']: msg[:500]}

    for r in recs:
        if needs:
            try:
                for fmt, b in needs:
                    cut(r, fmt, b)
            except (ValueError, OSError) as e:
                fail(r['id'], 'cutout: %s' % e)
                continue
        ok_recs.append(r)

    # ---- resume
    if opts.get('resume') and output and os.path.exists(output):
        h, old = records.parse_tsv_table(open(output, encoding='utf-8').read())
        n_re = 0
        for oid, row in old.items():
            if oid in ids and not row.get(pc['ERROR']) and row.get(pc['SERVICE']) == prof['name']:
                table[oid] = {k: v for k, v in row.items()}
                n_re += 1
        ok_recs = [r for r in ok_recs if r['id'] not in table]
        log('resume: %d objects already done in %s' % (n_re, output))

    stats = Stats()
    cch = cachemod.Cache(opts.get('cache_dir'), prof['name'],
                         enabled=not opts.get('no_cache') and not dry and prof['transport'] != 'mock')
    limiter = RateLimiter(float(prof['rate_limit_per_s']))
    mapped_cols = None
    bs = int(prof['batch_size'])
    batches = [ok_recs[i:i + bs] for i in range(0, len(ok_recs), bs)]
    models = set()
    reqids = []
    previews = []
    failures = []
    max_payload = int(prof['max_payload_bytes'])
    for bi, batch in enumerate(batches):
        bids = [r['id'] for r in batch]
        try:
            prep = ad.prepare(batch, cut)
        except (TemplateError, ValueError, OSError, KeyError) as e:
            for i in bids:
                fail(i, 'request: %s' % e)
            continue
        if prep.size_bytes > max_payload:
            for i in bids:
                fail(i, 'request payload %d bytes exceeds max_payload_bytes=%d' % (prep.size_bytes, max_payload))
            continue
        if dry:
            stats.dry_run_requests += 1
            previews.append(prep.preview)
            continue
        key = cachemod.make_key(phash, prep.key_material)
        hit = cch.get(key)
        when = None
        if hit is not None:
            obj, meta, when = hit['response']['obj'], hit['response']['meta'], hit.get('time')
            stats.cache_hits += 1
        else:
            try:
                obj, meta = send_with_retry(ad, prep, prof, stats, limiter, log)
            except (RequestFailed, AuthError) as e:
                for i in bids:
                    fail(i, ('auth: ' if isinstance(e, AuthError) else 'request: ') + str(e))
                continue
            when = adapters.MOCK_TIME if prof['transport'] == 'mock' else now_iso()
            res_try, _m, _r = ad.extract(obj, meta, batch)
            # only cache responses that mapped for at least one object (never cache pure garbage)
            if any('values' in v for v in res_try.values()):
                cch.put(key, {'obj': obj, 'meta': meta}, {}, when)
        res, model, rid = ad.extract(obj, meta, batch)
        if model:
            models.add(model)
        if rid:
            reqids.append(rid)
        for i in bids:
            e = res.get(i, {'error': 'no result'})
            if 'values' in e:
                row = dict(e['values'])
                row[pc['SERVICE']] = prof['name']
                row[pc['MODEL']] = model or ''
                row[pc['REQID']] = rid or ''
                row[pc['TIME']] = when or ''
                row[pc['ERROR']] = ''
                table[i] = row
            else:
                fail(i, e['error'])
        if (bi + 1) % 20 == 0 or bi == len(batches) - 1:
            log('  %d/%d requests done' % (bi + 1, len(batches)))

    summary = {'service': prof['name'], 'task': task, 'objects': len(recs), 'dry_run': dry}
    if dry:
        summary.update(requests=len(previews), previews=previews, exit_code=0)
        return {'summary': summary, 'exit_code': 0, 'table': {}, 'columns': [], 'ids': ids, 'previews': previews}

    # ---- columns: mapped names in a stable order (task contract order first, then profile order)
    fields = prof['response'].get('fields')
    if not fields:
        have = []
        for oid in ids:
            for k in table.get(oid, {}):
                if k not in have and k not in pc.values():
                    have.append(k)
        fields_order = have
        std = contracts.output_names(task) if task != 'generic' else []
        fields_order = [c for c in std if c in have] + sorted(c for c in have if c not in std)
        types = {c: (contracts.output_spec(task, c) or {}).get('type') or _guess_type(table, c) for c in fields_order}
    else:
        fields_order = list(fields)
        types = {c: fields[c].get('type', 'str') for c in fields_order}
    columns = [(c, types[c]) for c in fields_order] + [(pc[k], 'str') for k in contracts.PROVENANCE]
    nfail = sum(1 for oid in ids if table.get(oid, {}).get(pc['ERROR']) or oid not in table)
    for oid in ids:
        if oid not in table:
            fail(oid, 'no result')
    nok = len(ids) - nfail
    exit_code = 0 if nfail == 0 else (2 if nok == 0 else 0)
    summary.update(objects_ok=nok, objects_failed=nfail, requests_sent=stats.requests_sent, cache_hits=stats.cache_hits,
                   retries=stats.retries, request_failures=stats.request_failures, batches=len(batches),
                   models=sorted(models), request_ids=len(reqids))
    for oid in ids:
        e = table[oid].get(pc['ERROR'])
        if e:
            failures.append({'id': oid, 'error': e})
    return {'summary': summary, 'exit_code': exit_code, 'table': table, 'columns': columns, 'ids': ids,
            'id_col': id_col, 'failures': failures, 'profile': prof, 'profile_hash': phash, 'stats': stats,
            'cache': cch, 'models': sorted(models), 't_start': t_start, 'network': network or agent,
            'send_images': bool(opts.get('send_images')),
            'catalog_header': header}


def _guess_type(table, col):
    for row in table.values():
        v = row.get(col)
        if v is None:
            continue
        return 'float' if isinstance(v, float) else 'int' if isinstance(v, int) and not isinstance(v, bool) else \
            'bool' if isinstance(v, bool) else 'json' if isinstance(v, (dict, list)) else 'str'
    return 'str'


def apply_rename(res, mapping):
    """Rename output columns {old: new} (e.g. PHOTOZ -> PHOTO_Z for the local-step column names).  In place."""
    if not mapping:
        return res
    names = [c for c, _ in res['columns']]
    for old, new in mapping.items():
        if old not in names:
            raise BridgeError('--rename: column %s is not produced by this run (%s)' % (old, ', '.join(names)))
        if new in names and new != old:
            raise BridgeError('--rename: %s already exists among the output columns' % new)
    res['columns'] = [(mapping.get(c, c), t) for c, t in res['columns']]
    for oid, row in res['table'].items():
        for old, new in mapping.items():
            if old in row:
                row[new] = row.pop(old)
    return res


def write_outputs(res, out_path, args_record, catalog, images):
    """Write TSV (to out_path or stdout when None) and return the provenance dict (sidecar written by the caller)."""
    import io
    buf = io.StringIO()
    records.write_tsv(buf, res['id_col'], res['ids'], res['columns'], res['table'])
    text = buf.getvalue()
    if out_path:
        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        tmp = out_path + '.tmp'
        with open(tmp, 'w', encoding='utf-8', newline='') as f:
            f.write(text)
        os.replace(tmp, out_path)
    else:
        sys.stdout.write(text)
        sys.stdout.flush()
    prof = res['profile']
    env, st = profmod.env_status(prof)
    s = res['stats']
    prov = {
        'ai_bridge_version': __version__, 'contract': CONTRACT_VERSION, 'python': platform.python_version(),
        'platform': platform.platform(),
        'service': prof['name'], 'transport': prof['transport'], 'task': res['summary']['task'],
        'profile_sha256': res['profile_hash'],
        'profile': profmod.redacted({k: v for k, v in prof.items() if not str(k).startswith('_')}),
        'auth_env_var': env, 'auth_env_status_at_run': st,
        'service_reported_models': res['models'],
        'started': time.strftime('%Y-%m-%dT%H:%M:%S%z', time.localtime(res['t_start'])),
        'finished': now_iso(), 'seconds': round(time.time() - res['t_start'], 3),
        'catalog': os.path.abspath(catalog), 'catalog_sha256': sha256_file(catalog),
        'images': [{'path': os.path.abspath(p), 'sha256': sha256_file(p)} for p in images],
        'args': args_record,
        'output': os.path.abspath(out_path) if out_path else '<stdout>',
        'output_sha256': hashlib.sha256(text.encode('utf-8')).hexdigest(),
        'output_columns': [c for c, _ in res['columns']],
        'objects': res['summary']['objects'], 'objects_ok': res['summary']['objects_ok'],
        'objects_failed': res['summary']['objects_failed'],
        'requests': {'sent': s.requests_sent, 'cache_hits': s.cache_hits, 'retries': s.retries,
                     'failed_attempts': s.request_failures, 'batches': res['summary']['batches'],
                     'bytes_sent': s.bytes_sent},
        'cache': {'enabled': res['cache'].enabled, 'dir': res['cache'].dir, 'hits': res['cache'].hits,
                  'misses': res['cache'].misses, 'stores': res['cache'].stores},
        'failures': res['failures'][:1000],
    }
    if prof['transport'] == 'agent_cli':
        from . import agent_cli
        path, _t = agent_cli.find_executable(prof)
        prov['agent_cli'] = {'backend': prof.get('backend'), 'flavor': agent_cli.flavor_of(prof, path),
                             'executable': path, 'version': agent_cli.version_of(path) if path else None,
                             'images_sent': bool(res.get('send_images')),
                             'env_names_passed': [n for n, st in agent_cli.passthrough_names(prof) if st],
                             'extra_args': list(prof.get('extra_args') or []),
                             'note': 'no credentials are stored; the CLI used its own login'}
    return prov, text
