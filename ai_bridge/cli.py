"""Command line driver (also reachable as `python ds9_ai_bridge.py ...` and `python -m ai_bridge`)."""
import argparse
import json
import os
import shutil
import sys

from . import CONTRACT_VERSION, __version__, adapters, contracts, profile as profmod, runner
from .errors import BridgeError, ProfileError

EXIT_OK, EXIT_USAGE, EXIT_TOTAL_FAILURE, EXIT_PARTIAL_STRICT = 0, 1, 2, 3


def build_parser():
    ap = argparse.ArgumentParser(
        prog='ds9_ai_bridge.py',
        description='OGFinder AI bridge: send catalog objects / cutouts to a user-configured external service and '
                    'return new catalog columns (TSV on stdout or --output).  No models are included.')
    ap.add_argument('--mode', required=True, choices=['list-services', 'check-service', 'run', 'dry-run', 'validate-profile', 'set-enabled'])
    ap.add_argument('--service', help='service name in the profile file (or "mock")')
    ap.add_argument('--task', help='task type: ' + ', '.join(contracts.TASK_NAMES))
    ap.add_argument('--catalog', help='catalog TSV (OGFinder format; NUMBER column required)')
    ap.add_argument('--image', help='FITS image(s), comma separated; the first one is the pixel grid of X_IMAGE/Y_IMAGE')
    ap.add_argument('--bands', help='band names for --image files, comma separated (default main / band1,band2..)')
    ap.add_argument('--output', help='output TSV (NUMBER + new columns); default stdout')
    ap.add_argument('--provenance-output', '--provenance', dest='provenance', help='provenance JSON sidecar path (default: <output>.provenance.json if --output is given)')
    ap.add_argument('--services-file', help='profile JSON (default $OGF_AI_SERVICES_FILE or ~/.ds9/ai_services.json)')
    ap.add_argument('--numbers', help='only these NUMBER values (comma separated or @file with one per line)')
    ap.add_argument('--max-objects', type=int, help='process only the first N selected objects')
    ap.add_argument('--id-col', default='NUMBER')
    ap.add_argument('--mag-columns', help='comma separated MAG_* columns to send (default: all MAG_<band>, excluding APER/ISO..)')
    ap.add_argument('--size-pix', type=int, help='cutout size in pixels')
    ap.add_argument('--size-arcsec', type=float, help='cutout size in arcsec (needs a pixel scale)')
    ap.add_argument('--normalize', choices=['none', 'linear', 'zscale', 'asinh'])
    ap.add_argument('--cutout-format', choices=['png', 'fits', 'npy'])
    ap.add_argument('--cutout-dir', help='keep cutouts here (default: temporary directory, removed after the run)')
    ap.add_argument('--pixel-scale', type=float, help='arcsec / pixel when the WCS has none')
    ap.add_argument('--psf-fwhm', type=float, help='PSF FWHM in arcsec, passed to the service as psf_fwhm_arcsec')
    ap.add_argument('--param', action='append', default=[], metavar='KEY=VALUE',
                    help='service parameter available as {param_KEY} and in params (repeatable); numbers are parsed as numbers')
    ap.add_argument('--rename', action='append', default=[], metavar='OLD=NEW',
                    help='rename an output column (e.g. PHOTOZ=PHOTO_Z); repeatable')
    ap.add_argument('--no-cache', action='store_true', help='do not read or write ~/.ds9/ai_cache')
    ap.add_argument('--cache-dir', help='cache directory (default $OGF_AI_CACHE_DIR or ~/.ds9/ai_cache)')
    ap.add_argument('--resume', action='store_true', help='keep rows already completed (no error) in an existing --output')
    ap.add_argument('--allow-network', action='store_true',
                    help='required for services that send data over the network (like the session script)')
    ap.add_argument('--strict', action='store_true', help='exit 3 when some objects failed (default: 0 unless all failed)')
    ap.add_argument('--enabled', choices=['yes', 'no'], help='set-enabled: new state of --service in the profile file')
    ap.add_argument('--json', action='store_true', help='list-services / check-service: machine readable output')
    ap.add_argument('--mock-fail-ids', help=argparse.SUPPRESS)
    ap.add_argument('-q', '--quiet', action='store_true')
    ap.add_argument('--version', action='version', version='ai_bridge %s (%s)' % (__version__, CONTRACT_VERSION))
    return ap


def _parse_params(items):
    out = {}
    for it in items:
        if '=' not in it:
            raise BridgeError('--param %r: expected KEY=VALUE' % it)
        k, v = it.split('=', 1)
        try:
            out[k] = int(v)
        except ValueError:
            try:
                out[k] = float(v)
            except ValueError:
                out[k] = v
    return out


def _numbers(s):
    if not s:
        return None
    if s.startswith('@'):
        with open(s[1:], encoding='utf-8') as f:
            return [l.strip() for l in f if l.strip()]
    return [x.strip() for x in s.split(',') if x.strip()]


def cmd_list(a, log):
    lst, path, err = profmod.all_services(a.services_file)
    rows = []
    for p in lst:
        env, st = profmod.env_status(p)
        errs, warns = profmod.validate(profmod.with_defaults(p))
        tgt = p.get('base_url') or ((p.get('command') or [''])[0] if p.get('transport') == 'local_command' else '') or \
            p.get('callable') or ''
        rows.append({'target': tgt, 'profile_sha256': profmod.profile_hash(profmod.with_defaults(p)),
                     'name': p.get('name'), 'task': p.get('task'), 'transport': p.get('transport'),
                     'enabled': bool(p.get('enabled', True)), 'template': bool(p.get('template', False)),
                     'network': p.get('transport') in profmod.NETWORK_TRANSPORTS,
                     'auth_env': env, 'auth_env_status': st, 'valid': not errs,
                     'errors': errs, 'warnings': warns, 'description': p.get('description', '')})
    if a.json:
        print(json.dumps({'profile_file': path, 'file_error': err, 'services': rows}, indent=1))
    else:
        print('# profile file: %s%s' % (path, '' if os.path.exists(path) else '  (does not exist)'))
        if err:
            print('# ERROR: ' + err)
        print('name\ttask\ttransport\tenabled\tvalid\tauth_env\tenv_status\ttarget\tprofile_sha256\tnote')
        for r in rows:
            note = 'TEMPLATE' if r['template'] else ('MOCK (FAKE VALUES)' if r['transport'] == 'mock' else '')
            print('\t'.join([str(r['name']), str(r['task']), str(r['transport']), 'yes' if r['enabled'] else 'no',
                             'yes' if r['valid'] else 'NO', r['auth_env'] or '-', r['auth_env_status'],
                             r['target'] or '-', r['profile_sha256'], note]))
    return EXIT_USAGE if err else EXIT_OK


def cmd_set_enabled(a, log):
    """Rewrite the profile file with one service enabled/disabled (atomic; a .bak copy is kept)."""
    if not a.service or not a.enabled:
        raise BridgeError('--service and --enabled yes|no are required')
    path = profmod.profile_path(a.services_file)
    lst, p, err = profmod.load_file(a.services_file)
    if err:
        raise BridgeError(err)
    with open(path, 'r', encoding='utf-8') as f:
        doc = json.load(f)
    svcs = doc['services'] if isinstance(doc, dict) else doc
    if isinstance(svcs, dict):
        raise BridgeError('%s uses the {"services": {name: ...}} form; edit it by hand' % path)
    for x in svcs:
        if x.get('name') == a.service:
            x['enabled'] = (a.enabled == 'yes')
            break
    else:
        raise BridgeError('service %r is not defined in %s (the built-in mock cannot be disabled)' % (a.service, path))
    shutil.copyfile(path, path + '.bak')
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(doc, f, indent=2, ensure_ascii=False)
        f.write('\n')
    os.replace(tmp, path)
    print('%s: enabled=%s (%s)' % (a.service, a.enabled, path))
    return EXIT_OK


def cmd_validate(a, log):
    if a.service:
        lst = [profmod.find(a.service, a.services_file)]
    else:
        lst, path, err = profmod.load_file(a.services_file)
        if err:
            print('ERROR: ' + err)
            return EXIT_USAGE
        if not lst:
            print('no profiles in %s' % path)
            return EXIT_USAGE
    bad = 0
    for p in lst:
        errs, warns = profmod.validate(profmod.with_defaults(p), a.task)
        if p.get('transport') == 'mock':
            errs = [e for e in errs if True]
        print('%s: %s' % (p.get('name'), 'OK' if not errs else 'INVALID'))
        for e in errs:
            print('  error: ' + e)
        for w in warns:
            print('  warning: ' + w)
        bad += bool(errs)
    return EXIT_USAGE if bad else EXIT_OK


def cmd_check(a, log):
    if not a.service:
        raise BridgeError('--service is required')
    p = profmod.find(a.service, a.services_file)
    pr = runner.effective_profile(p, a.task or ('generic' if p.get('task') == 'any' else p.get('task')))
    errs, warns = profmod.validate(pr)
    env, st = profmod.env_status(pr)
    info = {'service': pr['name'], 'transport': pr['transport'], 'valid': not errs, 'errors': errs, 'warnings': warns,
            'auth_env': env, 'auth_env_status': st, 'network': pr['transport'] in profmod.NETWORK_TRANSPORTS}
    ok = not errs
    if ok and info['network'] and not a.allow_network:
        info.update(connection='not tested (needs --allow-network)', ok=False)
        ok = False
    elif ok:
        try:
            ad = adapters.make_adapter(pr, {'params': {}, 'bands': []})
            good, msg = ad.check()
        except Exception as e:                                          # noqa
            good, msg = False, '%s: %s' % (type(e).__name__, e)
        info.update(connection=msg, ok=good)
        ok = good
    else:
        info.update(connection='not tested (profile invalid)', ok=False)
    if a.json:
        print(json.dumps(info, indent=1))
    else:
        print('service %s (%s): %s' % (info['service'], info['transport'], 'OK' if ok else 'NOT OK'))
        for e in errs:
            print('  error: ' + e)
        for w in warns:
            print('  warning: ' + w)
        if env:
            print('  auth env var %s: %s' % (env, st))
        print('  ' + info['connection'])
    return EXIT_OK if ok else EXIT_USAGE


def cmd_run(a, log, dry):
    for k in ('service', 'task', 'catalog'):
        if not getattr(a, k):
            raise BridgeError('--%s is required for this mode' % k)
    if not os.path.exists(a.catalog):
        raise BridgeError('catalog not found: %s' % a.catalog)
    images = [x for x in (a.image or '').split(',') if x]
    for p in images:
        if not os.path.exists(p):
            raise BridgeError('image not found: %s' % p)
    prof = profmod.find(a.service, a.services_file)
    if not prof.get('enabled', True) and not dry:
        raise BridgeError('service %r is disabled in the profile file' % a.service)
    opts = dict(numbers=_numbers(a.numbers), id_col=a.id_col, mag_columns=(a.mag_columns.split(',') if a.mag_columns else None),
                size_pix=a.size_pix, size_arcsec=a.size_arcsec, normalize=a.normalize, cutout_format=a.cutout_format,
                cutout_dir=a.cutout_dir, bands=(a.bands.split(',') if a.bands else None), pixel_scale=a.pixel_scale,
                psf_fwhm=a.psf_fwhm, params=_parse_params(a.param), dry_run=dry, no_cache=a.no_cache, cache_dir=a.cache_dir,
                resume=a.resume, allow_network=a.allow_network, max_objects=a.max_objects,
                mock_fail_ids=(a.mock_fail_ids.split(',') if a.mock_fail_ids else None))
    try:
        res = runner.run(prof, a.task, a.catalog, images, a.output, opts, log)
        if not dry and a.rename:
            runner.apply_rename(res, dict(x.split('=', 1) for x in a.rename if '=' in x))
    finally:
        pass
    tmpd = opts.get('_tmp_cutout_dir')
    try:
        if dry:
            print(json.dumps({'service': res['summary']['service'], 'task': a.task, 'objects': res['summary']['objects'],
                              'requests_that_would_be_sent': res['summary']['requests'], 'NOTE': 'dry run: nothing was sent',
                              'requests': res['previews'][:3] if not a.json else res['previews']}, indent=1, default=str))
            return EXIT_OK
        args_record = {k: v for k, v in vars(a).items() if v not in (None, [], False)}
        prov, text = runner.write_outputs(res, a.output, args_record, a.catalog, images)
        pp = a.provenance or ((a.output + '.provenance.json') if a.output else None)
        if pp:
            os.makedirs(os.path.dirname(os.path.abspath(pp)), exist_ok=True)
            with open(pp, 'w', encoding='utf-8') as f:
                json.dump(prov, f, indent=1, sort_keys=True)
                f.write('\n')
        s = res['summary']
        log('%s/%s: %d objects, %d ok, %d failed; requests sent %d, cache hits %d, retries %d%s' % (
            s['service'], a.task, s['objects'], s['objects_ok'], s['objects_failed'], s['requests_sent'],
            s['cache_hits'], s['retries'], '' if not res['models'] else '; model(s) reported: ' + ', '.join(res['models'])))
        for f in res['failures'][:10]:
            log('  object %s failed: %s' % (f['id'], f['error']))
        if len(res['failures']) > 10:
            log('  ... %d more failures (see provenance)' % (len(res['failures']) - 10))
        if res['exit_code'] == 0 and s['objects_failed'] and a.strict:
            return EXIT_PARTIAL_STRICT
        if s['objects_ok'] == 0:
            log('ERROR: every object failed')
            return EXIT_TOTAL_FAILURE
        return EXIT_OK
    finally:
        if tmpd and os.path.isdir(tmpd) and not a.cutout_dir:
            shutil.rmtree(tmpd, ignore_errors=True)


def main(argv=None):
    a = build_parser().parse_args(argv)
    log = (lambda m: None) if a.quiet else (lambda m: print(m, file=sys.stderr, flush=True))
    try:
        if a.mode == 'list-services':
            return cmd_list(a, log)
        if a.mode == 'validate-profile':
            return cmd_validate(a, log)
        if a.mode == 'set-enabled':
            return cmd_set_enabled(a, log)
        if a.mode == 'check-service':
            return cmd_check(a, log)
        return cmd_run(a, log, dry=(a.mode == 'dry-run'))
    except ProfileError as e:
        print('ERROR: %s' % e, file=sys.stderr)
        return EXIT_USAGE
    except BridgeError as e:
        print('ERROR: %s' % e, file=sys.stderr)
        return EXIT_USAGE
    except (OSError, ValueError) as e:
        print('ERROR: %s: %s' % (type(e).__name__, e), file=sys.stderr)
        return EXIT_USAGE
    except KeyboardInterrupt:
        return 130
