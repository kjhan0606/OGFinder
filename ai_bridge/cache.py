"""On-disk response cache, keyed by a hash of (service profile, request) - never by credentials."""
import hashlib
import json
import os
import tempfile

DEFAULT_DIR = os.path.join('~', '.ds9', 'ai_cache')


def cache_dir(path=None):
    return os.path.expanduser(path or os.environ.get('OGF_AI_CACHE_DIR') or DEFAULT_DIR)


def make_key(profile_hash, key_material):
    h = hashlib.sha256()
    h.update(profile_hash.encode())
    h.update(b'\0')
    h.update(json.dumps(key_material, sort_keys=True, separators=(',', ':'), ensure_ascii=True, default=str).encode())
    return h.hexdigest()


class Cache:
    def __init__(self, base, service, enabled=True):
        self.dir = os.path.join(cache_dir(base), service)
        self.enabled = enabled
        self.hits = 0
        self.misses = 0
        self.stores = 0

    def _path(self, key):
        return os.path.join(self.dir, key[:2], key + '.json')

    def get(self, key):
        if not self.enabled:
            return None
        try:
            with open(self._path(key), 'r', encoding='utf-8') as f:
                d = json.load(f)
        except (OSError, ValueError):
            self.misses += 1
            return None
        if d.get('key') != key:
            self.misses += 1
            return None
        self.hits += 1
        return d

    def put(self, key, response, headers, when):
        if not self.enabled:
            return
        p = self._path(key)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        fd, tmp = tempfile.mkstemp(dir=os.path.dirname(p), suffix='.tmp')
        try:
            with os.fdopen(fd, 'w', encoding='utf-8') as f:
                json.dump({'key': key, 'time': when, 'headers': headers, 'response': response}, f)
            os.replace(tmp, p)
            self.stores += 1
        except OSError:
            try:
                os.unlink(tmp)
            except OSError:
                pass
