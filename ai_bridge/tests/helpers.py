import http.server
import json
import os
import sys
import threading
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

CATALOG = ("NUMBER\tX_IMAGE\tY_IMAGE\tALPHA_J2000\tDELTA_J2000\tMAG_AUTO\tMAGERR_AUTO\tMAG_F160W\tMAGERR_F160W\n"
           "1\t100.5\t200.5\t202.4696\t47.1952\t20.1\t0.05\t20.0\t0.04\n"
           "2\t300\t300\t202.45\t47.20\t22.5\t0.2\t99\t99\n"
           "3\t50\t60\t202.50\t47.18\t18.3\t0.01\t18.2\t0.01\n"
           "4\t10\t10\t202.51\t47.17\t21.0\t0.1\t21.1\t0.1\n")


def write_catalog(d, text=CATALOG, name='cat.tsv'):
    p = os.path.join(d, name)
    with open(p, 'w') as f:
        f.write(text)
    return p


class Server:
    """Local test HTTP server.  `behaviour(handler, body_json) -> (status, dict|bytes, headers)`; records every request."""

    def __init__(self, behaviour):
        self.requests = []
        outer = self

        class H(http.server.BaseHTTPRequestHandler):
            def log_message(self, *a):
                pass

            def _do(self):
                n = int(self.headers.get('Content-Length') or 0)
                raw = self.rfile.read(n) if n else b''
                try:
                    js = json.loads(raw.decode()) if raw and 'json' in (self.headers.get('Content-Type') or '') else None
                except ValueError:
                    js = None
                rec = {'method': self.command, 'path': self.path, 'headers': {k: v for k, v in self.headers.items()},
                       'raw': raw, 'json': js}
                outer.requests.append(rec)
                st, body, hd = behaviour(self, rec)
                data = body if isinstance(body, bytes) else json.dumps(body).encode()
                self.send_response(st)
                for k, v in (hd or {}).items():
                    self.send_header(k, v)
                if 'Content-Type' not in (hd or {}):
                    self.send_header('Content-Type', 'application/json')
                self.send_header('Content-Length', str(len(data)))
                self.end_headers()
                self.wfile.write(data)
            do_GET = do_POST = do_PUT = _do

        self.httpd = http.server.ThreadingHTTPServer(('127.0.0.1', 0), H)
        self.port = self.httpd.server_address[1]
        self.url = 'http://127.0.0.1:%d' % self.port
        self.t = threading.Thread(target=self.httpd.serve_forever, daemon=True)
        self.t.start()

    def close(self):
        self.httpd.shutdown()
        self.httpd.server_close()


def rest_profile(url, **kw):
    p = {'name': 'testrest', 'task': 'photoz', 'transport': 'https_json', 'base_url': url, 'endpoint': '/z',
         'method': 'POST',
         'request': {'template': {'oid': '{id}', 'pos': {'ra': '{ra}', 'dec': '{dec}'}, 'm': '{mags}'}},
         'response': {'model_path': 'model', 'request_id_path': 'rid',
                      'fields': {'PHOTOZ': {'path': 'out.z', 'type': 'float'},
                                 'PHOTOZ_ERR': {'path': 'out.e', 'type': 'float'}}},
         'retries': 3, 'backoff_s': 0.5, 'backoff_factor': 2.0, 'timeout_s': 5}
    p.update(kw)
    return p
