#!/usr/bin/env python3
"""fake_sky_tap.py CORE_DIR RA DEC PORTFILE - a local fake TAP service for the GUI checks (no network): answers the catalogue layer's ADQL box
queries from a synthetic catalogue around (RA, DEC) using tests/catmock.py of the Astrafex Web checkout (CORE_DIR = .../server).  Writes the port to
PORTFILE once it listens; stop it with SIGTERM.  Point the catalogue layer at it with ASTRAFEX_CATALOG_ENDPOINTS='{"*": "http://127.0.0.1:PORT/tap"}'."""
import os
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

core, ra, dec, portfile = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
sys.path.insert(0, os.path.join(core, "tests"))
from catmock import FakeTap, synthetic  # noqa: E402

tap = FakeTap(cat=synthetic(n=600, ra0=ra, dec0=dec, half=0.012, seed=3))


class H(BaseHTTPRequestHandler):
    def _go(self, data):
        status, body, _ = tap(f"http://{self.headers.get('host')}{self.path}", data, dict(self.headers), 30)
        self.send_response(status)
        self.send_header("content-type", "text/plain")
        self.send_header("content-length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        self._go(None)

    def do_POST(self):
        self._go(self.rfile.read(int(self.headers.get("content-length") or 0)))

    def log_message(self, *a):
        pass


srv = ThreadingHTTPServer(("127.0.0.1", 0), H)
with open(portfile + ".tmp", "w") as fh:
    fh.write(str(srv.server_address[1]))
os.replace(portfile + ".tmp", portfile)
srv.serve_forever()
