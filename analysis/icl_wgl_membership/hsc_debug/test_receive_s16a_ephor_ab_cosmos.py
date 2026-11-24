from http.server import BaseHTTPRequestHandler, HTTPServer
import io
import json
from threading import Thread
from urllib.request import build_opener

from astropy.io import fits
import pytest

from receive_s16a_ephor_ab_cosmos import CHECKPOINT, NoRedirect, receive_file, verify_state


def data():
    stream = io.BytesIO()
    fits.BinTableHDU.from_columns([fits.Column(name='ID', format='K', array=[2**53+1])]).writeto(stream)
    return stream.getvalue()


@pytest.mark.parametrize('mode', ['valid', 'changed_etag', 'truncated', 'redirect'])
def test_one_shot_transport_and_payload_cap(tmp_path, mode):
    payload = data()
    requests = []

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            requests.append((self.path, self.headers.get('If-Match')))
            if mode == 'redirect':
                self.send_response(302)
                self.send_header('Location', '/unexpected')
                self.end_headers()
                return
            self.send_response(200)
            self.send_header('ETag', '"changed"' if mode == 'changed_etag' else '"pinned"')
            self.send_header('Content-Length', str(len(payload)))
            self.end_headers()
            self.wfile.write(payload[:100] if mode == 'truncated' else payload)

    server = HTTPServer(('127.0.0.1', 0), Handler)
    background = Thread(target=server.serve_forever, daemon=True)
    background.start()
    entry = dict(filename='sample.fits', url=f'http://127.0.0.1:{server.server_port}/sample',
                 content_length=len(payload), etag='"pinned"')
    state = dict(plan_sha256='fixture', reserved_bytes=0, files={})
    plan = dict(plan_sha256='fixture', files=[entry], maximum_total_payload_bytes=len(payload))
    opener = build_opener(NoRedirect())
    try:
        if mode == 'valid':
            receive_file(entry, plan, tmp_path, state, opener.open, 'Basic fixture')
            verify_state(plan, tmp_path, state)
            assert state['files']['sample.fits']['tables'][0]['columns'] == ['ID']
            with fits.open(tmp_path/'sample.fits') as hdus:
                assert hdus[1].data['ID'][0] == 2**53+1
            with (tmp_path/'sample.fits').open('ab') as stream:
                stream.write(b'!')
            with pytest.raises(ValueError, match='changed'):
                verify_state(plan, tmp_path, state)
        else:
            with pytest.raises(ValueError):
                receive_file(entry, plan, tmp_path, state, opener.open, 'Basic fixture')
            reloaded = json.loads((tmp_path/CHECKPOINT).read_text())
            with pytest.raises(ValueError, match='did not complete'):
                verify_state(plan, tmp_path, reloaded)
            assert not (tmp_path/'sample.fits').exists()
            if mode == 'truncated':
                assert (tmp_path/'sample.fits.part').stat().st_size == 100
        with pytest.raises(ValueError, match='budget'):
            receive_file({**entry, 'filename':'extra.fits'}, plan, tmp_path, state, opener.open, 'Basic fixture')
        assert requests == [('/sample', '"pinned"')]
        assert 'Basic fixture' not in (tmp_path/CHECKPOINT).read_text()
        assert state['reserved_bytes'] == len(payload)
    finally:
        server.shutdown()
        background.join()
        server.server_close()
