"""Transfer-completeness, polling and redaction checks for the 2026-11-25 S16A receivers.

Kept separate so the executed receivers (and their pinned code hashes) stay unchanged.
budgeted_download refuses a short body only when Content-Length is sent, so a response must declare
either Content-Length or chunked transfer (http.client raises IncompleteRead on a short chunked body),
and the stored CSV must end with a line break.
"""
from __future__ import annotations

import html
import json
from pathlib import Path
import urllib.parse

from prepare_s23b_s16a_quality_supplement import require

MAX_POLL_SECONDS = 6*3600


def transport(response):
    """Which length mechanism http.client will actually use; refuses a body read until close.
    Real HTTPResponse objects expose .chunked and .length; test doubles fall back to the headers."""
    if hasattr(response, 'chunked') and hasattr(response, 'length'):
        if response.chunked:
            return 'chunked'
        require(response.length is not None,
                'CAS download declared neither Content-Length nor chunked transfer; completeness cannot be checked')
        return 'content-length'
    headers = response.headers
    if headers.get('Content-Length') is not None:
        return 'content-length'
    require((headers.get('Transfer-Encoding') or '').strip().lower() == 'chunked',
            'CAS download declared neither Content-Length nor chunked transfer; completeness cannot be checked')
    return 'chunked'


def require_line_end(path):
    with Path(path).open('rb') as stream:
        stream.seek(0, 2)
        require(stream.tell() > 0, 'empty CAS payload')
        stream.seek(-1, 2)
        require(stream.read(1) == b'\n', 'CAS payload does not end with a line break (possible truncation)')


def require_poll_budget(waited_seconds):
    require(waited_seconds < MAX_POLL_SECONDS, 'CAS job not done within the polling budget; inspect it before any retry')


def redact(message, client):
    """As the shared receiver redaction, plus quote-free HTML, decimal-entity and URL-encoded forms."""
    for secret in client.credential.values():
        if secret:
            forms = {secret, secret.replace('/', '\\/'), json.dumps(secret)[1:-1], json.dumps(secret, ensure_ascii=False)[1:-1],
                     html.escape(secret), html.escape(secret, quote=False), html.escape(secret).replace('&#x27;', '&#39;'),
                     html.escape(secret).replace('&#x27;', '&#039;'), repr(secret)[1:-1],
                     ''.join(f'&#{ord(ch)};' for ch in secret), ''.join(f'&#x{ord(ch):x};' for ch in secret),
                     urllib.parse.quote(secret, safe=''), urllib.parse.quote(secret), urllib.parse.quote_plus(secret),
                     urllib.parse.quote(secret, safe='').lower(), html.escape(secret).replace('&quot;', '&#34;')}
            forms |= {f.lower() for f in forms if '%' in f}
            for form in sorted(forms, key=len, reverse=True):
                message = message.replace(form, '[redacted]')
    return message[:2000]
