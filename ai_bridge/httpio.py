"""Minimal HTTPS client on urllib (stdlib only).  Redirects to another host are refused so that credentials can
never be forwarded to a third party.  Errors are classified as retryable / not retryable."""
import socket
import ssl
import urllib.error
import urllib.parse
import urllib.request

from .errors import RequestFailed

RETRYABLE_STATUS = {408, 425, 429, 500, 502, 503, 504}
MAX_RESPONSE_BYTES = 64 * 1024 * 1024


class _SameHostRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        a, b = urllib.parse.urlsplit(req.full_url), urllib.parse.urlsplit(newurl)
        if (a.scheme, a.hostname, a.port) != (b.scheme, b.hostname, b.port):
            raise urllib.error.HTTPError(req.full_url, code, 'redirect to a different host refused (%s)' % b.hostname,
                                         headers, fp)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _retry_after(headers):
    v = headers.get('Retry-After') if headers else None
    try:
        return float(v) if v is not None else None
    except ValueError:
        return None


def request(method, url, headers=None, body=None, timeout=30.0):
    """-> (status, response_headers(dict, lower-case keys), bytes). Raises RequestFailed."""
    req = urllib.request.Request(url, data=body, method=method, headers=dict(headers or {}))
    opener = urllib.request.build_opener(_SameHostRedirect)
    try:
        with opener.open(req, timeout=timeout) as r:
            data = r.read(MAX_RESPONSE_BYTES + 1)
            if len(data) > MAX_RESPONSE_BYTES:
                raise RequestFailed('response larger than %d bytes' % MAX_RESPONSE_BYTES, retryable=False)
            return r.status, {k.lower(): v for k, v in r.headers.items()}, data
    except urllib.error.HTTPError as e:
        try:
            full = e.read(4000).decode('utf-8', 'replace')
            snippet = full[:300].strip().replace('\n', ' ')
        except Exception:
            snippet = full = ''
        raise RequestFailed('HTTP %d %s%s' % (e.code, e.reason, (': ' + snippet) if snippet else ''),
                            retryable=e.code in RETRYABLE_STATUS, status=e.code, retry_after=_retry_after(e.headers), body=full)
    except (socket.timeout, TimeoutError):
        raise RequestFailed('timeout after %ss' % timeout, retryable=True)
    except urllib.error.URLError as e:
        r = e.reason
        if isinstance(r, ssl.SSLCertVerificationError) or 'CERTIFICATE_VERIFY_FAILED' in str(r):
            raise RequestFailed('TLS certificate verification failed: %s' % r, retryable=False)
        if isinstance(r, (socket.timeout, TimeoutError)):
            raise RequestFailed('timeout after %ss' % timeout, retryable=True)
        raise RequestFailed('connection error: %s' % r, retryable=True)
    except (ConnectionError, OSError) as e:
        raise RequestFailed('connection error: %s' % e, retryable=True)


def multipart(fields, files):
    """fields: {name: str}; files: [(name, filename, content_type, bytes)] -> (content_type_header, body)."""
    import uuid
    b = '----ogfbridge' + uuid.uuid4().hex
    out = []
    for k, v in fields.items():
        out.append(('--%s\r\nContent-Disposition: form-data; name="%s"\r\n\r\n%s\r\n' % (b, k, v)).encode('utf-8'))
    for name, fn, ct, data in files:
        out.append(('--%s\r\nContent-Disposition: form-data; name="%s"; filename="%s"\r\nContent-Type: %s\r\n\r\n'
                    % (b, name, fn, ct)).encode('utf-8') + data + b'\r\n')
    out.append(('--%s--\r\n' % b).encode())
    return 'multipart/form-data; boundary=' + b, b''.join(out)
