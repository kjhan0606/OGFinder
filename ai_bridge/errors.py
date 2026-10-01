"""Exception types.  CLI exit codes: see cli.py."""


class BridgeError(Exception):
    """Configuration / usage error (exit code 1)."""


class ProfileError(BridgeError):
    """Invalid service profile."""


class TemplateError(BridgeError):
    """Request template references an unknown placeholder."""


class RequestFailed(Exception):
    """A request to a service failed.  `retryable` tells the retry loop whether to try again."""

    def __init__(self, message, retryable=False, status=None, retry_after=None, body=''):
        super().__init__(message)
        self.body = body
        self.retryable = retryable
        self.status = status
        self.retry_after = retry_after


class AuthError(RequestFailed):
    """The environment variable named in the profile is not set (never retried)."""

    def __init__(self, message):
        super().__init__(message, retryable=False)
