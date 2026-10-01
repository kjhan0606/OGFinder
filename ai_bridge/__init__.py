"""ai_bridge - provider-agnostic connection points for external AI / astronomy services.

This package contains NO models.  It turns catalog rows (+ optional image cutouts) into
requests for a user-configured service (HTTPS REST, multipart upload, a local command,
a python callable, or an IVOA TAP service), maps the response onto named catalog
columns and records provenance.  See docs/ai_services.md.
"""
__version__ = "0.1.0"
CONTRACT_VERSION = "ogf-ai-bridge/1"
