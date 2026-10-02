"""ogfkit: small shared helpers of the OGFinder analysis plugins (TSV catalogs, FITS/mask I/O, image models).

Everything here is a pure function on numpy arrays / plain Python containers with JSON-serialisable results, so
that the same code can run in a Phase II web worker (docs/phase2_web_design.md).  No Tk, no global state."""
__version__ = "1.0"
