# Browser runtime provenance and licence review

The app and worker bundles are original Web code, with no runtime npm dependency.
The worker loads separately self-hosted **unmodified Pyodide 0.27.7** assets, CPython
3.12.7, NumPy 2.0.2, SciPy 1.14.1 and OpenBLAS 0.3.26. The original `ds10core`
subset is built from this repository; nothing from OGFinder/ogfkit is imported or
packaged. The core ZIP is checked by SHA-256 before loading.

`tools/browser_runtime_lock.json` pins every distributed runtime file, wheel and
notice by hash. `tools/install_browser_runtime.py` verifies downloads; the Web
build verifies installed assets again. Runtime downloads by users stay on the
application origin. No CDN or additional Python package is requested at execution.
The release lock contains metadata about other packages, but their binaries are
neither installed nor distributed. Runtime package dependency closure is exactly
NumPy, SciPy and OpenBLAS (verified against the upstream lock).

The reviewed component inventory is included by `tools/license_check.py` in
`THIRD_PARTY_NOTICES.md`. Pyodide uses MPL-2.0, permitted by project policy for
unmodified files with notices and source availability. Its source and the build
patches for this exact release are available from
[Pyodide 0.27.7](https://github.com/pyodide/pyodide/tree/0.27.7).
Our original application/core files retain their own licence.

CPython uses PSF and incorporated permissive notices; its licence documentation
and Expat/HACL notices accompany the assets. Emscripten, musl, libc++/libc++abi,
libffi and hiwire carry MIT/BSD/Apache licences and the LLVM exception. The runtime
links zlib and bzip2; the known permissive bzip2 licence is retained verbatim.
OpenBLAS uses BSD and includes LAPACK and f2c notices. Pyodide builds Fortran paths
via f2c; it does not distribute a gfortran shared runtime in the selected OpenBLAS
asset (that archive contains only `libopenblas.so`).

NumPy and SciPy wheel notices are retained, including random generator, pocketfft,
uarray, ARPACK and Qhull copyrights. Additional SciPy bundled licence files are
extracted from its hash-pinned upstream source archive **as notices only**, without
shipping the source archive or third-party implementation code. This covers
SuperLU, PROPACK, HiGHS, COBYQA, fast_float, ryu, TNC, Boost and other bundled
components. Qhull has a known permissive BSD-like licence with notice/source
attribution obligations; its full text accompanies the runtime. SciPy's UNU.RAN
notice explicitly records BSD permission and exclusion of the separately licensed
uniform generator; this permission, rather than the general upstream UNU.RAN
licence, is the basis for its inclusion.

The installer creates `/browser-runtime/THIRD_PARTY.txt`, containing full notices
and versioned source links. A visible link is provided in local mode. Deploy the
complete runtime directory and its notices together; do not publish only the WASM
binary. The installer may fetch the SciPy source archive when rebuilding missing
notices, but does not include that archive in the deployment.

Pinned assets are a reproducible baseline, not a claim that their older versions
are the latest. Updating a runtime requires reviewing the package closure, notices,
hashes and actual browser/native numerical tests again. Backend-only tools and
unreviewed binaries are never loaded into this worker.
