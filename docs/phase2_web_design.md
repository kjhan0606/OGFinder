# Phase II - OGFinder on the web (design)

> **Naming note (2026-10):** the product family is now called Astrafex ("ds10-web" = Astrafex Web, the stand-alone OGFinder = Astrafex Desktop, package `ds10core` = Astrafex Core). This document is kept as written at the time; see Astrafex Web `docs/naming.md`.

Status: **design proposal, nothing of it is implemented.**  Revision 2 (2026-10-03): the UI decision is taken - **shared UI components plus
user-selectable layout presets** (sections 6.1, 6.4-6.7, 9, 10).  The three mock-ups behind it are in `/workspace/webdesign/` (`README.md`,
`proposal_A.html`, `proposal_B.html`, `proposal_C.html`, PNG screenshots); they are not part of this repository and use example data only.  Revision 4 (2026-10-03): the **default layout preset is Guided** (proposal C); Classic (A) and Workspace (B) are user-selectable alternatives (6.5, 6.6, 9, 10).  Revision 3 (2026-10-03): the product is **commercial** with a user-chosen compute mode (local flat subscription / server with tokens and limits, chapter 13), and every feature is developed for **both the stand-alone and the web version** (chapter 12).  Target repository: `github.com/kjhan0606/ds10-web` (new, to be created by
the owner; this document does not push anywhere).  Numbers marked *measured* were taken on the development box (8 cores, 16 GB) with the
test images in `/workspace/fits`; everything else is an estimate and says so.

## 1. Goal and non-goals

Phase I is OGFinder as a customized SAOImageDS9 8.7: a Tk GUI around Python/C back ends.  Phase II is a **custom-built web version**
(not a port of ds9): a browser viewer with the same workflow, in which

* **small files** are read and rendered **in the browser** (WebGL): no upload, nothing leaves the machine;
* **large files** are **uploaded** to a server that builds a tile pyramid and runs the heavy computations (extraction, deblending, mask,
  ICL/LSBG, moving-object pipeline, AI services);
* **both paths drive one viewer interface**: the viewer asks an *image source* for pixels and never knows where they come from;
* **one set of independent UI components, several layouts**: viewer, catalog table, parameter panel, object inspector, command palette,
  light-curve/tracklet dock, job list and status bar - plus, for the commercial model, the cost-estimate dialog, the account/usage page and the usage meter (13.11) - are built once and arranged by a *layout preset* that the user chooses
  (**Guided** = pipeline stepper, the default; Classic = ds9-like; Workspace = icon rail + drawer + dock).  A preset arranges components and owns no analysis logic (section 6);
* **one product, two shells, developed together**: the stand-alone (ds9/Tk) version and the web version share one compute core per tool, one plugin manifest and one test suite; a tool is not done until both have an entry (chapter 12);
* **the user chooses where a job runs**: locally (flat subscription, stand-alone-like, no metering) or on the servers (tokens with daily/weekly/monthly limits) - chapter 13.

Non-goals for Phase II: full ds9 feature parity (regions language, XPA/SAMP, colourmap editor, 3-D cubes beyond a slider), a
multi-tenant public service on day one (the product is commercial since revision 3, chapter 13, but day one is still one server deployment), replacing the desktop version (it stays the reference implementation).  Team sharing is a
later decision (section 9, risks).

## 2. What exists today and what it gives the web version

The desktop code is already split in the way a web version needs: **UI and orchestration in Tcl, computation in stand-alone Python
CLIs, a declarative description of every step.**  That is the main asset.

| existing piece | where | reuse in the web version |
|---|---|---|
| computation CLIs (argv in, TSV/FITS/JSON out) | `ds9/library/ds9_*.py` (26 scripts: `ds9_sextract` binary in `bin/`, `ds9_mask.py`, `ds9_icl.py`, `ds9_lsbg.py`, `ds9_photo_z.py`, `ds9_sed_fit.py`, `ds9_sersic.py`, `ds9_bulge_disk.py`, `ds9_morphometry.py`, `ds9_psf_phot.py`, `ds9_crowded_phot.py`, `ds9_psf_deconv.py`, `ds9_separate.py`, `ds9_ai_merge.py`, `ds9_catalog_edit.py`, `ds9_fits_export.py`, ...) | run unchanged as **jobs** in worker processes; no rewrite |
| Moving-object / transient pipeline | `moving/` (`align.py`, `detect.py`, `zogy.py`, `tracklet.py`, `orbit*.py`, `pipeline.py`, `mast.py`, ...) driven by `ds9/library/ds9_moving.py --mode setup\|fetch\|align\|difference\|link\|identify\|orbit\|transients\|lightcurve\|export` and `scripts/run_moving_pipeline.py` | server-side job type `moving`; `mast.py` already caches downloads |
| AI service bridge | `ai_bridge/` (`runner.py`, `adapters.py`, `contracts.py`, `cache.py`, `auth.py`, `cutouts.py`) | server-side as is; keys stay on the server |
| plugin manifests | `plugins/*/plugin.json` (18 plugins, 131 top-level steps plus 39 menu variants, 148 parameters: name, label, type `int/float/string/bool/choice/file`, default, min/max, group, help, `expert`; steps with `cli` argv templates and `output.mode`) | **the API description of every analysis step and the source of every settings dialog** (section 6) |
| declarative dialog | `ds9/library/ogf_dialog.tcl` (`OGFParamDialog`: groups, Basic/Expert, range validation, presets, Reset) | re-implemented once in the browser from the same JSON; ~250 lines of Tcl, no per-plugin UI code |
| step runner / job runner | `ogf_core.tcl` (`::ogf::step::run`, `::ogf::job::run`, argv templates `{python} {image} {catalog} {param}`) | server job queue with the same templates |
| session recorder + exporter | `ogf_session.tcl`, `ogf_session_template.py` (the exported Python script with `--mode pipeline\|replay`), `docs/session_python_script.md` | **pipeline export on the server unchanged**: the web session records the same steps (`step`, `class`, templated `argv_t`, `post`, `payload`) |
| table / tile / tabs / link / click chooser | `ogf_td.tcl`, `ogf_tile.tcl`, `ogf_link.tcl`, `ogf_pick.tcl`, `ogf_ui.tcl` | UI behaviour specification (what is linked to what), not code |
| catalog service | `::ogf::cat` in `ogf_core.tcl` (`tsv`, `rows`, `values`, `add_columns`, `set_cells`, filters, `registry` of keys) | becomes the **catalog API** (section 5.4); the key registry is the state list |
| review + report | `plugins/report/` (`report.py` already makes the HTML/PDF; REVIEW columns) | works unchanged as a server job; the review columns become cells edited through the catalog API |

Consequences: (1) the server is mostly a thin job scheduler around existing CLIs; (2) the web UI is mostly a renderer of manifests;
(3) the new engineering is the **image path** (tiles, WebGL, source abstraction), the **upload/storage** layer, and **accounts**.
The desktop version's coupling to Tk globals (`catpanel(...)`, see `docs/architecture.md` section 7) is irrelevant for the web version, which
takes the *key registry* and the manifests as the specification and does not run any of the Tcl.

## 3. Architecture overview

```
 browser                                                       server (one host first, containers later)
 +-------------------------------------------+   HTTPS/WSS   +------------------------------------------------+
 | Viewer (TypeScript)                        |  REST + SSE   | API gateway (FastAPI)                           |
 |  ImageSource  <- LocalFileSource (FITS.js)  +-------------->|  auth/session, upload, tiles, catalog, jobs,    |
 |               <- RemoteTileSource (/tiles)  |              |  masks, sessions                                |
 |  WebGL renderer (scale, colormap, tiles)    |              +-------+--------------------+-------------------+
 |  UI components + store + presets (sect. 6)  |                      |                    |
 |  Dialog renderer (plugin.json -> form)      |              +-------v------+     +-------v-----------------+
 |  Local catalog (small files: WASM/Worker)   |              | Object store |     | Job queue + workers      |
 +-------------------------------------------+              | (originals,   |     | (the existing ds9_*.py,   |
                                                              |  tiles, masks,|     |  ds9_sextract, moving/,   |
                                                              |  catalogs)    |     |  ai_bridge/)              |
                                                              +---------------+     +---------------------------+
                                                              Postgres (users, files, jobs, sessions, quotas)
```

## 4. The image-source abstraction (one viewer, two paths)

The viewer depends only on this interface (TypeScript sketch):

```ts
interface ImageSource {
  readonly id: string;
  readonly width: number; readonly height: number;       // full-resolution pixels
  readonly bitpix: number; readonly nHdu: number;
  readonly levels: number;                                // 0 = full resolution
  header(hdu?: number): Promise<FitsHeader>;              // WCS keywords for coordinates
  // a tile of one level: Float32Array (raw) or Uint8 (pre-stretched) - the renderer accepts both
  tile(level: number, tx: number, ty: number, signal?: AbortSignal): Promise<Tile>;
  stats(region?: Rect): Promise<Stats>;                   // min/max/zscale/percentiles for the scale dialog
  pixel(x: number, y: number): Promise<number | null>;    // pixel table / cursor value
  mask?: MaskSource;                                      // optional bit-flag mask, same tiling
}
```

Two implementations:

* **`LocalFileSource`** - the file is read in the browser (`File`/`FileReader`, FITS parser in a Web Worker: header cards, BITPIX,
  big-endian conversion, BSCALE/BZERO, tile-compressed `.fits.fz` via WASM cfitsio or a small Rice/GZIP decoder).  It cuts the image into the
  *same* tile geometry in memory and builds the lower levels with the same 2x2 mean in a Worker.  Nothing is uploaded.
* **`RemoteTileSource`** - tiles come from `GET /api/v1/files/{id}/tiles/{level}/{x}/{y}` (raw float32 or 16-bit, see 5.2); header and stats
  from the API.

The renderer (WebGL2) holds a texture cache of tiles, applies **scale (linear/log/sqrt/asinh/zscale), limits, colormap and contrast/bias in
the fragment shader** from float tiles (`R32F`, or `R16F` where float32 textures are not filterable), so changing the scale never re-fetches.
Overlays (catalog markers, ellipses, tracks, masks, regions) are a 2-D canvas/SVG layer above it, in image coordinates mapped through the
WCS (the same fk5-based drawing rule as `::ogf::td::draw_markers`).

**Why raw float tiles and not pre-stretched PNG tiles?**  PNG/JPEG tiles (what Leaflet/OpenSeadragon map viewers use) fix the stretch on
the server: every scale change needs new tiles and the pixel value is not available.  Float tiles keep ds9's interactive scale semantics.
Cost: bytes.  *Measured* on HUDF F160W (3600x3600 float32, 52 MB): an 8-bit PNG tile of 256x256 is 1.7 kB (2.2 ms to encode), 512x512 is
21 kB (7 ms); a raw float32 256x256 tile is 256 kB (about 60-120 kB gzip, estimate).  Therefore the tile endpoint offers **three
encodings** chosen by the client: `f32` (exact; default for the visible level on a good connection), `u16` (scaled by a per-tile min/max sent in
a header; 128 kB, half the bytes, 16 bit is enough to display and to do zscale), and `png8` (pre-stretched, for thumbnails, overview and slow
links).  The viewer can mix them (overview from `png8`, current level from `f32`).

Pyramid (*measured*): 4 reduced levels for 3600x3600, +33 % pixels, 0.09 s with a 2x2 mean in numpy.  Large survey images
(e.g. 20000x20000 float32 = 1.6 GB) take seconds to a minute and dominate by disk I/O, not by the pyramid.

## 5. API (v1)

REST/JSON over HTTPS; long operations return a job id and report progress by Server-Sent Events.  `{id}` is an opaque file id.

### 5.1 Upload (resumable chunks)

| method, path | purpose |
|---|---|
| `POST /api/v1/uploads` `{name, size, sha256?, kind: "image"\|"catalog"\|"mask"\|"other", chunk_size?}` -> `{upload_id, chunk_size, chunks, already}` | open an upload; if the sha256 is already stored for this user the response says `already: <file_id>` (de-duplication, no transfer) |
| `PUT /api/v1/uploads/{upload_id}/chunks/{n}` (raw body, `Content-Length`, `X-Chunk-SHA256`) | one chunk (default 8 MB); idempotent: re-sending a stored chunk is a no-op |
| `GET /api/v1/uploads/{upload_id}` -> `{received:[...], size}` | resume after a broken connection: the client sends only the missing chunks |
| `POST /api/v1/uploads/{upload_id}/complete` -> `{file_id, job_id}` | assemble, verify the whole-file sha256, then start the **ingest job** (validate FITS, read header, build the pyramid, compute stats) |
| `DELETE /api/v1/uploads/{upload_id}` | abort and free the partial data |

The protocol is deliberately close to tus.io (so an existing client such as `tus-js-client` and a server library can be used instead of writing it);
section 8 lists that as an option.

### 5.2 Files and tiles

| method, path | purpose |
|---|---|
| `GET /api/v1/files` / `GET /api/v1/files/{id}` | list / metadata: size, sha256, state `uploading\|ingesting\|ready\|failed`, dimensions, `levels`, HDU list, `tile_size`, owner, expiry date (section 7) |
| `GET /api/v1/files/{id}/header?hdu=0` | the FITS header as JSON cards (WCS, instrument, filter, `EXPTIME`, ...) |
| `GET /api/v1/files/{id}/tiles/{level}/{x}/{y}?enc=f32\|u16\|png8&hdu=0` | one tile (256x256 default; 512 optional); `ETag` = file sha256 + level/x/y, `Cache-Control: immutable`, edge tiles are padded with NaN |
| `GET /api/v1/files/{id}/stats?region=x0,y0,x1,y1&level=` | min/max/mean/median/rms, percentiles, zscale limits, histogram (feeds the scale dialog) |
| `GET /api/v1/files/{id}/pixel?x=&y=` and `.../cutout?x=&y=&size=&format=fits\|png\|f32` | cursor value; cut-outs (the same job `plugins/report/report.py` and `ai_bridge/cutouts.py` do) |
| `GET /api/v1/files/{id}/download` / `DELETE /api/v1/files/{id}` | original file (range requests) / delete the file and everything derived from it |

### 5.3 Jobs (analysis steps)

A job is "run step S of plugin P with parameters X on these inputs".  The request is derived from the **plugin manifest**: the server loads
the same `plugin.json` files, validates parameters with the same rules as `::ogf::params::validate`, expands the `cli` template (`{python}`,
`{image}`, `{catalog}`, `{work}`, `{param}`, `{"if": ...}` as in `::ogf::step::build_argv`) and runs the resulting argv in a worker.

| method, path | purpose |
|---|---|
| `GET /api/v1/plugins` / `GET /api/v1/plugins/{id}` | manifests (JSON as they are, minus server-internal fields): tabs, steps, parameters, `needs`, `session` class, `output` mode.  Steps without a `cli` block are listed with `server: false` until converted (section 6.3) |
| `POST /api/v1/jobs` `{plugin, step, inputs:{image:file_id, catalog:cat_id, ...}, params:{...}, session_id}` -> `{job_id}` | queue a job; the response is 422 with the validation message of the failing parameter |
| `GET /api/v1/jobs/{id}` -> `{state: queued\|running\|done\|failed\|cancelled, progress, started, finished, argv, exit_code, stderr_tail, outputs:[...]}` | status; `argv` is the exact command (what the recorder stores) |
| `GET /api/v1/jobs/{id}/events` (SSE) | progress and log lines (the progress strip of `ogf_ui.tcl`: stage + elapsed) |
| `POST /api/v1/jobs/{id}/cancel` | the Stop button (`::ogf::job::cancel`) |
| `GET /api/v1/jobs/{id}/outputs/{name}` | result files (TSV, FITS, PNG, JSON, HTML report) |

The `output.mode` of a step decides what the server does with stdout, exactly as `::ogf::step::done` does today: `add_columns` merges the columns
into the catalog (5.4), `set` replaces it, `text` returns text to show.  Steps with `"proc"` (Tcl-only steps wrapping UI code) need a
decision per step (section 6.3).

### 5.4 Catalogs

The catalog stays what it is today - a TSV-shaped table keyed by `NUMBER` - but is a server resource, versioned (every change gets a
revision number and a recorded step), so that undo, session export and the review feature work.

| method, path | purpose |
|---|---|
| `POST /api/v1/catalogs` `{from_file?: file_id, image: file_id}` / `GET /api/v1/catalogs/{id}` | create / metadata: columns, row count, revision |
| `GET /api/v1/catalogs/{id}/rows?cols=&sort=&filter=&offset=&limit=` | the table view for the virtual-scrolling grid; `filter` is the same model as `::ogf::cat::filter_set` (column in values) plus text search |
| `GET /api/v1/catalogs/{id}/tsv`, `.../csv`, `.../fits` | exports (`ds9_fits_export.py`) |
| `PATCH /api/v1/catalogs/{id}/cells` `{changes: {NUMBER: {COLUMN: value}}}` | `::ogf::cat::set_cells`: review state, notes |
| `POST /api/v1/catalogs/{id}/edit` `{op: delete\|merge\|separate\|add_object\|sort\|trim, ...}` | the operations of `ds9_catalog_edit.py apply` / `plugins/objects` |
| `GET /api/v1/catalogs/{id}/revisions` | revision list (step, time, who) |

Markers are drawn client-side from the rows (X_IMAGE/Y_IMAGE, A/B/THETA); at about 10^5 sources the client asks for a spatial subset
(`bbox=`), the same "visible only" idea as `CatalogPanelShowVisible`.

### 5.5 Masks

One shared bit-flag mask per image, as in `ds9_mask.py` / `plugins/mask`.

| method, path | purpose |
|---|---|
| `GET /api/v1/files/{id}/mask` / `.../mask/tiles/{level}/{x}/{y}` | mask metadata and tiles (uint8 bit planes, same tiling as the image; PNG-compressed) |
| `POST /api/v1/files/{id}/mask/ops` `{op: add\|erase\|grow\|shrink\|invert\|clear\|undo\|redo\|auto, regions?, params?}` | edits; regions come from the browser as polygons/ellipses in image coordinates (replaces the ds9 region file); `auto` is a job |
| `GET /api/v1/files/{id}/mask/download?format=fits\|bool` | export |

### 5.6 Sessions and pipeline export

| method, path | purpose |
|---|---|
| `POST /api/v1/sessions` / `GET /api/v1/sessions/{id}` | a session = the recorded step list (the `ogfsess(steps)` records: `seq, step, class, argv, argv_t, params, payload, post, requires, ms, failed, cat_in, cat_after, note`) |
| `POST /api/v1/sessions/{id}/steps` | informational steps written by the client (review, sort, save: `class: note`) |
| `GET /api/v1/sessions/{id}/export.py` | the **Python pipeline script**, rendered by the existing `ogf_session_template.py` with `@@DATA@@` filled from JSON (the Tcl `OGFSessExport` logic - key assignment, argv templating - is ported once to Python) |
| `GET /api/v1/sessions/{id}/log` | the text of `OGFSessLogText` |
| `POST /api/v1/sessions/{id}/report` | the `report` plugin step (HTML, optional PDF) |

Because the recorder already stores *templated* argv (`@{PY} @{LIB}/ds9_x.py @{IMG:main} ...`), a web session exports to the same script that
runs on a cluster or the desktop; `scripts/verify_session_replay.py`-style tests apply unchanged.

### 5.7 Other

`GET /api/v1/me`, `GET /api/v1/quota`, `GET /api/v1/health`; WebSocket is not needed (SSE is one-way and passes proxies more easily).

## 6. UI

### 6.1 Components (independent, shared by all presets)

Every visible part of the application is one of ten components (#1-#8 below; #9 estimate dialog and #10 account page belong to the billing chapter, 13.11).  A component takes **props = selectors on the state store (6.4)** and
emits **actions** (commands); it never imports another component, never calls another component's methods and never reads the DOM of a
neighbour.  The behaviour the desktop UI already defines (what is linked to what) is the specification; the Tk code is not reused.

| # | component | what it does (desktop origin) | data it needs (API, section 5) | store slices it reads / writes |
|---|---|---|---|---|
| 1 | **Image viewer** | tile pyramid + WebGL2 layers (section 4): pan/zoom, scale/colormap, cursor readout, WCS, catalog markers/ellipses, mask and trail overlays, **Single / Tile / split-compare** (several canvases or one canvas with a split uniform, shared view state; `ogf_tile.tcl`) | `ImageSource` (local or `/tiles`), `/header`, `/stats`, `/pixel`, mask tiles, catalog `bbox` rows | reads `view`, `catalog` (markers), `selection`; writes `view`, `selection` (marker click; the `ogf_pick.tcl` chooser is a popup of this component) |
| 2 | **Catalog table** | virtual-scrolling grid, kind filter (Galaxy / Moving / Transient / All), search, column filters, sort, multi-select, review column (`ogf_td.tcl`, `ogf_link.tcl`) | `/catalogs/{id}/rows`, `PATCH .../cells`, `.../edit` | reads `catalog`, `selection`, `view` (visible only); writes `catalog.filter/sort`, `selection`, cell edits |
| 3 | **Parameter panel** | form generated from `plugin.json` `params`: groups, Basic/Expert, range checks, presets, Reset (`OGFParamDialog`, 6.2).  Used as a modal dialog, a docked drawer or an inline "cell" - the same component | `/plugins/{id}`, `POST /jobs` | reads `plugins`, `params`; writes `params`, dispatches `job.submit`; also hosts the per-job **Run on: local / server** selector with the cost estimate (13.2, 13.11) |
| 4 | **Object inspector** | selected object: cut-out, Sersic model, residual, radial profile, p(z), catalog row, review buttons + note, provenance | `/cutout`, job outputs (`multifit`, `isophote`, `photoz_sed`), `PATCH .../cells` | reads `selection`, `catalog`; writes review cells |
| 5 | **Command palette (Ctrl-K)** | fuzzy list of every manifest step, view action, session action and "go to object #n"; arrow keys / Enter; same entries as the menus because both read one **command registry** | `/plugins` (registry built from the manifests at start-up) | reads `commands`; acts by dispatching the chosen command; `ui.paletteOpen` |
| 6 | **Light-curve / tracklet dock** | light curve of a transient (`lightcurves` viewer), tracklet and orbit summary of a moving object (`moving`), chosen from the table | job outputs of `moving` / `lightcurves`, catalog rows | reads `selection`, `jobs` (outputs) |
| 7 | **Job list** | running/queued/finished jobs with progress, Stop, exact `argv`, log tail (progress strip of `ogf_ui.tcl`) | `/jobs`, SSE `/jobs/{id}/events` | reads `jobs`, `session`; dispatches `job.cancel` |
| 8 | **Status bar** | WCS and cursor value, local/server badge, zoom, running job, session-recording indicator ("REC n steps"), tile/cache state, **usage meter** (day/week/month tokens, compute-mode badge; 13.11) | `/usage` (+ SSE) | reads `view`, `source`, `jobs`, `session`, `usage`, `account` |
| 9 | **Estimate and approval dialog** | pre-run cost estimate, hold, remaining allowance, alternatives (run locally, buy tokens); the over-limit variant (13.5, 13.7) | `POST /estimates`, `POST /jobs` | reads `usage`, `jobs`; dispatches `job.submit` |
| 10 | **Account, usage and settings page** | plan and licence, usage windows and ledger, top-ups, default compute mode, retention (13.11) | `/me/plan`, `/usage`, `/ledger`, `/topups` | reads `account`, `usage`; writes `ui.computeDefault` and account settings |

Rules that keep them independent (they are also the review checklist):

* A component renders from its selectors alone; the same props give the same DOM (testable without the other components).
* Each component declares its **minimum size** and what it does when it gets less (collapse to an icon, scroll, or hide a column group).  A preset that
  cannot give the minimum must not include the component.
* High-rate state (cursor position, pan during a drag, slider drag) lives **inside** the viewer/slider and is committed to the store on
  pointer-up or throttled to 10 Hz; the cursor readout reaches the status bar through a separate lightweight event channel, not through the store.
* All user actions - menu item, ribbon button, palette entry, keyboard shortcut, table context menu - are **commands** in one registry
  (`id`, title, optional Korean title, manifest step or built-in, enabled predicate).  Presets only decide where a command is shown.
* UI strings are English with an optional Korean gloss (as in the mock-ups: *Detect / 탐지, Classify / 분류, Measure / 측정, Review / 검토, Report / 보고서*);
  strings live in one table per locale and widths are never hard-coded to English text (Korean-text check in 6.7).

### 6.2 Declarative dialogs

`plugin.json` `params` -> a form component: groups (fieldsets or tabs when `dialog_tabs`), Basic/Expert toggle, per-field help line, range checks,
presets (stored per user instead of `~/.ds9/ogf_params/<id>.presets`), Reset to defaults.  A choice with `choices_proc` (AI services) becomes
a request to `/api/v1/plugins/{id}/choices/{name}`.  Dialog values are sent as `params` of a job; they are validated on the client for
feedback and again on the server for safety.

### 6.3 Mapping plugins to the web

Classification of the manifest steps.  Counted on the current tree: of 131 top-level steps **only 1 has a `cli` template** (`example_hello`); **120 are Tcl `proc` steps** that
wrap the old GUI procedures (which build argv lists in Tcl and call the `ds9_*.py` scripts); the rest are variant cascades and toggles.  So the manifests describe *what* exists and its parameters,
but for almost every analysis step they do **not yet** contain the command line.  The rule per kind:

| kind of step | examples | web handling |
|---|---|---|
| `cli` step with `output.mode` add_columns / set / text | `example_hello`, the AI services, `extract` | generic job (5.3), no custom code |
| `proc` step that opens an `OGFParamDialog` and then calls a script | most of `icl`, `lsbg`, `star_psf`, `morphology`, `photoz_sed` | the dialog is generic; the script call needs a `cli` template - **the first task of Milestone 3 is to write the missing `cli` blocks** (the argv lists exist only as Tcl code in `plugins/*/*.tcl`, e.g. `CatalogPanelSEDFitRun`; `scripts/verify_cat_behavior.tcl` already captures the exec argv of ~50 of them as a golden file, which is the acceptance test for each conversion) |
| `proc` step that edits the table/markers interactively | `objects` (delete, merge, separate, add object, AI-merge accept/reject), `bands` registry, `mask` region editing | client behaviour + catalog `edit` API (5.4) + mask API (5.5) |
| display-only | Mark All, Show Visible Only, overlay toggle | client only |
| `moving` | fetch / align / difference / link / identify / orbit / lightcurve | `ds9_moving.py --mode ...` as jobs; needs network access to MAST/JPL from the **server** |

Modules that are Python already (`ds9/library/ds9_*.py`, `moving/`, `ai_bridge/`) need no change except a `--workdir`/path discipline (they read
`~/.ds9`; on the server each job gets a private work directory via `HOME`/`--workdir`).

### 6.4 State store (single source of truth)

One store, one direction of data flow: `component -> command/action -> reducer (or job/API call) -> store -> selectors -> components`.
The store is framework-agnostic TypeScript (a small vanilla store such as Zustand-vanilla/nanostores, or about 200 lines of our own; the choice
is closed in M0) holding plain serialisable data only (no DOM nodes, no WebGL objects), so a state snapshot can be logged, replayed in tests and
attached to a bug report.

```ts
interface AppState {
  source:   { id: string; kind: "local" | "remote"; bands: Band[]; detectionBand: string; hdu: number };   // section 4
  view:     { center: [number, number]; zoom: number; rotation: number;
              scale: { type: "linear"|"log"|"sqrt"|"asinh"|"zscale"; lo: number; hi: number }; colormap: string;
              layers: Record<string, { visible: boolean; opacity: number }>;                                // per band / overlay
              compare: { mode: "off"|"split"|"blink"|"tile"; split: number; left: string; right: string };
              lockWcs: boolean };
  catalog:  { id: string | null; revision: number; kind: "all"|"galaxy"|"moving"|"transient";
              filter: { text: string; magMax: number | null; review: string | null; visibleOnly: boolean };
              sort: { col: string; asc: boolean }; showMarkers: boolean };
  selection:{ primary: number | null; set: number[] };                                                    // NUMBER keys
  plugins:  Record<string, Manifest>;  params: Record<string, Record<string, unknown>>;
  jobs:     { byId: Record<string, Job>; order: string[] };  session: { steps: SessionStep[]; recording: boolean };   // Job gains mode: "local"|"server", estimate, hold (ch. 13)
  account:  { plan: Plan; entitlements: Entitlements; license: LicenseState };                              // written by server events only (ch. 13)
  usage:    { windows: Record<"day"|"week"|"month", { limit: string; used: string; held: string; resetsAt: string }>; topupBalance: string };
  ui:       { preset: "guided"|"classic"|"workspace";  /* default "guided" */ theme: "dark"|"light"; locale: "en"|"ko"; computeDefault: "local"|"server"|"ask";
              paletteOpen: boolean; panels: Record<string, PanelState> };                                   // per-preset layout, see 6.6
}
```

Ownership of each slice (who may write it) is fixed, so two components cannot fight over a value:

| slice | written by | notes |
|---|---|---|
| `source`, `view` | viewer, view commands, palette ("Compare: split") | the same view state drives N canvases in Tile / split mode; `source` may change local -> remote ("Analyse on server") without touching `view` (section 7.1) |
| `catalog`, `selection` | table, viewer (click), inspector (review), palette ("go to object") | server catalog revisions are the truth for cell values; the store keeps the revision number and re-fetches on change |
| `plugins`, `params` | parameter panel | defaults from the manifest; per-user presets in the account (section 6.2) |
| `jobs`, `session` | job events (SSE) | components only dispatch `job.submit` / `job.cancel` |
| `ui` | preset shell, palette, theme/locale controls, account settings (`computeDefault`) | the only slice persisted in the browser (6.6) |
| `account`, `usage` | server events (SSE `/usage/events`) only; no component writes them | amounts are decimal strings (micro-tokens), never floats; the UI shows, the server decides (13.6, 13.9) |

Switching a preset changes **only `ui.preset`**: the same store feeds a different arrangement, so selection, view, scale, catalog filter,
running jobs and unsaved parameter edits survive a switch (acceptance test in 6.7).

### 6.5 Layout presets

A preset is a shell: a CSS-grid template, the list of components it mounts with their slots and minimum sizes, and a few **shell widgets**
that exist only to navigate (built from manifests and commands, holding no analysis logic).  Presets are the only place where components
are arranged; none of them forks a component.  The three mock-ups in `/workspace/webdesign/` illustrate them.

| | **Guided** (`proposal_C`) - **default** | **Classic** (`proposal_A`) | **Workspace** (`proposal_B`) |
|---|---|---|---|
| idea | a pipeline you can follow and replay: *what do I do next?* | ds9 and OGFinder Phase I, almost no re-learning | modes + tool drawer + dock for multi-band, time-domain and tablets |
| for whom | first-time users, routine review and reporting, teaching (the default for everyone until they choose otherwise) | **ds9 veterans** and dense catalog work; one click away (Layout control, palette, first-run chooser 6.6) | multi-band / time-domain work, tablets |
| default theme | light | dark | dark, light toggle |
| shell widgets | left **stepper** (Load, Calibrate, Detect, Measure, Review, Report) with sub-steps, status and progress from `session.steps`/`jobs`; step **"cell" header** (title, parameter chips, Run); result strip (counts, histogram); recorded **session log**; split-compare toolbar | menu bar (File ... Analysis ... Window); workflow tabs *Detect / Classify / Measure / Low-SB / Time-domain / Results* with a **plugin ribbon** (one chip per plugin: primary step, step list, gear = settings); info panel (File, Object, Value, fk5 alpha/delta, x, y, Zoom); button bar; frame tabs (Single/Tile/Blink); colorbar | left **icon rail** of modes (Detect, Morph, LSB/ICL, Moving, Photo-z, Trails, AI, Stack); **collapsible drawer** (steps of the mode + parameter form); top bar (file chip, Local / Analyse-on-server switch, theme, language); floating viewer toolbar |
| components | viewer with split compare (1), inspector as right column (4), parameter panel inline in the cell (3), palette (5, always), job list inside the stepper (7), status/log (8); **catalog table (2) as a resizable bottom sheet** opened from the result strip or key `T` (open by default once a catalog exists, pinnable) - the table must never be more than one action away in the default preset | viewer (1), catalog table (2) with Galaxy/Moving/Transient/All tabs, parameter panel as a dialog (3), a reduced inspector as selected-object summary (4; the full inspector can be docked, since it already exists), palette (5, optional), job list + status in the status bar (7, 8) | viewer (1), parameter panel in the drawer (3), table (2) and light curve/tracklet/jobs (6, 7) in the bottom **dock**, status bar (8), palette (5) |
| manifest fields used | an optional **`ui.stage`** (`load\|calibrate\|detect\|measure\|review\|report`); default derived from `tab` (Detect -> detect; trails and mask -> calibrate; Classify, Measure, Low-SB, Time-domain -> measure; Results -> review/report) | `tab`, `order`, `short`, `primary`, `steps`, `params` (already present) | optional **`ui.mode`** and **`ui.icon`**; default: the plugin grouping of the mock-up, the rest under a "More..." flyout |
| minimum window | 1024 x 700 (panels collapse to drawers below 1280) | 1280 x 720 | 1024 x 700 (rail only below 1100; **best tablet landscape**) |
| strengths | review, reporting, teaching; the log is the exported script; works on tablets | cheapest to build, ds9-familiar, best catalog density | layers, time-domain and jobs visible at once; touch-friendly |
| delivery | **Phase 1** (M4; the default and primary preset) | **Phase 2** (M9) | **Phase 3** (M10), after the viewer is stable (9.1) |

`ui.stage`, `ui.mode` and `ui.icon` are additive manifest fields: the desktop loader ignores unknown keys, and a plugin without them still
appears (derived defaults).  Plugins that fit no stage (`lensmodel`, `spectra`, `cluster`, `xmatch`, ...) are listed under *Other tools* in
Guided (the default preset) and are reachable through the palette in every preset.

**Deferred in Workspace v1** (backlog, not in the first Workspace release): RGB composite and per-band blend layers, the overview minimap, and free
docking (drag/split/float).  v1 uses one fixed, drag-resizable dock below the viewer and one band visible at a time.  Reason: these are the expensive
parts of proposal B (9.2) and neither of the other two presets needs them.

**Why Guided is the default** (decision of 2026-10-03): it is the only preset that serves new users, tablets and the reproducibility story (the visible log is the exported script) at once.  Its price is that dense catalog work is slower than in Classic, hence the table sheet rule above, the one-click switch to Classic and the first-run chooser (6.6).  Guided is the **primary preset**: the reference for acceptance tests (6.7) and the first one built (M4).

Preset switch: a "Layout: Guided | Classic | Workspace" control and a palette command ("Layout: Classic").  Switching re-mounts the shell, not the
viewer: the viewer's canvas element is *moved* into the new shell rather than re-created, so tiles stay in the texture cache without a re-fetch
(risk: browsers that lose a WebGL context on re-parenting, section 10).

### 6.6 Persistence of the layout choice (localStorage + URL parameter)

* **Key and value**: `localStorage["ogf.ui.v1"] = { preset, firstRunDone, theme, locale, layouts: { guided: {...}, classic: {...}, workspace: {...} } }`.  **Default value of `preset` is `"guided"`** (and of `theme` the preset's own default, 6.5).
  `layouts[name]` holds only that preset's panel state (panel widths, collapsed drawers, dock height, open stepper steps).
  Nothing else - no file names, ids, paths, catalog data or tokens - is written to this key (a mirror of `computeDefault` is allowed; local-mode results live in OPFS/IndexedDB, 13.8, and are not part of this key).
* **URL parameter**: `?layout=guided|classic|workspace` (also `&theme=dark|light`, `&lang=en|ko`).  Precedence at start-up: **URL parameter >
  localStorage > default (`guided`)** (precedence unchanged; only the default changed).  A URL value applies to that page load and is written to localStorage only when the user then switches the
  layout explicitly, so following a link never silently changes the stored preference.  The URL is updated with `history.replaceState` when the user
  switches, so a copied address reopens the same layout.  There is no per-file or per-session link in Phase II (section 7.2); these are the only
  parameters the application reads.
* **First-run chooser (첫 실행 선택)**: when no stored object exists and no `?layout=` is given, Guided is applied immediately (no blocking dialog) and a dismissible, non-blocking card "Choose your layout / 레이아웃 선택" is shown once, with a thumbnail of Guided (selected), **Classic ("ds9 style", for ds9 users)** and, once M10 ships, Workspace.  Choosing writes `preset` and `firstRunDone: true`; dismissing keeps Guided and also sets `firstRunDone`; the card can be reopened from the Layout control.  It is not shown when `?layout=` is present (a link already chose a layout; `firstRunDone` stays unset).  Whether to show the card at all is an owner option (10, open decisions).
* **Robustness**: unknown or invalid values (URL, stored preset, missing `layouts` entry) fall back to the default **`guided`** with one console warning; the stored object is schema-versioned (`v1`) and
  validated field by field (a bad `layouts.guided` resets only that preset); storage that is disabled, full or blocked (private mode, enterprise
  policy) falls back to an in-memory store (preset `guided`) with a one-line notice; a `storage` event from another tab updates the *stored* preference but does not
  re-layout an open tab.  A preset whose minimum window size (6.5) is not met is offered but flagged, and a compatible one is suggested; it is never forced.
* **Reset**: "Reset layout" clears `layouts[current]`; "Reset all UI settings" removes the key.  A server-side copy per account is a later option (section 10).

### 6.7 Verification plan (DOM/layout checks per preset)

The mock-ups were checked with a Playwright/headless-Chrome script (`/workspace/webdesign/_src/shot.py`, `extra_checks.py`, `overlap.py`); it is the
starting point of the CI job.  The same checks run **for every preset x every tested state**, at 1440x900 (all presets), 1280x720 (all),
and 1024x768 / 1180x820 tablet landscape (Guided and Workspace).  **Guided is the primary preset**: the full matrix runs on every pull request and its golden screenshots are the reference; Classic and Workspace run the full matrix nightly and at release, plus a smoke subset (layout, images, hygiene, state) on every pull request that touches shared components or their shell.  Already measured on the mock-ups at 1440x900: 0 overflowing or out-of-bounds elements,
all images loaded, no console errors, no external requests, no overlapping floating pieces.

| layer | check | applies to |
|---|---|---|
| layout | `document` scroll size equals the viewport (no page scroll); **no element overflows its box or lies outside the viewport**, except containers marked `data-scroll` / `data-clip` (table body, viewer, log, step list); no clipped text | all presets, all states |
| layout | no overlap between floating pieces (viewer toolbar, layer/selection cards, toast, readout, split handle, palette) | Guided, Workspace |
| layout | minimum-size contract of each component: render it at its declared minimum and assert the collapse behaviour | component tests |
| images | every viewer/cut-out image or texture loads (`naturalWidth > 0`, WebGL context created, first tile drawn); the cursor readout over the viewer returns RA/Dec and a pixel value | all presets |
| hygiene | no console errors/warnings; **no request to an external host** (only the application's own API); no CSP violations | all presets |
| state | the same scripted scenario in each preset ends with the same store snapshot for `selection`, `view`, `catalog.filter`, `params` (select an object in the table -> marker highlighted, inspector/summary updated, dock shows the light curve for a transient and the tracklet for a moving object) | all presets |
| switching | Guided -> Classic -> Workspace -> Guided with a selected object, a changed scale, a running job and an edited, unsaved parameter: all survive; the viewer canvas is the same element (identity test) and no tile is re-fetched | all |
| persistence | empty storage, no parameter: **Guided** is shown and the first-run card appears once (dismiss and choose Classic both tested; not shown again after reload); choose a preset, reload: same preset (localStorage); open `?layout=classic` with a stored `guided`: Classic is shown and the stored value is unchanged (no first-run card); invalid parameter or stored value -> **Guided** + one warning; blocked storage -> in-memory Guided, no exception; corrupted JSON -> reset to Guided | all (default case: Guided) |
| i18n | long Korean strings (breadcrumb, step names, tooltips) and a 30-character English parameter label do not break the checks above; Noto Sans CJK / Nanum Gothic / system fallback present | all, `lang=ko` |
| theme | dark and light: contrast of text and of the marker colours on the image; the viewer background stays dark | Workspace, Guided |
| billing UI | usage meter, run-on selector, estimate dialog and account page in each state of 13.11 (normal, near limit, over limit, estimate failed, payment past due, tool server-only) render without overflow; the over-limit dialog offers local mode only for local-capable tools; no price or token number is hard-coded in the DOM (all from the API) | all presets |
| parity | the same parity cases (12.3) pass through the UI path of each preset (Run on local and on server give the same result envelope) | all presets |
| accessibility | every command reachable by keyboard; palette opens with Ctrl-K / Cmd-K and closes with Esc; focus returns to the previous element; roles/labels on rail buttons and tabs | all |
| visual | golden screenshots per preset/state with a tolerance, regenerated deliberately (the images in `/workspace/webdesign/` are the first references) | all |
| performance (informational) | store-update-to-DOM time for a selection change, table scroll at 100k rows, preset switch time | all |

A preset is *done* when this matrix is green, its user-acceptance scenarios (the M4/M7 scenarios re-written for the preset) pass, and its DOM metrics are
recorded in the release notes.

## 7. Size threshold, storage, accounts

### 7.1 Where the switch between "browser" and "server" is

Rule: **local if the file fits comfortably in the browser's memory budget and no server computation is needed; otherwise upload.**  Proposed
defaults (all configurable, to be tuned by measurement in Milestone 1):

* `size <= 200 MB` (decoded float32 image <= ~400 MB with the pyramid) **and** a single 2-D HDU -> *local*.  Reason: a browser tab can hold about 1-2 GB
  of typed arrays on a desktop (less on mobile, and limited by the 32-bit/ArrayBuffer limits of older engines); 200 MB keeps the Worker
  well inside that, and covers HST drizzled stamps, DES tiles, JWST NIRCam cutouts; the HUDF F160W test image is 52 MB.
* anything larger, multi-extension (e.g. HST FLC with SCI/ERR/DQ x 2 chips), `.fits.fz` that expands beyond the limit, or any file for which the user starts an analysis
  step that needs the server (extraction on a mosaic of 20k x 20k, LSBG, moving pipeline) -> *upload*.
* **"Analyse on server" for a small file** is allowed (one click): the file is uploaded and the same viewer switches from `LocalFileSource` to `RemoteTileSource`
  without losing view state (the view state lives outside the source).  Small files are never uploaded implicitly.
* **Revision 3**: this size rule is now only the *default suggestion* of the compute-mode selector; the user chooses local or server per job, and the stand-alone version is the unmetered local mode (13.2, 12.5).  WASM extraction is no longer an option but a milestone item (M14/M15, 13.10).
* Computations on small local files that are cheap enough for the browser (statistics, cut-outs, simple aperture photometry) can run in a
  Worker; **extraction with `ds9_sextract` (C) runs on the server only**, unless a WASM build is made (option in section 8; the HUDF F160W full
  extraction takes 1.9 s natively, *measured*, so the WASM path is plausible but unproven).

### 7.2 Login, storage, retention and deletion - options

This is the area where the owner's decision is needed; the design supports all of them through a per-deployment policy file
(`policy.yaml`: `auth`, `quota_bytes`, `retention_days`, `max_upload_bytes`, `allow_anonymous`).

**Revision 3**: with the commercial model the quota, retention and storage numbers below become plan parameters (`plans.yaml`), and disk usage is metered (13.4, 13.8); the table is kept as the option space.

| aspect | option A (recommended start) | option B | option C |
|---|---|---|---|
| login | institutional/OIDC (Google, ORCID, university SSO) via an existing library; no passwords stored | local accounts (email + password/argon2, TOTP) | anonymous sessions with a random token (no account), upload limit small |
| storage | one object store per deployment: local disk first (`/data/{user}/{file}`), S3-compatible (MinIO/S3) when more than one host | S3 only | disk quota per user, no object store |
| quota | 20 GB per user, 5 GB per upload, 3 concurrent jobs | per project/group quota | 2 GB, 1 job |
| retention | originals: 30 days after last access, derived data (tiles, catalogs, masks, sessions): same; user can pin ("keep") up to the quota | keep until deleted (+ yearly confirmation) | 24 h |
| deletion | user deletes any file -> immediate removal of file + tiles + masks + job outputs, catalogs/sessions keep their record but not the pixels; account deletion removes everything within 7 days; a nightly purge enforces retention; deletion is logged without file names | soft delete 14 days (undo) then purge | purge on session end |
| privacy statement | files are private to the owner, never used for training, not shared; admins can see metadata (size, time), not pixels, unless the owner opens a support ticket | same | same |

Data of other people's observations is common in astronomy (proprietary periods): the default must be **private**, **no public links** in
Phase II, and the retention default short.  Because of this the design contains **no sharing endpoint**; any sharing is a separate design.

## 8. Technology choices and trade-offs

| layer | choice | alternatives and trade-off |
|---|---|---|
| front-end language | **TypeScript**, framework **Svelte or React** (small component set: dialogs from JSON, grid, tabs; the eight components of 6.1 are the unit of work) | plain JS: no type safety for the tile/WCS code; Vue equal; the decision matters little - the viewer is framework-free |
| UI state | one **vanilla store** (Zustand-vanilla / nanostores / ~200 lines own code) with serialisable slices (6.4) and a command registry; only the `ui` slice is persisted (6.6) | Redux Toolkit: more boilerplate, good devtools; framework-local state only: preset switching and tests would need each component's internals |
| layout presets | CSS-grid shells + `ResizeObserver`-driven collapse rules; panel resize with small own code | docking library (golden-layout/dockview) - **only for the later Workspace docking**; not needed for Classic, Guided or Workspace v1 |
| renderer | **WebGL2** (R32F/R16F textures, fragment-shader scale and colormap), canvas/SVG overlay | **WebGPU**: better compute but not universal yet (as of this writing Safari/Firefox support is partial; verify before starting); Canvas2D only: no interactive scale on large images; OpenSeadragon/Leaflet: mature tiling but 8-bit stretched tiles only, no pixel values |
| FITS in the browser | own worker parser for the subset needed (image HDUs, BSCALE/BZERO, WCS) + WASM cfitsio only for compressed tiles | `astro.js`/`fitsjs` libraries: quicker start, uncertain maintenance; **decision after a one-week spike** |
| WCS | `wcslib`/AST compiled to WASM, or a small TAN/SIP + lookup in the worker | AST is already built in this repo (`ast/`) which makes a WASM build plausible; a hand-written TAN/SIP is enough for most HST/JWST/ground images; non-standard projections go to the server |
| server framework | **FastAPI** (Python, same language as all back ends, async upload, SSE) | Django: heavier, better admin/auth; Node: duplicates Python code; Go: fast but no reuse |
| job queue | **Postgres-based queue** (e.g. `procrastinate`/`pg_boss`-style) or Redis+RQ/Celery | Celery: most features, more moving parts; Slurm/K8s Jobs: when the workload becomes cluster-sized; a simple subprocess pool is enough for one host |
| execution isolation | each job in a **container / bubblewrap sandbox** with CPU, memory, time limits and no network (except `moving`/`ai` jobs through an egress proxy) | plain subprocess: easy, but jobs run user-provided FITS through C and Python code (see risks) |
| storage | local disk + SQLite/Postgres metadata first, S3-compatible later | S3 from day one: portable but adds ops |
| tile format | raw `f32` / `u16` plus `png8`, **256 px** default | Zarr/COG-like chunked formats: standard, tooling, but FITS tile-compression and WCS handling are not their strength; HiPS/HiPS-like for all-sky surveys (a later option for survey browsing, not for user images) |
| resumable upload | **tus** protocol (existing server and client libraries) | custom chunk API (5.1) as specified: fewer dependencies in the browser, more code to test |
| auth | OIDC (via `authlib`) | passwords: more security work for us |
| deployment | one VM with docker compose (API, worker, Postgres, reverse proxy with TLS), then split | Kubernetes at the start: too much for a first version |
| testing | contract tests of the API against the plugin manifests; the existing harness (`scripts/run_all_checks.sh`) runs against the server's job runner; Playwright for the browser; golden images for the renderer | - |

## 9. Milestones, phases and cost

Effort numbers are estimates in person-weeks for one experienced developer who knows the code base; they are ranges, not commitments.
M0-M8 deliver the shared components and the **Guided** preset (the default and primary preset, built first); **M9** adds Classic (the second preset), **M10** Workspace.  "(+x)" is the extra cost of building with independent
components, a state store and preset plumbing instead of one fixed layout; it is already included in the range shown.

| M | content | deliverable and acceptance test | estimate (extra) |
|---|---|---|---|
| **M0** spikes | (1) FITS parse + WebGL2 render of the 52 MB HUDF image with scale/colormap, (2) tile endpoint over the existing numpy code, (3) WCS in WASM vs hand-written, (4) the 200 MB threshold on three browsers, (5) **state store + command registry skeleton and two throw-away shells** (Guided grid first, Classic grid) around a stub viewer, to confirm the component contract and the canvas hand-over of 6.5 | go/no-go notes with measurements; decisions in section 8 closed (store library included) | 2 (0) |
| **M1** viewer + local files | `ImageSource`, `LocalFileSource`, WebGL renderer, pan/zoom/scale/colormap, pixel readout, WCS coordinates, image tabs, Tile mode; **viewer and status bar as stand-alone components** with the store contract (6.1, 6.4); split-compare uniform (small; used by Guided from M4) | open m51 and the HUDF crop offline; scale changes at 60 fps on a 4k tile set; test against ds9 pixel values (exact); component tests at minimum size | 6-8 (+1) |
| **M2** server core | auth (OIDC), resumable upload (5.1), ingest job, pyramid, tile API (5.2), `RemoteTileSource`, quotas, retention purge, delete | upload a 1.6 GB FITS over a throttled connection with a forced disconnect; resume; tiles identical to the local path for the same file | 5-6 (0) |
| **M3** jobs + plugins | job runner, manifests API, **parameter panel and job list components** (modal / drawer / inline hosting), `cli` blocks for all `proc` analysis steps, SSE progress, cancel; steps `extract`, `mask` auto, photo-z/SED, morphology, ICL, LSBG | each plugin's CLI step reproduces the desktop output byte-identically on m51/HUDF (re-use `scripts/verify_*`) | 6.5-8.5 (+0.5) |
| **M4** catalog + table + **Guided preset (default)** | catalog API (5.4), virtual grid, markers, selection link, kind filter, click chooser, review columns + filter, report export; **Guided shell**: stepper with the `ui.stage` mapping (status from jobs/session steps on the client until M6 makes it recorder-driven), step "cell" header, result strip, session log view, split-compare toolbar, catalog table as bottom sheet, collapse rules for 1024-1280 px; **object inspector** (cut-out, model/residual from multifit/isophote outputs, radial profile, p(z), review + note); **command palette** (registry from the manifests); **preset plumbing: layout switch, default `guided`, localStorage + URL parameter, first-run chooser (6.6)**; the 6.7 matrix for Guided | the `verify_review_gui`/`verify_click_*` scenarios re-written for the browser (Playwright); the scenario "accept 20 objects, export report"; 6.7 matrix green for Guided at 1440x900, 1280x720 and 1024x768 | 7-8 (+1) |
| **M5** mask + editing | mask API (5.5), polygon/ellipse editing, objects plugin (delete/merge/separate/add/AI-merge) | operations produce the same catalog as `ds9_catalog_edit.py` golden cases (`scripts/golden/`) | 4-5 (0) |
| **M6** sessions | server-side recorder, export of the Python script (`ogf_session_template.py`), replay check; the Guided stepper and session log switch from the client-side view to the recorder | exported script from a web session passes the same replay checks as the desktop (`verify_session_replay.py` equivalent) | 3-4 (0) |
| **M7** time-domain | `moving` plugin jobs, MAST cache, orbit/lightcurve views, Moving/Transient kinds in the table; **light-curve / tracklet dock component** (a tab beside the inspector or a bottom sheet in Guided; the dock of Workspace later) | the BB89 HST field gives the same tracklet as the desktop (`moving/tests`, `scripts/run_moving_pipeline.py`) | 4.5-5.5 (+0.5) |
| **M8** hardening | sandbox, rate limits, backups, monitoring, load test, security review, documentation, accessibility pass; preset-matrix CI job (6.7) | pen-test checklist, 20 concurrent users on the reference VM | 4.5 (+0.5) |
| **M9** Classic preset (second) | Classic shell: menu bar, workflow tabs + plugin ribbon from the manifests, info panel, button bar, frame tabs, colorbar, catalog panel with Galaxy/Moving/Transient/All tabs, parameter panel as a dialog; reuses the inspector (as an optional docked panel) and the palette that M4 built; the 6.7 matrix for Classic | 6.7 matrix green for Classic at 1440x900 and 1280x720; the scenario "accept 20 objects, export report" runs identically in Guided and Classic and exports the same script | 2-3 |
| **M10** Workspace preset (v1) | icon rail with the mode grouping, drawer (steps + parameter panel), top bar, fixed resizable dock (table, light curve/tracklet, jobs), light/dark theme, tablet layout; **no** RGB composite, minimap or free docking | 6.7 matrix green for Workspace at 1440x900, 1280x720, 1180x820; touch test on a tablet | 3-4 |

Totals: **M0-M8 (shared components + Guided, the default preset): 42.5-51.5** person-weeks (a single fixed layout was 38-48); with **M9** Classic (2-3) and **M10** (3-4):
**47.5-58.5** - unchanged by the revision of the default: the Guided shell, inspector and palette (3-4, formerly M9) moved into M4 while the Classic shell (2-3) moved out of it.  A useful first release (M0-M3 + extraction + table, **Guided only**) is about **21-27** (was 20-26 for a Classic-only release; +1 because the Guided shell, inspector and palette are bigger than the Classic shell).

Chapters 12 and 13 add the tracks **M3p** (parity harness) and **M11-M17** (metering/limits, payments, billing UI, local mode in waves): +39-57 person-weeks, total 86.5-115.5 (13.12).  They are not in the totals above.

### 9.1 Phased roadmap

1. **Phase 1 - Guided (M0-M8).**  The default and primary preset and the reference for all acceptance tests.  Built *as components* from the first day (viewer, table,
   parameter panel, inspector, palette, job list, status bar), so that the next presets only add shells.  Release 1 = Guided (Classic may be added to it as soon as M9 is done).
2. **Phase 2 - Classic (M9)**, can start as soon as M4 is done (it needs only the parameter panel, table/selection and the inspector/palette that M4 delivers) and is mostly shell work.
   Its purpose is familiarity for ds9 veterans and dense catalog work; reachable from the first-run chooser and the Layout control.  Release 1.x.
3. **Phase 3 - Workspace (M10)**, starts only when the **viewer-stable gate** is met: the `ImageSource`/renderer API unchanged for two releases, M1 acceptance
   tests green on the three target browsers, tile cache and WebGL context survive a canvas move (preset switch) without leaks, and the multi-band layer API
   specified (even if RGB/minimap are not built).  Backlog after M10: RGB composite with per-band blend, minimap, free docking.

### 9.2 Cost estimate: about +20-30 % over a single layout

Baseline (one fixed layout) 38-48 person-weeks.  Extra for the decision: component/store discipline and preset plumbing inside M1, M3, M4, M7, M8 = about
**2.5-3.5**; Guided shell, inspector and palette being larger than the baseline's single fixed layout (inside M4) **1**; Classic as the second preset (M9) **2-3**; Workspace v1 (M10) **3-4**; total roughly **9-11.5 person-weeks**, i.e. **about +20-30 %** of the baseline (about +21-27 % of its
mid-point, 43).  The components are built once; the extra is the shells, the test matrix (6.7) and the two other navigation widgets (menu bar + ribbon, rail/drawer) in addition to the stepper of the default preset.
What would raise it: Workspace with the deferred RGB/minimap/docking (+4-6, not included), a fourth preset, or per-preset forks of a component (excluded by
the rules in 6.1).  What would lower it: shipping Guided (M4) without the split-compare toolbar or the inspector's model/residual tabs in its first release.

## 10. Risks and open questions

* **Browser memory and speed** for 100-400 MB local files - threshold is a guess until M0; mitigation: lower the threshold, prefer upload.
* **Compute cost and abuse**: the pipelines are CPU- and memory-heavy (LSBG, ICL, moving pipelines on 4 x 1240² chips take minutes); without quotas and sandboxing a user can saturate the host.
  Mitigation: per-user job limits, time limits, container sandbox, queue priorities.
* **Security of file parsing**: the back ends parse user FITS with C code (`ds9_sextract`, `astropy`, `sep`); a malformed file is an attack vector.
  Mitigation: sandbox (no network, read-only root, seccomp), size/HDU limits at ingest, FITS validation before any job.
* **Data protection and proprietary data**: private by default; retention and deletion must be *enforced and tested* (purge test in CI); legal review of the privacy statement before inviting users.  This design deliberately contains **no sharing** feature.
* **Plugin manifests are not yet sufficient for the server**: 120 of 131 steps are Tcl `proc` steps; converting them to `cli` templates is real work (M3, probably the largest single item; the 6-8 week estimate assumes the golden argv file as test oracle) and may reveal Tcl logic that has to be reimplemented in Python (e.g. selection of stars for PSF building, trim, band registry, forced-photometry band loops).
* **State that lives only in Tcl globals** (`catpanel(...)`, `ed(...)`; `docs/architecture.md` section 7 lists what is still entangled) must be rebuilt as explicit session state; the key registry gives the list but not the semantics of the 52 table/widget references.
* **Reproducibility across the two paths**: local (browser) statistics and server statistics must agree (same zscale, same float handling); needs golden tests that run both.
* **Colour/scale parity with ds9** (zscale, asinh parameters, colormap tables): copy ds9's algorithms (`tksao`) and test numerically; display differences are a user-trust issue.
* **WebGL limits**: max texture size, float texture filtering (`OES_texture_float_linear`), mobile GPUs - fall back to `u16`/`png8`.
* **Licensing**: ds9 (GPL) code reused in the web viewer would make it GPL as well; clean-room the renderer or accept GPL; the same question applies to `ast` (LGPL) and `sep`/SExtractor-derived code.
* **Network-dependent steps** (MAST, JPL Horizons, AI services) need server egress and credentials; secrets stay in the server (`ai_bridge/auth.py` pattern); a user's own API key must be stored encrypted or not at all.
* **Cost** of storage and compute for a public instance is unknown; the retention policy is the main lever.
* **Preset sprawl and component coupling**: every extra preset multiplies the 6.7 test matrix, and a "small exception" inside a component for one layout breaks the independence rule.
  Mitigation: three presets only in Phase II; component-contract and minimum-size tests in CI; preset-specific behaviour must be a prop or a command, reviewed against 6.1.
* **State-store drift**: if components keep local copies of selection or view state, the preset switch loses them.  Mitigation: the ownership table (6.4), a lint rule against cross-component imports, the "switching" scenario of 6.7.
* **Viewer canvas hand-over on preset switch**: moving a WebGL canvas between shells may lose the context or flash on some browsers/GPUs.  Mitigation: tested in the M0 spike (5); fallback is re-creating the canvas and re-uploading from the tile cache (slower, still no network).
* **Manifest additions** (`ui.stage`, `ui.mode`, `ui.icon`) must stay optional and ignored by the desktop loader; wrong defaults misplace a step but never hide it (the palette always lists every step).
* **Layout preference storage**: localStorage can be blocked, cleared by the browser, or shared by two accounts on one machine; only layout, theme and language are stored there (6.6); an account-level copy is a later option.
* **Tablet and touch** are designed for Guided (the default, so part of release 1) and Workspace; Classic is desktop-only (1280 x 720 minimum); phone layouts are out of scope for Phase II.
* **Korean UI text**: glosses are partial (stage names, a few labels); a full translation, input-method (IME) handling in the palette/search and CJK line breaking in tables are not budgeted.
* **Mock-ups are not specifications**: they use example data and plain DOM (no WebGL); sizes, colours and wording will change.  What is decided is the component split and the state/preset rules above.
* **Dual-shell parity** (chapter 12): the rule that no tool is done without desktop and web entries slows early tools; legacy Tcl logic must move into Python cores first; native-vs-WASM float differences need tolerance classes (12.6).
* **Commercial model** (chapter 13): local-mode porting is the largest and least certain cost; distributing compiled GPL/LGPL-derived code needs a licence review; metering correctness, payment/tax law and abuse are new risk areas (13.13).
* **Default preset (resolved, 2026-10-03)**: Guided is the default for first-time users; ds9 veterans can switch to Classic with one action; a **first-run chooser** is offered (6.6) - the owner may still decide to drop it or to make it blocking.  Risk: veterans who never find the switch judge the product by Guided's slower catalog work; mitigation: the table-sheet rule (6.5), the visible Layout control and the first-run card.
* **Open decisions for the owner**: login method, quota/retention numbers, the first deployment target (institutional server vs cloud), whether the web version must support offline use (PWA), whether the layout preference also belongs to the account (server copy), licence, and the name of the repository (`ds10-web` is assumed here); the billing and compute-mode decisions are listed in 13.14.

## 11. First concrete steps (when the owner approves)

1. Create `ds10-web` (owner action); copy `plugins/*/plugin.json`, `ai_bridge/`, `moving/`, `ds9/library/ds9_*.py`, `ogf_session_template.py` as a **git submodule or a vendored directory with a sync script** so the desktop repo stays the single source of truth for algorithms; do not fork them.
2. M0 spikes (section 9), including the state-store/command-registry skeleton and the Guided/Classic throw-away shells.  Copy the three mock-ups from `/workspace/webdesign/` (and the DOM-check scripts in its `_src/`) into `ds10-web/docs/ui/` as the visual reference and the first golden screenshots (owner action; not done here).
3. Write the OpenAPI document for section 5 first; generate the TypeScript client and the FastAPI stubs from it, so that the browser and server work in parallel.
4. Add a `docs/` index in the new repository that points back to this file.
5. Before any billing code: settle the blocking items of 13.14 (plan structure, launch order, local runtime incl. the licence review) and fix the parity schemas of chapter 12 in M0 (result envelope, case format, `surface`/`compute` manifest fields).

## 12. 두 버전 동시 개발 원칙 (Dual-shell development principle)

**Principle.** Every feature and every analysis tool is developed **for both shells at the same time**: the stand-alone desktop version (SAOImageDS9/Tk, "Phase I") and the web version
(Phase II).  The two shells are two *front ends* of the same product, not two products.  The rules below are binding for all milestones (section 9) and for the billing
model (chapter 13).  The desktop version stays the reference implementation (section 1); the web version must not fork tool logic.

### 12.1 Rules

1. **One compute core per tool.**  The algorithm of a tool lives in exactly one place: `ogfkit/` (pure functions on numpy arrays and JSON-serialisable
   results, no Tk, no global state - already the stated contract of `ogfkit/__init__.py`), the tool's own package (`sersic_fit/`, `moving/`, `icl/`, ...) or
   a C library/binary (`ds9_sextract`, SEP).  Desktop and web call the same core through the same **CLI driver** (`ds9/library/ds9_*.py`: argv in, TSV/FITS/JSON out).
   No algorithm is re-implemented in Tcl or TypeScript; UI code may format and validate, never compute science results.
2. **One plugin manifest is the single source of truth for both shells.**  `plugins/<id>/plugin.json` defines tab, steps, parameters (type, default, range, group,
   help, `expert`), the `cli` argv template, the `output` mode, the session class - and, new, the additive blocks `ui` (`stage`, `mode`, `icon`; 6.5),
   `compute` (13.3) and `surface` (12.4).  The desktop dialog, menu, chip and the web parameter panel, ribbon, palette and API are generated from it.  A Tcl `proc` step may only open a
   dialog and then run the manifest step; it must not build its own argv (the 120 existing `proc` steps are converted to `cli` blocks in M3 - this principle is the reason that conversion is a
   precondition, not an option).
3. **Shared parameter and result schema.**  Parameters are the manifest `params` (JSON Schema is generated from them: `schemas/<plugin>.params.json`).  Results use one **result envelope**
   (v1): `{schema, core_version, step, status, outputs:[{name, kind: tsv|fits|png|json|html, ...}], columns:[...], metrics:{...}, messages:[...]}` wrapped around the
   existing stdout/TSV conventions (`output.mode` add_columns / set / text is kept).  Catalog column names come from the shared key registry.  Both shells read the same envelope; compute mode
   (local, server, stand-alone) never changes its content (12.5).
4. **One shared test suite, run against both shells.**  Parity cases (12.3) are executed on the desktop path (native CLI = reference) and on every web path (API/server mode, WASM/local mode).
5. **Definition of done: both entries or an explicit exemption.**  A new tool or feature is *not done* until it has a desktop entry **and** a web entry, or is marked desktop-only /
   web-only with a written reason in its manifest (12.4).  CI enforces this.
6. **Same release tag.**  Desktop build, server workers and the WASM bundle for local mode are built from the same git tag; the result envelope carries `core_version`; a mismatch is a visible warning, not silent.

### 12.2 Layers

```
 manifest (plugin.json: params, steps, cli, output, ui, compute, surface)         <- single source of truth
        |                                   |
  Desktop shell (Tk, ds9)             Web shell (components + presets, chapter 6)
  Tcl generates dialogs/menus         TS generates panel/ribbon/palette
        \                                   /
         \-- same step request: {plugin, step, params, inputs} --/
                          |
        execution backends (same CLI contract, same envelope)
          N  native process         desktop (stand-alone) and server workers (server mode)
          W  WASM/Pyodide worker    web local mode (13.3, 13.10)
          C  (option) local companion  native process on the user's PC behind the web UI (13.10, open decision)
                          |
        CLI drivers ds9_*.py  ->  ogfkit / tool packages / C (ds9_sextract, SEP)
```

### 12.3 Parity tests (CLI/API parity)

* **Case files** `tests/parity/<plugin>/<case>.json`: input fixture (m51, HUDF crop, the BB89 field; tiny versions for CI), step, parameters, expected outputs (golden) and a **tolerance class**:
  `exact` (byte-identical TSV/FITS: extraction, masks, catalog edits, stacking with fixed seed), `float` (relative 1e-6: fits, profiles, photometry), `statistical` (seeded random steps: summary statistics within bounds).
  Native-vs-WASM differences (libm, fused multiply-add, thread order) are expected at the 1e-7 level; a case declares which class applies and the class is reviewed, not loosened ad hoc.
* **Runners**: (a) *CLI runner* - native, the reference, also used by the desktop; (b) *API runner* - `POST /jobs` against a test server (server mode); (c) *WASM runner* - headless browser (Playwright)
  executing the same case in the local worker; (d) *desktop GUI smoke* - the existing `scripts/verify_*.tcl` scenarios drive the desktop entry of the step.  A case passes when (a) equals (b) equals (c) within its class.
* **Reuse**: `scripts/verify_cat_behavior.tcl` (golden argv of ~50 steps), `scripts/golden/`, `verify_session_replay.py`, `scripts/run_all_checks.sh` are the seed; the argv golden file becomes the manifest `cli` conformance test.
* **Gate**: a pull request that touches a core, a CLI driver or a manifest runs the parity cases of the affected plugins; releases run all of them on all backends that exist for the plugin (a tool without a local backend is simply not run there).

### 12.4 Definition of done and exemptions

A tool is done when all of the following exist (checklist in the PR template):

| item | desktop | web |
|---|---|---|
| core + unit tests | same | same |
| manifest block (`params`, `steps`, `cli`, `output`, `session`, `ui`, `compute`) | same | same |
| entry | menu/chip/dialog generated from the manifest (or a documented custom Tcl dialog) | generic parameter panel + ribbon/stepper/palette entry (or a documented custom component) |
| parity case(s) | CLI runner + GUI smoke | API runner (+ WASM runner if the tool is local-capable) |
| capability + cost calibration | n/a | row in the capability table (13.3), cost model calibrated on the fixtures (13.4) |
| documentation | `docs/<plugin>.md` | same document, web notes section |

Exemptions are declared in the manifest, never in a wiki: `"surface": {"desktop": true, "web": false, "reason": "uses Tk-only canvas plot; web replacement planned M9"}`.  An exemption needs a reason and a
review date; CI lists all exemptions at release time.  Known candidates today: Tk-only viewers (`ds9_analysis_gui.py`, "Interactive Plot...", "Analysis Viewer...") until the web plot components exist;
ds9 region language, XPA/SAMP (non-goals, section 1) - desktop-only by design; web-only: account, billing, upload and quota features (chapter 13), which have no desktop meaning.

### 12.5 Stand-alone vs web compute modes

The stand-alone desktop version **is the local mode without metering**: same cores, same CLI contract, executed as native processes on the user's computer, no tokens, no estimate dialog.
The web version offers the same step in three ways - stand-alone-like **local mode** (W: WASM in the browser, or C: companion process, 13.10), and **server mode** (N on our workers, metered, 13.4-13.6).
Consequences: a session recorded in any of them exports the same Python script and replays on any other (`--mode pipeline|replay`, section 5.6); the `mode` of each step (`desktop|local|server`) is recorded as
provenance but does not change the argv; a tool added for the desktop automatically has its web entry (generic panel) and only needs the capability row and cost calibration to be sold in server mode.
Whether the stand-alone version needs a licence/login is a business decision (13.14); technically it does not depend on the server.

### 12.6 Roadmap, cost and risks of the principle

* **Roadmap**: M0 fixes the schemas (result envelope, case format, `surface`/`compute` fields); **M3p** (with M3) builds the parity harness and converts `proc` steps; from M4 on every tool PR follows 12.4; the WASM runner appears with M14 (13.12).
* **Cost**: harness and conversion tooling 3-4 person-weeks (M3p, in addition to the `cli` conversion already counted in M3); afterwards a running **tax of roughly 10-15 % on each new tool** (parity cases, capability/cost rows, two entries).
  Cheap compared with the alternative (two diverging code bases), but it is not free and is counted in 13.12.
* **Risks**: (1) *schedule*: a rule that blocks merges slows the first tools; mitigation: generic panel/entry needs no custom code, exemptions are allowed with a reason.
  (2) *legacy logic in Tcl* (selection of PSF stars, trim, band loops, section 6.3) must move into Python cores before the parity cases can pass - this is the real cost driver and is already the largest M3 item.
  (3) *float differences* native vs WASM may force looser tolerances for iterative fits; tolerance classes are reviewed per case.
  (4) *desktop regressions*: changing a Tcl entry to a manifest-driven one can alter behaviour; the existing GUI scenarios (`verify_*.tcl`) stay mandatory.
  (5) *two release trains* (desktop installers vs web deploy) drift; mitigated by one tag (12.1 rule 6) and `core_version` checks.

## 13. 과금과 한도 (Billing and Limits)

Status: design proposal.  **This chapter contains no prices.**  Every price, token weight, allowance, limit, window length, retention period, refund percentage and threshold is a named
parameter (`<PARAM>` here, a key in `pricing.yaml` / `plans.yaml` in the implementation) to be set by the owner.  Effort numbers are estimates in person-weeks like in chapter 9.
Cross-references: chapter 12 (both shells are developed together; stand-alone = local mode without metering), section 4 and 7.1 (local vs upload), 5.3 (jobs), 6.1/6.4 (components, state).

### 13.1 Product model and vocabulary

OGFinder Web is a **commercial product** with two ways to compute, chosen by the **user per job** (and by a default in the settings):

| | **로컬 모드 (Local mode)** | **서버 모드 (Server mode)** |
|---|---|---|
| where the computation runs | the user's browser/computer: WebAssembly/JS (or the optional local companion, 13.10); data never leaves the machine | our servers (existing job runner, section 5.3, sandboxed workers) |
| what the server does | login, **licence** check, updates, billing - nothing else | everything: storage, tiles, jobs, outputs |
| charging | cheap **flat monthly subscription** `<PLAN_LOCAL_MONTHLY>`; no per-job charge; no token deduction | **tokens** are deducted for compute, disk and network (13.4); **daily / weekly / monthly limits** (13.7) like the plans of AI coding assistants |
| over the limit | n/a (no metering) | the user can **switch the job to local mode** (if the tool supports it, 13.3), wait for the reset, or **buy extra tokens** (13.7) |
| limits of the mode | tool coverage, browser memory, speed of the user's machine | cost, quota, queue |

Vocabulary used in the UI (English with Korean gloss): *token / 토큰*, *usage / 사용량*, *limit / 한도*, *estimate / 예상 비용*, *approve / 승인*, *top-up / 충전*, *refund / 환불*, *resets in / 재설정까지*, *local / 로컬*, *server / 서버*.
A **plan** bundles: the local subscription (optional), a server allowance per window (tokens), concurrency, storage allowance and retention.  A user may have both (local flat fee + server allowance) - open decision 13.14.
The stand-alone desktop version is outside this model: it is the local mode without metering (12.5).

### 13.2 Compute-mode selection (per job, default in settings)

* **Per-job selector**: every place that starts a job (parameter panel Run button, ribbon/stepper/drawer Run, palette command) shows **"Run on: 로컬 | 서버"** next to Run.  The control is part of the parameter panel
  component (6.1 #3), so all presets get it.  Server shows the estimate (`≈ x tokens`, 13.5) next to the label; Local shows "included in subscription".
* **Default** (settings, 13.11): `compute.default = local | server | ask`.  `local`: local whenever the tool is local-capable, otherwise the dialog asks to use server; `server`: server whenever allowance is available; `ask`: always show the choice.
  Per-tool overrides are allowed (`compute.override[plugin]`).
* **Suggestion rule** (replaces the fixed size rule of 7.1, which becomes only the *default suggestion*): suggest local if the tool is local-capable **and** the input fits the tool's `local.max_input_mb` **and** the browser reports enough memory;
  otherwise suggest server.  The user can always override; unsupported combinations are disabled with the reason ("Server only: ... / 이 도구는 서버 전용").
* **Recorded as provenance**: each session step stores `mode: local|server` and the `core_version`; the exported script is identical (12.5).
* **Switching a running job** is not supported; a queued/blocked server job can be re-submitted as local from the over-limit dialog (13.7) without re-entering parameters.
* **Files**: a local job needs the file in the browser.  A server job needs the data on the server; for a file that is only in the browser the approval dialog states the upload size and offers "upload only the cut-outs/region needed" for per-object jobs (13.8).

### 13.3 Per-tool capability table (which tools can run locally)

The classification is a **proposal built from the current dependencies of the CLI drivers** (`ds9/library/ds9_*.py`, imports of `ogfkit`/tool packages) and from the section 7.1 numbers; it is **confirmed or changed by measurement in M14** (13.10) and then stored as data.
Classes: **L** light, local-capable and suggested locally; **B** both - local possible (port needed, memory/time depend on input), server for large inputs; **S** server-only in v1 (local not offered).  Server mode is always available unless noted.

| plugin / steps (tab) | compute core and dependencies | class | local (WASM/JS) | server | remarks |
|---|---|---|---|---|---|
| `extract`: Extract, Dual-Image Extract, Trim catalog (Detect) | `ds9_sextract` (C), SEP | L / B | C -> WASM (Emscripten) plausible, unproven; native full HUDF F160W 1.9 s (section 7.1) | yes | mosaics beyond `local.max_input_mb` -> server |
| `objects`, `bands`, `catalog`, `session`, `repro` (Detect/Results) | catalog edit (pure Python), registry, session recorder | L | JS/Pyodide, small | not needed | multi-band loops inherit the cost class of the measuring step |
| `noisemodel`, `xmatch` (local catalogs), `report` (HTML) | numpy, astropy | L | Pyodide (numpy/astropy available there; versions to verify) | yes | VizieR cross-match and PDF report need network/tools -> S |
| `star_psf`: Find Stars, Build PSF (Classify) | SEP, `star_finder`, numpy | L / B | Pyodide + SEP build | yes | AI Star Classification (torch) = S; WebbPSF, TinyTim = external tools, S |
| `morphology`: Non-Parametric (CAS/Gini/M20); `morph_ext`; `isophote`; `stacking`; `cluster` (Measure) | `morphometry`, `ogfkit`, numpy/scipy | B | Pyodide; memory grows with number of cut-outs | yes | |
| `morphology`: **Sersic Fitting**, Bulge+Disk; `multifit` (**GALFIT-like fits**); `psf_deconv` Deconvolve | `sersic_fit`, `bulge_disk`, `ogfkit.multifit` (numpy, `parallel` multi-process) | **S** | not offered in v1 (iterative, many sources, multi-process); revisit after M14 measurements | yes (default) | GALFIT feedme import/export is just file conversion (L). The external GALFIT binary is closed source: not ported; server use only if its licence allows it |
| `photoz_sed`: Photo-z, SED fitting, calibration; `sedcodes` (EAZY) | `photo_z`, `sed_fit`, external codes, trained checkpoints | **S** | not offered | yes | model files and external codes stay on the server |
| `galaxy_model`: AI Morphology Classification; `ai_services` | PyTorch (`galaxy_morph`), `ai_bridge`, provider APIs | **S** | not offered; an ONNX-Runtime-Web (WebGPU) inference port is a possible later project | yes | keys and provider egress stay on the server |
| `trails`: Remove Trails, Detect Trails, Stack excluding trails | Radon + contrast scan on full frames | **S** | not offered in v1 | yes | cost measured in M11; may become B |
| `lsbg`, `icl`, `mask`: background model, auto mask, detection, SVM classify | SEP, scipy, `icl`, `lsbg` | B / S | small frames possible | yes | Run Full Pipeline on large mosaics = S |
| `moving`: Difference (image subtraction, ZOGY), Align, Detect, Link, Identify, Orbit; `lightcurves` | `moving/` (numpy/scipy), MAST/JPL/MPC access | **S** | not offered | yes | network access (MAST, JPL Horizons) only from the server; Light Curve/Export display = L |
| `photometry`: PSF photometry (B), Crowded field (S), Cross-Match VizieR (S), `completeness` (S), `daophot` (B/S), `lensmodel` (S), `spectra` kinematics (S), `batch` (S) | `psf_phot`, `crowded_phot`, `parallel`, injection loops | B / S | per step | yes | many-trial simulations are server work |

The decision is stored in the manifest and the server (additive `compute` block, same single-source-of-truth rule as 12.1):

```json
"compute": { "modes": ["local","server"], "default": "server",
             "local": { "runtime": "wasm|pyodide|companion", "max_input_mb": "<N>", "min_memory_mb": "<N>" },
             "cost_model": "<id of the calibrated cost model, 13.4>" }
```

`GET /api/v1/plugins` returns it; the UI disables the unsupported mode with the reason.  A step without a `compute` block is server-only until classified.

### 13.4 Metering model (server mode)

Three resource families are metered, each in physical units, then converted to **tokens** by weights that live in a versioned `pricing.yaml` (`pricing_version` is stored on every ledger entry, so past usage is never re-priced).

| family | measured quantity | how it is measured | token weight (placeholder) |
|---|---|---|---|
| **Compute** | CPU core-seconds, GPU-seconds (per GPU class), memory GB-seconds (resident set above a free baseline) | worker cgroup counters sampled every `<METER_INTERVAL_S>` and at exit; wall time x allocated cores for reserved-core jobs | `w_cpu`, `w_gpu[class]`, `w_mem` |
| **Disk** | GB x retention days: originals, tile pyramids, masks, catalogs, job outputs, sessions | daily snapshot of bytes per account from the object store (`bytes_day` accrues once per day); pinned data keeps accruing | `w_disk` |
| **Network** | egress GB (downloads of originals/outputs, exports, tile traffic beyond a fair-use allowance); ingress (uploads) not charged | reverse-proxy byte counters per account | `w_net`, `NET_FAIR_USE_GB_PER_DAY` |

```
tokens(job)  = ceil_micro( w_cpu*cpu_core_s + w_gpu[c]*gpu_s + w_mem*mem_gb_s            // compute (measured during the job)
                         + w_net*egress_gb )                                             // outputs fetched by the user (later, to the same job id)
tokens(day)  = ceil_micro( w_disk * stored_gb )                                          // per day, per account (not per job)
```

* Token amounts are integers in micro-tokens (no floating point in the ledger).  A minimum charge per server job `MIN_JOB_TOKENS` covers scheduling overhead; interactive operations (tiles, catalog rows, cut-outs of already uploaded data) are **not charged per request** but count against the network fair-use and rate limits.
* Each plugin step has a **cost model** (`cost_model` in the manifest): resources as a function of input size and parameters, calibrated on the fixtures (m51, HUDF crop, BB89 field) and refined from production history (median and p90 per size bucket).  Calibration is part of the definition of done (12.4).
* **Local mode is not metered**; the only server cost is login/licence/update traffic, covered by the subscription.
* Transparency: every job shows its resources (core-seconds, memory peak, GB, egress) and the resulting tokens next to the log (job list component).

### 13.5 Pre-run estimate and approval dialog

1. The user presses Run with **Server** selected -> `POST /api/v1/estimates` with the same body as a job (plugin, step, inputs, params).  The server returns an `estimate_id` valid for `<ESTIMATE_TTL_MIN>` minutes and bound to the input hashes and parameters.
2. The **estimate and approval dialog** (component, 13.11) shows: tool and step; mode; expected tokens as **p50 and upper bound (p90)**; the **hold** that will be reserved (the cap, default `p90 x ESTIMATE_SAFETY`, editable by the user within `[p50, MAX_JOB_TOKENS]`);
   resources behind it (core-hours, memory, disk GB x days, egress); data to be uploaded (bytes, or "only cut-outs: x MB"); retention of the outputs and the auto-delete date; remaining allowance per window (day / week / month) and top-up balance **after** the hold;
   and the alternatives: **Run locally instead** (if capable, with the local limits: input size, memory), **Lower cost options** (e.g. "max-sources 100 -> 20"), **Cancel**.
3. **Approve** sends `POST /api/v1/jobs` with `estimate_id` and `approve: {max_tokens}`.  The server re-checks the estimate (not expired, inputs unchanged), creates the hold (13.6) and queues the job.
4. **Auto-approve** (settings): jobs whose upper bound is below `<AUTO_APPROVE_TOKENS>` start without the dialog (the status bar shows the deduction); anything above, any upload, and any job that would cross 80 % (`<WARN_FRACTION>`) of a window always ask.
5. If the estimate is **unavailable** (new tool without calibration) the dialog says so and requires an explicit cap; the job is killed at the cap.
6. If the user is **over a limit** the dialog does not offer Approve; it shows which window blocks, when it resets, and the buttons *Run locally*, *Buy tokens*, *Wait* (13.7).

### 13.6 Holds, settlement and refunds

* **Hold at start** (reserve): the cap is reserved against the allowance of **every** window (day, week, month) and, if the allowance does not cover it, against the top-up pool.  No hold, no start.  Holds appear in the ledger (`kind=hold`) and reduce "available", not "used".
* **During the job**: usage is appended every `<METER_INTERVAL_S>` (`kind=usage`, running total ≤ cap).  When the running total reaches the cap the worker is stopped (SIGTERM, then kill); partial outputs are kept for `<PARTIAL_KEEP_DAYS>` and marked `capped`.
* **Settlement** (post-run): at exit the actual tokens are computed from the final counters; `settle` converts the used part of the hold into usage and `release` returns the rest to the windows.  The user sees estimate vs. actual in the job list.  Late network egress (downloading outputs) is added to the same job id.
* **Refund policy** (percentages are parameters):

| outcome | charge |
|---|---|
| success | actual usage |
| platform failure (node lost, our bug, timeout caused by the platform) | refund **100 %** (`REFUND_PLATFORM`) |
| input validation failure before compute starts | no charge (`MIN_JOB_TOKENS` waived) |
| tool exit with error after compute (user data/parameters) | actual usage minus `REFUND_USER_ERROR` % (could be 0) |
| cancelled by the user | actual usage up to the stop |
| killed at the cap | the cap |

* Refunds restore the window allowance of the **window in which the usage happened** only if that window is still open; otherwise they are credited to the top-up pool as `refund` tokens (rule is a parameter, avoids reopening closed windows).
* Everything is idempotent (`idempotency_key` = job id + phase), so a retried webhook or worker restart cannot double-charge.

### 13.7 Limit manager, top-ups and payments

**Windows.**  Server allowance is enforced in three nested windows: **daily** `<LIMIT_DAILY_TOKENS>`, **weekly** `<LIMIT_WEEKLY_TOKENS>`, **monthly** `<LIMIT_MONTHLY_TOKENS>` per plan (tokens that are available = min over windows of `limit - used - holds`).
Window kind is a parameter: *fixed calendar windows in the account time zone* (default Asia/Seoul, reset 00:00; weeks from Monday) or *rolling windows* (e.g. `<ROLLING_HOURS>` as in some assistant plans) - open decision 13.14.
The limit manager exposes, per window: `limit`, `used`, `held`, `resets_at`; the UI shows "재설정까지 3시간 20분 / resets in 3 h 20 min".  Plan allowance does not roll over unless `<ROLLOVER_FRACTION>` > 0.

**Order of consumption**: plan allowance (windows) first, then the **top-up pool** (not subject to windows, except the per-job cap and the rate limits).

**Over the limit.**  A job that cannot be held opens the **limit dialog** listing the blocking window(s) and the exact missing tokens, with: (1) *Run locally* (only if the tool is local-capable; re-submits the same step, 13.2), (2) *Buy tokens* (13.7 top-ups),
(3) *Schedule at reset* (queue the job to start automatically at `resets_at` if the user approves the cap now; the hold is taken then), (4) *Change plan*.  A job already running is never stopped by a window; its hold already covers it.

**Top-ups.**  Extra tokens bought in packs `<TOPUP_PACKS>`; purchased tokens expire after `<TOPUP_EXPIRY_DAYS>` (0 = never; decide with legal advice).  Auto top-up (optional, user-set monthly maximum `<AUTO_TOPUP_MAX>`) is off by default.
Top-ups appear in the ledger as `topup`; unused expired tokens as `expire`.  A top-up is available as soon as the payment provider confirms (webhook), the UI shows "pending" until then.

**Payments and billing integration (placeholders).**  The application talks to a `PaymentProvider` interface and never to a vendor SDK directly in business logic:
`create_customer`, `create_subscription(plan)`, `change_plan`, `cancel_subscription`, `create_checkout(topup_pack)`, `portal_url`, `refund(payment_id)`; incoming webhooks `payment.succeeded`, `payment.failed`, `subscription.renewed|canceled|past_due`,
`refund.succeeded`.  Provider candidates to evaluate (no choice made here): an international card processor and a Korean payment aggregator (card, bank transfer, domestic invoicing).  Placeholders to be filled by the owner and legal: currency and VAT handling,
tax invoice (세금계산서) / receipt rules, invoicing for institutions (purchase orders, lab/group accounts with a shared pool), dunning (failed payment -> grace `<PAST_DUE_GRACE_DAYS>` -> downgrade to read-only), terms of service and refund terms.
Entitlements are derived from the provider state via the webhook only (`entitlements` table), so a UI bug cannot grant tokens.

### 13.8 Data retention, privacy and abuse prevention

**Where data lives, by mode**
* **Local mode**: files are read with the browser `File` API and never uploaded; results (catalogs, masks, sessions) are kept in the browser (OPFS/IndexedDB, within the browser quota) and exported by the user; the server sees only account, licence and version information (no file names, no telemetry about images; any usage analytics is opt-in).
* **Small files can stay in the browser even when the user works in server mode for other things**: a file is uploaded only when a server job on it is approved (the dialog states the size) or when the user chooses "Analyse on server" (7.1).  For per-object jobs (photo-z, fits on selected rows) the option **"upload cut-outs only"** sends the needed regions instead of the whole image: smaller disk/network bill and less data leaving the machine.
* **Server mode**: private to the owner, no public links, no use for training (7.2).  The disk meter (13.4) and the retention rule below apply.

**Retention and auto-delete rule** (all numbers are parameters)
1. Default retention `<RETENTION_DAYS>` days after last access for originals and derived data; the deletion date is shown in the file list and in the approval dialog.
2. **Pin** ("keep") exempts data from auto-delete; pinned data keeps accruing disk tokens daily.  `<FREE_STORAGE_GB>` of pinned data may be included in the plan.
3. Notices: `<NOTICE_DAYS>` and 1 day before deletion (in-app and email).  After the deadline files, tiles, masks and job outputs are deleted at once; catalogs/sessions keep their record but not the pixels (7.2); the ledger keeps usage without file names.
4. **Out of tokens while data is pinned**: disk accrual that cannot be covered moves the data to `expiring` (read-only, no new jobs) for `<DISK_GRACE_DAYS>`; then it is deleted by the same rule.  Never silently: notices at the start of the grace period and 1 day before.
5. Account deletion removes everything within 7 days (7.2); the immutable ledger is retained only as long as accounting law requires (legal item).

**Abuse prevention**
* **Rate limits** per account and per IP: API requests, uploads (bytes/hour), estimate requests, login attempts, webhook endpoints (signed); parameters `<RATE_*>`; 429 with `Retry-After`.
* **Concurrency**: maximum concurrent server jobs per plan `<MAX_CONCURRENT_JOBS>`, queue length per account, per-job caps (cores, memory, wall time) from the plan, fair queueing so one account cannot starve others.
* **Holds** make spam expensive: a job cannot start without a hold, and a stalled queue does not accrue usage.  **Estimate flooding** is limited by the estimate rate limit.
* **Trial and free tier**: a one-time token grant `<TRIAL_TOKENS>` per verified account (email + payment method or institutional login), device/IP heuristics against multiple trials; grants are ledger entries (`grant`) with expiry.
* **Licence enforcement for local mode**: signed licence token with expiry, renewed online; offline grace `<OFFLINE_GRACE_DAYS>`; maximum devices `<MAX_DEVICES>`.  Client-side code can always be copied, so licence checks deter casual sharing only (13.13).
* **Payment fraud and chargebacks**: provider risk signals; chargeback -> account freeze and ledger adjustment.
* **Sandboxing and file parsing** as in section 10 (no network except egress proxy for `moving`/AI jobs).  **Anomaly alerts** to operators: consumption x times the user's median, many accounts from one payment fingerprint.

### 13.9 API and data schema (usage ledger)

New or changed endpoints (v1, JSON, same auth as section 5):

| method, path | purpose |
|---|---|
| `GET /api/v1/me/plan` | plan, entitlements (local licence, server allowance), concurrency, retention, payment status |
| `GET /api/v1/usage?window=day\|week\|month` | per window: `limit`, `used`, `held`, `resets_at`; top-up balance; a short history; SSE `GET /api/v1/usage/events` for live updates (status-bar meter) |
| `POST /api/v1/estimates` -> `{estimate_id, expires_at, resources:{cpu_core_s,gpu_s,mem_gb_s,disk_gb_day,egress_gb}, tokens:{p50,p90,cap_default,min_job}, upload_bytes, retention_days, delete_on, windows_after:[...], can_run_local, local_reason}` | pre-run estimate (13.5) |
| `POST /api/v1/jobs` (extended) `{plugin, step, inputs, params, session_id, compute_mode:"server", estimate_id, approve:{max_tokens}}` -> `{job_id, hold_id}` or `402 limit_exceeded {blocking_windows:[...], missing_tokens, can_run_local}` | local jobs create no server job; the client only appends the session step (`mode:"local"`) |
| `GET /api/v1/jobs/{id}` (extended) | adds `mode`, `estimate`, `hold`, `usage_so_far`, `settled`, `refund` |
| `GET /api/v1/ledger?from=&to=&kind=&job_id=` and `.../ledger.csv` | the user's usage ledger (read-only) |
| `POST /api/v1/topups` `{pack}` -> `{checkout_url}`; `GET /api/v1/topups` | buy / list top-ups |
| `GET /api/v1/billing/portal` | provider portal (invoices, card) |
| `POST /api/v1/webhooks/payments` (provider signature) | entitlement and top-up updates |
| `GET /api/v1/license` / `POST /api/v1/license/refresh` | signed licence token for local mode (device id, expiry, `core_version` range) |
| `GET /api/v1/files/{id}` (extended) | adds `retention: {delete_on, pinned}`; `PUT .../pin` |

```ts
// append-only; balances are derived, never stored as the only truth
interface LedgerEntry {
  id: string;                    // ULID
  account_id: string;
  ts: string;                    // UTC
  kind: "hold"|"usage"|"settle"|"release"|"refund"|"topup"|"grant"|"expire"|"adjust"|"disk_day"|"egress";
  job_id?: string; hold_id?: string; file_id_hash?: string;          // no file names
  resource?: "cpu_core_s"|"gpu_s"|"mem_gb_s"|"disk_gb_day"|"egress_gb";
  quantity?: string;             // decimal string, physical unit
  pricing_version: string;       // weights used (13.4)
  tokens_micro: string;          // signed integer; hold/release are reservations, others consume or credit
  pool: "allowance"|"topup";
  window_effects: { window: "day"|"week"|"month"; window_start: string; delta_micro: string }[];
  idempotency_key: string;       // unique per account
  note?: string;
}
```

Tables (server DB): `plans`, `subscriptions`, `entitlements`, `limit_windows(account_id, kind, window_start, limit_micro, used_micro, held_micro)`, `holds(id, job_id, cap_micro, state: active|settled|released)`, `estimates(id, inputs_hash, params_hash, resources, expires_at)`,
`usage_ledger` (above), `topups`, `tool_capabilities` (from the manifest `compute` blocks + measured limits), `cost_models`, `licenses(device_id, expires_at)`, `rate_counters`.  `window_*` rows are caches rebuilt from the ledger (a nightly check recomputes and alerts on drift).

### 13.10 Local mode: runtime, licence, updates and the porting cost (honest assessment)

**What "local mode" needs in the browser**
* a **worker pool** (Web Workers) executing the same CLI contract: argv in, files in a virtual file system (OPFS/MEMFS), envelope out (12.2 backend W);
* a **runtime**: C tools compiled with Emscripten to WebAssembly; Python/numpy tools run in **Pyodide** (CPython in WASM).  Pyodide ships numpy, scipy, astropy and matplotlib; the exact set and versions, and whether **SEP** (C extension, used by most `ds9_*.py` drivers) is available or must be built by us, are to be verified in the M14 spike;
* **threads**: `SharedArrayBuffer` needs cross-origin isolation (COOP/COEP headers) on our pages; Python `multiprocessing` (our `parallel` helper) does not work in Pyodide and must be replaced by a worker pool or run serially;
* **licence and updates**: the signed licence token (13.9) is checked at start and renewed online; WASM bundles are versioned with the release tag (12.1 rule 6), integrity-checked (hash pinning) and cached for offline use within `<OFFLINE_GRACE_DAYS>`.

**Porting cost - what is easy and what is not**

| area | assessment |
|---|---|
| `ds9_sextract` and SEP (C) | the easiest: Emscripten builds of C libraries are routine; the extraction of a 52 MB image takes 1.9 s natively (section 7.1), a 1.5-3x slowdown in WASM is a common expectation (**unmeasured**); I/O adapters and FITS reading still need work |
| pure numpy/scipy/astropy drivers (morphometry, isophote, stacking, noise, xmatch, report HTML) | moderate: run in Pyodide with adapters; cost is in file/IO shims, startup time (tens of MB of runtime to download and cache), memory, and removing `multiprocessing`, `psutil`, `tkinter` imports |
| iterative fitters (`sersic_fit`, `multifit`, `bulge_disk`) | hard to make *good*: they run, but a many-source fit that takes minutes natively on 8 cores becomes much slower on one or a few browser threads; usable only for a handful of sources -> kept **S** in v1 |
| PyTorch models (`galaxy_morph`, `ai_merge`, star classifier) | **no PyTorch for Pyodide**; needs an export to ONNX and ONNX Runtime Web (WebGPU/WASM) per model, inference only, plus validation against the native outputs; separate project |
| network-dependent tools (MAST, JPL, VizieR, AI providers) | cannot be local (keys, CORS, egress) |
| external closed or licensed tools (GALFIT binary, EAZY, WebbPSF, TinyTim) | not ours to compile or redistribute |
| memory | wasm32 addresses at most 4 GiB (browsers often allow less) and a tab competes with the viewer's textures; the 200 MB local threshold (7.1) is a guess until measured |
| determinism | native vs WASM float results can differ at ~1e-7 (12.3 tolerance classes) |
| correctness cost | every ported tool needs the parity cases (12.3) on the WASM runner: this is a large share of the work |

**Option C - local companion (decision needed)**: the stand-alone desktop application (or a small headless build of it) already runs every tool natively.  It can expose the *same job API on localhost* so that the web UI uses **the user's own computer at native speed**
without porting to WASM.  This is still "local mode" (the user's computer, flat subscription, no metering, no upload) and costs far less than porting wave 3, but the user must install it, and browsers restrict calls from a public https page to localhost
(mixed content, CORS, Private Network Access prompts) - to be tested in M14.  Recommended plan: WASM for the light class L/B tools (waves 1-2) for the zero-install experience, and the companion for the heavy tools instead of WASM wave 3, subject to the M14 spike.
The stand-alone version is then not a competitor but the engine of the heavy local mode.

**Waves** (M15-M17, 13.12): wave 1 = L-class tools (extract, objects/catalog/session/repro, noise, xmatch local, report HTML, star finder, Find Stars); wave 2 = B-class tools (non-parametric morphology, isophote, stacking, cluster, mask/background for small frames); wave 3 = optional heavy tools (only if the companion route is rejected).

### 13.11 UI components for billing, limits and compute mode (shared by all presets)

The shared component list of 6.1 gains two components (#9, #10) and extends two others (#3 parameter panel: "Run on" selector; #8 status bar: usage meter and mode badge).

| # | component | content | states to design and test |
|---|---|---|---|
| 8+ | **Usage meter** (in the status bar component) | compact bars for day/week/month (`used + held` of `limit`), top-up balance, compute-mode badge (로컬 / 서버), click -> usage popover (resets in ..., last jobs, "Buy tokens", "Open account") | normal, ≥ `WARN_FRACTION`, at limit, no plan/trial, offline (local only), payment past due |
| 3+ | **Run-on selector** (in the parameter panel) | 로컬 / 서버 segmented control, estimate text, disabled-with-reason states | tool local-only impossible, server unavailable, estimate loading/failed |
| 9 | **Estimate and approval dialog** (13.5) | tokens p50/p90 and cap, resources, upload size, retention/delete date, windows after hold, alternatives, Approve / Run locally / Cancel; the over-limit variant (13.7) | estimate ok, no calibration, over limit (which window, reset time), upload required, auto-approved toast |
| 10 | **Account, usage and settings page** | route `/account` (full page, not a modal): plan and licence (devices), usage windows and history chart, **ledger table** (filters, CSV), top-up and payment (provider portal), **default compute mode** and per-tool overrides, auto-approve threshold, retention default and pinned files, notifications | loading, empty ledger, payment failed, licence expired, after top-up |

Placement per preset (all components read the same `account`/`usage` state, see the state store additions below):

| | **Guided** (default) | **Classic** | **Workspace** |
|---|---|---|---|
| usage meter | pill in the header (next to the avatar) and in the cell header; status bar shows the running job's tokens | right side of the status bar; menu `File > Account...` | pill in the top bar; run bar of the drawer shows the estimate |
| run-on selector | in the step cell next to *Run* | in the parameter dialog next to *Run*; ribbon chips run with the default mode (right-click: "Run on ...") | in the drawer run bar |
| estimate dialog | modal (the cell shows the estimate inline before the dialog) | modal | modal (popover from the run bar for the quick view) |
| settings page | avatar menu -> Account | `Edit > Preferences...` -> Account tab; full page route | avatar menu -> Account |

State store (6.4) additions: `account: { plan, entitlements, license }`, `usage: { windows: { day, week, month: { limit, used, held, resetsAt } }, topupBalance }` (written only by usage events from the server), `jobs.byId[id].{mode, estimate, hold}`,
`ui.computeDefault: "local"|"server"|"ask"` (server copy in the account, mirrored in the browser, 6.6).  Verification (6.7) gets the rows for these states.
Korean UI strings (examples): 사용량, 한도, 오늘/이번 주/이번 달, 예상 비용, 승인, 서버 전용 도구, 로컬에서 실행, 토큰 충전, 재설정까지.

### 13.12 Roadmap and added cost

Dependencies: the server track needs M2 (storage, auth) and M3 (jobs); the UI components follow the presets they live in; local mode needs M14 first.  New milestones (person-weeks; **estimates**, not commitments; legal, tax, pricing research, customer support and operations are **not** included):

| M | content | deliverable and acceptance test | estimate |
|---|---|---|---|
| **M3p** (with M3) | parity harness (12.3): result envelope, case format, CLI/API runners, `surface` and `compute` fields, PR checklist; conversion tooling for `proc` -> `cli` | the first 10 plugins' cases pass on CLI and API runners; CI blocks a plugin without entries/exemption | 3-4 |
| **M11** | metering and ledger: resource counters, cost-model framework + calibration of the first tools, estimate API, hold/settle/release/refund, limit manager (3 windows, resets), rate limits and concurrency, usage events, disk-day accrual and auto-delete rule | simulated accounts: holds never exceed limits, settlement reconciles to the byte counters, failure refund, nightly ledger-vs-window check, retention purge test | 5-7 |
| **M12** | billing integration: `PaymentProvider` interface + one provider adapter (placeholder), subscriptions, entitlements via webhooks, top-ups, invoices/portal, dunning, trial grant, licence token service | test-mode purchase end to end; webhook replay is idempotent; entitlement changes only through webhooks | 4-6 |
| **M13** | UI components of 13.11 in **Guided** (the default preset), then Classic (with M9) and Workspace (with M10): +0.5 each; account page; 6.7 matrix rows | all billing states pass layout/DOM checks in each preset at 1440x900 and 1280x720 | 3-4 (+1 over the presets) |
| **M14** | local runtime: worker pool, Pyodide/WASM loader, licence check and offline grace, signed/cached bundles, capability negotiation (`compute` blocks), WASM parity runner; spikes: SEP in WASM, ds9_sextract in WASM, COOP/COEP, companion-over-localhost, ONNX feasibility, measured native-vs-WASM timings that decide the capability table | spike report with numbers; `extract` parity case green on the WASM runner | 4-6 |
| **M15** | local wave 1 (L tools, 13.10) | parity cases on the WASM runner for each tool | 6-10 |
| **M16** | local wave 2 (B tools) | same | 6-8 |
| **M17** | local wave 3 (heavy tools) - **only if the companion route is not chosen** | same; or instead: companion packaging and the localhost job API (estimate 3-5) | 8-12 |

**Added cost** over the 47.5-58.5 person-weeks of chapter 9 (single layout x 3 presets): M3p 3-4 + M11 5-7 + M12 4-6 + M13 3-4 + M14 4-6 + M15 6-10 + M16 6-8 + M17 8-12 = **39-57 person-weeks**, giving **86.5-115.5** in total (roughly +70 to +120 %).
Breakdown: parity principle 3-4; server billing (M11-M13) 12-17; local mode 24-36 (**the dominant and the most uncertain part**).  Recurring: about 10-15 % extra on every new tool (12.6) and operations/support.
Suggested order: M3p -> M11 -> M12 -> M13 (Guided) = a sellable **server-mode** product (about 15-21 person-weeks on top of the base release); in parallel M14 and wave 1 so that the flat local plan can be sold with a meaningful tool set; waves 2-3 or the companion afterwards.
If local mode is cut to wave 1 only (no M16, M17), the added cost is about 25-37.

### 13.13 Risks

* **Local-mode porting is the largest cost and schedule risk** (13.10): speed, memory, SEP/Pyodide availability, multiprocessing, determinism.  Mitigation: M14 spike with numbers before promising tools; capability table is data; companion option; sell local mode with the tool list that exists.
* **Licensing of distributed code**: local mode *distributes* compiled code (WASM bundles, companion) to users.  ds9 (GPL), AST (LGPL), SEP and SExtractor-derived code have terms that apply to distribution and may conflict with a closed commercial product; server-only use is a different situation.  Needs a licence review **before** M14 (not legal advice; section 10 already notes the GPL question for the viewer).
* **Metering correctness and trust**: wrong estimates, double charges or lost refunds destroy trust.  Mitigation: append-only ledger with idempotency, nightly reconciliation, user-visible estimate vs actual, caps, simulated-account tests, a manual `adjust` path with audit.
* **Estimate quality**: cost models for new tools are uncalibrated at first; cap-and-kill and the "no calibration" flow (13.5) limit damage; p90 buffers cost users money - tune `ESTIMATE_SAFETY`.
* **Limit-window design** (rolling vs fixed, nested windows) is easy to get confusing; show the blocking window and reset time everywhere; test the over-limit dialog with real users.
* **Subscription gaming and copying**: client-side code cannot be made un-copyable; licence tokens, device limits and updates only deter casual sharing.  Keep tool value (server models, data, updates) on the server.
* **Payments and legal**: tax/VAT, invoicing, refund law, chargebacks, data-protection terms for a Korean/international user base; not budgeted in 13.12.
* **Cost of the free path**: trial grants, platform-failure refunds and cached tile traffic have real server cost; set `TRIAL_TOKENS`, `NET_FAIR_USE_GB_PER_DAY` conservatively and monitor.
* **Two modes double the test matrix** (12.3, 6.7); local mode may drift from server mode if parity gates are skipped.
* **Data loss expectations**: auto-delete and `expiring` states can surprise users; notices (13.8), explicit pin and export prompts; deletion is tested in CI (section 10).
* **Over-engineering for the first release**: the first commercial release can ship fixed windows, one provider, no auto top-up, no rollover, no group accounts.

### 13.14 Open decisions for the owner

1. Plan structure: local-only plan, server-only plan, or bundle; is the stand-alone desktop licensed/priced separately, free, or included (technically independent of the server, 12.5)?
2. Launch order: server mode first, local mode first, or both; which tool list is the minimum for the flat local plan?
3. Local runtime: WASM only, WASM + companion (13.10), or companion for heavy tools; acceptance of the licence review result (13.13).
4. Window kind (fixed calendar vs rolling), window lengths, whether allowance rolls over, time zone default (Asia/Seoul proposed).
5. Token definition and all weights (`pricing.yaml`), allowance per plan, `MIN_JOB_TOKENS`, `ESTIMATE_SAFETY`, auto-approve threshold; whether disk and network are charged separately or only shown.
6. Top-up pack sizes, expiry of purchased tokens (legal), auto top-up default.
7. Retention default, free pinned storage, disk-grace period, notice schedule.
8. Refund policy percentages for platform failure / user error / cancel.
9. Payment provider(s), currency, tax-invoice handling, institutional/group accounts and invoicing.
10. Trial/free tier size and verification method; device limits for the local licence; offline grace.
11. Whether usage analytics for local mode exist at all (default proposed: none, opt-in only).
12. Whether tiles/viewing traffic is metered (proposed: fair-use allowance, not per request).


---

## Memo (ds10-web M12, 2026-10-04): the plan field that exists today, and what chapter 13 can build on

Status: **implemented in the web repo as a plain account field, no billing logic.**  Tiers are *ordered*: 1 = free, 2 = pro, 3 = max, 4 = ultimate (higher = higher tier); they replace the earlier two-step
free/paid idea.  Facts for the billing chapter (13.4-13.11):

* **Where it lives.** `users.plan` (INTEGER 1..4, default 1) next to the hash-chained usage ledger; nothing in the ledger refers to it and the ledger still writes `pricing_version = "none"`, `tokens_micro = "0"`.
  A future billing layer decides *which plan a payment buys* and calls the same administrator route (`PATCH /api/v1/admin/users/{id}` `{"plan": n}`); only administrators (or a billing service acting as one) may change it,
  and every change is an audit-log entry `admin.user.plan` ("1 free -> 2 pro").
* **Per-plan settings table = the `plans.yaml` of 13.8, as JSON/env today** (`server/ds10web/plans.py`, `DS10_PLANS_FILE`, `DS10_PLAN_<NAME>_<FIELD>`): storage capacity, concurrent jobs, queue length, finished-job
  retention days, keep-original flag, original retention days, local-file-session auto-restore flag.  **No prices anywhere.**  All numbers are provisional placeholders (free = long-standing defaults 20 GiB / 3 / 20; pro
  100 GiB / 4 / 50; max 500 GiB / 8 / 100; ultimate 2 TiB / 16 / 200; retention 30 / 90 / 180 / forever days; originals kept from pro up for 180 / 365 / forever days) and must be decided by the owner (open decision 10, 13.14).
* **Behaviour already wired to it.** account limits are copied from the plan row when the plan is set; retention (job outputs and original uploads) follows the plan; originals are kept from pro up; sessions recorded on
  files opened from the user's own computer restore automatically from a copy kept with the session on pro and up, a free account gets a notice ("유료 계정에서 원본 보관 시 자동 복원").
* **Not done, deliberately:** windows/allowances per plan (13.6), holds and settlement, top-ups, taxes, invoices, trial logic, plan expiry or proration, a user-facing plan page (the account menu only shows the plan),
  per-plan rate limits.  The first administrator's starting plan is an administrator setting (default free).
* **Compatibility note for 13.9:** because plan and ledger are separate, changing a plan never rewrites ledger entries; a billing layer that wants plan history can read the audit log.
