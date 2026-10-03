# Phase II - OGFinder on the web (design)

Status: **design proposal, nothing of it is implemented.**  Revision 2 (2026-10-03): the UI decision is taken - **shared UI components plus
user-selectable layout presets** (sections 6.1, 6.4-6.7, 9, 10).  The three mock-ups behind it are in `/workspace/webdesign/` (`README.md`,
`proposal_A.html`, `proposal_B.html`, `proposal_C.html`, PNG screenshots); they are not part of this repository and use example data only.  Target repository: `github.com/kjhan0606/ds10-web` (new, to be created by
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
  light-curve/tracklet dock, job list and status bar are built once and arranged by a *layout preset* that the user chooses
  (Classic = ds9-like, Guided = pipeline stepper, Workspace = icon rail + drawer + dock).  A preset arranges components and owns no analysis logic (section 6).

Non-goals for Phase II: full ds9 feature parity (regions language, XPA/SAMP, colourmap editor, 3-D cubes beyond a slider), a
multi-tenant public service on day one, replacing the desktop version (it stays the reference implementation).  Team sharing is a
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

Every visible part of the application is one of eight components.  A component takes **props = selectors on the state store (6.4)** and
emits **actions** (commands); it never imports another component, never calls another component's methods and never reads the DOM of a
neighbour.  The behaviour the desktop UI already defines (what is linked to what) is the specification; the Tk code is not reused.

| # | component | what it does (desktop origin) | data it needs (API, section 5) | store slices it reads / writes |
|---|---|---|---|---|
| 1 | **Image viewer** | tile pyramid + WebGL2 layers (section 4): pan/zoom, scale/colormap, cursor readout, WCS, catalog markers/ellipses, mask and trail overlays, **Single / Tile / split-compare** (several canvases or one canvas with a split uniform, shared view state; `ogf_tile.tcl`) | `ImageSource` (local or `/tiles`), `/header`, `/stats`, `/pixel`, mask tiles, catalog `bbox` rows | reads `view`, `catalog` (markers), `selection`; writes `view`, `selection` (marker click; the `ogf_pick.tcl` chooser is a popup of this component) |
| 2 | **Catalog table** | virtual-scrolling grid, kind filter (Galaxy / Moving / Transient / All), search, column filters, sort, multi-select, review column (`ogf_td.tcl`, `ogf_link.tcl`) | `/catalogs/{id}/rows`, `PATCH .../cells`, `.../edit` | reads `catalog`, `selection`, `view` (visible only); writes `catalog.filter/sort`, `selection`, cell edits |
| 3 | **Parameter panel** | form generated from `plugin.json` `params`: groups, Basic/Expert, range checks, presets, Reset (`OGFParamDialog`, 6.2).  Used as a modal dialog, a docked drawer or an inline "cell" - the same component | `/plugins/{id}`, `POST /jobs` | reads `plugins`, `params`; writes `params`, dispatches `job.submit` |
| 4 | **Object inspector** | selected object: cut-out, Sersic model, residual, radial profile, p(z), catalog row, review buttons + note, provenance | `/cutout`, job outputs (`multifit`, `isophote`, `photoz_sed`), `PATCH .../cells` | reads `selection`, `catalog`; writes review cells |
| 5 | **Command palette (Ctrl-K)** | fuzzy list of every manifest step, view action, session action and "go to object #n"; arrow keys / Enter; same entries as the menus because both read one **command registry** | `/plugins` (registry built from the manifests at start-up) | reads `commands`; acts by dispatching the chosen command; `ui.paletteOpen` |
| 6 | **Light-curve / tracklet dock** | light curve of a transient (`lightcurves` viewer), tracklet and orbit summary of a moving object (`moving`), chosen from the table | job outputs of `moving` / `lightcurves`, catalog rows | reads `selection`, `jobs` (outputs) |
| 7 | **Job list** | running/queued/finished jobs with progress, Stop, exact `argv`, log tail (progress strip of `ogf_ui.tcl`) | `/jobs`, SSE `/jobs/{id}/events` | reads `jobs`, `session`; dispatches `job.cancel` |
| 8 | **Status bar** | WCS and cursor value, local/server badge, zoom, running job, session-recording indicator ("REC n steps"), tile/cache state | none of its own | reads `view`, `source`, `jobs`, `session` |

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
  jobs:     { byId: Record<string, Job>; order: string[] };  session: { steps: SessionStep[]; recording: boolean };
  ui:       { preset: "classic"|"guided"|"workspace"; theme: "dark"|"light"; locale: "en"|"ko";
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
| `ui` | preset shell, palette, theme/locale controls | the only slice persisted in the browser (6.6) |

Switching a preset changes **only `ui.preset`**: the same store feeds a different arrangement, so selection, view, scale, catalog filter,
running jobs and unsaved parameter edits survive a switch (acceptance test in 6.7).

### 6.5 Layout presets

A preset is a shell: a CSS-grid template, the list of components it mounts with their slots and minimum sizes, and a few **shell widgets**
that exist only to navigate (built from manifests and commands, holding no analysis logic).  Presets are the only place where components
are arranged; none of them forks a component.  The three mock-ups in `/workspace/webdesign/` illustrate them.

| | **Classic** (`proposal_A`) | **Guided** (`proposal_C`) | **Workspace** (`proposal_B`) |
|---|---|---|---|
| idea | ds9 and OGFinder Phase I, almost no re-learning | a pipeline you can follow and replay: *what do I do next?* | modes + tool drawer + dock for multi-band, time-domain and tablets |
| default theme | dark | light | dark, light toggle |
| shell widgets | menu bar (File ... Analysis ... Window); workflow tabs *Detect / Classify / Measure / Low-SB / Time-domain / Results* with a **plugin ribbon** (one chip per plugin: primary step, step list, gear = settings); info panel (File, Object, Value, fk5 alpha/delta, x, y, Zoom); button bar; frame tabs (Single/Tile/Blink); colorbar | left **stepper** (Load, Calibrate, Detect, Measure, Review, Report) with sub-steps, status and progress from `session.steps`/`jobs`; step **"cell" header** (title, parameter chips, Run); result strip (counts, histogram); recorded **session log**; split-compare toolbar | left **icon rail** of modes (Detect, Morph, LSB/ICL, Moving, Photo-z, Trails, AI, Stack); **collapsible drawer** (steps of the mode + parameter form); top bar (file chip, Local / Analyse-on-server switch, theme, language); floating viewer toolbar |
| components | viewer (1), catalog table (2) with Galaxy/Moving/Transient/All tabs, parameter panel as a dialog (3), a reduced inspector as selected-object summary (4), palette (5, optional), job list + status in the status bar (7, 8) | viewer with split compare (1), inspector as right column (4), parameter panel inline in the cell (3), palette (5, always), job list inside the stepper (7), status/log (8); catalog table as a drawer/modal (2) | viewer (1), parameter panel in the drawer (3), table (2) and light curve/tracklet/jobs (6, 7) in the bottom **dock**, status bar (8), palette (5) |
| manifest fields used | `tab`, `order`, `short`, `primary`, `steps`, `params` (already present) | plus an optional **`ui.stage`** (`load\|calibrate\|detect\|measure\|review\|report`); default derived from `tab` (Detect -> detect; trails and mask -> calibrate; Classify, Measure, Low-SB, Time-domain -> measure; Results -> review/report) | plus optional **`ui.mode`** and **`ui.icon`**; default: the plugin grouping of the mock-up, the rest under a "More..." flyout |
| minimum window | 1280 x 720 | 1024 x 700 (panels collapse to drawers below 1280) | 1024 x 700 (rail only below 1100; **best tablet landscape**) |
| strengths | cheapest, ds9-familiar, best catalog density | review, reporting, teaching; the log is the exported script | layers, time-domain and jobs visible at once; touch-friendly |
| delivery | **Phase 1** (default preset) | **Phase 2** | **Phase 3**, after the viewer is stable (9.1) |

`ui.stage`, `ui.mode` and `ui.icon` are additive manifest fields: the desktop loader ignores unknown keys, and a plugin without them still
appears (derived defaults).  Plugins that fit no stage (`lensmodel`, `spectra`, `cluster`, `xmatch`, ...) are listed under *Other tools* in
Guided and are reachable through the palette in every preset.

**Deferred in Workspace v1** (backlog, not in the first Workspace release): RGB composite and per-band blend layers, the overview minimap, and free
docking (drag/split/float).  v1 uses one fixed, drag-resizable dock below the viewer and one band visible at a time.  Reason: these are the expensive
parts of proposal B (9.2) and neither of the other two presets needs them.

Preset switch: a "Layout: Classic | Guided | Workspace" control and a palette command ("Layout: Guided").  Switching re-mounts the shell, not the
viewer: the viewer's canvas element is *moved* into the new shell rather than re-created, so tiles stay in the texture cache without a re-fetch
(risk: browsers that lose a WebGL context on re-parenting, section 10).

### 6.6 Persistence of the layout choice (localStorage + URL parameter)

* **Key and value**: `localStorage["ogf.ui.v1"] = { preset, theme, locale, layouts: { classic: {...}, guided: {...}, workspace: {...} } }`.
  `layouts[name]` holds only that preset's panel state (panel widths, collapsed drawers, dock height, open stepper steps).
  Nothing else - no file names, ids, paths, catalog data or tokens - is written to browser storage.
* **URL parameter**: `?layout=classic|guided|workspace` (also `&theme=dark|light`, `&lang=en|ko`).  Precedence at start-up: **URL parameter >
  localStorage > default (`classic`)**.  A URL value applies to that page load and is written to localStorage only when the user then switches the
  layout explicitly, so following a link never silently changes the stored preference.  The URL is updated with `history.replaceState` when the user
  switches, so a copied address reopens the same layout.  There is no per-file or per-session link in Phase II (section 7.2); these are the only
  parameters the application reads.
* **Robustness**: unknown or invalid values fall back to the default with one console warning; the stored object is schema-versioned (`v1`) and
  validated field by field (a bad `layouts.guided` resets only that preset); storage that is disabled, full or blocked (private mode, enterprise
  policy) falls back to an in-memory store with a one-line notice; a `storage` event from another tab updates the *stored* preference but does not
  re-layout an open tab.  A preset whose minimum window size (6.5) is not met is offered but flagged, and a compatible one is suggested; it is never forced.
* **Reset**: "Reset layout" clears `layouts[current]`; "Reset all UI settings" removes the key.  A server-side copy per account is a later option (section 10).

### 6.7 Verification plan (DOM/layout checks per preset)

The mock-ups were checked with a Playwright/headless-Chrome script (`/workspace/webdesign/_src/shot.py`, `extra_checks.py`, `overlap.py`); it is the
starting point of the CI job.  The same checks run **for every preset x every tested state**, at 1440x900 (all presets), 1280x720 (all),
and 1024x768 / 1180x820 tablet landscape (Guided and Workspace).  Already measured on the mock-ups at 1440x900: 0 overflowing or out-of-bounds elements,
all images loaded, no console errors, no external requests, no overlapping floating pieces.

| layer | check | applies to |
|---|---|---|
| layout | `document` scroll size equals the viewport (no page scroll); **no element overflows its box or lies outside the viewport**, except containers marked `data-scroll` / `data-clip` (table body, viewer, log, step list); no clipped text | all presets, all states |
| layout | no overlap between floating pieces (viewer toolbar, layer/selection cards, toast, readout, split handle, palette) | Guided, Workspace |
| layout | minimum-size contract of each component: render it at its declared minimum and assert the collapse behaviour | component tests |
| images | every viewer/cut-out image or texture loads (`naturalWidth > 0`, WebGL context created, first tile drawn); the cursor readout over the viewer returns RA/Dec and a pixel value | all presets |
| hygiene | no console errors/warnings; **no request to an external host** (only the application's own API); no CSP violations | all presets |
| state | the same scripted scenario in each preset ends with the same store snapshot for `selection`, `view`, `catalog.filter`, `params` (select an object in the table -> marker highlighted, inspector/summary updated, dock shows the light curve for a transient and the tracklet for a moving object) | all presets |
| switching | Classic -> Guided -> Workspace -> Classic with a selected object, a changed scale, a running job and an edited, unsaved parameter: all survive; the viewer canvas is the same element (identity test) and no tile is re-fetched | all |
| persistence | choose a preset, reload: same preset (localStorage); open `?layout=guided` with a stored `classic`: Guided is shown and the stored value is unchanged; invalid parameter -> default + one warning; blocked storage -> in-memory fallback, no exception; corrupted JSON -> reset | all |
| i18n | long Korean strings (breadcrumb, step names, tooltips) and a 30-character English parameter label do not break the checks above; Noto Sans CJK / Nanum Gothic / system fallback present | all, `lang=ko` |
| theme | dark and light: contrast of text and of the marker colours on the image; the viewer background stays dark | Workspace, Guided |
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
* Computations on small local files that are cheap enough for the browser (statistics, cut-outs, simple aperture photometry) can run in a
  Worker; **extraction with `ds9_sextract` (C) runs on the server only**, unless a WASM build is made (option in section 8; the HUDF F160W full
  extraction takes 1.9 s natively, *measured*, so the WASM path is plausible but unproven).

### 7.2 Login, storage, retention and deletion - options

This is the area where the owner's decision is needed; the design supports all of them through a per-deployment policy file
(`policy.yaml`: `auth`, `quota_bytes`, `retention_days`, `max_upload_bytes`, `allow_anonymous`).

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
M0-M8 deliver the shared components and the **Classic** preset; **M9** adds Guided, **M10** Workspace.  "(+x)" is the extra cost of building with independent
components, a state store and preset plumbing instead of one fixed layout; it is already included in the range shown.

| M | content | deliverable and acceptance test | estimate (extra) |
|---|---|---|---|
| **M0** spikes | (1) FITS parse + WebGL2 render of the 52 MB HUDF image with scale/colormap, (2) tile endpoint over the existing numpy code, (3) WCS in WASM vs hand-written, (4) the 200 MB threshold on three browsers, (5) **state store + command registry skeleton and two throw-away shells** (Classic grid, Guided grid) around a stub viewer, to confirm the component contract and the canvas hand-over of 6.5 | go/no-go notes with measurements; decisions in section 8 closed (store library included) | 2 (0) |
| **M1** viewer + local files | `ImageSource`, `LocalFileSource`, WebGL renderer, pan/zoom/scale/colormap, pixel readout, WCS coordinates, image tabs, Tile mode; **viewer and status bar as stand-alone components** with the store contract (6.1, 6.4); split-compare uniform (small; used from M9) | open m51 and the HUDF crop offline; scale changes at 60 fps on a 4k tile set; test against ds9 pixel values (exact); component tests at minimum size | 6-8 (+1) |
| **M2** server core | auth (OIDC), resumable upload (5.1), ingest job, pyramid, tile API (5.2), `RemoteTileSource`, quotas, retention purge, delete | upload a 1.6 GB FITS over a throttled connection with a forced disconnect; resume; tiles identical to the local path for the same file | 5-6 (0) |
| **M3** jobs + plugins | job runner, manifests API, **parameter panel and job list components** (modal / drawer / inline hosting), `cli` blocks for all `proc` analysis steps, SSE progress, cancel; steps `extract`, `mask` auto, photo-z/SED, morphology, ICL, LSBG | each plugin's CLI step reproduces the desktop output byte-identically on m51/HUDF (re-use `scripts/verify_*`) | 6.5-8.5 (+0.5) |
| **M4** catalog + table + **Classic preset** | catalog API (5.4), virtual grid, markers, selection link, kind filter, click chooser, review columns + filter, report export; **Classic shell** (menu bar, workflow tabs + plugin ribbon from manifests, info panel, button bar, frame tabs, colorbar); **preset plumbing: layout switch, localStorage + URL parameter (6.6)**; the 6.7 matrix for Classic | the `verify_review_gui`/`verify_click_*` scenarios re-written for the browser (Playwright); 6.7 matrix green for Classic at 1440x900 and 1280x720 | 6-7 (+1) |
| **M5** mask + editing | mask API (5.5), polygon/ellipse editing, objects plugin (delete/merge/separate/add/AI-merge) | operations produce the same catalog as `ds9_catalog_edit.py` golden cases (`scripts/golden/`) | 4-5 (0) |
| **M6** sessions | server-side recorder, export of the Python script (`ogf_session_template.py`), replay check; the session log later feeds the Guided stepper | exported script from a web session passes the same replay checks as the desktop (`verify_session_replay.py` equivalent) | 3-4 (0) |
| **M7** time-domain | `moving` plugin jobs, MAST cache, orbit/lightcurve views, Moving/Transient kinds in the table; **light-curve / tracklet dock component** (a bottom/side tab in Classic) | the BB89 HST field gives the same tracklet as the desktop (`moving/tests`, `scripts/run_moving_pipeline.py`) | 4.5-5.5 (+0.5) |
| **M8** hardening | sandbox, rate limits, backups, monitoring, load test, security review, documentation, accessibility pass; preset-matrix CI job (6.7) | pen-test checklist, 20 concurrent users on the reference VM | 4.5 (+0.5) |
| **M9** Guided preset | stepper with the `ui.stage` mapping and recorder-driven status, step "cell" header, result strip, session log, split-compare toolbar, **object inspector** (cut-out, model/residual from multifit/isophote outputs, radial profile, p(z), review + note), **command palette** (registry from the manifests), collapse rules for 1024-1280 px | 6.7 matrix green for Guided; the scenario "accept 20 objects, export report" runs identically in Classic and Guided and exports the same script | 3-4 |
| **M10** Workspace preset (v1) | icon rail with the mode grouping, drawer (steps + parameter panel), top bar, fixed resizable dock (table, light curve/tracklet, jobs), light/dark theme, tablet layout; **no** RGB composite, minimap or free docking | 6.7 matrix green for Workspace at 1440x900, 1280x720, 1180x820; touch test on a tablet | 3-4 |

Totals: **M0-M8 (Classic, shared components): 41.5-50.5** person-weeks (a single fixed layout was 38-48); with **M9** (3-4) and **M10** (3-4):
**47.5-58.5**.  A useful first release (M0-M3 + extraction + table, Classic only) is about **20-26** (was 18-24).

### 9.1 Phased roadmap

1. **Phase 1 - Classic (M0-M8).**  Default preset and the reference for all acceptance tests.  Built *as components* from the first day, so Phase 2
   needs no refactoring.  Release 1 = Classic.
2. **Phase 2 - Guided (M9)**, starts when M3, M4 and M6 are done (it needs the parameter panel, table/selection and the recorder).  It brings the inspector and
   the palette, which are then also offered in Classic (Ctrl-K, and the inspector as an optional right-hand panel).  Release 1.x.
3. **Phase 3 - Workspace (M10)**, starts only when the **viewer-stable gate** is met: the `ImageSource`/renderer API unchanged for two releases, M1 acceptance
   tests green on the three target browsers, tile cache and WebGL context survive a canvas move (preset switch) without leaks, and the multi-band layer API
   specified (even if RGB/minimap are not built).  Backlog after M10: RGB composite with per-band blend, minimap, free docking.

### 9.2 Cost estimate: about +20-30 % over a single layout

Baseline (one fixed layout) 38-48 person-weeks.  Extra for the decision: component/store discipline and preset plumbing inside M1, M3, M4, M7, M8 = about
**2.5-3.5**; Guided (M9) **3-4**; Workspace v1 (M10) **3-4**; total roughly **9-11.5 person-weeks**, i.e. **about +20-30 %** of the baseline (about +21-27 % of its
mid-point, 43).  The components are built once; the extra is the shells, the test matrix (6.7) and the two navigation widgets (stepper, rail/drawer).
What would raise it: Workspace with the deferred RGB/minimap/docking (+4-6, not included), a fourth preset, or per-preset forks of a component (excluded by
the rules in 6.1).  What would lower it: shipping Guided without the split-compare toolbar or the inspector's model/residual tabs in its first release.

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
* **Tablet and touch** are designed for Guided and Workspace only; Classic is desktop-only (1280 x 720 minimum); phone layouts are out of scope for Phase II.
* **Korean UI text**: glosses are partial (stage names, a few labels); a full translation, input-method (IME) handling in the palette/search and CJK line breaking in tables are not budgeted.
* **Mock-ups are not specifications**: they use example data and plain DOM (no WebGL); sizes, colours and wording will change.  What is decided is the component split and the state/preset rules above.
* **Open decisions for the owner**: login method, quota/retention numbers, the first deployment target (institutional server vs cloud), whether the web version must support offline use (PWA), whether the layout preference also belongs to the account (server copy), the default preset for first-time users (Classic is assumed), licence, and the name of the repository (`ds10-web` is assumed here).

## 11. First concrete steps (when the owner approves)

1. Create `ds10-web` (owner action); copy `plugins/*/plugin.json`, `ai_bridge/`, `moving/`, `ds9/library/ds9_*.py`, `ogf_session_template.py` as a **git submodule or a vendored directory with a sync script** so the desktop repo stays the single source of truth for algorithms; do not fork them.
2. M0 spikes (section 9), including the state-store/command-registry skeleton and the Classic/Guided throw-away shells.  Copy the three mock-ups from `/workspace/webdesign/` (and the DOM-check scripts in its `_src/`) into `ds10-web/docs/ui/` as the visual reference and the first golden screenshots (owner action; not done here).
3. Write the OpenAPI document for section 5 first; generate the TypeScript client and the FastAPI stubs from it, so that the browser and server work in parallel.
4. Add a `docs/` index in the new repository that points back to this file.
