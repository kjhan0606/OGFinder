#!/usr/bin/env python3
"""Generate docs/architecture.md: handwritten text + tables derived from the code (audit.json, plugin.json files, menu dump)."""
import os as _os
_HERE=_os.path.dirname(_os.path.abspath(__file__)); _ROOT=_os.environ.get('OGF_ROOT',_os.path.dirname(_HERE))
import json, os, glob
R=_ROOT
audit=json.load(open(_HERE+'/audit.json'))
base_layout_lines=12958
cur_layout_lines=sum(1 for _ in open(R+'/ds9/library/layout.tcl'))
feat_table=open(_HERE+'/feature_table.md').read()
verif=open(_HERE+'/verification.md').read()

def table_inventory():
    # menu inventory from the baseline dump
    rows=[l.rstrip('\n').split('\t') for l in open(_HERE+'/data/menu_baseline.tsv',encoding='utf-8')]
    top={}
    for path,t,lab,cmd,var in rows:
        m=path.split(' > ')[0].strip()
        top.setdefault(m,[]).append((path,t,lab,cmd))
    out=['| Old top-level menu | leaf entries | procs behind them (first ones) |','|---|---|---|']
    for m,es in top.items():
        leaves=[e for e in es if e[1]!='cascade']
        procs=[]
        for e in leaves:
            c=e[3].split()[0] if e[3] else ''
            if c and c not in procs: procs.append(c)
        out.append(f'| {m} | {len(leaves)} | {", ".join(procs[:6])}{" ..." if len(procs)>6 else ""} |')
    return '\n'.join(out)

def table_var_usage():
    import collections
    W=collections.defaultdict(set);Rd=collections.defaultdict(set)
    for f,v in audit.items():
        for k in v['write']: W[k].add(f)
        for k in v['read']: Rd[k].add(f)
    keep=['alldata','tbl','tbldb','sel','hover','cache','status','param','psf','icl','lsbg','markall','merge','ai','trim','visible_mode','add_objects_mode','sort','filename','detached','photoz','sed','bd','plot','morph']
    out=['| `catpanel(<prefix>,...)` | written by | read by |','|---|---|---|']
    for k in keep:
        out.append(f'| `{k}` | {", ".join(sorted(W[k])) or "-"} | {", ".join(sorted(Rd[k])) or "-"} |')
    return '\n'.join(out)

def table_classes():
    out=['| Plugin | Step | Recorder step name(s) | Class | Proc / CLI |','|---|---|---|---|---|']
    for f in sorted(glob.glob(R+'/plugins/*/plugin.json')):
        m=json.load(open(f))
        if m.get('example'): continue
        for s in m['steps']:
            cls=s.get('session','NONE')
            if s.get('variants'):
                pr=', '.join(v['proc'] for v in s['variants'][:2])+(' ...' if len(s['variants'])>2 else '')
            else:
                pr=s.get('proc','(cli)')
            out.append(f"| {m['id']} | {s['label']} | {', '.join(s.get('records',[])) or '-'} | {cls} | `{pr}` |")
    return '\n'.join(out)

def table_calls():
    core={'ds9-layout','catalog-core','sessionlog','cli-script'}
    out=['| Feature | reaches into | via |','|---|---|---|']
    for f,v in audit.items():
        c={}
        for k,n in v['calls'].items():
            ft,p=k.split(':')
            c.setdefault(ft,set()).add(p)
        nc={k:sorted(x) for k,x in c.items() if k not in core and k!='ds9-layout'}
        for ft,ps in nc.items():
            out.append(f"| {f} | {ft} | {', '.join('`%s`'%p for p in ps[:4])}{' ...' if len(ps)>4 else ''} |")
    return '\n'.join(out)

doc=f'''# OGFinder architecture

Status: stages (a) and (b) of the restructuring are implemented (core layer, plugin registry, declarative dialog,
workflow tabs, image tab strip with Tile/Single, one time-domain table); stage (c) (moving feature code out of
`layout.tcl` into `plugins/`) is done for Moving, AI services, Bands and Mask and **not** for the rest (section 6).
Sections 1-4 are an audit of the code *before* the restructuring (all numbers measured on the tree at commit `d6d9dff82`,
`layout.tcl` = {base_layout_lines} lines); section 5 describes the target architecture and what is implemented;
section 6 is the remaining migration plan.  Tables marked "generated" come from `tools/` scripts that parse the Tcl
sources and the `plugin.json` manifests, not from memory.

## 1. Files and sizes (before)

| File | lines | role |
|---|---|---|
| `ds9/library/layout.tcl` | {base_layout_lines} | DS9 window layout **and** the whole catalog panel: 250 procs, 11 menus, extraction, table, markers, PSF, ICL, LSBG, ... |
| `ogf_session.tcl` (+ `ogf_session_template.py`) | 594 (+1409) | session recorder `OGFSessLog`, export of the GUI session as a Python pipeline script |
| `ogf_moving.tcl` (now `plugins/moving/moving.tcl`) | 582 | Moving Objects menu, dialogs, tables, orbit window |
| `ogf_ai.tcl` (now `plugins/ai_services/ai.tcl`) | 566 | AI service bridge GUI |
| `ogf_mask.tcl` (now `plugins/mask/mask.tcl`) | 539 | shared mask manager |
| `ogf_bands.tcl` (now `plugins/bands/bands.tcl`) | 475 | multi-band registry, forced photometry, tile view |
| `ogf_link.tcl` | 398 | image <-> table link (click, hover, multi-select, keys) |
| `ogf_util.tcl` | 103 | `OGFPython`, `OGFForm`, `OGFTextWindow`, ... |

Python back ends are `ds9/library/ds9_*.py` (one per computational feature) plus the packages `moving/`,
`ai_bridge/`, `icl/`, `lsbg/`, `sersic_fit/`, ... at the repository root; the C binary `bin/ds9_sextract`.

## 2. Menus (before)

The catalog panel had a flat menubar of **11 top-level menus** (`SExtractor, Objects, Galaxy Model, Star(PSF),
Deconvolution, Bands, Mask, ICL, LSBG, Analysis, Moving Objects`) which requested 620 px and so fixed the width of the
right pane; next to the 12 standard DS9 menus of the left pane.  Generated inventory (measured by walking the real
menu widgets in a running ds9, `fits/menu_baseline.tsv`):

{table_inventory()}

## 3. Procs per feature (before; line ranges in `layout.tcl` unless another file is named)

Generated by `tools/audit.py` (proc-name groups; `ds9-layout` is the rest of the file: canvas, view/tile layout,
theme glue).  Ranges are inclusive line numbers of `proc` bodies in the baseline tree.

{feat_table}

Note the interleaving: ICL and LSBG procs are scattered over 7755-11946 because the CLI export/import code
(`CatalogPanelExportCLIScript`, `CatalogPanelImportCLIScript`, 8744-9593) and the shared helpers sit between them.

## 4. Coupling

### 4.1 Global `catpanel(...)` variables per feature

One array (`::catpanel`) holds *everything*: table state, per-feature parameters, status text, selection.  Generated
from the proc bodies (`set catpanel(x,...)` = write, any `catpanel(x,...)` = read):

{table_var_usage()}

Observations (all from the table above):
* `catpanel(status)` is written by every feature (23 groups) - there is no status service.
* `catpanel(alldata)` (the catalog as one TSV string) is read by 20 groups and written by 7: catalog access is
  direct string parsing in each feature.  `ogf_link.tcl` additionally keeps a parsed cache keyed on a write trace.
* `catpanel(icl,*)` / `catpanel(lsbg,*)` are read by the **mask** manager (`OGFMaskInit` copies defaults from
  `icl,param,*`; `OGFMaskAfterEdit` sets `icl,has_mask` / `lsbg,has_mask`), i.e. the supposedly shared mask service is
  coupled to two pipelines.
* `catpanel(psf,*)` is read by Deconvolution, Bulge+Disk, PSF photometry (state of the Star/PSF feature).
* Parameters of each feature live in `catpanel(<f>,param,<name>)` with a hand-written `*ParamLoad/*ParamSave`
  pair and a hand-written settings dialog (9 dialogs, 3 shapes: `ed()` copy + Apply, `ttk::spinbox` + OK, `OGFForm`).

### 4.2 Calls between features (excluding the catalog core, session recorder and the menu code in `layout.tcl`)

{table_calls()}

Plus: `ds9.tcl` creates `ds9(catalog_frame)`; `frame.tcl` calls `CatalogPanelHover/LinkPress/LinkShiftClick/
LinkRelease/StepKey` (the image <-> table link entry points).  The four AI-capable steps (`CatalogPanelGalaxyMorphology`,
`CatalogPanelStarFinder`, `CatalogPanelPhotoZ`, `CatalogPanelSEDFit`) start with
`if {{[OGFAIBackendHook ...]}} return` - the AI bridge hooks into other features by name.

### 4.3 Session recorder classification (AUTO / CONFIG / MANUAL)

`OGFSessLog STEP CLASS ARGV ...` is called from 41 places.  AUTO steps are replayed by the exported script in
`--mode pipeline` on new data, CONFIG steps (band registration) only in replay mode, MANUAL steps (hand edits, BCG
centre pick, ...) only with `--include-manual` / replay.  Generated from the manifests (`session` field was
transcribed from the `OGFSessLog` call sites above, see the `records` column for the recorder step names):

{table_classes()}

## 5. Target architecture

```mermaid
flowchart TB
  subgraph UI["Workflow UI (ogf_ui.tcl) - built from the registry"]
    MB["menubar: Workflow | Tools + progress strip + Stop"]
    TABS["tabs: Detect | Classify | Measure | Low-SB | Time-domain | Results"]
    CHIP["one chip per plugin: [Name v] [>] [gear]"]
    DLG["declarative parameter dialog (ogf_dialog.tcl)"]
  end
  subgraph CORE["core layer (ogf_core.tcl)"]
    REG["plugin registry\\nOGFRegisterPlugin / plugin.json loader"]
    CAT["::ogf::cat  catalog table"]
    MASK["::ogf::mask  mask service"]
    BANDS["::ogf::bands  band manager"]
    SESS["::ogf::session  recorder (OGFSessLog)"]
    AI["::ogf::ai  AI bridge"]
    JOB["::ogf::job  runner: progress, Stop, recorder"]
    PAR["::ogf::params  parameter store, presets, JSON"]
    LOG["::ogf::log"]
  end
  subgraph PLUG["plugins/<id>/  (plugin.json + optional .tcl + python driver)"]
    P1[extract] --- P2[star_psf] --- P3[morphology] --- P4[icl / lsbg] --- P5[moving] --- P6[example_hello]
  end
  subgraph LEG["legacy (to be migrated)"]
    L1["layout.tcl procs CatalogPanel*"]
  end
  PLUG -->|register| REG
  REG --> TABS
  REG --> CHIP
  PLUG -->|params| PAR --> DLG
  PLUG -->|steps| JOB --> SESS
  JOB --> CAT
  PLUG --> CAT & MASK & BANDS & AI
  CHIP -->|proc steps| L1
  L1 --> SESS
  SESS --> EXPORT["exported Python pipeline script"]
```

### 5.1 What exists now

* `ogf_json.tcl` - JSON reader (the embedded Tcl has no tcllib json; also accepts the `NaN`/`Infinity` Python writes).
* `ogf_core.tcl` - logging, the services `::ogf::cat/mask/bands/session/ai`, `::ogf::job`, `::ogf::params`,
  `::ogf::reg` (registry: discovery of `plugins/*/plugin.json`, validation, enable/disable, `OGFRegisterPlugin`,
  `"required": 1` plugins cannot be disabled), `::ogf::step` (argv templating incl. `{{plugin_dir}}`, running a plugin
  step through the job runner and the recorder, stage progress).
* `ogf_dialog.tcl` - `OGFParamDialog`: one dialog builder for all parameter specs (groups, Basic/Expert, validation,
  presets, Reset, Apply/OK).  Used by extract, deconv, morphology (bulge+disk), example_hello.
* `ogf_ui.tcl` - menubar (`Workflow`, `Tools`), tab bar, plugin chips, progress strip + Stop.
* `ogf_tile.tcl` - the image tab strip with `Single` / `Tile` / `Layout` (section 5.4) and marker drawing in every tile.
* `ogf_td.tcl` - the time-domain table: kind registry, filter, sort, save, markers (section 5.3).
* `plugins/*/plugin.json` - 16 built-in plugins that describe every existing step (+ `example_hello`, disabled by
  default).  For the plugins that are not migrated the steps call the *existing* procs (thin wrappers): computation,
  CLI argv and recorder calls are untouched.
* `plugins/moving/moving.tcl`, `plugins/ai_services/ai.tcl`, `plugins/bands/bands.tcl`, `plugins/mask/mask.tcl` - the
  former `ogf_moving.tcl`, `ogf_ai.tcl`, `ogf_bands.tcl`, `ogf_mask.tcl`, moved with `git mv` (history kept) and sourced
  through the manifest's `"tcl"` field.  Their proc names are unchanged, so `layout.tcl`, the recorder and
  `scripts/verify_ai_gui.tcl` keep working.

### 5.2 Dependency rules

1. A plugin may call `::ogf::*` and its own procs.  It must not read `catpanel(...)`, `ogfmask(...)`, `ogfband(...)`
   directly (legacy procs still do; they are the migration backlog, section 6).
2. Core services may call legacy procs (they are adapters over `layout.tcl`/plugin files today), never plugins.
3. Everything the recorder must replay goes through `OGFSessLog` (directly or via `::ogf::job::run` /
   `::ogf::step::run`), so exported scripts stay valid.
4. The UI is generated from the registry; no UI code names a feature.

### 5.3 One table for galaxies and time-domain objects

`catpanel(alldata)` (the galaxy catalogue, used by every galaxy feature and by the recorder) is **never modified** by the
time-domain kinds.  `::ogf::td::register_kind` (called by the `moving` plugin's init for `moving`, `transient`,
`detection`) gives a kind a prefix, marker colour and column list; `::ogf::td::show KIND` fills the same Tktable with the
rows of that kind (`All` = galaxies + movers + transients, common columns plus `kind`).  Row selection, sort, filter and
Save use the same widgets and callbacks as for galaxies (non-galaxy sort/save are logged as informational `note` steps).
Extra columns: `host_id`/`host_sep` (transients, nearest galaxy within 5 arcsec), `overlap`/`overlap_id` (movers).
Markers are fk5 `point(...)` regions tagged `ogf_td`/`ogf_td.<KEY>` in every frame; selection reuses the `sextract_sel`
tag.  The table geometry (y=181, h=769, info 154) is the same in every view.

### 5.4 Image tabs and Tile (ds9's mosaic is kept)

`ogf_tile.tcl` adds a 24 px strip under the image canvas: one tab per frame (band name or file name; click = go to the
frame) and on the right `Single`, `Tile`, `Layout` (grid / column / row, lock pan-zoom by WCS, share scale limits) and a
frame list.  `Tile` is ds9's own one-canvas `tile` mode (`OGFUIDisplay tile`), `Single` returns to tab view
(`OGFUIDisplay single`); the right-click item *Tile all* calls the same function.  Catalog markers are drawn into every
registered band frame (`OGFMarkSend`/`OGFMarkDelete`, positions mapped with `OGFBandPos`) and a click on a marker in any
tile (`Button1Frame` -> `OGFTileClick`) selects the same table row.

### 5.5 Menus after

Top level of the catalogue panel: **Workflow | Tools** (+ the progress strip and Stop), then the 6 tabs
**Detect | Classify | Measure | Low-SB | Time-domain | Results**, each showing one chip `[Name v][>][gear]` per plugin.
Measured by walking the real menu widgets (`fits/menu_new3.tsv`, compared by `tools/compare_menus.py`): all 145 leaf
entries of the old 11 menus are reachable (144 as the same command, the old `Show Results Table` is replaced by the kind
filter of the shared table, the three old settings dialogs by `OGFParamDialog`); 23 entries are new (tabs, tile options,
plugin dialogs, details window).  Panel width is 559 px instead of 620 px.

## 6. Remaining migration

Done in this pass: **moving**, **ai_services**, **bands**, **mask** (own file under `plugins/<id>/`, loaded through the
manifest), **example_hello** (new, with a Python step), and removal of the superseded settings dialogs from `layout.tcl`
(`CatalogPanelSettingsDialog`, `CatalogPanelDeconvSettings*`, `CatalogPanelBulgeDiskSettings`, `CatalogPanelBDSettingsApply`
and the PSF compatibility aliases).  `layout.tcl` went from {base_layout_lines} to {cur_layout_lines} lines.

Not migrated (the manifests wrap the legacy procs which still live in `layout.tcl`; line ranges in section 3):
`extract` (+ dual image, trim), `merge`, `ai-merge`, `galaxy-model`, `star-psf`, `deconv` (compute part), `separate/edit`,
`icl`, `lsbg`, `measure-misc` (morphometry, Sersic, PSF photometry, ...), `photo-z`, `sed`, `bulge-disk` (compute part),
`cli-script` (CLI export/import), `plot/viewer`, and the catalogue core (table, markers, filter, sort, save).  The remaining
settings dialogs (`CatalogPanelStarPSFSettings`, `...SeparateSettings` (objects), `...ICLSettings`, `...LSBGSettings`,
`...MaskOverlaySettings`, `OGFAIRegistry`) are still hand-written.  All of them still read `catpanel(...)`.
Suggested order: photo-z/sed, bulge-disk and deconv compute, star-psf, galaxy-model, then icl/lsbg (largest, interleaved
with the CLI export code), finally extract + catalogue core behind `::ogf::cat`.

## 7. Verification (measured on the final build, see the final report)

{verif}
'''
open(R+'/docs/architecture.md','w').write(doc)
print(len(doc.splitlines()),'lines')
