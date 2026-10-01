# Writing an OGFinder plugin

A plugin is a directory `plugins/<id>/` with a `plugin.json` manifest, optionally a Tcl file and a Python driver.
The workflow tabs, chips, menus, parameter dialogs, progress strip, session-recorder entries and the
"Pipeline overview" are all generated from the manifests; no UI code names a feature.  The worked example is
`plugins/example_hello/` (a parameter dialog, a CLI step that adds a catalog column, and a Tcl step).

## Five steps

1. **Create** `plugins/<id>/plugin.json` (reference below).  `id` is the directory name; `tab` is one of
   `Detect | Classify | Measure | Low-SB | Time-domain | Results`.
2. **Declare the parameters** in `"params"` (name, label, type `int|float|string|bool|choice|file`, default, group,
   help, optional min/max/unit/choices/`expert`).  One declarative dialog (`OGFParamDialog`, the gear chip /
   Tools > Plugin settings) edits them: groups, Basic/Expert toggle, help line, range validation, Reset, presets.
   Values persist in `~/.ds9/ogf_params/<id>.json` (or in a legacy array, see `store`).
3. **Declare the steps** in `"steps"`.  A step is either a *CLI step* (`"cli": [argv template]`; the job runner
   runs it asynchronously, the session recorder logs the exact argv, `output` says how stdout is applied) or a
   *Tcl step* (`"proc": "MyProc"`; use this to wrap existing code).  Put the driver script next to the manifest;
   `{plugin_dir}` expands to the directory on disk.
4. **(Optional) Tcl** in `"tcl": "x.tcl"`: sourced once at start-up (only when the plugin is enabled).  Talk to the
   application through `::ogf::*` (core API below), not through `catpanel(...)`.
5. **Test**: start with `OGFINDER_PLUGINS=<id>` for plugins that are disabled by default (`"example": true`), look
   at Tools > Plugins..., run the step, check Workflow > Session recorder > Export, and `scripts/verify_*`.

Plugin search path (later directories override an equal `id`): `<ds9 zip>/plugins` (built in), `<root>/plugins`
(next to `bin/`), `~/.ds9/plugins`, every directory in `$OGFINDER_PLUGIN_PATH` (`:` separated).
Enable/disable: `$OGFINDER_PLUGINS=a,b,-c` (comma list; a leading `-` forces a plugin off), or
`~/.ds9/ogf_plugins.json` `{"enabled":[...],"disabled":[...]}` (also written by Tools > Plugins...).
A broken manifest is rejected with a message in Tools > Plugin log; other plugins are unaffected.

## `plugin.json` reference

| key | meaning |
|---|---|
| `schema` | always `1` |
| `id`, `name`, `short` | identifier, long name (Plugins window), chip label |
| `tab`, `order` | workflow tab, sort key inside the tab (chips left to right) |
| `description` | tooltip and "Pipeline overview" |
| `example` / `enabled` | `"example": true` = disabled unless enabled explicitly; `"enabled": false` likewise |
| `menu` | `"Workflow"` puts the plugin's entries into that main menu instead of a chip |
| `requires` | `{"binaries":[...], "python":[modules]}`; missing items disable the step with a message |
| `network` | `1` when the plugin contacts a remote service (recorded; the exported script asks before running it) |
| `primary` | step id run by the chip's ▶ button |
| `settings` | `"params"` (declarative dialog) or the name of a Tcl proc to call from the gear button |
| `store` | `{"array":"catpanel","key":"bd,param,%s","save":"ProcName"}`: bind the parameters to an existing Tcl array instead of the JSON file (legacy code keeps reading `catpanel(...)`); `save` is called after the dialog applies |
| `tcl`, `init` | Tcl file to source (a string **or a list** of files, relative to the plugin dir); `init` proc called once the panel exists |
| `tcl_always` | `1`: the `tcl` files are sourced even when the plugin is disabled (the procs are still called by `layout.tcl`, the recorder or other plugins; only the chip/menu entries disappear) |
| `params` | list of parameter specs (above) |
| `dialog_tabs` | `1`: the settings dialog is a notebook with one tab per parameter group (long dialogs: star-psf, ICL, LSBG) |
| `on_open` | Tcl proc called before the dialog reads its values (lazy initialisation of the legacy store) |
| `on_apply` | Tcl proc called after Apply/OK has validated, stored and saved the values |
| `choices_proc` | on a parameter: Tcl proc returning the list of choices of a `choice` parameter at open time (e.g. AI backends) |
| `steps` | list of step objects (below) |

Step object: `id`, `label`, `title` (recorder title), `proc` **or** `cli`, `session` (`AUTO|CONFIG|MANUAL|NONE`),
`stage` (`detect|classify|measure`, drives the progress strip), `records` (recorder step names this step produces,
used for the strip and the overview), `needs` (`image`, `catalog`), `network`, `legacy` (the old menu path, for the
menu-parity table), `settings_only`, `variants` (a cascade of `{id,label,proc}`), `output`
(`{"mode":"add_columns","columns":[..]}` | `{"mode":"set"}` replaces the catalog | `{"mode":"text"}` shows stdout),
`after` (Tcl called when the job is done).

`cli` templates: `{python}` `{plugin_dir}` `{work}` (`~/.ds9`) `{root}` `{image}` (current FITS) `{catalog}`
(temporary TSV of the current catalog) and `{param-name}`; an element `{"if":"bool-param","argv":[...]}` is included
only when the parameter is true.  The argv is what the recorder stores: paths under `{root}` become `@{ROOT}/...`
in the exported script, so drivers must live on disk (not in the embedded zip).

Session classes: **AUTO** replayed by the exported script in pipeline mode on new data; **CONFIG** (band registry)
only in replay mode; **MANUAL** (hand edits) only with `--include-manual` or replay; **NONE** not recorded.
CLI steps go through the recorder automatically (`class` from `session`).  Tcl steps record themselves via
`OGFSessLog`/`::ogf::session::log` exactly as before.

## Core API (`ogf_core.tcl`)

| namespace | procs |
|---|---|
| `::ogf` | `log LEVEL MSG`, `status MSG` (status line) |
| `::ogf::cat` | `tsv`, `has`, `columns`, `nrows`, `values COL`, `rows`, `selection`, `select NUMS ?mode? ?pan?`, `clear_selection`, `add_columns TSV NAMES`, `load_tsv TSV NAME`, `temp_file SUFFIX`, `image_file` |
| `::ogf::mask` | `paths`, `exists`, `bool_path`, `run MODE ...`, `overlay` |
| `::ogf::bands` | `names`, `detect`, `path BAND`, `frames`, `sorted`, `register`, `set_detect`, `pos FRAME NUM X Y` |
| `::ogf::session` | `log STEP CLASS ARGV ...`, `set_field SEQ KEY VAL`, `steps`, `crc TEXT` |
| `::ogf::ai` | `run`, `services` |
| `::ogf::job` | `run ARGV ?-step -class -title -network -outputs -payload -post -requires -done -plugin?` (async, recorder entry written first, elapsed time, `cancel`, `busy`) |
| `::ogf::params` | `get PLUGIN NAME`, `put`, `specs`, `defaults`, `save`, `restore`, `preset_list/save/load`, `validate` |
| `::ogf::reg` | `get ID`, `ids`, `step ID SID`, `on_tab TAB`, `register DICT`, `missing ID`; Tcl plugins may call `OGFRegisterPlugin DICT` to register without a JSON file |
| `::ogf::td` | time-domain table (next section) |

`OGFParamDialog PLUGIN ?-step S?` opens the dialog from code.

## Time-domain table (`ogf_td.tcl`)

One table shows rows of several **kinds** (`galaxy | moving | transient | detection`, plugin-defined) with a kind
filter (Galaxies / Moving / Transients / All) on the Time-domain tab.  The galaxy catalog (`catpanel(alldata)`) is
never modified; the other kinds live in `ogftd(raw,KIND)`.

```tcl
::ogf::td::register_kind mykind -label "My objects" -prefix X -color orange -order 40 \
    -columns {n score host_id} -decorate MyDecorate -point {point=cross 12}
::ogf::td::set_rows mykind {id ra dec mag n score} $rows 1     ;# 1 = switch the table to this kind
proc MyDecorate {cols rows} {...; return [list $newcols $newrows]}   ;# computed columns, e.g. host_id
::ogf::td::galaxy_near $ra $dec $radius_arcsec     ;# {{number sep_arcsec extent_arcsec} ...} from the CURRENT galaxy catalog
::ogf::td::on_select MyProc                          ;# MyProc KEY RAWROW when a row of any non-galaxy kind is selected
```

Common columns (always): `NUMBER` (key; prefix + id), `kind`, `X_IMAGE`, `Y_IMAGE` (detection grid, from the
WCS), `ALPHA_J2000`, `DELTA_J2000`, `MAG_AUTO`; with one kind selected its own columns follow.  A `_pos` column
(`ra1 dec1 ra2 dec2 ...`) draws several markers per object (a moving object's track).  Markers are drawn in sky
coordinates into every frame/tile, with the kind's colour; a click on one selects the same row.  Sort, filter,
save, hover, multi-select and the session recorder behave as for galaxies (non-galaxy sort/save are recorded as
informational `note` steps; the exported script's galaxy catalog is not touched).

## Image tabs and tile

The strip under the image has one tab per frame plus `Single`, `Tile` and `Layout ▾` (grid / column / row, WCS lock,
share scale).  `Tile` is ds9's own one-canvas mosaic; the tabs never replace it.  Right-click on the strip:
"Tile all" / "Single", the same functions (`OGFUIDisplay tile|single`).  Catalog markers are drawn in every tile
(registered bands: positions mapped through the WCS); a click on a marker in any tile selects its table row.

## Migration status

See `docs/architecture.md` section 6.  All feature procs now live in `plugins/<id>/*.tcl` (loaded through `"tcl"`);
`layout.tcl` is 1873 lines (was 12958) and only holds the panel construction, detach/layout code and ds9 core.
`bands`, `mask` are `"required": 1`; `catalog`, `deconv`, `extract`, `galaxy_model`, `icl`, `lsbg`, `morphology`, `objects`, `photometry`, `photoz_sed`, `star_psf` use `tcl_always` because `layout.tcl`
and the recorder call their procs unconditionally.  Parameter dialogs generated from the manifest
(`"settings": "params"`): extract, deconv, morphology (bulge+disk), star-psf, objects, ICL, LSBG, mask overlay,
AI services, example_hello.  The plugins still read and write `catpanel(...)` / `ed()` globals (not decoupled).
