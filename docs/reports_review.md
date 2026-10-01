# Candidate review and reports (plugin `report`)

Two features in one plugin (`plugins/report/`): a per-candidate **review state** and an **Export Report...** step that writes
a self-contained HTML file (optionally a PDF).  Nothing is sent anywhere; there is no sharing feature.

## 1. Review state

Each candidate (a row of the galaxy catalog) can be **accepted**, **rejected** or marked **uncertain**, with a free-text note.
The state is stored as three ordinary catalog columns, so it travels with the catalog:

| column | values |
|---|---|
| `REVIEW` | `accept`, `reject`, `uncertain`, or empty (not reviewed) |
| `REVIEW_NOTE` | free text; tabs and line breaks are folded to a space |
| `REVIEW_TIME` | `YYYY-MM-DD HH:MM:SS` (local time) of the last decision |

They are appended at the right end of the table the first time a decision is made.  Because they are catalog columns they are
saved and loaded with *Save/Load Catalog* (TSV and CSV; a reloaded catalog keeps its decisions), exported by *Export FITS*,
and follow the catalog when the frame is switched.  The `Review` chip (tab *Results*) has:

* **Accept / Reject / Mark Uncertain / Clear Review of Selected** - apply to the selection (click on the image, in the table,
  or Ctrl-click for several).  **Set Note on Selected...** edits the note.
* **Review All Shown Rows** (accept / reject / uncertain) - applies to every row the table shows now (after the text search and
  the review filter), after a confirmation.
* **Show by Review** - table filter: All, Accepted, Rejected, Uncertain, Unreviewed, Not Rejected.  The filter is ANDed with the
  search box and stays in force across sort, reload and further decisions.  It is dropped automatically when a catalog without a
  `REVIEW` column is loaded.  Rows are tinted (green / red / amber) in the `REVIEW` and `NUMBER` cells.
* The same operations exist as Tcl: `::ogf::review::set_status accept|reject|uncertain|clear ?NUMBERS? ?note?`,
  `::ogf::review::set_note TEXT ?NUMBERS?`, `::ogf::review::show all|accept|reject|uncertain|none|notrej`,
  `::ogf::review::counts`.

Review applies to the **galaxy** catalog.  With the time-domain kind filter on Moving / Transients / All the commands refuse with
a message ("switch the kind filter to Galaxies first"): those rows are generated views, not stored cells.  This is a limit, not a
design choice.

### Session recorder

Every decision is recorded as an **informational (`note`) step**, so the session log and the exported script list it, and the
script skips it (like the sort/save steps of the time-domain table):

| step | payload |
|---|---|
| `review.set` | `numbers`, `status` (`accept` / `reject` / `uncertain` / `clear`), `note` |
| `review.note` | `numbers`, `note` |
| `report.export` | `file`, `rows`, `thumbnails` |

No existing step signature changes (`verify_session_replay.py`: 78 checks, 0 failed).  The catalog fingerprint the recorder stores after
each step (`OGFSessCatInfo`) ignores the three review columns, because the exported script never reproduces hand annotations;
without that, replay mode would report "catalog differs from the GUI session" after any review.  If the catalog has none of the
columns the text is untouched.

### Core hooks added for this feature

Small, general additions to `ogf_core.tcl` (all in the `::ogf::cat` namespace; documented in `docs/architecture.md` section 7):
`filter_set / filter_clear / filter_active / filters / filter_text` (column-value filters applied by `CatalogPanelFilter`),
`shown_numbers` (NUMBERs of the rows shown), `set_cells {NUMBER {COLUMN VALUE}}` (write cells, appending missing columns, then reload
the table with the selection kept), `on_table_filled CMD` (hook after the table is filled; the tint uses it).

## 2. Export Report...

*Review > Export Report...* (or `OGFReportExport ?file.html?`) runs `plugins/report/report.py` as a job (the Stop button cancels it)
and writes one HTML file.  Dialog: *Report Options...* (declarative, from `plugin.json`).

Contents of the report:

1. **Summary** - catalog rows, rows in the report, counts accepted / rejected / uncertain / unreviewed, which rows were selected
   (the rows the table shows, or by review state), table filter at export time.
2. **Candidates** - one row per candidate: PNG **cut-out** (`cutout_size` px square around `X_IMAGE`,`Y_IMAGE`, one global stretch
   zscale / asinh / linear for all cut-outs, optional yellow circle), the key measurements (NUMBER, X/Y, RA/Dec, MAG_AUTO and error,
   FLUX_RADIUS, FWHM, ellipticity, CLASS_STAR, plus analysis columns such as photo-z or morphology if present), the review state with
   time and note, and a collapsible line with every other column.  Links at the top filter the list by review state (a little
   JavaScript; the file works without it).
3. **Review history** - the `review.*` steps of the session log.
4. **Processing provenance** - input files with size, SHA-256 and modification time (image and catalog); software (OGFinder git HEAD and
   whether sources were modified, ds9 version, the Python and library versions - numpy, astropy, Pillow, sep, scipy - the report was made
   with); the **session log summary** (every recorded step: number, class, step, title, seconds, command line); the stored
   parameters of every plugin; the extraction parameters of the panel; the report options.

The HTML has no external resources (cut-outs are `data:` URIs), so it can be mailed or archived.  All catalog and provenance text is
HTML-escaped.  A cut-out costs about 3 kB (64 px, zoom 2); `max_rows` (default 500) limits the file size.

**PDF** (`pdf` option or `--pdf`): converted by a headless Chrome/Chromium/Edge when one is on `PATH`; otherwise only the HTML is written and a warning is
logged (print the HTML to PDF in a browser).  No other converter is attempted.

CLI (the same code, e.g. for a batch run on a saved catalog):

    python3 plugins/report/report.py --catalog cat.tsv --image image.fits --output report.html \
        [--provenance prov.json] [--include all|accepted|accepted_uncertain|not_rejected|uncertain|rejected|unreviewed] \
        [--sort-by -MAG_AUTO] [--max-rows 500] [--cutout-size 64] [--thumb-zoom 2] [--stretch zscale|asinh|linear] [--no-thumbs] [--pdf]

stdout: `REPORT<TAB>html<TAB>listed<TAB>selected<TAB>thumbnails<TAB>pdf-or-"-"`.

## Tests

* `plugins/report/tests` (pytest, `scripts/run_all_checks.sh --only report_tests`): the report is generated from a real extraction of
  `m51.fits` and of a 600x600 crop of `hudf_f160w.fits` and the HTML is parsed from Python: counts, review classes, escaped notes, one
  valid PNG per candidate with the expected size, centring of the cut-out on the real image, row selection modes, sorting, column
  errors, SHA-256 of the inputs, provenance steps, no external URL, XSS escaping, the CLI stdout contract, the optional PDF.
* `scripts/verify_review_gui.tcl` (`--only review_gui`): in a real ds9 on m51 - columns appended in order, values, time format, selection
  kept, counts, session steps and payloads, the five table filters, filter + search, filter survives sort and reload, tint cells, review of
  all shown rows, clear, save/load TSV and CSV round trip, filter dropped on a catalog without `REVIEW`, the export job (real `report.py`
  subprocess), provenance JSON, session export lists the notes, layout 181/769/154 at 1300x950.

## Limits

* Galaxy catalog only (see above).  Decisions are keyed by `NUMBER`: after a new extraction the numbers change and the old decisions do not
  apply (load the saved catalog instead).
* A review is not replayed by the exported script: it is a hand annotation.  Re-running the pipeline on new data produces a catalog
  without decisions.
* Cut-outs are of the *displayed image file* (first 2-D HDU), not of the current ds9 scale or other bands; one stretch for all.
* The report shows the table as it was at export; it is a snapshot, not a live view.
* Real mouse clicks on the new menu entries were not tested (the entries call the same procs the test calls); the tint colours were
  verified by cell tags, not by looking at the screen.
