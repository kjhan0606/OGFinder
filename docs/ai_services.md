# AI Services: connection points for external AI / astronomy services

OGFinder does **not** ship or call any foundation model. *Analysis > AI Services* is a small adapter layer: it turns
catalogue rows (and optional image cutouts) into requests for a service **you** configure, maps the answer onto named
catalogue columns, and records everything so the session script can replay it. Nothing in this repository claims that a
particular commercial or research model or endpoint is supported; all example profiles that point at a "service" use the
reserved host `example.invalid` and are **templates**.

* GUI: `plugins/ai_services/ai.tcl`; CLI: `ds9/library/ds9_ai_bridge.py`; library: `ai_bridge/` (Python 3 standard library;
  `numpy` + `astropy` only for cutouts); example profiles: `ai_services.example.json`; tests: `ai_bridge/tests/`,
  `scripts/verify_ai_gui.py`.

## Menu

*Analysis > AI Services* (three entries, nothing else):

| label | what it does |
|---|---|
| `Service Registry...` | lists the profiles (task, transport, enabled, auth env-var name and **set/unset**, target); **Enable / Disable**, **Test Connection**, **Detect Agent CLIs**, profile-file path (browse), **Create from Example**; "backend" selectors for the local Photo-z / SED / morphology / star steps |
| `Run Task on Catalog...` | choose task, **Backend** (Codex CLI / Claude Code / agy-Gemini CLI / Grok / generic services - see *Agent CLIs* below) and service, rows (all / selected in the table), cutout size (pix or arcsec) / normalisation / format, **Dry run** |
| `Show Last Run Log` | command line, exit status, bridge messages, provenance path of the last run |

Results become new catalogue columns (same machinery as the other analysis steps, `CatalogPanelAddColumnsFromTSV`), so
sort / filter / save work. Editing a profile is done in a text editor (JSON); the registry shows the path.

## Security notes (read this)

1. **Keys only through environment variables.** A profile stores only the *name* of the variable (`"auth": {"scheme":
   "bearer_env", "env": "MY_TOKEN"}`). The value is read when a request is sent, is never written to a profile, cache,
   log, provenance file or session script, is replaced by `***` in error text and appears as `<value of $MY_TOKEN: set>`
   in dry-run output. The GUI shows only *set* / *unset*; it never asks for a key. Profiles containing something that
   looks like a literal secret (`Bearer <long token>`, unknown `auth.*` keys) are rejected by `validate-profile`.
   Start ds9 from a shell in which the variable is exported (GUI-launched apps do not see variables set only in another terminal).
2. **Payloads leave your machine.** Positions, magnitudes, the catalogue-row text and image cutouts are sent to whatever
   `base_url` says. Check the service's terms and whether your data may be sent. The GUI asks for confirmation **once per
   service and session** before the first transfer (and before *Test Connection*); the dry-run shows exactly what would be sent.
3. Network services need `--allow-network` on the command line (the GUI adds it after you confirm). The exported session
   script treats `ai.run` like cross-match: in **pipeline** mode it is skipped with a message unless the script is run with
   `--allow-network`; in **replay** mode it is run.
4. `http://` is refused except for `localhost`/`127.0.0.1`; TLS certificates are verified; redirects to another host are
   refused (credentials are never forwarded). Responses are limited to 64 MB. `local_command` runs an argv list without a
   shell, with a minimal environment (`PATH HOME LANG TMPDIR VIRTUAL_ENV PYTHONPATH` + `env_passthrough`).
5. Responses are cached under `~/.ds9/ai_cache/<service>/` keyed by a hash of the *profile and the request* (never of
   credentials). Cached files contain the service responses, so treat that directory like the data. `--no-cache` bypasses it.
6. Results from a service are only as good as that service. Columns carry provenance (`AI_<P>_SERVICE/_MODEL/_REQID/_TIME/_ERROR`);
   the model string is whatever the **service returns** (OGFinder never invents one).

## Agent CLIs: Codex CLI, Claude Code, agy / Gemini CLI, Grok

Besides REST / local command / python callable, the external AI connection can be one of four **agent command-line
tools that are installed on your machine**. This is the transport `agent_cli` (`ai_bridge/agent_cli.py`); four ready-made
profiles are built in (no profile file needed) and appear in the registry and in the *Run Task* **Backend** dropdown:

| ready-made service | backend | executable looked up on `PATH` | command line built by the bridge (argv list, never a shell) |
|---|---|---|---|
| `agent_codex` | Codex CLI | `codex` | `codex exec --skip-git-repo-check --ephemeral --ignore-user-config --ignore-rules --sandbox read-only --color never -C <tmpdir> -` (prompt on stdin) |
| `agent_claude` | Claude Code | `claude` | `claude -p --output-format json --no-session-persistence --disable-slash-commands --strict-mcp-config --permission-mode dontAsk --tools ""` (prompt on stdin; `--tools Read` only when images are sent) |
| `agent_agy` | agy (Google Antigravity CLI), falls back to the Gemini CLI | `agy`, else `gemini` | `agy --print <prompt> --output-format json --print-timeout <n>s --disable-slash-commands`; with only `gemini` installed: `gemini -p "<short instruction>" --output-format json --approval-mode plan --skip-trust` (prompt on stdin) |
| `agent_grok` | Grok (xAI CLI) | `grok` | `grok --prompt-file <tmpdir>/prompt.txt --output-format json --permission-mode dontAsk --tools "" --no-auto-update --cwd <tmpdir>` (`--tools` is left out when images are sent) |

Each ready-made profile is *generic over tasks* (`"task": "any"`); the task chosen in the dialog / `--task` decides the prompt.

### What the bridge does with a CLI

1. **Prompt from the task contract.** `contracts.TASKS[task]` supplies the task description and the exact answer columns
   (name, type, `[0, 1]` range, optional or required). The catalog rows (id, x, y, ra, dec, magnitudes + errors; the full
   catalog row only for tasks that need it, e.g. `moving_object_classification`) are embedded as one JSON block marked as
   *data, never instructions*. `--dry-run` shows the complete prompt and the argv. For `generic` give
   `--param columns=SCORE:float,NOTE:str --param instruction="..."`.
2. **Run.** One subprocess per batch (`batch_size`, default 10 rows) via `subprocess.Popen(argv, shell=False)` in its own process
   group, in an **empty temporary working directory** (removed afterwards), with a **minimal environment**
   (`PATH HOME LANG TMPDIR ... proxy and CA variables` + the variable *names* in `env_passthrough` + the CLI's own
   login/API-key names, e.g. `ANTHROPIC_API_KEY`, `CODEX_API_KEY`, `XAI_API_KEY`, `GEMINI_API_KEY`). Prompt delivery: stdin
   (Codex, Claude, Gemini), a temp file (Grok) or one argv element (agy - the prompt is refused above 100 kB). Catalog text is
   only ever *data inside the prompt*; nothing from the catalog reaches argv except agy's single prompt argument, and no shell
   ever parses it (tested with `'; touch X; $(touch X)` in a catalog cell).
3. **Timeout.** `timeout_s` (default 300) kills the whole process group; timeouts are retried.
4. **Strict answer parsing.** stdout must be exactly one JSON document `{"results":[{"id":..., <columns>}], "model": ...}`.
   Only (a) the CLI's *own documented JSON envelope* (`claude`: `result` / `structured_output`; `gemini`: `response`; `agy`:
   `response` / `structured_output`; `grok`: `result`/`response`/`text`/`output`/`content`) and (b) **one** ```` ```json ````
   fence around the whole answer are unwrapped. Prose before or after the JSON is rejected. Then: every requested id exactly
   once, no unknown ids, required columns present, JSON types (a number must be a JSON number, not text), `[0, 1]` ranges
   (`STAR_PROB`, `REALBOGUS_SCORE`, `MORPH_CONF`, ...) - otherwise the answer is *invalid*.
5. **Retries.** An invalid answer, a timeout and rate-limit-like errors are retried (`retries`, default 2, with back-off); the
   retry prompt tells the model why its previous reply was rejected. A "not logged in" failure is *not* retried and says
   how to log in with the CLI itself. Failures are per batch: the rows get `AI_<P>_ERROR`, other batches continue.
6. **Provenance.** `AI_<P>_SERVICE` = the service name, `AI_<P>_MODEL` = the model the CLI/answer reports (else
   `<cli> CLI (model not reported)` - never invented), `AI_<P>_REQID` = the CLI's session/conversation id when it prints one.
   The provenance JSON gets an `agent_cli` block: backend, executable, `--version` string, whether images were sent, the
   names of the environment variables handed over, `extra_args`.

### Detection, executable path, extra args

* Detection is `shutil.which` (no model is contacted; the version comes from `<cli> --version`):
  `ds9_ai_bridge.py --mode detect-agents [--json]`, registry button **Detect Agent CLIs**, and the Backend dropdown shows
  *installed* / *NOT FOUND on PATH* for every backend.
* A CLI that is not on `PATH` is configured with `"executable": "/full/path/to/claude"` and optionally `"extra_args":
  ["--model", "sonnet"]` (appended to the argv). In the GUI: Run Task > Backend > *Executable* / *Extra args (JSON list)* >
  **Save to profile** (writes only these two keys into your profile file, via `--mode set-agent`; the file is copied to
  `.bak` first and a ready-made profile is copied into the file on first save). `"args"` replaces the built-in argv completely
  (and `"backend": "custom"` + `"executable"` + `"args"` + `"prompt_via": "stdin|file|argv"` define any other CLI; `{prompt_file}`
  `{workdir}` `{prompt}` `{timeout_s}` are replaced; see the template `template_custom_agent_cli` in `ai_services.example.json`).
* `check-service` / **Test Connection** on an agent service only checks the executable and prints its version:
  **no prompt, no data, no cloud call.**

### Security model of the agent backends

* **No secrets in profiles.** `validate-profile` refuses: `auth` blocks, anything that looks like a credential in `extra_args`
  (`--api-key=sk-...`, long random tokens), and arguments that switch off the CLI's own safety (`--dangerously-*`,
  `--yolo`, `--always-approve`, `bypassPermissions`, `danger-full-access`, `--approve-for-me`). Authentication is the CLI's own login
  (`codex login`, `claude` + `/login`, `agy`, `gemini`, `grok login`) or an API-key environment variable that you export
  yourself; the profile may list only the variable **names** (`env_passthrough`). Values are never printed, cached or recorded
  and are replaced by `***` in error text.
* **Tools are off.** The ready-made argv gives the model no shell / file / web tool where the CLI allows it (`--tools ""` for Claude
  and Grok, `--sandbox read-only` + empty working directory for Codex, `--approval-mode plan` for Gemini; agy has no switch for
  this in print mode: it runs with its default permission policy, where shell commands are soft-denied, and its working directory is the empty temp folder).
  These are *mitigations against prompt injection from catalog text*, not a guarantee.
* **Images only if you tick it.** Without `--send-images` (GUI: *Send image cutouts*, default off) no cutout is even made
  and the CLI's working directory stays empty; with it, the cutouts are copied into the temp directory, listed in the prompt,
  and (Claude, Grok) the file-reading tool is enabled / (Codex) attached with `--image`. The exported script records the
  choice; pipeline mode needs **`--allow-agent-images`** in addition to `--allow-agent-cli` to repeat such a step.
* **Confirmation before the first run (GUI).** The dialog is built from the bridge's own dry run (`--mode dry-run
  --summary-only`, so it cannot disagree with what is executed) and states: the CLI and its path, the provider whose cloud
  model receives the prompt (OpenAI / Anthropic / Google / xAI), the number of objects and **which fields per object**, whether the
  full catalog row text is included, **whether image data is included (and how many files)**, the prompt size, the environment
  variable names passed, and that OGFinder stores no credentials. It is asked once per service and session, and again when
  *Send image cutouts* is ticked for the first time. *Cancel* sends nothing; *Dry run* shows the exact prompt first.
* **Command-line consent.** `run` without `--allow-agent-cli` fails before anything is executed (`--allow-network` does not
  count). The GUI adds the flag after you confirm; it is stored in the recorded step.
* **The CLI is a third-party program** with its own telemetry, auto-update and data policy; read the provider's terms. The prompt (with
  your catalog rows) is processed by that provider under *your* account. Ready-made argv uses `--no-auto-update` (Grok) and
  `--ephemeral` / `--no-session-persistence` (Codex, Claude) so that the CLI does not keep the session on disk; the others may keep a history.

### Session recorder and the exported Python pipeline

An agent run is the existing step `ai.run` (same argv/step signature as before) with the extra payload keys `backend` and
`send_images`, `network=1`. Whole-catalogue runs are **AUTO**, row-subset runs are **MANUAL** (as for every `ai.run`).
In the exported script: **replay** mode runs it (you confirmed it in the GUI); **pipeline** mode skips it with
`use --allow-agent-cli` (and with `use --allow-agent-images as well` for a step recorded with *Send image cutouts*), exactly like
network services need `--allow-network`. `--allow-network` does not enable agent steps.

### Verification status of the flags (read this before relying on a CLI)

| CLI | version checked | what was verified against the **real binary** | what is only from the vendor documentation |
|---|---|---|---|
| Codex CLI | `codex-cli 0.160.0` | `codex exec --help` lists `--skip-git-repo-check --ephemeral --ignore-user-config --ignore-rules --sandbox read-only --color --image -C -o --output-schema --json`, prompt from stdin with `-`; the real binary accepted our whole argv and failed with `401 Unauthorized` (no login) which the bridge reports with the login hint | the plain-text final message on stdout (we never saw a model answer); `--image=<file>` form (only the flag is verified, not that it is accepted attached to a text prompt); `--output-schema` is *not* used |
| Claude Code | `2.1.287` | `claude --help` lists `-p --output-format json --no-session-persistence --disable-slash-commands --strict-mcp-config --permission-mode dontAsk --tools --json-schema --max-turns --bare`; real run printed the error envelope `{"type":"result","is_error":true,"result":"Not logged in..."}` which is parsed | the success envelope (`result`, `structured_output`, `modelUsage`, `session_id`) - from the docs; `--tools Read` for images |
| agy | `1.2.14` | `agy --help` lists `--print/-p --output-format --print-timeout --disable-slash-commands --json-schema --model --effort --sandbox`; the real binary accepted our argv and answered `authentication required` with the documented error envelope `{"status":"ERROR","error":...}` | the success envelope (`status` `SUCCESS`, `response`); **upstream issue [antigravity-cli#408](https://github.com/google-antigravity/antigravity-cli/issues/408): `agy --print` writes nothing to stdout when stdout is a pipe** - the bridge reports "empty stdout" if that happens to your version and you should then use the `gemini` executable (set `"executable": "gemini"` in the profile) |
| Gemini CLI (as the agy fallback) | `0.62.0` | `gemini --help` lists `-p --output-format -o --approval-mode plan --skip-trust`; real run failed with exit 41 and the JSON error object on **stderr** (parsed) | the success envelope (`response`, `stats.models`) |
| Grok (xAI) | `grok 1.0.46` | `grok --help` lists `-p/--single --prompt-file --output-format {plain,json,...} --permission-mode dontAsk --tools --no-auto-update --cwd --json-schema --max-turns`; real run printed `{"type":"error","message":"Not signed in..."}` which is parsed | **the success JSON envelope was never seen** - the parser accepts `{"results": ...}` directly or the first of `result/response/text/output/content`, otherwise fails with a message naming the keys; treat Grok as the least verified |

**Grok authentication (round 2, item 6; grok 1.0.46, checked on this box without reading any credential).** `grok login --help` offers `--oauth` (browser, auth.x.ai) and
`--device-auth` / `--device-code` (headless: the CLI prints a URL and a one-time code, a human approves it in a browser on any machine; **needs a person, cannot be done
non-interactively**). `~/.grok` holds only `config.toml` and `bin/` here, i.e. no stored login; `grok models` answers "You are not authenticated". The non-interactive route
is the environment variable **`XAI_API_KEY`** (a key from console.x.ai; the binary also knows `GROK_CODE_XAI_API_KEY`, `GROK_DISABLE_API_KEY_AUTH`, `GROK_HOME`,
`GROK_CONFIG_PATH`, `GROK_AUTH_PROVIDER_COMMAND|ACCESS_TOKEN|REFRESH_TOKEN|EXPIRES_AT` for an external token provider): with a dummy value `grok models` prints "You are
using XAI_API_KEY", but `grok -p` / `--prompt-file` rejects an invalid key with the same `Not signed in` error envelope, so the key is validated at the first request.
No key exists on the box and the `SAND_*` variables are the sandbox's own (unrelated to Grok), so no model answer was obtained and the success envelope is still unverified.
`ai_bridge/tests/test_grok_live.py` runs the real round trip on a 2-object star/galaxy task and is SKIPPED unless `grok` is installed and authenticated (it also
skips on an authentication error at run time); the first successful run will show which envelope keys the parser actually needs.
Steps for the user: `grok login --device-code` (approve the code in a browser) or `export XAI_API_KEY=...`, then `python -m pytest ai_bridge/tests/test_grok_live.py -v -rs`.

**No real agent CLI has produced a model answer in any test**: the box has no login for any of them and no API key was
used or exists, so nothing was sent to any provider. All five binaries were installed only to read their `--help` and to confirm
that our argv is accepted and that their *not-logged-in* failure is understood (`scripts/verify_agent_cli_real.py`).
Everything about a *successful* answer (envelope shape, JSON-only compliance by the model, timing, cost) is verified only
against the fake executables of `ai_bridge/tests/fake_agents.py`, which mimic the documented formats.
Vendor flags change quickly (these CLIs were at versions 0.160 / 2.1 / 1.2 / 0.62 / 1.0 in October 2026); if a flag is rejected,
use `"args"` in the profile to override the argv. Vendor documentation used:
[Codex non-interactive](https://developers.openai.com/codex/noninteractive),
[Claude Code headless](https://code.claude.com/docs/en/headless),
[Antigravity CLI headless](https://antigravity.google/docs/cli/headless/),
[Gemini CLI headless](https://google-gemini.github.io/gemini-cli/docs/cli/headless.html),
[Grok CLI headless](https://docs.x.ai/build/cli/headless-scripting).

Model answers are *not deterministic*: a re-run (replay) may give different numbers; the response cache (key = profile +
exact prompt) makes an identical repeated request return the first answer. Treat the results as the opinion of a language
model, not a measurement.


## Offline testing with the mock service

`mock` is built in (no profile needed): deterministic **fake** values that depend only on (object id, column), so replays
are bit-identical. Every output row has `AI_<P>_SERVICE = mock` and `AI_<P>_MODEL = mock-model (FAKE VALUES)`.

    python ds9/library/ds9_ai_bridge.py --mode run --service mock --task photoz --catalog cat.tsv --output out.tsv

## Command line

    python ds9_ai_bridge.py --mode {list-services,check-service,run,dry-run,validate-profile,set-enabled,detect-agents,set-agent} --service NAME \
        --task TASK --catalog CAT.tsv [--image A.fits[,B.fits] --bands F105W,F160W] --output OUT.tsv [options]

| option | meaning |
|---|---|
| `--services-file F` | profile file (default `$OGF_AI_SERVICES_FILE`, else `~/.ds9/ai_services.json`) |
| `--numbers 3,17,@file` | only these NUMBER values (output stays in catalogue order) |
| `--size-pix N` / `--size-arcsec X` | cutout size; `--normalize none|linear|zscale|asinh`; `--cutout-format png|fits|npy`; `--cutout-dir` |
| `--pixel-scale`, `--psf-fwhm`, `--mag-columns`, `--id-col` | inputs for the standard record |
| `--param KEY=VALUE` | service parameter, available as `{param_KEY}` and in `params` |
| `--rename OLD=NEW` | rename output columns (the GUI uses it to give local-step column names, e.g. `PHOTOZ=PHOTO_Z`) |
| `--no-cache`, `--cache-dir` | response cache (`~/.ds9/ai_cache`, or `$OGF_AI_CACHE_DIR`) |
| `--resume` | keep already finished rows (no error, same service) of an existing `--output` and request only the rest |
| `--allow-network` | required for `https_*` / `tap_query` services in `run` mode |
| `--allow-agent-cli` | required for `agent_cli` services in `run` mode (see *Agent CLIs*); `--allow-network` does **not** substitute |
| `--send-images` | `agent_cli` only: also hand the image cutouts to the CLI (default: no image data) |
| `--strict` | exit 3 if any object failed |
| `--provenance-output P` | provenance JSON (default `<output>.provenance.json`) |

Exit status: `0` success (also when *some* objects failed - their row has `AI_<P>_ERROR` set; use `--strict` to change),
`1` usage / configuration error (invalid profile, missing `--allow-network`, ...), `2` **every** object failed, `3` partial
failure with `--strict`. Output rows are always in catalogue order (deterministic); the output TSV has `NUMBER` plus the new
columns and is merged into the table by NUMBER.

`dry-run` prints the exact request(s) (URL, headers with credentials redacted, body; base64 blobs are shown as length + hash)
and sends nothing; `validate-profile` checks a profile offline; `check-service` tests the connection (a GET of
`base_url` + `health.path`, or the VOSI `/availability` of a TAP service; local commands: executable exists).

**Provenance sidecar** (`OUT.tsv.provenance.json`): service, transport, task, profile SHA-256 and a redacted profile copy, auth env
var *name* and set/unset, service-reported model strings, catalogue / image / output SHA-256, command-line arguments, object
counts (ok / failed), request counts (sent / cache hits / retries / failed attempts / bytes), cache statistics, per-object
failures, ai_bridge version, python version, timings.

## Profile format

A profile file is JSON: `{"schema": "ogf-ai-bridge/1", "services": [ {...}, ... ]}`. JSON has no comments, so every key that
starts with `_` is a note and is ignored (and not part of the profile hash). The shipped `ai_services.example.json` contains
five **template** profiles (disabled, placeholder URLs) and one real, verified non-ML service (`simbad_cone`).

Annotated example (`//` comments are explanatory only; they are not valid JSON):

    {
      "name": "template_photoz_rest",          // unique; letters, digits, _ . -
      "template": true,                        // informational: marks a TEMPLATE
      "enabled": false,                        // disabled profiles cannot run (Registry / set-enabled toggles it)
      "task": "photoz",                        // photoz sed_fit morphology star_galaxy real_bogus transient_classification
                                               // moving_object_classification anomaly_detection embedding_similarity
                                               // captioning image_retrieval generic        ("any" only for mock)
      "transport": "https_json",               // https_json | https_multipart | local_command | python_callable | tap_query | mock
      "base_url": "https://example.invalid/api",   // https:// (http:// only for localhost)
      "endpoint": "/v1/photoz",                // appended to base_url
      "method": "POST",                        // GET | POST | PUT
      "headers": {"Accept": "application/json"},   // static headers (no secrets here)
      "auth": {"scheme": "bearer_env",         // none | bearer_env | header_env (+"header") | query_env (+"param")
               "env": "MY_PHOTOZ_TOKEN"},      // NAME of the environment variable; "prefix": "Bearer " is the default
      "request": {
        "template": {                          // JSON request body; placeholders in {braces}
          "object_id": "{id}",                 //   "{x}" alone -> typed value (number/object); "a {x} b" -> text
          "coordinates": {"ra": "{ra}", "dec": "{dec}"},
          "photometry": "{mags}"               //   {{ and }} are literal braces; unknown placeholder = error, nothing sent
        },
        "query": {"band": "{param_band}"}      // optional query string (also templated)
      },
      "response": {
        "model_path": "model.version",         // where the SERVICE reports its model name/version (-> AI_*_MODEL)
        "request_id_path": "request_id",       // or "request_id_header": "X-Request-Id"  (-> AI_*_REQID)
        "error_path": "error",                 // if present and non-empty the object is marked failed
        "items_path": "results",               // batch responses: list of per-object items ...
        "id_path": "id",                       // ... matched back to catalogue rows by this id (order may differ)
        "fields": {                            // COLUMN -> where to find it, with type conversion and units
          "PHOTOZ":     {"path": "result.z", "type": "float"},
          "PHOTOZ_P84": {"path": "result.pdf_pct[2]", "type": "float",
                         "scale": 1.0, "offset": 0.0,                   // value*scale+offset
                         "missing_values": [-99],                       // -> empty cell
                         "unit": "arcsec", "to_unit": "arcmin"}         // declared-unit conversion (angle/flux/length/time)
        }
      },
      "timeout_s": 30, "retries": 3, "backoff_s": 1.0, "backoff_factor": 2.0, "backoff_max_s": 30,  // retry 408 425 429 5xx, timeouts, connection errors
      "rate_limit_per_s": 5,                   // minimum spacing between requests (0 = unlimited)
      "batch_size": 1,                         // objects per request; >1 needs "{items}" in request.template (+ "item_template")
      "max_payload_bytes": 2000000,            // larger requests are not sent (objects marked failed)
      "params": {"band": "r"},                 // default {param_*} values (override with --param)
      "cutouts": {"size_arcsec": 12, "normalize": "asinh", "format": "png"},   // see below
      "health": {"path": "/health"}            // used by Test Connection
    }

**Response mapping.** Paths are JSONPath-like: `a.b`, `a.items[0].c`, `a.items[-1]`, `a.items[*].c` (wildcard -> list), optional
leading `$.`; `$root.x` in a batched response reads from the whole response. Types: `float int str bool json`. Strings like
`nan`/`null` and `missing_values` become empty cells; NaN/inf are empty. A probability column (`STAR_PROB`, `REALBOGUS_SCORE`,
`MORPH_CONF`, `CLASS_PROB`, `MOVING_SCORE`, `ANOMALY_SCORE`) outside [0, 1] makes that object fail (a wrong mapping is
reported, not stored). If none of the mapped paths exist in a response, the object fails; the failure is shown in
`AI_<P>_ERROR` and not cached.

**Placeholders** (request templates): `{id} {x} {y} {ra} {dec} {mags} {mag_errs} {catalog_row} {mags_json} {mag_errs_json}
{catalog_row_json} {pixel_scale_arcsec} {psf_fwhm_arcsec} {task} {service} {bands} {radius_deg} {param_NAME}` and, for
cutouts, `{cutout_png_base64}`, `{cutout_fits_base64}`, `{cutout_npy_base64}`, `{cutout_<fmt>_path}`, each with an optional
band suffix (`{cutout_png_base64_F160W}`; without suffix the first image), `{cutouts_json}` (all bands as base64 JSON) and
`{cutout_fits_url}` (built from `cutouts.url_template`, for services that fetch the image themselves; nothing is uploaded).

**Transports**

| transport | what it does |
|---|---|
| `https_json` (`generic_rest`) | one request per object or per batch; JSON body (or query string for GET) |
| `https_multipart` | image upload: `request.multipart: {"fields": {...}, "files": {"image": "{cutout_png_path}"}}`, one object per request |
| `local_command` | runs `"command": [argv...]` (`{python}` = the interpreter running the bridge), JSON on stdin / stdout (below) |
| `python_callable` | `"callable": "module:function"` (+ optional `"sys_path": [dirs]`): `function(records, params, context) -> list of dicts` |
| `mock` | built-in FAKE values (above) |
| `tap_query` | IVOA TAP `POST <base_url>/sync` (`REQUEST=doQuery LANG=ADQL FORMAT=csv`), `request.adql` is an ADQL template (`{ra} {dec} {radius_deg} {id} {param_*}`); the CSV is exposed to the mapping as `rows[i].column`, `n_rows`. A cone search without a match yields empty cells (not an error). Non-ML example; use task `generic` |

**local_command / python_callable contract** (`ai_bridge/examples/echo_local_command.py` is a toy that shows the plumbing):

    stdin  : {"contract": "ogf-ai-bridge/1", "task": "star_galaxy", "service": "NAME", "params": {...},
              "records": [ <input record>, ... ]}
    stdout : {"results": [{"id": "7", "STAR_PROB": 0.93, "CLASS_LABEL": "STAR"}, {"id": "8", "error": "why"}],
              "model": "<name+version of your model>", "request_id": "<optional>"}

Column names in `results` are the task contract names (or define `response.fields` with paths into each result).

**Cutouts.** `cutouts` (profile) and `--size-pix/--size-arcsec/--normalize/--cutout-format` (CLI/GUI override) control
`ai_bridge/cutouts.py`: per object and band, centred on the catalogue position (X_IMAGE/Y_IMAGE refer to the first image; other
bands are located through their WCS), edge objects are NaN-padded, size in pixels or arcsec (arcsec needs a WCS or
`--pixel-scale`), normalisation `asinh` (zscale limits, a=0.1) / `zscale` / `linear` (1-99.5 %) / `none` (raw values),
format PNG (8-bit grey, lossy by construction, normalised) / FITS (raw pixels with a shifted WCS) / `.npy` (float32).
Images are opened only when a template or the task needs cutouts. Objects outside an image fail individually.

## Task contracts

Standard **input record** (what every adapter gets for each object; built from the catalogue TSV and the images):

| field | type | meaning |
|---|---|---|
| `id` | str | catalog object identifier (NUMBER column by default) |
| `x, y` | float | pixel position in the first image (1-based, SExtractor convention: X_IMAGE, Y_IMAGE) |
| `ra, dec` | float | degrees, ICRS/J2000 (ALPHA_J2000, DELTA_J2000) or null |
| `mags` | dict band->float | from MAG_<band> columns (MAG_AUTO -> band 'AUTO'); null when 99/blank |
| `mag_errs` | dict band->float | from MAGERR_<band> columns |
| `cutouts` | dict band->file | per-band cutout, FITS / PNG / npy, produced lazily (see cutouts.py) |
| `wcs` | dict | FITS WCS keywords of the cutout (CTYPE, CRVAL, CRPIX, CD/CDELT, ...) |
| `pixel_scale_arcsec` | float | arcsec/pixel (WCS, or --pixel-scale) |
| `psf_fwhm_arcsec` | float | from --psf-fwhm, else null |
| `catalog_row` | dict col->str | the full catalog row, as text |

Standard **output columns** per task (a profile maps a service response onto these names; extra columns are allowed
with a warning). Provenance columns (all text) are added to every run: `AI_<P>_SERVICE`, `AI_<P>_MODEL`, `AI_<P>_REQID`,
`AI_<P>_TIME`, `AI_<P>_ERROR` (empty when the object succeeded), `<P>` = the prefix in the table (`provenance_prefix` in
the profile overrides it).

| task | provenance prefix | output columns (type) |
|---|---|---|
| `photoz` | `AI_PHOTOZ_*` | `PHOTOZ` (float), `PHOTOZ_ERR` (float), `PHOTOZ_P16` (float), `PHOTOZ_P50` (float), `PHOTOZ_P84` (float) |
| `sed_fit` | `AI_SED_*` | `LOG_MASS` (float), `LOG_MASS_ERR` (float), `LOG_AGE` (float), `LOG_AGE_ERR` (float), `LOG_Z` (float), `AV` (float), `SFR` (float), `SED_CHI2` (float) |
| `morphology` | `AI_MORPH_*` | `MORPH_TYPE` (str), `MORPH_CONF` (float), `MORPH_DESC` (str) |
| `star_galaxy` | `AI_SG_*` | `STAR_PROB` (float), `CLASS_LABEL` (str) |
| `real_bogus` | `AI_RB_*` | `REALBOGUS_SCORE` (float), `CLASS_LABEL` (str) |
| `transient_classification` | `AI_TRANS_*` | `CLASS_LABEL` (str), `CLASS_PROB` (float), `CLASS_PROBS` (json) |
| `moving_object_classification` | `AI_MOV_*` | `CLASS_LABEL` (str), `CLASS_PROB` (float), `MOVING_SCORE` (float) |
| `anomaly_detection` | `AI_ANOM_*` | `ANOMALY_SCORE` (float), `ANOMALY_FLAG` (bool) |
| `embedding_similarity` | `AI_EMB_*` | `EMBEDDING_ID` (str), `SIMILAR_IDS` (str), `SIMILARITY_TOP` (float) |
| `captioning` | `AI_CAP_*` | `CAPTION` (str) |
| `image_retrieval` | `AI_RET_*` | `RETRIEVED_IDS` (str), `RETRIEVAL_SCORE` (float) |
| `generic` | `AI_GEN_*` | defined by the profile's response mapping |

`needs` of each task (used by the mock to decide whether to produce cutouts): photoz/sed_fit need magnitudes;
morphology/real_bogus/captioning/image_retrieval need cutouts; the others take magnitudes or cutouts.
The columns `PHOTOZ/PHOTOZ_ERR` are renamed to the local step's `PHOTO_Z/PHOTO_Z_ERR` by the GUI backend hook;
`LOG_MASS ... SFR` (SED) and `MORPH_TYPE/MORPH_CONF/MORPH_DESC` (morphology) already use the local names.

## Session recorder

Each run from the GUI is recorded as step **`ai.run`**: class **AUTO** when it ran on the whole catalogue, **MANUAL** when it
ran on a hand-selected row subset (`--numbers`, replayable only in `--mode replay`). The recorded argv is exactly the
`ds9_ai_bridge.py` command line (python, script, catalogue and image paths are templated like the other steps); the step
declares the new columns (`post.kind = "ai"`) and `network = 1` for network services. Dry runs are not recorded.
`ogf_session_template.py` therefore replays it as `python ds9_ai_bridge.py --mode run ...` and merges the columns exactly like the GUI;
in pipeline mode the step needs `--allow-network` when the service is a network service.
The service profile is read from the profile file on the machine that replays the script (the recorded profile hash is in
the step payload: *the script does not embed the profile and does not check the hash*), and the provenance file is volatile
in the manifest. Timestamps (`AI_*_TIME`) of non-mock services are the time of the request, so a replay that makes new requests
differs from the GUI catalogue in that column (and in anything the service returns differently). Running the replay with
`OGF_AI_CACHE_DIR` pointing at the GUI's cache directory re-uses the stored responses (same request, same profile), which gives the
same values and timestamps; the mock service has no such difference (fixed fake timestamp).

## Backend selection of existing local steps

Photo-z (AI), SED Fitting (AI), Galaxy morphology and AI Star Classification keep their local implementation by default.
In *Service Registry...* the "backend" selectors can switch each to `local` or an enabled service of the matching task.
Implementation: one line at the top of each entry proc in `layout.tcl`
(`if {[info commands OGFAIBackendHook] ne {} && [OGFAIBackendHook photoz]} return`); the hook is a no-op when the backend is `local`.
Limits of this hook (honest list): the external path runs the service on the **whole catalogue with default options of the
profile** (no per-step dialog: bands/mag columns are all `MAG_<band>` columns, not the local dialog's choice);
local-only side effects are not reproduced - photo-z `PHOTO_Z_Q68/PHOTO_Z_OUTLIER`, the star finder's `AI_STAR/AI_STAR_CONF` columns and
marker recolouring, the morphology colour markers - the external step writes the contract columns (`STAR_PROB`, `CLASS_LABEL`, ...)
instead. The choice is stored in `~/.ds9/ai_services.prf`.

## Add a new service in 5 steps

1. Read the service's own API documentation; note URL, method, authentication, request and response JSON.
2. `export MY_KEY=...` in the shell that starts ds9 / the script (never in a file you commit).
3. Registry > *Create from Example* (or copy `ai_services.example.json` to `~/.ds9/ai_services.json`) and edit a copy of the closest
   template: `base_url`, `endpoint`, `auth.env`, `request.template`, `response.fields`; set `"enabled": true`.
4. `python ds9_ai_bridge.py --mode validate-profile --service NAME`, then `--mode dry-run ... --max-objects 3` and read the request.
5. `--mode check-service --allow-network`, then a `run` on 3 objects (`--max-objects 3`), compare with what the service's
   documentation says, then run the whole catalogue (GUI: *Run Task on Catalog...*).

If the service has no HTTP API you can use (or you want a model of your own), wrap it in a `local_command` or a `python_callable`.

## Tests and what they cover

* `python -m pytest ai_bridge/tests -v` (needs numpy + astropy for the cutout tests; `OGF_AI_OFFLINE=1` skips the live TAP tests):
  templating, response mapping, units, auth-from-env (bearer/header/query, unset variable, no secret in provenance/cache),
  retry/backoff against a local `http.server` (503, `Retry-After`, 4xx, timeout, giving up), cache hit/miss and key sensitivity,
  partial failure + `--resume`, batching + id matching, payload limit, rate limit, redirect refusal, dry run, multipart upload,
  local command, python callable, cutouts (FITS/PNG/npy, arcsec size, edges, determinism), CLI modes and exit codes, profile validation.
* `ai_bridge/tests/test_tap_live.py`: **live** SIMBAD TAP and VizieR TAP (Gaia DR3 table `I/355/gaiadr3`) cone searches and a real TAP error.
  Verified on 2026-10-01. These services are public and may change; they are not an ML test.
* `ai_bridge/tests/test_agent_cli.py` (+ `fake_agents.py`): 41 tests with **fake `codex` / `claude` / `agy` / `gemini` / `grok`
  executables on `PATH`** that mimic each CLI's documented output: detection (found / not found / agy->gemini fallback / configured path),
  argv + prompt delivery per CLI, envelope unwrapping, strict parsing (prose rejected, one fence tolerated, bad range / type /
  missing id / extra id / garbage), retry with the rejection reason, timeout + process-group kill, auth failure not retried, rate limit
  retried, consent flag missing -> nothing executed, dry run executes nothing, shell-injection text in the catalog, minimal
  environment and secret scrubbing, empty temp dir, images off by default / on with `--send-images`, `set-agent`, cache.
* `scripts/verify_agent_gui.py` (+ `.tcl`): 34 checks under Xvfb with fake CLIs on `PATH` (agy and grok absent): Backend dropdown,
  installed / not found, confirmation text, declined / dry-run / not-installed paths, image tick, *Save to profile*, recorder,
  replay in a plain shell and the four pipeline-mode flag combinations.
* `scripts/verify_agent_cli_real.py`: runs whichever **real** CLIs are installed without a login (see the verification table).
* `scripts/verify_ai_gui.py`: GUI smoke test under Xvfb (registry dialog, dry run, mock runs with/without cutouts, selected rows,
  confirmation dialog, backend hook, recorder) + replay of the exported script in a plain shell + pipeline mode with/without `--allow-network`.

## Not verified / limits

* **No real foundation-model endpoint has been tested** - only the mock, local command / callable examples, a local HTTP test
  server and the public TAP services. Templates may need adjusting to a service's actual API (authentication styles beyond
  bearer / header / query parameter, asynchronous job APIs, streaming, OAuth flows are not implemented).
* TAP: only synchronous queries; no authentication flows beyond the three schemes; ADQL is sent as written in the profile.
* The registry dialog does not edit profiles (edit the JSON file); the profile is re-read on every action.
* Rows selected in the GUI are replayed only in `--mode replay`.
* **Agent CLIs: see the verification table above - no real model answer was ever obtained.** Not implemented: streaming output,
  `--json-schema` / `--output-schema` server-side schema enforcement (we validate ourselves), per-call model selection other than
  `extra_args`, sending FITS (only PNG/npy/FITS cutouts as files, the CLI decides what it can read), Windows-specific process
  handling (the process-group kill is POSIX; on Windows only the CLI process is killed).
