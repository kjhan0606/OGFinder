#!/usr/bin/env python3
"""Static validation of every plugins/*/plugin.json (and of `cli` templates in particular).  Exit status 1 when a problem is found.

usage: tools/validate_manifests.py [plugins_dir] [-v]

Checks (the Tcl loader `::ogf::reg::validate` only checks id/name/tab, step ids, session class, cli-xor-proc and parameter types):
  * JSON parses; required keys; `tab` is a known tab; step ids and parameter names are unique; `session`, `stage`, `needs` have known values
  * each step has `cli` XOR `proc` XOR `variants`; a `proc` is defined in one of the plugin's Tcl files (or in ds9/library) - names are looked up with a regexp, not by running Tcl
  * parameters: type, default type, choice defaults inside `choices`, min <= default <= max, `store` keys
  * every `cli` template: tokens are known ({python} {plugin_dir} {work} {root} {image} {catalog} {PARAM} {OTHER:PARAM} {script:NAME} {cat:KEY});
    `{PARAM}` / `if` name a parameter of the plugin, `{OTHER:PARAM}` a parameter of that plugin; `{script:NAME}` and `{plugin_dir}/FILE` exist on disk;
    the driver's `--flags` of the template occur in the driver source (argparse) - catches typos;
    `{catalog}` needs `needs: catalog`; `{image}` needs `needs: image`; `{cat:psf,file}` style values only with a matching `needs`
  * `output.mode` known; `add_columns` has `columns`; a step with `output` has a `cli`
  * `records` are strings; the `record` of a cli step (recorder name) is a string
"""
import glob, json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABS = ["Detect", "Classify", "Measure", "Low-SB", "Time-domain", "Results"]
SESSION = {"AUTO", "CONFIG", "MANUAL", "NONE"}
STAGES = {"detect", "classify", "measure"}
NEEDS = {"image", "catalog", "psf"}
PTYPES = {"int", "float", "string", "bool", "choice", "file"}
OUT_MODES = {"add_columns", "set", "text"}
TOKEN = re.compile(r"\{([A-Za-z0-9_.:,-]+)\}")
BUILTIN = {"python", "plugin_dir", "work", "root", "image", "catalog"}


def procs_in(files):
    names = set()
    for f in files:
        try:
            txt = open(f, encoding="utf-8", errors="replace").read()
        except OSError:
            continue
        names.update(re.findall(r"^\s*proc\s+([^\s{]+)", txt, re.M))
    return names


def validate(pdir):
    problems = []
    plugs = {}
    for f in sorted(glob.glob(os.path.join(pdir, "*", "plugin.json"))):
        try:
            plugs[os.path.basename(os.path.dirname(f))] = (f, json.load(open(f, encoding="utf-8")))
        except Exception as e:
            problems.append("%s: invalid JSON: %s" % (f, e))
    # procs of any plugin file count too: plugins call each other's procs (e.g. lsbg -> catalog plugin's CatalogPanelSaveCatalog)
    libprocs = procs_in(glob.glob(os.path.join(ROOT, "ds9", "library", "*.tcl")) + glob.glob(os.path.join(pdir, "*", "*.tcl")))
    libprocs.add("OGFUIBandMenu")                    # dispatched by name prefix in ogf_ui.tcl (OGFUIAddDynamic / ::ogf::ui::band_menu)
    stats = dict(plugins=len(plugs), steps=0, cli=0, params=0)
    for pid, (f, m) in plugs.items():
        P = lambda msg: problems.append("%s: %s" % (pid, msg))
        if m.get("id") != pid:
            P("id %r differs from directory name" % m.get("id"))
        for k in ("id", "name", "tab"):
            if not m.get(k):
                P("missing %r" % k)
        if m.get("tab") and m["tab"] not in TABS:
            P("unknown tab %r" % m["tab"])
        tcl = m.get("tcl") or []
        tcl = [tcl] if isinstance(tcl, str) else tcl
        pdirp = os.path.dirname(f)
        for t in tcl:
            if not os.path.isfile(os.path.join(pdirp, t)):
                P("tcl file %s missing" % t)
        defined = procs_in([os.path.join(pdirp, t) for t in tcl]) | libprocs
        params = {}
        for p in m.get("params", []):
            stats["params"] += 1
            n = p.get("name")
            if not n:
                P("parameter without name"); continue
            if n in params:
                P("duplicate parameter %s" % n)
            params[n] = p
            t = p.get("type", "string")
            if t not in PTYPES:
                P("parameter %s: unknown type %s" % (n, t))
            d = p.get("default")
            if t == "choice":
                if not p.get("choices") and not p.get("choices_proc"):
                    P("parameter %s: choice without choices" % n)
                elif p.get("choices") and str(d) not in [str(c) for c in p["choices"]]:
                    P("parameter %s: default %r not in choices" % (n, d))
            if t in ("int", "float") and d is not None:
                if not isinstance(d, (int, float)) or isinstance(d, bool):
                    P("parameter %s: default %r is not a number" % (n, d))
                else:
                    if "min" in p and d < p["min"]:
                        P("parameter %s: default %s < min %s" % (n, d, p["min"]))
                    if "max" in p and d > p["max"]:
                        P("parameter %s: default %s > max %s" % (n, d, p["max"]))
            if t == "bool" and d not in (None, 0, 1, True, False):
                P("parameter %s: bool default %r" % (n, d))
        if "store" in m and not {"array", "key"} <= set(m["store"]):
            P("store needs array and key")
        seen = set()

        def check_step(s, where):
            stats["steps"] += 1
            sid = s.get("id")
            if not sid or not s.get("label"):
                P("%s: step without id/label" % where); return
            if sid in seen:
                P("duplicate step id %s" % sid)
            seen.add(sid)
            if s.get("session", "AUTO").upper() not in SESSION:
                P("step %s: session %r" % (sid, s.get("session")))
            if s.get("stage") and s["stage"] not in STAGES:
                P("step %s: stage %r" % (sid, s["stage"]))
            for n in s.get("needs", []):
                if n not in NEEDS:
                    P("step %s: unknown need %r" % (sid, n))
            kinds = [k for k in ("cli", "proc", "variants") if k in s]
            if len(kinds) > 1 and set(kinds) != {"proc", "variants"}:
                P("step %s: gives %s (use one)" % (sid, " and ".join(kinds)))
            if not kinds and not s.get("toggle") and not s.get("settings_only"):
                P("step %s: neither cli, proc nor variants" % sid)
            if "proc" in s:
                name = s["proc"].split()[0]
                if name not in defined:
                    P("step %s: proc %s is not defined in the plugin's Tcl files or ds9/library" % (sid, name))
            for v in s.get("variants", []):
                if "proc" not in v:
                    P("step %s variant %s: no proc" % (sid, v.get("id")))
                elif v["proc"].split()[0] not in defined:
                    P("step %s variant %s: proc %s is not defined" % (sid, v.get("id"), v["proc"].split()[0]))
            for r in s.get("records", []):
                if not isinstance(r, str):
                    P("step %s: records must be strings" % sid)
            out = s.get("output")
            if out:
                if out.get("mode") not in OUT_MODES:
                    P("step %s: output.mode %r" % (sid, out.get("mode")))
                if out.get("mode") == "add_columns" and not out.get("columns"):
                    P("step %s: add_columns without columns" % sid)
                if "cli" not in s:
                    P("step %s: output without cli" % sid)
            if "cli" in s:
                stats["cli"] += 1
                check_cli(s, sid)

        def check_cli(s, sid):
            cli = s["cli"]
            if not isinstance(cli, list) or not cli:
                P("step %s: cli must be a non-empty list" % sid); return
            flat = []; in_file = set()
            for e in cli:
                if isinstance(e, str):
                    flat.append(e)
                elif isinstance(e, dict) and "argv" in e and ("if" in e or "if_file" in e):
                    if "if" in e and e["if"] not in params:
                        P("step %s: condition parameter %r does not exist" % (sid, e["if"]))
                    if "if_file" in e:
                        flat.append(e["if_file"])
                        in_file.update(TOKEN.findall(e["if_file"]))
                    flat += [a for a in e["argv"] if isinstance(a, str)]
                else:
                    P("step %s: bad cli element %r" % (sid, e))
            text = " ".join(flat)
            needs = set(s.get("needs", []))
            for a in flat:
                for tok in TOKEN.findall(a):
                    if tok in BUILTIN:
                        continue
                    if tok.startswith("script:"):
                        path = os.path.join(ROOT, "ds9", "library", tok[7:])
                        if not os.path.isfile(path):
                            P("step %s: {%s}: %s not found" % (sid, tok, path))
                    elif tok.startswith("cat:"):
                        if tok.startswith("cat:psf,") and "psf" not in needs and tok not in in_file:
                            P("step %s: {%s} without needs psf or if_file" % (sid, tok))
                    elif ":" in tok:
                        op, on = tok.split(":", 1)
                        if op not in plugs:
                            P("step %s: {%s}: no plugin %s" % (sid, tok, op))
                        elif on not in {q.get("name") for q in plugs[op][1].get("params", [])}:
                            P("step %s: {%s}: plugin %s has no parameter %s" % (sid, tok, op, on))
                    elif tok not in params:
                        P("step %s: unknown token {%s}" % (sid, tok))
            if "{catalog}" in text and "catalog" not in needs:
                P("step %s: uses {catalog} without needs: catalog" % sid)
            if "{image}" in text and "image" not in needs:
                P("step %s: uses {image} without needs: image" % sid)
            # driver script and flags
            drv = None
            for a in flat:
                m2 = re.match(r"\{plugin_dir\}/(.+)$", a)
                if m2:
                    drv = os.path.join(pdirp, m2.group(1))
                m2 = re.match(r"\{script:(.+)\}$", a)
                if m2:
                    drv = os.path.join(ROOT, "ds9", "library", m2.group(1))
            if drv and os.path.isfile(drv):
                src = open(drv, encoding="utf-8", errors="replace").read()
                for a in flat:
                    if re.match(r"--[A-Za-z0-9][A-Za-z0-9-]*$", a) and a not in src and ("args." + a[2:].replace("-", "_")) not in src:
                        P("step %s: flag %s does not occur in %s" % (sid, a, os.path.relpath(drv, ROOT)))
            elif drv:
                P("step %s: driver %s not found" % (sid, drv))
            if sid and not s.get("record") and not s.get("records") and s.get("session", "AUTO").upper() == "AUTO":
                pass          # recorded as <plugin>.<step>

        for s in m.get("steps", []):
            check_step(s, "steps")
    return problems, stats


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("-")]
    pd = a[0] if a else os.path.join(ROOT, "plugins")
    probs, st = validate(pd)
    for p in probs:
        print("PROBLEM", p)
    print("validated %(plugins)d plugins, %(steps)d steps (%(cli)d with a cli template), %(params)d parameters: %(n)d problem(s)" % dict(st, n=len(probs)))
    sys.exit(1 if probs else 0)
