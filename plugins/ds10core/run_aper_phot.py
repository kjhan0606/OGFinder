#!/usr/bin/env python3
"""Launch ``python -m ds10core.tools.aper_phot`` (star aperture photometry of the shared core) for the Desktop shell.

Driver options:  --core-dir DIR  --stars regions|catalog|detect  --picks FILE (ds9 regions of the clicked stars, image coordinates)
--frames "a.fits,b.fits" (time series)  --lc-png FILE (light-curve plot, matplotlib if available).
Flags of the aperture photometry tool that the plugin manifest passes through unchanged (listed here so that tools/validate_manifests.py can
check the manifest against this driver; the real definitions are the argparse options of ``ds10core/tools/aper_phot.py``):
--out --summary --regions --cog-out --series-out --lc-out --aavso-out --catalog --find-thresh --max-stars --centroid --recenter-max
--radii --unit --fwhm --sky-in --sky-out --sky-method --gain --saturate --apcor --zeropoint --ref-mag --ref-colour --target --comps
--check --track --starid --obscode --use-aperture
"""
from __future__ import annotations

import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ds10 import alias_env, find_core  # noqa: E402  (reuse core discovery)

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


def picks_from_regions(text: str) -> list[tuple[float, float, str]]:
    """Centres of the point / circle regions of a ds9 region file in image coordinates (1-based), with their text labels."""
    out = []
    image = True
    for ln in text.splitlines():
        s = ln.strip()
        low = s.lower()
        if low in ("image", "physical"):
            image = True
            continue
        if low in ("fk5", "icrs", "fk4", "galactic", "ecliptic", "wcs"):
            image = False
            continue
        m = re.match(r"^\s*(?:image;\s*)?(point|circle)\s*\(\s*(" + NUM + r")\s*,\s*(" + NUM + r")", s, re.I)
        if m and image:
            lab = re.search(r"text=\{([^}]*)\}", s)
            out.append((float(m.group(2)), float(m.group(3)), lab.group(1).strip() if lab else ""))
    return out


def main(argv: list[str]) -> int:
    alias_env()
    args = list(argv)
    opt = {"--core-dir": None, "--stars": "regions", "--picks": "", "--frames": "", "--lc-png": "", "--catalog": ""}
    rest = []
    i = 0
    while i < len(args):
        if args[i] in opt and i + 1 < len(args):
            opt[args[i]] = args[i + 1]
            i += 2
        else:
            rest.append(args[i])
            i += 1
    d, why = find_core(opt["--core-dir"] or None)
    if not d:
        sys.stderr.write(why + "\n")
        return 2
    mode = opt["--stars"]
    tool = []
    if mode == "regions":
        picks = []
        if opt["--picks"] and os.path.exists(opt["--picks"]):
            with open(opt["--picks"], encoding="utf-8", errors="replace") as fh:
                picks = picks_from_regions(fh.read())
        if not picks:
            sys.stderr.write("No point/circle regions on the image: click on the stars first (ds9 Region menu, shape circle or point); detecting stars instead.\n")
            tool += ["--find"]
        else:
            tool += ["--positions", ";".join(f"{x:.3f} {y:.3f} {lab}".strip() for x, y, lab in picks), "--coords", "pixel"]
            sys.stderr.write(f"{len(picks)} picked star(s) from the image regions\n")
    elif mode == "catalog":
        cat = opt["--catalog"]
        if not (cat and os.path.exists(cat) and os.path.getsize(cat) > 0):
            sys.stderr.write("Stars = catalog, but no catalogue is loaded in the table\n")
            return 2
        tool += ["--catalog", cat]
    else:
        tool += ["--find"]
    for f in [x.strip() for x in opt["--frames"].split(",") if x.strip()]:
        tool += ["--frame", f]
    env = os.environ.copy()
    env["PYTHONPATH"] = d + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    py = env.get("DS10_PYTHON") or env.get("ASTRAFEX_PYTHON") or sys.executable
    # output mode "set": the catalogue table is filled from this process's stdout -> run the tool (it writes --out), then echo the catalogue
    rc = subprocess.call([py, "-m", "ds10core.tools.aper_phot", *rest, *tool], env=env, stdout=sys.stderr)
    out = rest[rest.index("--out") + 1] if "--out" in rest[:-1] else None
    lc = rest[rest.index("--lc-out") + 1] if "--lc-out" in rest[:-1] else None
    if rc == 0 and lc and opt["--lc-png"] and os.path.exists(lc):
        lc_plot(lc, opt["--lc-png"])
    if rc == 0 and out and os.path.exists(out):
        with open(out, encoding="utf-8") as fh:
            sys.stdout.write(fh.read())
    return rc


def lc_plot(lc_tsv: str, png: str) -> None:
    """Light-curve PNG (target and check star, differential magnitudes) when matplotlib is available."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        sys.stderr.write("matplotlib not available: no light-curve plot\n")
        return
    import csv
    with open(lc_tsv, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    def col(k):
        return [float(r[k]) if r.get(k) not in (None, "", "na") else float("nan") for r in rows]
    t = col("BJD_TDB") if all(r.get("BJD_TDB") not in ("", "na", None) for r in rows) else col("JD_UTC")
    t0 = int(min(t))
    fig, ax = plt.subplots(figsize=(6.4, 3.4), dpi=100)
    ax.errorbar([x - t0 for x in t], col("DMAG"), yerr=col("DMAG_ERR"), fmt="o", ms=3, label="target - ensemble")
    if rows and "DMAG_CHECK" in rows[0]:
        k = col("DMAG_CHECK")
        import statistics
        km = statistics.median([v for v in k if v == v])
        dm = statistics.median([v for v in col("DMAG") if v == v])
        ax.errorbar([x - t0 for x in t], [v - km + dm + 0.1 for v in k], yerr=col("DMAG_CHECK_ERR"), fmt="s", ms=3, color="tab:green", label="check (offset)")
    ax.invert_yaxis()
    ax.set_xlabel(f"{'BJD_TDB' if 'BJD' in str(rows[0].get('BJD_TDB', '')) or rows[0].get('BJD_TDB') else 'JD'} - {t0}")
    ax.set_ylabel("delta mag")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(png)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
