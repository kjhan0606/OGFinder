#!/usr/bin/env python3
"""Run the OGFinder Moving Objects CLI steps on a set of FLC/FLT exposures and write a manifest.

Standard library only; it calls ds9/library/ds9_moving.py with the python given by --python.
Steps: align, difference, link, identify (each is also a menu entry of Moving Objects).
The manifest (manifest.json) records the command lines, the versions of the python packages used by
the driver, and the sha256 of every input and output file.  Nothing is uploaded anywhere.

Example:
  scripts/run_moving_pipeline.py --workdir wd --ra 150.1375 --dec 2.3610 --files a_flc.fits b_flc.fits c_flc.fits d_flc.fits
"""
import argparse, hashlib, json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
DRIVER = os.path.join(HERE, "..", "ds9", "library", "ds9_moving.py")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--files", nargs="+", required=True)
    ap.add_argument("--ra", type=float, required=True)
    ap.add_argument("--dec", type=float, required=True)
    ap.add_argument("--half-pix", type=int, default=700)
    ap.add_argument("--snr", type=float, default=8.0)
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--skip-identify", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    base = [a.python, os.path.normpath(DRIVER), "--workdir", a.workdir]
    steps = [
        ("align", ["--mode", "align", "--files"] + a.files),
        ("difference", ["--mode", "difference", "--ra", str(a.ra), "--dec", str(a.dec), "--half-pix", str(a.half_pix), "--snr", "5", "--files"] + a.files),
        ("link", ["--mode", "link", "--snr", str(a.snr), "--files"] + a.files),
    ]
    if not a.skip_identify:
        steps.append(("identify", ["--mode", "identify", "--files"] + a.files))
    man = dict(started=time.strftime("%Y-%m-%dT%H:%M:%S%z"), python=sys.version.split()[0], steps=[], inputs={}, outputs={}, packages={})
    for f in a.files:
        man["inputs"][f] = sha(f)
    for st, args in steps:
        t = time.time()
        r = subprocess.run(base + args, capture_output=True, text=True)
        man["steps"].append(dict(step=st, argv=base + args, returncode=r.returncode, seconds=round(time.time() - t, 1)))
        sys.stderr.write("[%s] rc=%d %.1fs\n" % (st, r.returncode, time.time() - t))
        if r.returncode:
            sys.stderr.write(r.stderr[-2000:])
            break
    try:
        out = subprocess.run([a.python, "-c", "import numpy,scipy,astropy;print(numpy.__version__,scipy.__version__,astropy.__version__)"],
                             capture_output=True, text=True).stdout.split()
        man["packages"] = dict(zip(["numpy", "scipy", "astropy"], out))
    except Exception:
        pass
    for root, _, fs in os.walk(a.workdir):
        for f in fs:
            p = os.path.join(root, f)
            if f != "manifest.json" and os.path.getsize(p) < 500e6:
                man["outputs"][os.path.relpath(p, a.workdir)] = sha(p)
    with open(os.path.join(a.workdir, "manifest.json"), "w") as f:
        json.dump(man, f, indent=1)
    print(os.path.join(a.workdir, "manifest.json"))


if __name__ == "__main__":
    main()
