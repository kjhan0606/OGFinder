#!/usr/bin/env python3
"""OGFinder "Moving Objects" CLI driver (thin wrapper; the code lives in the `moving` package next to the tree).

Every step reads and writes explicit files inside a work directory (--workdir), so it can be recorded by the session
recorder as one argv and replayed:

  --mode setup        report versions, ephemeris files and provider reachability (JSON on stdout)
  --mode fetch        query the reference/epoch providers for imaging covering --ra/--dec; optionally download
                        --provider mast|lsst|both   (recorded in the step args)
                        --lsst-local-dir DIR        user-supplied Rubin FITS (no credentials are ever used)
  --mode align        align the exposures (--files) to Gaia DR3 (else relative chain); writes align/*.json
  --mode difference   template + ZOGY difference + detection + single-exposure classification -> detections.tsv, diff/*.fits
  --mode link         link detections into tracklets -> tracklets.json, movers.tsv, movers.reg
  --mode nightlink    link tracklets of SEVERAL nights/visits (--tracklet-files a/tracklets.json b/tracklets.json ...) by 2-body orbit fits -> nightlinks.json, nightlinks.tsv
                        --nl-nbody   re-fit each linked group in a separate process (ASSIST, else CODES)
                        --nl-ranging sample the orbital-element distribution of each linked group and pair (simplified statistical ranging)
  --mode identify     known-object identification (SkyBoT + Horizons) -> identified.tsv
  --mode orbit        orbit determination for one tracklet (--tracklet N) or a known object (--designation D with MPC obs)
  --mode transients   static transient candidates + host association (--catalog TSV) -> transients.tsv, transients.reg
  --mode lightcurve   forced photometry on the difference images at the transient positions
  --mode export       MPC 80-column / JSON elements for the orbit result (never submitted anywhere)
  --mode run          fetch(optional)+align+difference+link+identify+transients in one go (scripts/run_moving_pipeline.py)

Machine-readable status lines on stdout start with '#MOVING' (parsed by the Tcl GUI).  No credentials are read, asked for or stored.
"""
import os, sys, json, argparse, time


def _find_root():
    here = os.path.dirname(os.path.abspath(__file__))
    cands = []
    if os.environ.get("OGF_MOVING_DIR"):
        cands.append(os.path.dirname(os.path.abspath(os.environ["OGF_MOVING_DIR"])))
    cands += [os.path.dirname(os.path.dirname(here)), os.path.dirname(here), here,
              os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))]
    for c in cands:
        if c and os.path.isfile(os.path.join(c, "moving", "__init__.py")):
            return c
    return None


_r = _find_root()
if _r and _r not in sys.path:
    sys.path.insert(0, _r)
try:
    import numpy as np
    from moving import util
except Exception as e:                                    # pragma: no cover
    sys.stderr.write("ERROR: cannot import the moving package (%s).\nSet OGF_MOVING_DIR=<OGFinder>/moving and OGFINDER_PYTHON to a Python "
                     "with numpy/astropy/sep/reproject/rebound/assist (see docs/moving_objects.md).\n" % e)
    sys.exit(2)


def status(**kw):
    print("#MOVING " + json.dumps(kw, default=str), flush=True)


def _wd(a):
    return util.ensure_dir(os.path.abspath(os.path.expanduser(a.workdir)))


def _files(a):
    fs = []
    for f in (a.files or []):
        for g in f.split(","):
            if g:
                fs.append(os.path.abspath(os.path.expanduser(g)))
    return fs


def _load_all(a, apply_align=True):
    from moving import imaging as I, align as A
    wd = _wd(a)
    chips = []
    for f in _files(a):
        chips += I.load_chips(f)
    if apply_align:
        for c in chips:
            A.load_sidecar_into(c, os.path.join(wd, "align"))
    return chips


# -------------------------------------------------------------------------------------------------- modes
def m_setup(a):
    from moving import orbit as O, lsst
    out = dict(versions=util.versions(), ephemeris=dict(zip(("planets", "asteroids"), O.ephem_paths())) if hasattr(O, "ephem_paths") else {},
               assist=O.assist_available())
    if a.probe:
        out["lsst_probe"] = lsst.probe()
    print(json.dumps(out, indent=1, default=str))


def m_fetch(a):
    from moving import mast, lsst
    wd = _wd(a)
    prov = a.provider
    res = dict(provider=prov, ra=a.ra, dec=a.dec, radius_arcmin=a.radius, candidates=[], status={})
    cands = []
    if prov in ("mast", "both"):
        cc = mast.query_candidates(a.ra, a.dec, a.radius, collections=tuple(a.collections.split(",")),
                                   exposures_only=a.exposures_only)
        for c in cc:
            c["provider"] = "MAST"
        cands += cc
        res["status"]["mast"] = dict(available=True, n=len(cc))
    if prov in ("lsst", "both"):
        lc, st = lsst.query_candidates(a.ra, a.dec, a.radius, local_dir=a.lsst_local_dir)
        for c in lc:
            c["provider"] = "LSST"
        cands += lc
        res["status"]["lsst"] = st
        if not st.get("available"):
            status(kind="provider_unavailable", provider="LSST", reason=st.get("reason"))
            if prov == "lsst" and not any(c["provider"] == "MAST" for c in cands):
                # documented fallback: MAST HST/JWST + Pan-STARRS
                cc = mast.query_candidates(a.ra, a.dec, a.radius, collections=("HST", "JWST", "PS1"))
                for c in cc:
                    c["provider"] = "MAST(fallback)"
                cands += cc
                res["status"]["fallback"] = "MAST/PS1 (LSST unavailable without Rubin Science Platform login)"
    mast_c = [c for c in cands if c["provider"].startswith("MAST")]
    mast.group_epochs(mast_c, gap_days=1.0)
    res["candidates"] = cands
    util.write_json(os.path.join(wd, "candidates.json"), res)
    cols = ["provider", "obs_id", "obsid", "collection", "instrument", "filter", "texp", "mjd", "epoch", "pixscale", "depth", "proposal", "kind", "target"]
    util.write_tsv(os.path.join(wd, "candidates.tsv"), cands, cols)
    paths = []
    if a.download:
        sel = [c for c in mast_c if c.get("kind") == "exposure" and (not a.filter or a.filter.lower() in c["filter"].lower())]
        if a.obsids:
            want = set(a.obsids.split(","))
            sel = [c for c in mast_c if c["obsid"] in want or c["obs_id"] in want]
        for c in sel[:a.max_download]:
            try:
                paths += mast.download(c["obsid"], collection=c["collection"])
            except Exception as e:
                util.log("download failed %s: %s" % (c["obs_id"], e))
        if any(c["provider"] == "LSST" for c in cands):
            paths += [c["path"] for c in cands if c["provider"] == "LSST" and c.get("path")]
        util.write_json(os.path.join(wd, "files.json"), dict(files=[dict(path=p, sha256=util.sha256_file(p), bytes=os.path.getsize(p)) for p in paths]))
    status(kind="fetch", n_candidates=len(cands), n_downloaded=len(paths), provider=prov, outputs=[os.path.join(wd, "candidates.tsv")] + paths)
    print("provider=%s candidates=%d downloaded=%d" % (prov, len(cands), len(paths)))


def m_align(a):
    from moving import align as A, imaging as I
    wd = _wd(a)
    chips = _load_all(a, apply_align=False)
    st = A.align_field(chips)
    ad = os.path.join(wd, "align")
    rep = []
    for c in chips:
        if getattr(c, "align", None):
            A.save_sidecar(c, ad)
        s = st.get(c.name, {})
        g = s.get("gaia", {}); r = s.get("relative", {})
        rep.append(dict(chip=c.name, gaia_ok=bool(g.get("ok")), gaia_n=g.get("n_match"), gaia_rms_mas=g.get("rms_mas"),
                        rel_ok=r.get("ok"), rel_n=r.get("n_match"), rel_rms_mas=r.get("rms_mas")))
    util.write_json(os.path.join(wd, "align_report.json"), rep)
    status(kind="align", n_chips=len(chips), outputs=[ad, os.path.join(wd, "align_report.json")])
    for r in rep:
        print("%-34s gaia n=%s rms=%s mas | rel n=%s rms=%s mas" % (r["chip"], r["gaia_n"], None if r["gaia_rms_mas"] is None else round(r["gaia_rms_mas"], 1), r["rel_n"], None if r["rel_rms_mas"] is None else round(r["rel_rms_mas"], 1)))


def _detect_opts(a):
    """Detection-stage options (item 1).  The defaults are ON and need no argv; the flags below turn single parts off (recorded in the session argv only when used)."""
    o = {}
    if getattr(a, "no_cr_reject", False): o["cr_reject"] = False
    if getattr(a, "no_trail_fit", False): o["trail_fit"] = False
    if getattr(a, "no_realbogus", False): o["realbogus"] = False
    # ZOGY options (item 2): all default OFF, so the default argv and results are unchanged
    if getattr(a, "source_noise", False): o["source_noise"] = True
    am = getattr(a, "astrom", None)
    if am not in (None, "", "off"):
        o["astrom_sigma"] = True if am == "sidecar" else "measure" if am == "measure" else float(am)
    if getattr(a, "psf_tile", None): o["psf_tile"] = int(a.psf_tile)
    if getattr(a, "psf_source", None) not in (None, "tiles"): o["psf_source"] = a.psf_source
    if getattr(a, "template_psf", None) not in (None, "target"): o["template_psf"] = a.template_psf
    return o


def m_difference(a):
    from moving import pipeline as P
    wd = _wd(a)
    chips = _load_all(a)
    center = (a.ra, a.dec) if a.ra is not None and a.dec is not None else None
    dets, infos = P.detect_in_region(chips, center=center, half_pix=a.half_pix, snr_det=a.snr, save_diff_dir=os.path.join(wd, "diff"),
                                     progress=lambda m: print(m, flush=True), **_detect_opts(a))
    rows = []
    for i, d in enumerate(dets):
        r = {k: d[k] for k in ("ex", "chip", "file", "t", "ra", "dec", "x_chip", "y_chip", "sign", "snr", "flux_e_s", "flux_err", "a_pix", "b_pix",
                               "elong", "pa_deg", "sharp", "tpl_snr", "neg_frac", "on_cr", "near_bad", "channel", "cls", "zp_ab", "pixscale", "filter", "texp", "rb", "lac3", "lac_n7", "trail_fit", "trail_len_pix", "theta_pix") if k in d}
        r["id"] = i; r["why"] = "; ".join(d.get("why", [])); r["trail_len_arcsec"] = d.get("trail_len_arcsec", "")
        r["sig_pos_arcsec"] = d.get("sig_pos_arcsec", "")
        rows.append(r)
    util.write_tsv(os.path.join(wd, "detections.tsv"), rows)
    util.write_json(os.path.join(wd, "difference_info.json"), infos)
    status(kind="difference", n_detections=len(rows), outputs=[os.path.join(wd, "detections.tsv"), os.path.join(wd, "diff")])
    print("detections=%d" % len(rows))


def _read_dets(wd):
    rows = util.read_tsv(os.path.join(wd, "detections.tsv"))
    for r in rows:
        r.setdefault("why", "")
        r["why"] = [r["why"]] if isinstance(r["why"], str) else r["why"]
        r["on_cr"] = str(r.get("on_cr")) in ("True", "1")
        r["near_bad"] = str(r.get("near_bad")) in ("True", "1")
    return rows


def m_link(a):
    from moving import pipeline as P
    wd = _wd(a)
    chips = _load_all(a, apply_align=False)
    offs, tm = P.hst_offsets(chips) if all(c.tel == "HST" for c in chips) else ({}, None)
    dets = _read_dets(wd)
    for d in dets:
        d["ex"] = int(d["ex"])
    obs_off = offs if offs else None
    # exposures index in detections.tsv follows the sorted file order used by detect_in_region
    shapes = {c.name: tuple(c.shape) for c in chips}
    vst = {}
    trs = P.link_detections(dets, obs_off, snr_min=a.snr, tol_arcsec=a.tol, min_exposures=a.min_exposures, max_per_exposure=a.max_per_exposure,
                            chip_shapes=shapes, veto_stats=vst, orbit_prior=a.orbit_prior, orbit_prior_fraction=a.orbit_prior_fraction)
    trs = trs[:a.max_tracklets]
    if vst:
        util.log("veto: " + ", ".join("%s=%d" % kv for kv in sorted(vst.items())))
    for k, t in enumerate(trs):
        t["id"] = k
    util.write_json(os.path.join(wd, "tracklets.json"), trs)
    rows = [dict(id=t["id"], ra=t["ra"][len(t["ra"]) // 2], dec=t["dec"][len(t["dec"]) // 2], n=t["n"], rate_arcsec_h=t["rate_ash"], sig_rate=t["sig_rate_ash"],
                 pa_deg=t["pa_deg"], rms_arcsec=t["rms_resid_arcsec"], score=t["score"], inv_delta=t["inv_delta"], status="candidate") for t in trs]
    util.write_tsv(os.path.join(wd, "movers.tsv"), rows)
    write_regions(os.path.join(wd, "movers.reg"), trs)
    status(kind="link", n_tracklets=len(trs), outputs=[os.path.join(wd, "tracklets.json"), os.path.join(wd, "movers.tsv"), os.path.join(wd, "movers.reg")])
    print("tracklets=%d" % len(trs))


def write_regions(path, trs, known=None):
    """DS9 region file: for each tracklet a vector (arrow, length ~ 30 s of motion scaled for visibility) at the first position, circles at
    every member, text label with rate / PA."""
    import math
    L = ["# Region file format: DS9 version 4.1", "global color=cyan width=1", "fk5"]
    for t in trs:
        col = "green" if known and known.get(t["id"]) else "cyan"
        ra0, de0 = t["ra"][0], t["dec"][0]
        length = max(2.0, min(8.0, t["rate_ash"] * 0.1))
        # DS9 stores vector regions as "# vector(...)" comment-style lines (a bare "vector(" is a syntax error)
        L.append("# vector(%.6f,%.6f,%.2f\",%.1f) vector=1 color=%s text={M%d %.1f as/h PA%.0f}" % (ra0, de0, length, (90 - t["pa_deg"]) % 360, col, t["id"], t["rate_ash"], t["pa_deg"]))
        for r, d in zip(t["ra"], t["dec"]):
            L.append("circle(%.6f,%.6f,0.35\") # color=%s" % (r, d, col))
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


def m_identify(a):
    from moving import identify
    wd = _wd(a)
    trs = util.read_json(os.path.join(wd, "tracklets.json"))
    ids = {}
    rows = []
    for t in trs[:a.max_tracklets]:
        r = identify.identify_tracklet(t, radius_arcmin=a.radius, tol_arcsec=a.tol)
        t["identification"] = r
        best = r["matches"][0] if r["matches"] else None
        known = r["status"] == "known"
        ids[t["id"]] = known
        rows.append(dict(id=t["id"], status=r["status"], name=best["name"] if best else "", number=best["num"] if best else "",
                         vmag=best["vmag"] if best else "", max_sep_arcsec=best["max_sep"] if best else "", n_skybot_candidates=len(r["matches"])))
    util.write_json(os.path.join(wd, "tracklets.json"), trs)
    util.write_tsv(os.path.join(wd, "identified.tsv"), rows)
    mv = util.read_tsv(os.path.join(wd, "movers.tsv"))
    st = {r["id"]: r for r in rows}
    for m in mv:
        s = st.get(int(m["id"]))
        if s:
            m["status"] = s["status"] if s["status"] == "unknown" else "known: " + (s["number"] + " " + s["name"]).strip()
    util.write_tsv(os.path.join(wd, "movers.tsv"), mv)
    write_regions(os.path.join(wd, "movers.reg"), trs, known=ids)
    status(kind="identify", n_known=sum(ids.values()), n_unknown=len(ids) - sum(ids.values()), outputs=[os.path.join(wd, "identified.tsv")])
    print("known=%d unknown=%d" % (sum(ids.values()), len(ids) - sum(ids.values())))


def _obs_from_tracklet(t, sig_floor=0.05, hst=True):
    from moving import orbit as O, horizons, obs as OBS
    n = len(t["t"])
    sig = np.maximum(np.array(t["sig"]), sig_floor)
    stn = "250"
    jd = util.utc_mjd_to_tdb_jd(np.array(t["t"]))
    sat = horizons.vectors("-48", jd, center="500@399")[:, :3] * util.AU_KM
    return O.Obs(t["t"], t["ra"], t["dec"], sig, sig, [stn] * n, sat_xyz_km=[s for s in sat])


def m_orbit(a):
    from moving import orbit as O, orbitfit as OF, photometry as PH, obs as OBS, mpc, kepler as K
    wd = _wd(a)
    out = dict(mode="orbit")
    if a.designation:
        j = mpc.get_obs(a.designation)
        allo = [o for o in OBS.parse_obs80(j["OBS80"]) if o["stn"] != "500" and
                (OBS.obscodes().get(o["stn"], (None,))[0] is not None or o.get("sat_xyz_km") is not None)]
        sel = [o for o in allo if (a.mjd_min is None or o["mjd_utc"] >= a.mjd_min) and (a.mjd_max is None or o["mjd_utc"] <= a.mjd_max)]
        ra = []; de = []
        for o in sel:
            dr, dd = OBS.debias(o["ra"], o["dec"], o["mjd_utc"], o["cat"]) if not a.no_debias else (0.0, 0.0)
            ra.append(o["ra"] - dr / np.cos(np.radians(o["dec"])) / 3600); de.append(o["dec"] - dd / 3600)
        yr = int(2000 + (np.mean([o["mjd_utc"] for o in sel]) - 51544) / 365.25)
        sig = [max(OBS.station_sigma(o["stn"], yr), 0.1) for o in sel]
        obs = O.Obs([o["mjd_utc"] for o in sel], ra, de, sig, sig, [o["stn"] for o in sel], sat_xyz_km=[o["sat_xyz_km"] for o in sel],
                    mag=[o.get("mag") for o in sel], band=[o.get("band", "") for o in sel])
        label = a.designation
    else:
        trs = util.read_json(os.path.join(wd, "tracklets.json"))
        t = trs[a.tracklet]
        obs = _obs_from_tracklet(t)
        label = "tracklet%d" % t["id"]
    prop = O.Propagator()
    res = OF.fit_orbit(obs, prop, ranging_samples=a.samples)
    el, sg = res["elements"], res["sigmas"]
    out.update(label=label, n_obs=int(obs.n), n_used=res["n_used"], rms_arcsec=res["rms_arcsec"], chi2_red=res["chi2_red"],
               arc_days=float(np.ptp(obs.mjd_utc)), epoch_jd=res["epoch_jd"], orbit_class=res["orbit_class"],
               class_status=res.get("class_status", ""), orbit_class_probs=res.get("orbit_class_probs", {}),
               class_note=res.get("class_note", ""), orbit_class_bestfit=res.get("orbit_class_bestfit", ""),
               elements={k: dict(value=el[k], sigma=sg[k]) for k in ("a", "e", "i", "om", "w", "ma", "q")},
               force_model=res["force_model"], residuals=dict(mjd=obs.mjd_utc.tolist(), dra=res["residual_ra"].tolist(), ddec=res["residual_dec"].tolist(),
                                                              active=res["active"].tolist(), stn=obs.stn))
    out["determined"] = res.get("determined", True)
    if not out["determined"]:
        out["note"] = res.get("status_note", "")
        out["ranging_quantiles_16_50_84"] = res["info"].get("ranging_quantiles_16_50_84")
        out["elements_unconstrained"] = out.pop("elements")
        out["elements"] = {k: dict(value=float("nan"), sigma=float("nan")) for k in ("a", "e", "i", "om", "w", "ma", "q")}
    rg = res["info"].get("ranging")
    if rg:
        out["ranging"] = dict(neff=rg["neff"], class_probs=rg["class_probs"])
        if not out["determined"]:
            out["ranging_note"] = "arc too short for a unique orbit: classification probabilities from statistical ranging (simplified); the elements are NOT constrained"
    util.write_json(os.path.join(wd, "orbit_%s.json" % label), out)
    with open(os.path.join(wd, "orbit_%s.kv" % label), "w") as f:         # flat file for the Tcl dialog
        for k in ("label", "n_obs", "n_used", "rms_arcsec", "chi2_red", "arc_days", "epoch_jd", "orbit_class", "force_model"):
            f.write("info\t%s\t%s\n" % (k, out[k]))
        f.write("info\tdetermined\t%s\n" % out["determined"])
        if not out["determined"]:
            f.write("info\tnote\t%s\n" % out["note"])
            for k, v in (out.get("ranging_quantiles_16_50_84") or {}).items():
                f.write("rq\t%s\t%.4g\t%.4g\t%.4g\n" % ((k,) + tuple(v)))
        for k, v in out["elements"].items():
            f.write("elem\t%s\t%.10g\t%.4g\n" % (k, v["value"], v["sigma"]))
        f.write("info\tclass_status\t%s\n" % out["class_status"])
        for k, v in (out.get("orbit_class_probs") or out.get("ranging", {}).get("class_probs", {}) or {}).items():
            f.write("class\t%s\t%.4f\n" % (k, v))
        r = out["residuals"]
        for i in range(len(r["mjd"])):
            f.write("res\t%.6f\t%.4f\t%.4f\t%d\t%s\n" % (r["mjd"][i], r["dra"][i], r["ddec"][i], int(r["active"][i]), r["stn"][i]))
    status(kind="orbit", label=label, outputs=[os.path.join(wd, "orbit_%s.json" % label)])
    print(json.dumps({k: v for k, v in out.items() if k != "residuals"}, indent=1, default=str))


def m_nightlink(a):
    """Cross-night linking of tracklets.json files (moving.nightlink): groups of tracklets from different nights consistent with one Sun-bound orbit."""
    from moving import nightlink as NL
    wd = _wd(a)
    files = a.tracklet_files or [os.path.join(wd, "tracklets.json")]
    trks = []
    src = []
    for fi, fn in enumerate(files):
        for t in util.read_json(fn):
            mjd = t.get("t") or t.get("mjd")
            if not mjd or len(mjd) < 2:
                continue
            sig = np.maximum(np.asarray(t.get("sig", [a.nl_sigma] * len(mjd)), float), 0.05)
            trks.append(NL.Tracklet(mjd, t["ra"], t["dec"], sig, code=a.obs_code, tid=len(trks), night=None))
            src.append(dict(file=fn, tracklet=t.get("id")))
    res = NL.link_nights(trks, max_gap_days=a.nl_max_gap, chi2_max=a.nl_chi2, floor0=a.nl_floor, min_tracklets=a.nl_min)
    if a.nl_nbody:
        from moving import nbody_refine as NR
        res = NR.refine_result(trks, res, chi2_max=a.nl_chi2)
    if a.nl_ranging:
        from moving import arc_ranging as AR
        res = AR.attach(trks, res, n_samples=a.nl_ranging_samples, n_propagate=a.nl_ranging_propagate)
    rows = []
    for gi, g in enumerate(res["groups"] + res["pairs"]):
        row = dict(group=gi, kind="linked" if gi < len(res["groups"]) else "pair", n_tracklets=len(g["ids"]), n_obs=g["n_obs"], nights=",".join(str(x) for x in g["nights"]),
                   members=";".join("%s#%s" % (os.path.basename(os.path.dirname(src[i]["file"])) or src[i]["file"], src[i]["tracklet"]) for i in g["ids"]),
                   a_au=round(g["elements"].get("a", float("nan")), 3), e=round(g["elements"].get("e", float("nan")), 3),
                   inc_deg=round(g["elements"].get("inc_deg", float("nan")), 2), chi2_red=round(g["chi2_red"], 3), rms_arcsec=round(g["rms_arcsec"], 3))
        if g.get("nbody"):
            row["nbody_status"] = g["nbody"].get("status")
            row["nbody_accepted"] = bool(g["nbody"].get("accepted"))
            row["two_body_raw_rms_arcsec"] = g.get("two_body_raw_rms_arcsec")
        if g.get("ranging"):
            row["ranging"] = g["ranging"]
        rows.append(row)
    rejected_rows = []
    for g in res.get("nbody_rejected", []):
        rejected_rows.append(dict(ids=list(g["ids"]), nights=list(g["nights"]), nbody=g.get("nbody"),
                                  two_body_elements=g.get("two_body_elements"), two_body_raw_rms_arcsec=g.get("two_body_raw_rms_arcsec")))
    util.write_json(os.path.join(wd, "nightlinks.json"), dict(groups=rows, unlinked=[src[i] for i in res["unlinked"]], stats=res["stats"],
                                                               nbody=res.get("nbody"), nbody_rejected=rejected_rows,
                                                               ranging=res.get("ranging"),
                                                               params=dict(obs_code=a.obs_code, min_tracklets=a.nl_min, chi2_max=a.nl_chi2, floor_arcsec=a.nl_floor, max_gap_days=a.nl_max_gap, nbody=bool(a.nl_nbody), ranging=bool(a.nl_ranging))))
    with open(os.path.join(wd, "nightlinks.tsv"), "w") as f:
        cols = ["group", "kind", "n_tracklets", "n_obs", "nights", "members", "a_au", "e", "inc_deg", "chi2_red", "rms_arcsec"]
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    status(kind="nightlink", n_tracklets=len(trks), n_groups=len(res['groups']), n_pairs=len(res['pairs']), n_unlinked=len(res["unlinked"]),
           outputs=[os.path.join(wd, "nightlinks.json"), os.path.join(wd, "nightlinks.tsv")])


def m_transients(a):
    from moving import transients as TR
    wd = _wd(a)
    dets = _read_dets(wd)
    trs = util.read_json(os.path.join(wd, "tracklets.json")) if os.path.exists(os.path.join(wd, "tracklets.json")) else []
    excl = set()
    for t in trs:
        if t.get("n", 0) >= 3:
            excl.update(int(m) for m in t["members"])
    for d in dets:
        d["ra"] = float(d["ra"]); d["dec"] = float(d["dec"]); d["snr"] = float(d["snr"]); d["sign"] = int(float(d["sign"]))
    cls_ok = ("point", "faint") + (("artefact_cr",) if a.include_cr_flagged else ())
    groups = TR.group_static(dets, tol_arcsec=a.tol, min_epochs=a.min_epochs, snr_min=a.snr, exclude=excl, classes=cls_ok)
    rows_cat, cmap = TR.read_catalog(a.catalog) if a.catalog else ([], {})
    rows = []
    for gi, g in enumerate(groups):
        ra = float(np.mean([dets[i]["ra"] for i in g])); de = float(np.mean([dets[i]["dec"] for i in g]))
        px = float(np.median([float(dets[i].get("pixscale") or 0.05) for i in g]))
        host = TR.associate_host(ra, de, rows_cat, cmap, px) if rows_cat else None
        cls, why = TR.heuristic_class(host, None)
        r = dict(id=gi, ra=ra, dec=de, n_det=len(g), files=len({dets[i]["file"] for i in g}), snr_max=max(float(dets[i]["snr"]) for i in g),
                 host_sep_arcsec=host["sep_arcsec"] if host else "", host_id=host.get("id") if host else "", host_z=host.get("z", "") if host else "",
                 offset_re=host.get("offset_in_re", "") if host else "", heuristic=cls, why="; ".join(why))
        rows.append(r)
    util.write_tsv(os.path.join(wd, "transients.tsv"), rows)
    with open(os.path.join(wd, "transients.reg"), "w") as f:
        f.write("# Region file format: DS9 version 4.1\nglobal color=yellow\nfk5\n")
        for r in rows:
            f.write("point(%.6f,%.6f) # point=cross 12 color=yellow text={T%d}\n" % (r["ra"], r["dec"], r["id"]))
    status(kind="transients", n_candidates=len(rows), outputs=[os.path.join(wd, "transients.tsv"), os.path.join(wd, "transients.reg")])
    print("transient candidates=%d (static, >=%d files, not in mover tracklets)" % (len(rows), a.min_epochs))


def m_lightcurve(a):
    from moving import transients as TR
    wd = _wd(a)
    rows = util.read_tsv(os.path.join(wd, "transients.tsv"))
    chips = _load_all(a)
    diffs = []
    for fn in sorted(os.listdir(os.path.join(wd, "diff"))) if os.path.isdir(os.path.join(wd, "diff")) else []:
        from astropy.io import fits
        from astropy.wcs import WCS
        with fits.open(os.path.join(wd, "diff", fn)) as h:
            diffs.append((None, WCS(h[0].header), np.array(h[0].data, float), None, fn))
    if not diffs:
        print("no difference images: run Difference first"); return
    # time and zero point from the chips by name
    cm = {c.name.replace("[", "_").replace("]", "").replace(",", "_"): c for c in chips}
    dd = []
    for t, w, im, zp, fn in diffs:
        key = fn[len("diff_"):-len(".fits")]
        c = cm.get(key) or next((v for k, v in cm.items() if key.startswith(k.split("_SCI")[0]) and key.endswith(k.split("_SCI")[1] if "_SCI" in k else "")), None)
        if c is not None:
            dd.append((float(c.mjd), w, im, c.zp_ab, c.texp))
    lcs = TR.lightcurve([(float(r["ra"]), float(r["dec"])) for r in rows], dd, r_pix=a.aperture)
    util.write_json(os.path.join(wd, "lightcurves.json"), [dict(id=int(r["id"]), **lc) for r, lc in zip(rows, lcs)])
    status(kind="lightcurve", n=len(rows), outputs=[os.path.join(wd, "lightcurves.json")])
    for r, lc in zip(rows, lcs):
        print("T%s: %d epochs, mags %s" % (r["id"], len(lc["t"]), [None if not np.isfinite(m) else round(m, 2) for m in lc["mag"]]))


def m_export(a):
    from moving import obs as OBS
    wd = _wd(a)
    oj = a.orbit_json
    if not oj:
        cands = ["orbit_tracklet%d.json" % a.tracklet] if a.tracklet is not None else []
        if a.designation:
            cands.append("orbit_%s.json" % a.designation)
        oj = next((c for c in cands if os.path.exists(os.path.join(wd, c))), None)
        if oj is None:
            import glob
            g = sorted(glob.glob(os.path.join(wd, "orbit_*.json")), key=os.path.getmtime)
            oj = os.path.basename(g[-1]) if g else None
    if not oj:
        print("no orbit_*.json in work directory: run Orbit Fit first"); return
    j = util.read_json(os.path.join(wd, oj))
    trs = util.read_json(os.path.join(wd, "tracklets.json")) if os.path.exists(os.path.join(wd, "tracklets.json")) else []
    out = {}
    if a.tracklet is not None and trs:
        t = trs[a.tracklet]
        o = [dict(mjd_utc=tt, ra=r, dec=d, stn="250", sat_xyz_km=None) for tt, r, d in zip(t["t"], t["ra"], t["dec"])]
        # HST positions as 's' lines need km; they are re-derived through Horizons
        from moving import horizons
        jd = util.utc_mjd_to_tdb_jd(np.array(t["t"]))
        sat = horizons.vectors("-48", jd, center="500@399")[:, :3] * util.AU_KM
        for rec, s in zip(o, sat):
            rec["sat_xyz_km"] = s
        txt = OBS.write_obs80(o, desig=a.designation or "", stn="250")
        fn = os.path.join(wd, "tracklet%d.obs80" % a.tracklet)
        open(fn, "w").write(txt)
        out["obs80"] = fn
    fn = os.path.join(wd, "elements_export.json")
    util.write_json(fn, dict(orbit=j, note="local export only; nothing is submitted to the MPC"))
    out["json"] = fn
    status(kind="export", outputs=list(out.values()))
    print(json.dumps(out))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", required=True, choices="setup fetch align difference link nightlink identify orbit transients lightcurve export".split())
    ap.add_argument("--workdir", default="moving_work")
    ap.add_argument("--files", nargs="*")
    ap.add_argument("--ra", type=float); ap.add_argument("--dec", type=float)
    ap.add_argument("--radius", type=float, default=1.0, help="arcmin (fetch/identify)")
    ap.add_argument("--provider", default="mast", choices=["mast", "lsst", "both"])
    ap.add_argument("--lsst-local-dir"); ap.add_argument("--collections", default="HST,JWST,PS1")
    ap.add_argument("--exposures-only", action="store_true")
    ap.add_argument("--download", action="store_true"); ap.add_argument("--max-download", type=int, default=6)
    ap.add_argument("--filter"); ap.add_argument("--obsids")
    ap.add_argument("--half-pix", type=int, default=700)
    ap.add_argument("--snr", type=float, default=8.0)
    ap.add_argument("--no-cr-reject", action="store_true", help="difference: skip the L.A.Cosmic features")
    ap.add_argument("--no-trail-fit", action="store_true", help="difference: skip the trailed-PSF centroid refit")
    ap.add_argument("--no-realbogus", action="store_true", help="difference: skip the real/bogus score (the linker then orders its pool by S/N)")
    ap.add_argument("--source-noise", action="store_true", help="difference: ZOGY source-noise variance terms (Vn, Vr; fewer negative/static residuals)")
    ap.add_argument("--astrom", default="off", help="difference: ZOGY astrometric-registration term: off (default) | sidecar (alignment rms per chip) | "
                    "measure (estimated from the star offsets of each target/template pair) | a number = sigma in pixels")
    ap.add_argument("--psf-tile", type=int, default=0, help="difference: spatially varying PSF on tiles of N pixels (tiled ZOGY); 0 = one PSF per chip")
    ap.add_argument("--psf-source", default="tiles", choices=["tiles", "model"], help="difference with --psf-tile: PSF per tile from stacked stars (default) or one spatially varying polynomial PSF model per chip (ogfkit.psfmodel)")
    ap.add_argument("--template-psf", default="target", choices=["target", "measure"], help="difference: PSF of the template (default: the target PSF)")
    ap.add_argument("--tol", type=float, default=0.5)
    ap.add_argument("--min-exposures", type=int, default=3)
    ap.add_argument("--orbit-prior", default=None, help="link: vet tracklets against the apparent-motion distribution of a bound-orbit population: element table .npz or 'auto' (SBDB, cached)")
    ap.add_argument("--orbit-prior-fraction", type=float, default=0.99, help="link: keep tracklets inside this highest-density fraction of the population (default 0.99)")
    ap.add_argument("--max-per-exposure", type=int, default=900)
    ap.add_argument("--max-tracklets", type=int, default=400)
    ap.add_argument("--tracklet", type=int)
    ap.add_argument("--tracklet-files", nargs="*", dest="tracklet_files", help="nightlink: tracklets.json files of the different nights/visits")
    ap.add_argument("--obs-code", default="500", help="nightlink: MPC observatory code of the tracklets (500 geocentre, 250 HST, ...)")
    ap.add_argument("--nl-chi2", type=float, default=4.0); ap.add_argument("--nl-floor", type=float, default=0.3, help="nightlink model-error floor [arcsec]")
    ap.add_argument("--nl-nbody", action="store_true", help="nightlink: re-fit each linked group in a separate process (ASSIST differential correction, else the CODES integrator)")
    ap.add_argument("--nl-ranging", action="store_true", help="nightlink: sample the orbital-element distribution of each linked group and pair")
    ap.add_argument("--nl-ranging-samples", type=int, default=400, help="nightlink --nl-ranging: admissible-region samples per group")
    ap.add_argument("--nl-ranging-propagate", type=int, default=32, help="nightlink --nl-ranging: accepted samples to propagate with CODES (0 skips)")
    ap.add_argument("--nl-max-gap", type=float, default=30.0); ap.add_argument("--nl-min", type=int, default=3, help="nightlink: tracklets needed for a confirmed group (2-tracklet links are reported as pairs)"); ap.add_argument("--nl-sigma", type=float, default=0.3)
    ap.add_argument("--designation"); ap.add_argument("--mjd-min", type=float); ap.add_argument("--mjd-max", type=float)
    ap.add_argument("--no-debias", action="store_true"); ap.add_argument("--samples", type=int, default=400)
    ap.add_argument("--catalog"); ap.add_argument("--min-epochs", type=int, default=2)
    ap.add_argument("--aperture", type=float, default=3.0)
    ap.add_argument("--orbit-json", default="")
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--include-cr-flagged", action="store_true",
                    help="transients: also group detections flagged by the archive CR mask (a CR cannot repeat at the same place in >=2 exposures)")
    a = ap.parse_args(argv)
    t0 = time.time()
    globals()["m_" + a.mode](a)
    util.log("[%s] done in %.1f s" % (a.mode, time.time() - t0))


if __name__ == "__main__":
    main()
