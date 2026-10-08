"""Likelihood contours of the Apollo Laplace sample, with the JPL state at that epoch.

The draws are the local Gaussian approximation of the likelihood.  In each
plane a Gaussian KDE of those draws is the density, and the contour lines
enclose 68% and 95% of the draws.  The comparison point is the Horizons
heliocentric ICRF state at the sample's UTC epoch, converted to ecliptic
elements with the same ``kepler.state_to_elements`` the sample uses.
"""
import json
import os
import sys

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "2"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))

import run
from moving import arc_ranging as AR
from moving import horizons
from moving import kepler as K
from moving import obs as O

PNG_PATH = os.path.join(HERE, "likelihood_contours.png")
JSON_PATH = os.path.join(HERE, "likelihood_comparison.json")
FRACTIONS = (0.68, 0.95)
PLANES = (
    ("a", "e", "a (AU)", "e"),
    ("a", "inc_deg", "a (AU)", "i (deg)"),
    ("q", "e", "q (AU)", "e"),
)


def _elements_of(states):
    rows = []
    for state in states:
        el = K.state_to_elements(state[:3], state[3:])
        rows.append({
            "a": float(el["a"]),
            "e": float(el["e"]),
            "inc_deg": float(el["i"]),
            "q": float(el["q"]),
            "om_deg": float(el["om"]),
            "w_deg": float(el["w"]),
            "ma_deg": float(el["ma"]),
        })
    return rows


def _hpd_levels(density):
    """Density thresholds that keep the densest ``fraction`` of the draws."""
    order = np.argsort(density)[::-1]
    ranked = density[order]
    cum = np.arange(1, len(ranked) + 1) / float(len(ranked))
    levels = []
    for fraction in FRACTIONS:
        k = int(np.searchsorted(cum, fraction, side="left"))
        k = min(k, len(ranked) - 1)
        levels.append(float(ranked[k]))
    return levels


def _plane(ax, x, y, xlabel, ylabel, actual_xy, mode_xy):
    samples = np.vstack([x, y])
    kde = gaussian_kde(samples)
    sample_density = kde(samples)
    levels = _hpd_levels(sample_density)
    pad_x = 0.35 * (np.ptp(x) + 1e-6)
    pad_y = 0.35 * (np.ptp(y) + 1e-6)
    lo_x = min(np.min(x), actual_xy[0], mode_xy[0]) - pad_x
    hi_x = max(np.max(x), actual_xy[0], mode_xy[0]) + pad_x
    lo_y = min(np.min(y), actual_xy[1], mode_xy[1]) - pad_y
    hi_y = max(np.max(y), actual_xy[1], mode_xy[1]) + pad_y
    gx = np.linspace(lo_x, hi_x, 160)
    gy = np.linspace(lo_y, hi_y, 160)
    mesh_x, mesh_y = np.meshgrid(gx, gy)
    grid = kde(np.vstack([mesh_x.ravel(), mesh_y.ravel()])).reshape(mesh_x.shape)
    peak = float(max(np.max(grid), np.max(sample_density)))
    ax.contourf(mesh_x, mesh_y, grid, levels=np.linspace(0.0, peak, 8)[1:], cmap="Blues", alpha=0.85)
    # sorted: lower density is the outer 95% line, higher density is the inner 68% line
    ax.contour(mesh_x, mesh_y, grid, levels=sorted(levels),
               colors=["#7aa0c4", "#16324f"], linewidths=[1.15, 1.55], linestyles=["--", "-"])
    ax.scatter(x, y, s=11, c="#1f4e79", alpha=0.45, linewidths=0, zorder=3)
    ax.scatter([mode_xy[0]], [mode_xy[1]], marker="+", s=90, c="black", linewidths=1.4, zorder=4)
    ax.scatter([actual_xy[0]], [actual_xy[1]], marker="*", s=130, c="#D97757", zorder=5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(lo_x, hi_x)
    ax.set_ylim(lo_y, hi_y)
    actual_density = float(kde(np.array([[actual_xy[0]], [actual_xy[1]]]))[0])
    return {
        "n": int(len(x)),
        "kde_covariance": [[float(v) for v in row] for row in np.asarray(kde.covariance)],
        "peak_density": peak,
        "level_68": levels[0],
        "level_95": levels[1],
        "level_68_over_peak": levels[0] / peak,
        "level_95_over_peak": levels[1] / peak,
        "density_at_horizons": actual_density,
        "horizons_inside_68": bool(actual_density + 1e-15 >= levels[0]),
        "horizons_inside_95": bool(actual_density + 1e-15 >= levels[1]),
    }


def main():
    parsed = [o for o in O.parse_obs80(open(run.OBS_PATH).read()) if o["stn"] == run.STATION and o["note2"] not in ("S", "s")]
    trks, _made = run._tracklets(parsed)
    summary, pack = AR.distribution(trks, n_samples=run.N_SAMPLES, seed=run.SEED)
    if summary.get("sampler") != "laplace":
        raise SystemExit("expected the laplace sampler, got %s" % summary.get("sampler"))
    saved = json.load(open(os.path.join(HERE, "result.json")))
    saved_a = saved["ranging"]["quantiles_16_50_84"]["a"]
    fresh_a = summary["quantiles_16_50_84"]["a"]
    if max(abs(float(a) - float(b)) for a, b in zip(saved_a, fresh_a)) > 1e-9:
        raise SystemExit("Laplace quantiles do not match result.json")
    rows = _elements_of(pack["state"])
    from astropy.time import Time
    epoch = Time(float(summary["mjd_epoch"]), format="mjd", scale="utc")
    jd_tdb = float(epoch.tdb.jd)
    state = horizons.vectors("1862;", [jd_tdb], center="500@10")[0]
    actual = _elements_of(state.reshape(1, 6))[0]
    mode = {
        "a": float(summary["best"]["a"]),
        "e": float(summary["best"]["e"]),
        "inc_deg": float(summary["best"]["inc_deg"]),
        "q": float(summary["best"]["q"]),
    }
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.55))
    planes = {}
    for ax, (xk, yk, xlabel, ylabel) in zip(axes, PLANES):
        planes["%s_%s" % (xk, yk)] = _plane(
            ax,
            np.array([r[xk] for r in rows]),
            np.array([r[yk] for r in rows]),
            xlabel, ylabel,
            (actual[xk], actual[yk]),
            (mode[xk], mode[yk]),
        )
    from matplotlib.lines import Line2D
    fig.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#1f4e79", markersize=5, label="Laplace draws"),
            Line2D([0], [0], color="#16324f", lw=1.55, label="68% of draws"),
            Line2D([0], [0], color="#7aa0c4", lw=1.15, ls="--", label="95% of draws"),
            Line2D([0], [0], marker="+", color="black", lw=0, markersize=8, label="Laplace mode"),
            Line2D([0], [0], marker="*", color="none", markerfacecolor="#D97757", markersize=11, label="Horizons"),
        ],
        loc="lower center", ncol=5, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.02),
    )
    fig.suptitle("Apollo 1862, local likelihood sample at MJD %.6f UTC" % float(summary["mjd_epoch"]), fontsize=11)
    fig.tight_layout(rect=(0.0, 0.06, 1.0, 0.96))
    fig.savefig(PNG_PATH, dpi=160)
    plt.close(fig)
    offset = {key: float(actual[key] - mode[key]) for key in ("a", "e", "inc_deg", "q")}
    doc = {
        "observations_file": "observations.obs80",
        "result_file": "result.json",
        "figure": "likelihood_contours.png",
        "n_samples": int(summary["n_samples"]),
        "seed": int(run.SEED),
        "sampler": summary["sampler"],
        "sample_mjd_utc": float(summary["mjd_epoch"]),
        "sample_utc": epoch.utc.isot,
        "horizons": {
            "command": "1862;",
            "center": "500@10",
            "jd_tdb": jd_tdb,
            "tdb": epoch.tdb.isot,
            "frame": (
                "Horizons VECTORS, heliocentric ICRF equatorial. Elements from "
                "moving.kepler.state_to_elements, which rotates equatorial vectors "
                "onto the ecliptic with the IAU 1976 obliquity."
            ),
            "elements": actual,
        },
        "laplace_mode": mode,
        "offset_horizons_minus_mode": offset,
        "planes": planes,
        "note": (
            "Contours are a Gaussian KDE of the Laplace draws in that plane. "
            "The 68% line is the density of the draw at which 68% of the draws, "
            "counted from the densest, are inside. Same rule for 95%. "
            "The Horizons elements are osculating at the sample UTC instant, "
            "expressed as a TDB Julian date for the vector request. "
            "They are not the SBDB elements at JD 2461200.5."
        ),
    }
    with open(JSON_PATH, "w") as handle:
        json.dump(doc, handle, indent=2)
        handle.write("\n")
    print(
        "CONTOURS a=%.5f e=%.5f i=%.5f q=%.5f offset_a=%+.5f offset_e=%+.5f offset_i=%+.5f offset_q=%+.5f"
        % (actual["a"], actual["e"], actual["inc_deg"], actual["q"],
           offset["a"], offset["e"], offset["inc_deg"], offset["q"])
    )
    for name, plane in planes.items():
        print("PLANE %s inside68=%s inside95=%s" % (name, plane["horizons_inside_68"], plane["horizons_inside_95"]))


if __name__ == "__main__":
    main()
