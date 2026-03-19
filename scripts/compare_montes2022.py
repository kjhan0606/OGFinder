#!/usr/bin/env python3
"""
Compare OGFinder ICL pipeline results with Montes & Trujillo (2022, ApJ 940 L51).

Generates publication-quality figures overlaying our measurements with
the Montes+2022 reference model (BCG + ICL power-law):
  - fig_sb_profile.png   : SB profile overlaid with M+T22 model
  - fig_icl_fraction.png : f_ICL multi-threshold sweep
  - fig_decompose.png    : BCG+ICL double Sersic decomposition
  - fig_color_profile.png: Color gradient (if available)

Usage: cd OGFinder && python3 scripts/compare_montes2022.py
"""

import os
import sys
import numpy as np

RESULTS_DIR = "results/montes2022"

# ---- Reference values from Montes & Trujillo 2022 (ApJ 940 L51) ----

# Key measured quantities
PAPER_ALPHA_3D = -2.47       # stellar mass density 3D slope
PAPER_ALPHA_3D_ERR = 0.13
PAPER_ALPHA_2D_E = -1.48     # projected 2D slope (east)
PAPER_ALPHA_2D_E_ERR = 0.13
PAPER_ALPHA_2D_W = -1.46     # projected 2D slope (west)
PAPER_ALPHA_2D_W_ERR = 0.08
PAPER_RADIAL_EXTENT_KPC = 400  # ICL traced to ~400 kpc
PAPER_SB_LIMIT = 28.0         # mag/arcsec^2 conservative limit
PAPER_PIXEL_SCALE = 0.063     # arcsec/pixel (NIRCam LW)

# SB limits (3-sigma, 10x10 arcsec^2) from the paper
PAPER_SB_3SIGMA = {'F277W': 31.28, 'F356W': 31.32, 'F444W': 31.10}

# Ellien 2025 (A&A 698 A45) - same cluster, DAWIS method
ELLIEN_F_ICL = 0.28           # within 400 kpc
ELLIEN_F_ICL_ERR = 0.02
ELLIEN_F_ICL_BCG = 0.34
ELLIEN_F_ICL_BCG_ERR = 0.02

# Cosmology for SMACS J0723.3-7327 (z=0.39)
CLUSTER_Z = 0.39
KPC_PER_ARCSEC = 5.277        # at z=0.39 (flat LCDM, H0=70)

# Tidal feature locations (approximate, from paper text)
TIDAL_FEATURES = {
    'East stream 1': 150,  # kpc
    'East stream 2': 243,  # kpc
    'West loop': 200,      # kpc (center of ~150-250 kpc feature)
}


def montes2022_model_profile(r_arcsec):
    """Generate approximate Montes+2022 SB profile model.

    Constructs a two-component (BCG + ICL) model based on:
    - BCG: Sersic n=4 (de Vaucouleurs), re~3"
    - ICL: power-law with alpha_proj = -1.47, reaching mu~28 at r=72"

    These components are added in flux space to produce the combined profile.
    The model approximates the data shown in Montes+2022 Figure 2.
    """
    r = np.asarray(r_arcsec, dtype=float)
    r = np.clip(r, 0.1, None)

    # BCG component: Sersic n=4
    n_bcg = 4.0
    re_bcg = 3.0  # arcsec
    mu_e_bcg = 21.5  # mag/arcsec^2 at re
    bn = 2.0 * n_bcg - 1.0 / 3.0 + 4.0 / (405.0 * n_bcg)
    mu_bcg = mu_e_bcg + 2.5 * bn / np.log(10.0) * ((r / re_bcg) ** (1.0 / n_bcg) - 1.0)

    # ICL component: power-law alpha_proj = alpha_3D + 1 = -1.47
    # Anchored at r=72" -> mu=28
    alpha_proj = PAPER_ALPHA_3D + 1.0  # -1.47
    r_anchor = 72.0  # arcsec
    mu_anchor = 28.0  # mag/arcsec^2
    # mu(r) = mu_anchor - 2.5 * alpha_proj * log10(r / r_anchor)
    mu_icl = mu_anchor - 2.5 * alpha_proj * np.log10(r / r_anchor)

    # Combine in flux space
    f_bcg = 10.0 ** (-0.4 * mu_bcg)
    f_icl = 10.0 ** (-0.4 * mu_icl)
    mu_total = -2.5 * np.log10(f_bcg + f_icl)

    return mu_total, mu_bcg, mu_icl


def load_profile(path):
    """Load SB profile TSV (R_PIX R_ARCSEC MU MU_ERR FLUX AREA NPIX)."""
    r_pix, r_arcsec, mu, mu_err = [], [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('R_PIX'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            try:
                r_pix.append(float(parts[0]))
                r_arcsec.append(float(parts[1]))
                mu.append(float(parts[2]))
                mu_err.append(float(parts[3]))
            except ValueError:
                continue
    return {
        'r_pix': np.array(r_pix),
        'r_arcsec': np.array(r_arcsec),
        'r_kpc': np.array(r_arcsec) * KPC_PER_ARCSEC,
        'mu': np.array(mu),
        'mu_err': np.array(mu_err),
    }


def load_icl_fraction(path):
    """Load multi-threshold ICL fraction TSV."""
    mu_thresh, f_icl, f_icl_err = [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('MU'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            try:
                mu_thresh.append(float(parts[0]))
                f_icl.append(float(parts[1]))
                f_icl_err.append(float(parts[2]))
            except ValueError:
                continue
    return {
        'mu_thresh': np.array(mu_thresh),
        'f_icl': np.array(f_icl),
        'f_icl_err': np.array(f_icl_err),
    }


def load_decompose(path):
    """Load decomposition result TSV."""
    result = {}
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
    if len(lines) < 2:
        return None
    keys = lines[0].split('\t')
    vals = lines[1].split('\t')
    for k, v in zip(keys, vals):
        try:
            result[k] = float(v)
        except ValueError:
            result[k] = v
    return result


def sersic_sb(r, mu_e, r_e, n):
    """Sersic surface brightness profile."""
    bn = 2.0 * n - 1.0 / 3.0 + 4.0 / (405.0 * n)
    return mu_e + 2.5 * bn / np.log(10.0) * ((r / r_e) ** (1.0 / n) - 1.0)


def plot_sb_profile(profile, profile_east=None, profile_west=None):
    """Figure 1: SB profile overlaid with Montes+2022 model (cf. Fig 2)."""
    import matplotlib.pyplot as plt

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # ---- Left panel: R in kpc ----
    valid = profile['mu'] < 99.0
    r_kpc = profile['r_kpc'][valid]
    r_arcsec = profile['r_arcsec'][valid]
    mu = profile['mu'][valid]
    mu_err = profile['mu_err'][valid]

    ax.errorbar(r_kpc, mu, yerr=mu_err, fmt='o', ms=4, color='navy',
                ecolor='lightblue', elinewidth=1, capsize=2,
                label='This work (circular)', zorder=3)

    # Sector profiles
    if profile_east is not None:
        ve = profile_east['mu'] < 99.0
        ax.plot(profile_east['r_kpc'][ve], profile_east['mu'][ve],
                's', ms=3, color='crimson', alpha=0.7, label='This work (East)')
    if profile_west is not None:
        vw = profile_west['mu'] < 99.0
        ax.plot(profile_west['r_kpc'][vw], profile_west['mu'][vw],
                '^', ms=3, color='forestgreen', alpha=0.7, label='This work (West)')

    # Montes+2022 model profile (BCG + ICL power-law)
    r_model_arcsec = np.logspace(np.log10(0.2), np.log10(80), 300)
    r_model_kpc = r_model_arcsec * KPC_PER_ARCSEC
    mu_model, mu_bcg, mu_icl = montes2022_model_profile(r_model_arcsec)

    ax.plot(r_model_kpc, mu_model, '-', color='red', linewidth=2, alpha=0.8,
            label='M+T22 model (BCG+ICL)', zorder=2)
    ax.plot(r_model_kpc, mu_bcg, ':', color='blue', linewidth=1, alpha=0.5,
            label=r'M+T22 BCG (S\'ersic n=4)')
    ax.plot(r_model_kpc, mu_icl, ':', color='orange', linewidth=1, alpha=0.5,
            label=r'M+T22 ICL ($\alpha_{3D}$=' + f'{PAPER_ALPHA_3D})')

    # Power-law slope reference anchored to data
    r_ref_kpc = 50.0
    idx_ref = np.argmin(np.abs(r_kpc - r_ref_kpc))
    mu_ref = mu[idx_ref]
    r_pw = np.linspace(20, 400, 200)
    alpha_proj = PAPER_ALPHA_3D + 1.0  # -1.47
    # FIXED: mu gets fainter (larger) at larger r
    mu_pw = mu_ref - 2.5 * alpha_proj * np.log10(r_pw / r_ref_kpc)
    ax.plot(r_pw, mu_pw, '--', color='red', linewidth=1.2, alpha=0.6,
            label=fr'$\alpha_{{3D}}={PAPER_ALPHA_3D}\pm{PAPER_ALPHA_3D_ERR}$ (anchored)')

    # Error band for alpha_3D uncertainty
    alpha_lo = (PAPER_ALPHA_3D - PAPER_ALPHA_3D_ERR) + 1.0
    alpha_hi = (PAPER_ALPHA_3D + PAPER_ALPHA_3D_ERR) + 1.0
    mu_pw_lo = mu_ref - 2.5 * alpha_lo * np.log10(r_pw / r_ref_kpc)
    mu_pw_hi = mu_ref - 2.5 * alpha_hi * np.log10(r_pw / r_ref_kpc)
    ax.fill_between(r_pw, mu_pw_lo, mu_pw_hi, color='red', alpha=0.08)

    # Tidal feature markers
    for name, r_feat in TIDAL_FEATURES.items():
        ax.axvline(r_feat, color='purple', linestyle=':', linewidth=0.7,
                   alpha=0.4)
        ax.text(r_feat, 17.5, name, fontsize=6, rotation=90,
                va='bottom', ha='right', color='purple', alpha=0.6)

    # Limits
    ax.axhline(PAPER_SB_LIMIT, color='gray', linestyle=':', linewidth=1,
               label=f'Conservative limit ({PAPER_SB_LIMIT} mag/"$^2$)')
    ax.axvline(PAPER_RADIAL_EXTENT_KPC, color='gray', linestyle='--',
               linewidth=0.8, alpha=0.5)
    ax.text(PAPER_RADIAL_EXTENT_KPC, 28.5, f'{PAPER_RADIAL_EXTENT_KPC} kpc',
            fontsize=8, ha='center', color='gray')

    ax.set_xlabel('R (kpc)', fontsize=13)
    ax.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=13)
    ax.set_title('SMACS J0723 ICL - SB Profile (F277W)\n'
                 'cf. Montes & Trujillo 2022 Fig.2', fontsize=12)
    ax.set_xscale('log')
    ax.invert_yaxis()
    ax.set_xlim(0.5, 500)
    ax.set_ylim(29, 17)
    ax.legend(fontsize=7, loc='upper right', ncol=1)
    ax.grid(True, alpha=0.3)

    # ---- Right panel: R in arcsec (same as paper's x-axis) ----
    ax2.errorbar(r_arcsec, mu, yerr=mu_err, fmt='o', ms=4, color='navy',
                 ecolor='lightblue', elinewidth=1, capsize=2,
                 label='This work', zorder=3)

    ax2.plot(r_model_arcsec, mu_model, '-', color='red', linewidth=2,
             alpha=0.8, label='M+T22 model (BCG+ICL)', zorder=2)
    ax2.plot(r_model_arcsec, mu_bcg, ':', color='blue', linewidth=1,
             alpha=0.5, label=r'BCG (S\'ersic n=4, $r_e$=3")')
    ax2.plot(r_model_arcsec, mu_icl, ':', color='orange', linewidth=1,
             alpha=0.5, label=r'ICL ($\alpha_{proj}=-1.47$)')

    ax2.axhline(PAPER_SB_LIMIT, color='gray', linestyle=':', linewidth=1)
    ax2.axvline(72, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax2.text(72, 28.5, '72" (380 kpc)', fontsize=8, ha='center', color='gray')

    ax2.set_xlabel('R (arcsec)', fontsize=13)
    ax2.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=13)
    ax2.set_title('Same data, arcsec scale\n'
                  '(Montes+2022 Fig.2 x-axis)', fontsize=12)
    ax2.set_xscale('log')
    ax2.invert_yaxis()
    ax2.set_xlim(0.05, 100)
    ax2.set_ylim(29, 17)
    ax2.legend(fontsize=8, loc='upper right')
    ax2.grid(True, alpha=0.3)

    out = os.path.join(RESULTS_DIR, 'fig_sb_profile.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

    # Fit and report slope
    mask_fit = (r_kpc > 30) & (r_kpc < 300) & (mu < 99)
    if np.sum(mask_fit) >= 3:
        log_r = np.log10(r_kpc[mask_fit])
        mu_fit = mu[mask_fit]
        coeffs = np.polyfit(log_r, mu_fit, 1)
        alpha_measured = coeffs[0] / 2.5
        print(f"  Measured projected slope (30-300 kpc): "
              f"alpha_proj = {alpha_measured:.3f}")
        print(f"  => alpha_3D ~ {alpha_measured - 1.0:.3f} "
              f"(paper: {PAPER_ALPHA_3D} +/- {PAPER_ALPHA_3D_ERR})")


def plot_icl_fraction(frac_data):
    """Figure 2: ICL fraction vs SB threshold."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))

    ax.errorbar(frac_data['mu_thresh'], frac_data['f_icl'] * 100,
                yerr=frac_data['f_icl_err'] * 100,
                fmt='o-', color='navy', ms=6, capsize=3, linewidth=1.5,
                label=r'This work: $f_{ICL}(\mu > \mu_{thresh})$')

    # Ellien 2025 reference (aperture-based, ~26.5 threshold equivalent)
    ax.axhline(ELLIEN_F_ICL * 100, color='red', linestyle='--', linewidth=1.5,
               label=f'Ellien+25 (DAWIS): {ELLIEN_F_ICL*100:.0f}%'
                     f' $\\pm$ {ELLIEN_F_ICL_ERR*100:.0f}%')
    ax.axhspan((ELLIEN_F_ICL - ELLIEN_F_ICL_ERR) * 100,
               (ELLIEN_F_ICL + ELLIEN_F_ICL_ERR) * 100,
               color='red', alpha=0.08)

    # Ellien+25 BCG+ICL
    ax.axhline(ELLIEN_F_ICL_BCG * 100, color='orange', linestyle='--',
               linewidth=1, alpha=0.7,
               label=f'Ellien+25 (BCG+ICL): {ELLIEN_F_ICL_BCG*100:.0f}%')

    ax.set_xlabel(r'$\mu_{threshold}$ (mag arcsec$^{-2}$)', fontsize=13)
    ax.set_ylabel(r'$f_{ICL}$ (%)', fontsize=13)
    ax.set_title('SMACS J0723 - ICL Fraction vs SB Threshold\n'
                 r'$f_{ICL} = L(\mu > \mu_{thresh}) / L_{total}$',
                 fontsize=12)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 100)

    # Annotate convention
    ax.text(0.03, 0.03,
            'Convention: ICL = light fainter than threshold\n'
            r'$f_{ICL}$ decreases as threshold gets fainter'
            '\n(fewer annuli qualify)',
            transform=ax.transAxes, fontsize=8, va='bottom',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow',
                      edgecolor='gray', alpha=0.7))

    out = os.path.join(RESULTS_DIR, 'fig_icl_fraction.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

    for mu_t in [25.0, 26.0, 26.5, 27.0]:
        idx = np.argmin(np.abs(frac_data['mu_thresh'] - mu_t))
        if abs(frac_data['mu_thresh'][idx] - mu_t) < 0.3:
            print(f"  f_ICL(mu>{mu_t}) = "
                  f"{frac_data['f_icl'][idx]*100:.1f} +/- "
                  f"{frac_data['f_icl_err'][idx]*100:.1f} %")


def plot_decompose(profile, decomp):
    """Figure 3: BCG+ICL double Sersic decomposition + M+T22 model."""
    import matplotlib.pyplot as plt

    if decomp is None:
        print("  Decomposition data not available, skipping figure.")
        return

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    # Data
    valid = profile['mu'] < 99.0
    r = profile['r_arcsec'][valid]
    mu = profile['mu'][valid]
    mu_err = profile['mu_err'][valid]

    ax.errorbar(r, mu, yerr=mu_err, fmt='o', ms=4, color='black',
                ecolor='gray', elinewidth=1, capsize=2, label='Data', zorder=5)

    # Our Sersic components
    r_model = np.logspace(np.log10(0.1), np.log10(80), 300)

    n_bcg = decomp.get('n_bcg', 4.0)
    re_bcg = decomp.get('re_bcg_arcsec', 5.0)
    n_icl = decomp.get('n_icl', 1.0)
    re_icl = decomp.get('re_icl_arcsec', 30.0)

    Ie_bcg = decomp.get('Ie_bcg', 0.1)
    Ie_icl = decomp.get('Ie_icl', 0.01)
    zp_sb = 20.47
    mu_e_bcg = -2.5 * np.log10(max(Ie_bcg, 1e-30)) + zp_sb
    mu_e_icl = -2.5 * np.log10(max(Ie_icl, 1e-30)) + zp_sb

    mu_bcg = sersic_sb(r_model, mu_e_bcg, re_bcg, n_bcg)
    mu_icl = sersic_sb(r_model, mu_e_icl, re_icl, n_icl)

    f_bcg = 10 ** (-0.4 * mu_bcg)
    f_icl_fit = 10 ** (-0.4 * mu_icl)
    mu_total = -2.5 * np.log10(f_bcg + f_icl_fit)

    ax.plot(r_model, mu_bcg, '--', color='blue', linewidth=1.5,
            label=f'BCG (n={n_bcg:.1f}, $r_e$={re_bcg:.1f}")')
    ax.plot(r_model, mu_icl, '--', color='red', linewidth=1.5,
            label=f'ICL (n={n_icl:.1f}, $r_e$={re_icl:.1f}")')
    ax.plot(r_model, mu_total, '-', color='orange', linewidth=2,
            label='BCG+ICL total')

    # Montes+2022 model for comparison
    mu_mt22, _, _ = montes2022_model_profile(r_model)
    ax.plot(r_model, mu_mt22, '-', color='green', linewidth=1.5, alpha=0.5,
            label='M+T22 reference model')

    # Transition radius
    r_trans = decomp.get('R_transition', 0)
    if r_trans > 0:
        r_trans_arcsec = r_trans * PAPER_PIXEL_SCALE
        ax.axvline(r_trans_arcsec, color='gray', linestyle=':', linewidth=1,
                   label=f'$R_{{trans}}$ = {r_trans_arcsec:.1f}"')

    chi2 = decomp.get('chi2', float('nan'))
    f_icl_sersic = decomp.get('f_ICL_sersic', float('nan'))

    ax.set_xlabel('R (arcsec)', fontsize=13)
    ax.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=13)
    ax.set_title(f'BCG+ICL Decomposition '
                 f'($f_{{ICL}}$={f_icl_sersic:.1%}, $\\chi^2$={chi2:.2f})',
                 fontsize=12)
    ax.set_xscale('log')
    ax.invert_yaxis()
    ax.set_xlim(0.05, 100)
    ax.set_ylim(30, 16)
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3)

    out = os.path.join(RESULTS_DIR, 'fig_decompose.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

    print(f"  BCG: n={n_bcg:.2f}, re={re_bcg:.2f}\"")
    print(f"  ICL: n={n_icl:.2f}, re={re_icl:.2f}\"")
    print(f"  f_ICL(Sersic) = {f_icl_sersic:.4f}")


def plot_color_profile():
    """Figure 4: Color gradient (F090W - F277W) if available."""
    import matplotlib.pyplot as plt

    color_file = os.path.join(RESULTS_DIR, 'color_profile.tsv')
    if not os.path.exists(color_file):
        print("  Color profile not available (Step 7 was skipped).")
        return

    r_arcsec, color_vals = [], []
    with open(color_file) as f:
        header = None
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if header is None:
                header = line.split('\t')
                continue
            parts = line.split('\t')
            if len(parts) < len(header):
                continue
            try:
                r_arcsec.append(float(parts[1]))
                color_vals.append(float(parts[-1]))
            except (ValueError, IndexError):
                continue

    if not r_arcsec:
        print("  No valid color data found.")
        return

    r_arcsec = np.array(r_arcsec)
    color_vals = np.array(color_vals)
    r_kpc = r_arcsec * KPC_PER_ARCSEC
    valid = np.isfinite(color_vals)

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    ax.plot(r_kpc[valid], color_vals[valid], 'o-', color='purple', ms=4,
            label='This work')

    # Montes+2022: negative color gradient at R > 100 kpc
    ax.axvline(100, color='gray', linestyle=':', linewidth=0.8, alpha=0.5)
    ax.text(105, ax.get_ylim()[1] * 0.9 if ax.get_ylim()[1] > 0 else 1,
            'M+T22: gradient\nchange ~100 kpc',
            fontsize=8, color='gray')

    ax.set_xlabel('R (kpc)', fontsize=13)
    ax.set_ylabel('F090W $-$ F277W (mag)', fontsize=13)
    ax.set_title('SMACS J0723 - ICL Color Gradient\n'
                 'cf. Montes+2022 Fig.3', fontsize=12)
    ax.set_xscale('log')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    out = os.path.join(RESULTS_DIR, 'fig_color_profile.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


def main():
    if not os.path.isdir(RESULTS_DIR):
        print(f"ERROR: Results directory not found: {RESULTS_DIR}")
        print("Run scripts/reproduce_montes2022.sh first.")
        sys.exit(1)

    print("=" * 60)
    print("Montes & Trujillo (2022) ICL Comparison")
    print("  Data: SMACS J0723.3-7327 (z=0.39)")
    print("  Filter: F277W (NIRCam LW, Module B)")
    print("  Pixel scale: 0.063\"/pixel")
    print("=" * 60)

    # Load data
    sb_file = os.path.join(RESULTS_DIR, 'sb_profile.tsv')
    if not os.path.exists(sb_file):
        print(f"ERROR: SB profile not found: {sb_file}")
        sys.exit(1)

    profile = load_profile(sb_file)
    print(f"\nLoaded SB profile: {len(profile['r_pix'])} points, "
          f"R = {profile['r_arcsec'][0]:.2f} - "
          f"{profile['r_arcsec'][-1]:.2f} arcsec "
          f"({profile['r_kpc'][0]:.1f} - {profile['r_kpc'][-1]:.1f} kpc)")

    # Sector profiles
    east_file = os.path.join(RESULTS_DIR, 'sb_profile_east.tsv')
    west_file = os.path.join(RESULTS_DIR, 'sb_profile_west.tsv')
    profile_east = load_profile(east_file) if os.path.exists(east_file) else None
    profile_west = load_profile(west_file) if os.path.exists(west_file) else None

    # Figure 1: SB profile with Montes+2022 overlay
    print("\n--- SB Profile (overlaid with M+T22 model) ---")
    plot_sb_profile(profile, profile_east, profile_west)

    # Figure 2: ICL fraction
    frac_file = os.path.join(RESULTS_DIR, 'icl_fraction_multi.tsv')
    if os.path.exists(frac_file):
        print("\n--- ICL Fraction ---")
        frac_data = load_icl_fraction(frac_file)
        if len(frac_data['mu_thresh']) > 0:
            plot_icl_fraction(frac_data)
        else:
            print("  No valid f_ICL data.")
    else:
        print(f"\n  f_ICL file not found: {frac_file}")

    # Figure 3: Decomposition
    decomp_file = os.path.join(RESULTS_DIR, 'decompose_result.tsv')
    if os.path.exists(decomp_file):
        print("\n--- BCG+ICL Decomposition ---")
        decomp = load_decompose(decomp_file)
        plot_decompose(profile, decomp)
    else:
        print(f"\n  Decomposition file not found: {decomp_file}")

    # Figure 4: Color profile
    print("\n--- Color Profile ---")
    plot_color_profile()

    print("\n" + "=" * 60)
    print("Done. Figures saved to", RESULTS_DIR)


if __name__ == '__main__':
    main()
