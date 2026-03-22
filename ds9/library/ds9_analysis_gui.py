#!/usr/bin/env python3
"""DS9 Analysis Viewer - Standalone multi-analysis GUI.

Provides an integrated analysis environment with:
- Left panel: matplotlib plots (SED, photo-z, morphology diagrams)
  + source cutout (FITS image)
- Right panel: sortable catalog table
- Analysis runners: Photo-z, SED Fitting, Morphology, Star/Galaxy,
  CAS/Gini/M20, Sérsic, Bulge+Disk, PSF Phot, Multi-Band,
  Cross-Match, Completeness
- Settings save/load, region export, FITS table export

Usage:
    ds9_analysis_gui.py [--fits image.fits] [--catalog catalog.tsv]
"""

import sys
import os
import subprocess
import argparse
import threading
import numpy as np

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.patches import Ellipse
from matplotlib.gridspec import GridSpec

# Try astropy for FITS I/O
try:
    from astropy.io import fits as pyfits
    HAS_ASTROPY = True
except ImportError:
    HAS_ASTROPY = False

# Project root
_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..', '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

# Filter wavelengths (Angstrom) for SED plotting
FILTER_WAVELENGTHS = {
    "F435W": 4328, "F475W": 4746, "F555W": 5346, "F606W": 5907,
    "F625W": 6311, "F775W": 7693, "F814W": 8057, "F850LP": 9033,
    "F098M": 9864, "F105W": 10552, "F110W": 11534, "F125W": 12486,
    "F140W": 13923, "F160W": 15369,
    "F225W": 2372, "F275W": 2710, "F336W": 3355, "F390W": 3923,
    "F070W": 7088, "F090W": 9010, "F115W": 11540, "F150W": 15010,
    "F200W": 19890, "F277W": 27620, "F356W": 35680, "F444W": 44050,
    "u": 3551, "g": 4686, "r": 6166, "i": 7480, "z": 8932,
}

MORPH_NAMES = [
    "Disturbed", "Merging", "Round Smooth", "In-between Smooth",
    "Cigar Smooth", "Barred Spiral", "Unbarred Tight Spiral",
    "Unbarred Loose Spiral", "Edge-on (no bulge)", "Edge-on (bulge)",
]

MODES = [
    ("overview",  "Overview"),
    ("photoz",    "Photo-z"),
    ("sed",       "SED Fitting"),
    ("classify",  "Classification"),
    ("structure", "Structure"),
    ("plot",      "Custom Plot"),
]

PLOT_PRESETS = {
    "(custom)":    ("", ""),
    "SFR vs Mass": ("LOG_MASS", "SFR"),
    "Mass vs z":   ("PHOTO_Z", "LOG_MASS"),
    "C-A":         ("ASYM", "CONC"),
    "Gini-M20":    ("M20", "GINI"),
    "n vs Re":     ("SERSIC_RE", "SERSIC_N"),
    "Mag vs z":    ("PHOTO_Z", "MAG_AUTO"),
    "B/T vs Mass": ("LOG_MASS", "BD_BULGE_FRAC"),
}

DEFAULT_SETTINGS = {
    'n_workers': 1,
    'cutout_size': 80,
    'pixel_scale': 0.0,
    'mag_zeropoint': 25.0,
    'sed_backend': 'auto',
    'sed_bands': 'g,r,i,z',
    'sed_mag_columns': '',
    'sed_photoz_column': 'PHOTO_Z',
    'sed_ckpt_emulator': '',
    'sed_ckpt_inverse': '',
    'photoz_bands': '',
    'photoz_mag_columns': '',
    'photoz_checkpoint': '',
}

_SETTINGS_FILE = os.path.join(os.path.expanduser("~"),
                               ".ds9", "analysis_gui.prf")


class AnalysisGUI:
    """Standalone analysis viewer with plots + cutout + catalog table."""

    def __init__(self, root, fits_path=None, catalog_path=None):
        self.root = root
        self.root.title("DS9 Analysis Viewer")
        self.root.geometry("1400x850")
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

        self.fits_path = fits_path
        self.catalog_path = catalog_path

        # Data
        self.header = []
        self.rows = []
        self.filtered_indices = []
        self.selected_idx = -1
        self.sort_col = None
        self.sort_asc = True
        self.fits_data = None       # cached 2D FITS array

        # State
        self.mode = "overview"

        # Settings
        self.settings = dict(DEFAULT_SETTINGS)
        self._load_settings()

        # Tk variables
        self.sed_backend = tk.StringVar(value=self.settings['sed_backend'])
        self.plot_xcol = tk.StringVar(value="MAG_AUTO")
        self.plot_ycol = tk.StringVar(value="FLUX_AUTO")
        self.plot_type = tk.StringVar(value="scatter")
        self.plot_preset = tk.StringVar(value="(custom)")

        self._build_ui()

        if catalog_path:
            self.load_catalog(catalog_path)

    # ================================================================
    # UI Construction
    # ================================================================

    def _build_ui(self):
        self._build_menu()
        self._build_toolbar()

        self.paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        self.paned.pack(fill=tk.BOTH, expand=True)

        self._build_plot_area()
        self._build_table_area()

        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(self.root, textvariable=self.status_var,
                  relief=tk.SUNKEN, anchor=tk.W).pack(
                      side=tk.BOTTOM, fill=tk.X)

    def _build_menu(self):
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        # File
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Open Catalog...",
                              command=self._open_catalog)
        file_menu.add_command(label="Open FITS...",
                              command=self._open_fits)
        file_menu.add_separator()
        file_menu.add_command(label="Export Catalog (TSV/CSV)...",
                              command=self._export_catalog)
        file_menu.add_command(label="Export Regions (.reg)...",
                              command=self._export_regions)
        if HAS_ASTROPY:
            file_menu.add_command(label="Export FITS Table...",
                                  command=self._export_fits_table)
        file_menu.add_separator()
        file_menu.add_command(label="Quit", command=self._on_close)
        menubar.add_cascade(label="File", menu=file_menu)

        # Analysis
        analysis_menu = tk.Menu(menubar, tearoff=0)
        analysis_menu.add_command(label="Photo-z (AI)",
                                  command=self.run_photoz)
        analysis_menu.add_command(label="SED Fitting",
                                  command=self.run_sed_fit)
        analysis_menu.add_separator()
        analysis_menu.add_command(label="Galaxy Morphology (CNN)",
                                  command=self.run_morphology)
        analysis_menu.add_command(label="Star/Galaxy (CNN)",
                                  command=self.run_starfinder)
        analysis_menu.add_separator()
        analysis_menu.add_command(label="CAS/Gini/M20",
                                  command=self.run_morphometry)
        analysis_menu.add_command(label="Sérsic Fitting",
                                  command=self.run_sersic)
        analysis_menu.add_command(label="Bulge+Disk Decomposition",
                                  command=self.run_bulge_disk)
        analysis_menu.add_separator()
        analysis_menu.add_command(label="PSF Photometry",
                                  command=self.run_psf_phot)
        analysis_menu.add_command(label="Multi-Band Photometry",
                                  command=self.run_multiband)
        analysis_menu.add_command(label="Crowded Field Photometry",
                                  command=self.run_crowded_phot)
        analysis_menu.add_separator()
        analysis_menu.add_command(label="Cross-Match (VizieR)",
                                  command=self.run_crossmatch)
        analysis_menu.add_command(label="Completeness Simulation",
                                  command=self.run_completeness)
        analysis_menu.add_separator()
        analysis_menu.add_command(label="Compute Color Columns",
                                  command=self._compute_colors)
        menubar.add_cascade(label="Analysis", menu=analysis_menu)

        # Settings
        settings_menu = tk.Menu(menubar, tearoff=0)
        settings_menu.add_command(label="General...",
                                  command=self._settings_general)
        settings_menu.add_command(label="SED Fitting...",
                                  command=self._settings_sed)
        settings_menu.add_command(label="Photo-z...",
                                  command=self._settings_photoz)
        menubar.add_cascade(label="Settings", menu=settings_menu)

        # View
        view_menu = tk.Menu(menubar, tearoff=0)
        for mode_id, mode_label in MODES:
            view_menu.add_command(
                label=mode_label,
                command=lambda m=mode_id: self._set_mode(m))
        menubar.add_cascade(label="View", menu=view_menu)

    def _build_toolbar(self):
        toolbar = ttk.Frame(self.root)
        toolbar.pack(side=tk.TOP, fill=tk.X, padx=2, pady=2)

        # Mode selector
        ttk.Label(toolbar, text="Mode:").pack(side=tk.LEFT, padx=(5, 2))
        self.mode_var = tk.StringVar(value="overview")
        mode_cb = ttk.Combobox(toolbar, textvariable=self.mode_var,
                               values=[m[0] for m in MODES],
                               state="readonly", width=12)
        mode_cb.pack(side=tk.LEFT, padx=2)
        mode_cb.bind("<<ComboboxSelected>>",
                      lambda e: self._set_mode(self.mode_var.get()))

        ttk.Separator(toolbar, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=5)

        # SED backend
        ttk.Label(toolbar, text="SPS:").pack(side=tk.LEFT, padx=2)
        ttk.Combobox(
            toolbar, textvariable=self.sed_backend, width=10,
            values=["auto", "analytic", "fsps", "bagpipes",
                    "prospector", "cigale", "dense_basis"],
            state="readonly",
        ).pack(side=tk.LEFT, padx=2)

        ttk.Separator(toolbar, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=5)

        # FITS path
        ttk.Label(toolbar, text="FITS:").pack(side=tk.LEFT, padx=2)
        self.fits_var = tk.StringVar(
            value=os.path.basename(self.fits_path) if self.fits_path else "")
        ttk.Label(toolbar, textvariable=self.fits_var, width=25,
                  relief=tk.SUNKEN).pack(side=tk.LEFT, padx=2)

        # Source count
        self.count_var = tk.StringVar(value="0 sources")
        ttk.Label(toolbar, textvariable=self.count_var).pack(
            side=tk.RIGHT, padx=10)

    def _build_plot_area(self):
        plot_frame = ttk.Frame(self.paned)
        self.paned.add(plot_frame, weight=1)

        self.fig = Figure(figsize=(7, 7), dpi=96)
        gs = GridSpec(2, 2, figure=self.fig, hspace=0.35, wspace=0.30)
        self.ax1 = self.fig.add_subplot(gs[0, 0])        # top-left
        self.ax_cutout = self.fig.add_subplot(gs[0, 1])   # top-right
        self.ax2 = self.fig.add_subplot(gs[1, :])          # bottom

        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        nav_frame = ttk.Frame(plot_frame)
        nav_frame.pack(side=tk.BOTTOM, fill=tk.X)
        self.nav_toolbar = NavigationToolbar2Tk(self.canvas, nav_frame)

    def _build_table_area(self):
        table_frame = ttk.Frame(self.paned)
        self.paned.add(table_frame, weight=1)

        # Search bar
        search_frame = ttk.Frame(table_frame)
        search_frame.pack(side=tk.TOP, fill=tk.X, padx=2, pady=2)
        ttk.Label(search_frame, text="Filter:").pack(side=tk.LEFT, padx=2)
        self.filter_var = tk.StringVar()
        self.filter_var.trace_add("write", lambda *a: self._apply_filter())
        ttk.Entry(search_frame, textvariable=self.filter_var,
                  width=20).pack(side=tk.LEFT, padx=2)
        ttk.Button(search_frame, text="Clear",
                   command=lambda: self.filter_var.set("")).pack(
                       side=tk.LEFT, padx=2)

        # Custom plot controls (shown in plot mode)
        self.plot_controls = ttk.Frame(table_frame)
        ttk.Label(self.plot_controls, text="Preset:").pack(side=tk.LEFT)
        preset_cb = ttk.Combobox(
            self.plot_controls, textvariable=self.plot_preset,
            values=list(PLOT_PRESETS.keys()), width=12, state="readonly")
        preset_cb.pack(side=tk.LEFT, padx=2)
        preset_cb.bind("<<ComboboxSelected>>", self._on_preset_select)

        ttk.Label(self.plot_controls, text="X:").pack(side=tk.LEFT)
        self.xcol_cb = ttk.Combobox(self.plot_controls,
                                     textvariable=self.plot_xcol,
                                     width=14, state="readonly")
        self.xcol_cb.pack(side=tk.LEFT, padx=2)
        ttk.Label(self.plot_controls, text="Y:").pack(side=tk.LEFT)
        self.ycol_cb = ttk.Combobox(self.plot_controls,
                                     textvariable=self.plot_ycol,
                                     width=14, state="readonly")
        self.ycol_cb.pack(side=tk.LEFT, padx=2)
        ttk.Combobox(self.plot_controls, textvariable=self.plot_type,
                     values=["scatter", "histogram"], width=10,
                     state="readonly").pack(side=tk.LEFT, padx=2)
        ttk.Button(self.plot_controls, text="Plot",
                   command=self._update_plots).pack(side=tk.LEFT, padx=5)

        # Treeview
        tree_container = ttk.Frame(table_frame)
        tree_container.pack(fill=tk.BOTH, expand=True)

        self.tree = ttk.Treeview(tree_container, show="headings",
                                  selectmode="browse")

        vsb = ttk.Scrollbar(tree_container, orient=tk.VERTICAL,
                             command=self.tree.yview)
        hsb = ttk.Scrollbar(tree_container, orient=tk.HORIZONTAL,
                             command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set,
                            xscrollcommand=hsb.set)

        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        tree_container.grid_rowconfigure(0, weight=1)
        tree_container.grid_columnconfigure(0, weight=1)

        self.tree.bind("<<TreeviewSelect>>", self._on_select)

    # ================================================================
    # Data Loading & Table
    # ================================================================

    def load_catalog(self, path):
        try:
            with open(path, 'r') as f:
                lines = f.readlines()
        except Exception as e:
            messagebox.showerror("Error", f"Cannot load catalog:\n{e}")
            return

        if not lines:
            return

        start = 0
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped and not stripped.startswith('#'):
                start = i
                break

        self.header = lines[start].strip().split('\t')
        self.rows = []
        for line in lines[start + 1:]:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            vals = line.split('\t')
            while len(vals) < len(self.header):
                vals.append("")
            self.rows.append(vals)

        self.catalog_path = path
        self.filtered_indices = list(range(len(self.rows)))
        self.selected_idx = -1
        self.fits_data = None  # reset cache

        self._rebuild_table()
        self._update_plot_col_lists()
        self.count_var.set(f"{len(self.rows)} sources")
        self.status_var.set(f"Loaded {len(self.rows)} sources from "
                            f"{os.path.basename(path)}")
        self._update_plots()

    def _rebuild_table(self):
        self.tree.delete(*self.tree.get_children())
        self.tree["columns"] = self.header
        for col in self.header:
            w = 90 if col != "NUMBER" else 60
            self.tree.column(col, width=w, minwidth=50, anchor=tk.E)
            self.tree.heading(
                col, text=col,
                command=lambda c=col: self._sort_by_column(c))

        for i in self.filtered_indices:
            self.tree.insert("", tk.END, iid=str(i),
                             values=self.rows[i])

    def _apply_filter(self):
        text = self.filter_var.get().strip().lower()
        if not text:
            self.filtered_indices = list(range(len(self.rows)))
        else:
            self.filtered_indices = [
                i for i, row in enumerate(self.rows)
                if any(text in v.lower() for v in row)]

        self.tree.delete(*self.tree.get_children())
        for i in self.filtered_indices:
            self.tree.insert("", tk.END, iid=str(i),
                             values=self.rows[i])
        self.count_var.set(
            f"{len(self.filtered_indices)}/{len(self.rows)} sources")

    def _sort_by_column(self, col):
        if self.sort_col == col:
            self.sort_asc = not self.sort_asc
        else:
            self.sort_col = col
            self.sort_asc = True

        col_idx = self.header.index(col)

        def sort_key(i):
            try:
                return float(self.rows[i][col_idx])
            except (ValueError, TypeError):
                return self.rows[i][col_idx]

        self.filtered_indices.sort(key=sort_key, reverse=not self.sort_asc)

        self.tree.delete(*self.tree.get_children())
        for i in self.filtered_indices:
            self.tree.insert("", tk.END, iid=str(i),
                             values=self.rows[i])

    def _on_select(self, event):
        sel = self.tree.selection()
        if sel:
            self.selected_idx = int(sel[0])
            self._update_plots()

    def _col_idx(self, name):
        try:
            return self.header.index(name)
        except ValueError:
            return -1

    def _col_values(self, name, as_float=True):
        idx = self._col_idx(name)
        if idx < 0:
            return None
        if as_float:
            vals = []
            for row in self.rows:
                try:
                    vals.append(float(row[idx]))
                except (ValueError, IndexError):
                    vals.append(np.nan)
            return np.array(vals)
        return [row[idx] for row in self.rows]

    def _has_col(self, name):
        return name in self.header

    def _update_plot_col_lists(self):
        self.xcol_cb['values'] = self.header
        self.ycol_cb['values'] = self.header

    def _on_preset_select(self, event=None):
        preset = self.plot_preset.get()
        xcol, ycol = PLOT_PRESETS.get(preset, ("", ""))
        if xcol and xcol in self.header:
            self.plot_xcol.set(xcol)
        if ycol and ycol in self.header:
            self.plot_ycol.set(ycol)
        if xcol:
            self._update_plots()

    # ================================================================
    # FITS Data & Cutout
    # ================================================================

    def _load_fits_data(self):
        """Load and cache the FITS image data."""
        if self.fits_data is not None:
            return True
        if not self.fits_path or not os.path.exists(self.fits_path):
            return False

        if not HAS_ASTROPY:
            return False

        try:
            with pyfits.open(self.fits_path) as hdul:
                for hdu in hdul:
                    if hdu.data is not None and hdu.data.ndim == 2:
                        self.fits_data = hdu.data.astype(np.float32)
                        return True
            return False
        except Exception:
            return False

    def _plot_cutout(self):
        """Plot FITS cutout for the selected source."""
        self.ax_cutout.clear()

        if self.selected_idx < 0:
            self.ax_cutout.text(
                0.5, 0.5, "Select a source\nto view cutout",
                transform=self.ax_cutout.transAxes,
                ha='center', va='center', color='gray', fontsize=9)
            self.ax_cutout.set_xticks([])
            self.ax_cutout.set_yticks([])
            return

        if not self._load_fits_data():
            self.ax_cutout.text(
                0.5, 0.5, "No FITS data" +
                ("\n(install astropy)" if not HAS_ASTROPY else ""),
                transform=self.ax_cutout.transAxes,
                ha='center', va='center', color='gray', fontsize=9)
            self.ax_cutout.set_xticks([])
            self.ax_cutout.set_yticks([])
            return

        xi = self._col_idx("X_IMAGE")
        yi = self._col_idx("Y_IMAGE")
        if xi < 0 or yi < 0:
            return

        try:
            x = float(self.rows[self.selected_idx][xi])
            y = float(self.rows[self.selected_idx][yi])
        except (ValueError, IndexError):
            return

        half = self.settings.get('cutout_size', 80) // 2
        ix, iy = int(round(x)), int(round(y))
        ny, nx = self.fits_data.shape

        y0 = max(0, iy - half)
        y1 = min(ny, iy + half)
        x0 = max(0, ix - half)
        x1 = min(nx, ix + half)

        cutout = self.fits_data[y0:y1, x0:x1]
        if cutout.size == 0:
            return

        # Log stretch
        vmin = np.nanpercentile(cutout, 1)
        vmax = np.nanpercentile(cutout, 99.5)
        self.ax_cutout.imshow(
            cutout, origin='lower', cmap='gray',
            vmin=vmin, vmax=vmax,
            extent=[x0, x1, y0, y1])

        # Ellipse overlay (3x Kron-like)
        ai = self._col_idx("A_IMAGE")
        bi = self._col_idx("B_IMAGE")
        ti = self._col_idx("THETA_IMAGE")
        if ai >= 0 and bi >= 0:
            try:
                a = float(self.rows[self.selected_idx][ai])
                b = float(self.rows[self.selected_idx][bi])
                theta = float(
                    self.rows[self.selected_idx][ti]) if ti >= 0 else 0
                ell = Ellipse(
                    (x, y), width=2 * a * 3, height=2 * b * 3,
                    angle=theta, fill=False, color='lime', lw=1.5)
                self.ax_cutout.add_patch(ell)
            except (ValueError, IndexError):
                pass

        # Crosshair
        self.ax_cutout.axhline(y, color='cyan', lw=0.5, alpha=0.5)
        self.ax_cutout.axvline(x, color='cyan', lw=0.5, alpha=0.5)

        num_idx = self._col_idx("NUMBER")
        num = self.rows[self.selected_idx][num_idx] if num_idx >= 0 else "?"
        self.ax_cutout.set_title(f"#{num}", fontsize=9)

    # ================================================================
    # Settings Save/Load
    # ================================================================

    def _load_settings(self):
        if not os.path.exists(_SETTINGS_FILE):
            return
        try:
            with open(_SETTINGS_FILE, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    eq = line.find('=')
                    if eq < 0:
                        continue
                    key = line[:eq].strip()
                    val = line[eq + 1:].strip()
                    if key in self.settings:
                        # Preserve type
                        orig = DEFAULT_SETTINGS.get(key, '')
                        if isinstance(orig, int):
                            try:
                                val = int(val)
                            except ValueError:
                                continue
                        elif isinstance(orig, float):
                            try:
                                val = float(val)
                            except ValueError:
                                continue
                        self.settings[key] = val
        except Exception:
            pass

    def _save_settings(self):
        os.makedirs(os.path.dirname(_SETTINGS_FILE), exist_ok=True)
        try:
            with open(_SETTINGS_FILE, 'w') as f:
                for key, val in self.settings.items():
                    f.write(f"{key}={val}\n")
        except Exception:
            pass

    def _settings_general(self):
        """General settings dialog."""
        w = tk.Toplevel(self.root)
        w.title("General Settings")
        w.geometry("350x200")
        w.transient(self.root)

        row = 0
        ttk.Label(w, text="Workers:").grid(
            row=row, column=0, padx=8, pady=5, sticky='w')
        nw_var = tk.IntVar(value=self.settings['n_workers'])
        ttk.Spinbox(w, from_=1, to=32, textvariable=nw_var,
                     width=6).grid(row=row, column=1, padx=8, pady=5)

        row += 1
        ttk.Label(w, text="Cutout size (pix):").grid(
            row=row, column=0, padx=8, pady=5, sticky='w')
        cs_var = tk.IntVar(value=self.settings['cutout_size'])
        ttk.Spinbox(w, from_=20, to=500, textvariable=cs_var,
                     width=6).grid(row=row, column=1, padx=8, pady=5)

        row += 1
        ttk.Label(w, text="Pixel scale (arcsec):").grid(
            row=row, column=0, padx=8, pady=5, sticky='w')
        ps_var = tk.DoubleVar(value=self.settings['pixel_scale'])
        ttk.Entry(w, textvariable=ps_var, width=8).grid(
            row=row, column=1, padx=8, pady=5)

        row += 1
        ttk.Label(w, text="Mag zeropoint:").grid(
            row=row, column=0, padx=8, pady=5, sticky='w')
        zp_var = tk.DoubleVar(value=self.settings['mag_zeropoint'])
        ttk.Entry(w, textvariable=zp_var, width=8).grid(
            row=row, column=1, padx=8, pady=5)

        def apply():
            self.settings['n_workers'] = nw_var.get()
            self.settings['cutout_size'] = cs_var.get()
            self.settings['pixel_scale'] = ps_var.get()
            self.settings['mag_zeropoint'] = zp_var.get()
            self._save_settings()
            self.fits_data = None  # reset cutout on size change
            w.destroy()

        row += 1
        ttk.Button(w, text="OK", command=apply).grid(
            row=row, column=0, columnspan=2, pady=10)

    def _settings_sed(self):
        """SED Fitting settings dialog."""
        w = tk.Toplevel(self.root)
        w.title("SED Fitting Settings")
        w.geometry("420x300")
        w.transient(self.root)

        entries = {}
        fields = [
            ('sed_backend',       'Backend:',           12),
            ('sed_bands',         'Bands:',             30),
            ('sed_mag_columns',   'Mag columns:',       30),
            ('sed_photoz_column', 'Photo-z column:',    15),
            ('sed_ckpt_emulator', 'Emulator ckpt:',     30),
            ('sed_ckpt_inverse',  'Inverse ckpt:',      30),
        ]

        for row, (key, label, width) in enumerate(fields):
            ttk.Label(w, text=label).grid(
                row=row, column=0, padx=8, pady=4, sticky='w')
            if key == 'sed_backend':
                var = tk.StringVar(value=self.settings[key])
                cb = ttk.Combobox(
                    w, textvariable=var, width=width, state="readonly",
                    values=["auto", "analytic", "fsps", "bagpipes",
                            "prospector", "cigale", "dense_basis"])
                cb.grid(row=row, column=1, padx=8, pady=4, sticky='w')
                entries[key] = var
            else:
                ent = ttk.Entry(w, width=width)
                ent.insert(0, str(self.settings[key]))
                ent.grid(row=row, column=1, padx=8, pady=4, sticky='w')
                entries[key] = ent

        def apply():
            for key, widget in entries.items():
                if isinstance(widget, tk.StringVar):
                    self.settings[key] = widget.get()
                else:
                    self.settings[key] = widget.get()
            self.sed_backend.set(self.settings['sed_backend'])
            self._save_settings()
            w.destroy()

        ttk.Button(w, text="OK", command=apply).grid(
            row=len(fields), column=0, columnspan=2, pady=10)

    def _settings_photoz(self):
        """Photo-z settings dialog."""
        w = tk.Toplevel(self.root)
        w.title("Photo-z Settings")
        w.geometry("420x180")
        w.transient(self.root)

        entries = {}
        fields = [
            ('photoz_bands',       'Bands:',      30),
            ('photoz_mag_columns', 'Mag columns:', 30),
            ('photoz_checkpoint',  'Checkpoint:',  30),
        ]

        for row, (key, label, width) in enumerate(fields):
            ttk.Label(w, text=label).grid(
                row=row, column=0, padx=8, pady=4, sticky='w')
            ent = ttk.Entry(w, width=width)
            ent.insert(0, str(self.settings[key]))
            ent.grid(row=row, column=1, padx=8, pady=4, sticky='w')
            entries[key] = ent

        def apply():
            for key, ent in entries.items():
                self.settings[key] = ent.get()
            self._save_settings()
            w.destroy()

        ttk.Button(w, text="OK", command=apply).grid(
            row=len(fields), column=0, columnspan=2, pady=10)

    # ================================================================
    # Compute Color Columns
    # ================================================================

    def _compute_colors(self):
        """Auto-detect magnitude columns and create color columns."""
        if not self.rows:
            messagebox.showwarning("Warning", "Load a catalog first.")
            return

        mag_cols = [c for c in self.header
                    if c.startswith("MAG_") and "ERR" not in c]
        if len(mag_cols) < 2:
            messagebox.showinfo("Info",
                                "Need at least 2 MAG_ columns for colors.")
            return

        n_added = 0
        for i in range(len(mag_cols) - 1):
            c1, c2 = mag_cols[i], mag_cols[i + 1]
            color_name = f"COLOR_{c1.replace('MAG_', '')}_{c2.replace('MAG_', '')}"
            if color_name in self.header:
                continue

            idx1 = self._col_idx(c1)
            idx2 = self._col_idx(c2)
            self.header.append(color_name)
            for row in self.rows:
                try:
                    v1 = float(row[idx1])
                    v2 = float(row[idx2])
                    if 0 < v1 < 90 and 0 < v2 < 90:
                        row.append(f"{v1 - v2:.4f}")
                    else:
                        row.append("-99")
                except (ValueError, IndexError):
                    row.append("-99")
            n_added += 1

        if n_added > 0:
            self._rebuild_table()
            self._update_plot_col_lists()
            self.status_var.set(f"Added {n_added} color columns")
        else:
            self.status_var.set("Color columns already exist")

    # ================================================================
    # Analysis Runners
    # ================================================================

    def _check_ready(self):
        if not self.rows:
            messagebox.showwarning("Warning", "Load a catalog first.")
            return False
        return True

    def _check_fits(self):
        if not self.fits_path or not os.path.exists(self.fits_path):
            messagebox.showwarning("Warning",
                                    "FITS file required. Use File > Open FITS.")
            return False
        return True

    def _save_temp_catalog(self, suffix="analysis"):
        tmpdir = os.path.join(os.path.expanduser("~"), ".ds9")
        os.makedirs(tmpdir, exist_ok=True)
        path = os.path.join(tmpdir, f"{suffix}_gui_catalog.tsv")
        with open(path, 'w') as f:
            f.write('\t'.join(self.header) + '\n')
            for row in self.rows:
                f.write('\t'.join(row) + '\n')
        return path

    def _find_script(self, name):
        path = os.path.join(_script_dir, name)
        return path if os.path.exists(path) else None

    def _run_script(self, script_name, extra_args, output_columns,
                    status_msg):
        if not self._check_ready():
            return

        script = self._find_script(script_name)
        if not script:
            messagebox.showerror("Error", f"Script not found: {script_name}")
            return

        tmpcat = self._save_temp_catalog(
            script_name.replace('.py', '').replace('ds9_', ''))

        self.status_var.set(status_msg)
        self.root.update_idletasks()

        def run():
            try:
                cmd = [sys.executable, script]
                if self.fits_path:
                    cmd.append(self.fits_path)
                cmd.extend(["--catalog", tmpcat])
                cmd.extend(extra_args)
                cmd.extend(["--n-workers",
                             str(self.settings['n_workers'])])

                result = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=600)

                if result.returncode != 0:
                    self.root.after(0, lambda: self._show_error(
                        f"{script_name} failed:\n{result.stderr[:500]}"))
                    return

                stdout = result.stdout.strip()
                if not stdout:
                    self.root.after(0, lambda: self.status_var.set(
                        f"{script_name}: no output"))
                    return

                self.root.after(0, lambda: self._merge_results(
                    stdout, output_columns,
                    status_msg.replace("...", " done")))

            except subprocess.TimeoutExpired:
                self.root.after(0, lambda: self._show_error(
                    f"{script_name} timed out"))
            except Exception as e:
                self.root.after(0, lambda: self._show_error(str(e)))

        threading.Thread(target=run, daemon=True).start()

    def _merge_results(self, tsv_data, col_names, done_msg):
        lines = tsv_data.strip().split('\n')

        header_idx = 0
        for i, line in enumerate(lines):
            if not line.startswith('#'):
                header_idx = i
                break

        res_header = lines[header_idx].split('\t')
        res_rows = {}
        num_col = -1
        for i, h in enumerate(res_header):
            if h.strip() == "NUMBER":
                num_col = i
                break

        if num_col < 0:
            self.status_var.set("Error: NUMBER column not found in output")
            return

        for line in lines[header_idx + 1:]:
            vals = line.strip().split('\t')
            if len(vals) > num_col:
                res_rows[vals[num_col]] = vals

        col_indices = {}
        for cname in col_names:
            for i, h in enumerate(res_header):
                if h.strip() == cname:
                    col_indices[cname] = i
                    break

        if not col_indices:
            self.status_var.set("No matching columns in output")
            return

        my_num_col = self._col_idx("NUMBER")
        for cname in col_names:
            if cname not in self.header:
                self.header.append(cname)
                for row in self.rows:
                    row.append("-99")

        for row in self.rows:
            num = row[my_num_col] if my_num_col >= 0 else ""
            if num in res_rows:
                rvals = res_rows[num]
                for cname, ridx in col_indices.items():
                    cidx = self._col_idx(cname)
                    if cidx >= 0 and ridx < len(rvals):
                        row[cidx] = rvals[ridx]

        self._rebuild_table()
        self._update_plot_col_lists()
        self._set_mode(self.mode)
        self.status_var.set(done_msg)

    def _show_error(self, msg):
        self.status_var.set("Error")
        messagebox.showerror("Analysis Error", msg)

    # --- Individual runners ---

    def run_photoz(self):
        if not self._check_fits():
            return
        args = []
        b = self.settings['photoz_bands']
        if b:
            args.extend(["--bands", b])
        m = self.settings['photoz_mag_columns']
        if m:
            args.extend(["--mag-columns", m])
        ckpt = self.settings['photoz_checkpoint']
        if ckpt and os.path.exists(ckpt):
            args.extend(["--checkpoint", ckpt])
        self._run_script(
            "ds9_photo_z.py", args,
            ["PHOTO_Z", "PHOTO_Z_ERR", "PHOTO_Z_Q68", "PHOTO_Z_OUTLIER"],
            "Running Photo-z (AI)...")
        self._set_mode("photoz")

    def run_sed_fit(self):
        if not self._check_fits():
            return
        s = self.settings
        args = ["--backend", self.sed_backend.get()]
        if s['sed_bands']:
            args.extend(["--bands", s['sed_bands']])
        if s['sed_mag_columns']:
            args.extend(["--mag-columns", s['sed_mag_columns']])
        if s['sed_photoz_column']:
            args.extend(["--photo-z-column", s['sed_photoz_column']])
        if s['sed_ckpt_emulator'] and os.path.exists(s['sed_ckpt_emulator']):
            args.extend(["--checkpoint-emulator", s['sed_ckpt_emulator']])
        if s['sed_ckpt_inverse'] and os.path.exists(s['sed_ckpt_inverse']):
            args.extend(["--checkpoint-inverse", s['sed_ckpt_inverse']])
        self._run_script(
            "ds9_sed_fit.py", args,
            ["LOG_MASS", "LOG_MASS_ERR", "LOG_AGE", "LOG_AGE_ERR",
             "LOG_Z", "AV", "SFR"],
            f"Running SED Fitting ({self.sed_backend.get()})...")
        self._set_mode("sed")

    def run_morphology(self):
        if not self._check_fits():
            return
        self._run_script("ds9_galaxy_morph.py", [],
                          ["MORPH_TYPE", "MORPH_CONF"],
                          "Running Galaxy Morphology (CNN)...")
        self._set_mode("classify")

    def run_starfinder(self):
        if not self._check_fits():
            return
        self._run_script("ds9_star_finder.py", [],
                          ["AI_STAR", "AI_STAR_CONF"],
                          "Running Star/Galaxy (CNN)...")
        self._set_mode("classify")

    def run_morphometry(self):
        if not self._check_fits():
            return
        self._run_script("ds9_morphometry.py", [],
                          ["CONC", "ASYM", "GINI", "M20", "R_PETRO"],
                          "Running CAS/Gini/M20...")
        self._set_mode("structure")

    def run_sersic(self):
        if not self._check_fits():
            return
        self._run_script("ds9_sersic.py", [],
                          ["SERSIC_N", "SERSIC_RE", "SERSIC_MU_E",
                           "SERSIC_CHI2"],
                          "Running Sérsic Fitting...")
        self._set_mode("structure")

    def run_bulge_disk(self):
        if not self._check_fits():
            return
        self._run_script(
            "ds9_bulge_disk.py", [],
            ["BD_BULGE_FRAC", "BD_BULGE_RE", "BD_BULGE_N",
             "BD_DISK_RE", "BD_DISK_N", "BD_CHI2"],
            "Running Bulge+Disk Decomposition...")

    def run_psf_phot(self):
        if not self._check_fits():
            return
        self._run_script("ds9_psf_phot.py", [],
                          ["MAG_PSF", "MAGERR_PSF", "PSF_CHI2"],
                          "Running PSF Photometry...")

    def run_multiband(self):
        if not self._check_fits():
            return
        self._run_script("ds9_multiband.py", [],
                          [],  # columns depend on bands
                          "Running Multi-Band Photometry...")

    def run_crowded_phot(self):
        if not self._check_fits():
            return
        self._run_script("ds9_crowded_phot.py", [],
                          ["MAG_CROWD", "MAGERR_CROWD"],
                          "Running Crowded Field Photometry...")

    def run_crossmatch(self):
        if not self._check_ready():
            return
        self._run_script("ds9_crossmatch.py", [],
                          [],  # columns depend on external catalog
                          "Running Cross-Match (VizieR)...")

    def run_completeness(self):
        if not self._check_fits():
            return
        self._run_script("ds9_completeness.py", [],
                          [],  # outputs statistics
                          "Running Completeness Simulation...")

    # ================================================================
    # Mode & Plot Management
    # ================================================================

    def _set_mode(self, mode):
        self.mode = mode
        self.mode_var.set(mode)

        if mode == "plot":
            self.plot_controls.pack(side=tk.TOP, fill=tk.X, padx=2, pady=2)
        else:
            self.plot_controls.pack_forget()

        self._update_plots()

    def _update_plots(self):
        self.ax1.clear()
        self.ax_cutout.clear()
        self.ax2.clear()

        # Remove leftover colorbars
        while len(self.fig.axes) > 3:
            self.fig.axes[-1].remove()

        if not self.rows:
            self.ax1.text(0.5, 0.5, "Load a catalog to begin",
                          transform=self.ax1.transAxes, ha='center',
                          fontsize=12, color='gray')
            self.ax_cutout.set_xticks([])
            self.ax_cutout.set_yticks([])
            self.canvas.draw()
            return

        # Always draw cutout
        self._plot_cutout()

        dispatch = {
            "overview":  self._plot_overview,
            "photoz":    self._plot_photoz,
            "sed":       self._plot_sed,
            "classify":  self._plot_classification,
            "structure": self._plot_structure,
            "plot":      self._plot_custom,
        }
        try:
            dispatch.get(self.mode, self._plot_overview)()
        except Exception as e:
            self.ax1.text(0.5, 0.5, f"Plot error: {e}",
                          transform=self.ax1.transAxes, ha='center',
                          fontsize=10, color='red', wrap=True)

        self.fig.tight_layout(pad=2.0)
        self.canvas.draw()

    # ================================================================
    # Plot Functions
    # ================================================================

    def _sel_val(self, col):
        """Get selected source's value for a column (float or None)."""
        if self.selected_idx < 0:
            return None
        idx = self._col_idx(col)
        if idx < 0:
            return None
        try:
            v = float(self.rows[self.selected_idx][idx])
            return v if np.isfinite(v) else None
        except (ValueError, IndexError):
            return None

    def _highlight_selected(self, ax, xcol, ycol, marker='*', s=100,
                             color='red'):
        """Plot a highlight marker for the selected source."""
        sx = self._sel_val(xcol)
        sy = self._sel_val(ycol)
        if sx is not None and sy is not None:
            ax.scatter([sx], [sy], s=s, c=color, marker=marker,
                       zorder=10, edgecolors='black', linewidths=0.5)

    def _plot_overview(self):
        # Magnitude histogram
        mag = self._col_values("MAG_AUTO")
        if mag is not None:
            valid = mag[np.isfinite(mag) & (mag > 0) & (mag < 90)]
            if len(valid) > 0:
                self.ax1.hist(valid, bins=40, color='steelblue',
                              edgecolor='white', alpha=0.8)
                self.ax1.set_xlabel("MAG_AUTO")
                self.ax1.set_ylabel("Count")
                self.ax1.set_title("Magnitude Distribution")
                sv = self._sel_val("MAG_AUTO")
                if sv is not None and 0 < sv < 90:
                    self.ax1.axvline(sv, color='red', lw=2, ls='--',
                                      label=f"Sel: {sv:.1f}")
                    self.ax1.legend(fontsize=8)

        # Source positions
        x = self._col_values("X_IMAGE")
        y = self._col_values("Y_IMAGE")
        if x is not None and y is not None:
            v = np.isfinite(x) & np.isfinite(y)
            self.ax2.scatter(x[v], y[v], s=2, alpha=0.5, c='steelblue')
            self.ax2.set_xlabel("X_IMAGE")
            self.ax2.set_ylabel("Y_IMAGE")
            self.ax2.set_title("Source Positions")
            self.ax2.set_aspect('equal')
            self._highlight_selected(self.ax2, "X_IMAGE", "Y_IMAGE",
                                      marker='+', s=80, color='red')

    def _plot_photoz(self):
        pz = self._col_values("PHOTO_Z")
        if pz is None:
            self.ax1.text(0.5, 0.5,
                          "Run Photo-z first\n(Analysis > Photo-z)",
                          transform=self.ax1.transAxes, ha='center',
                          fontsize=12, color='gray')
            return

        valid = np.isfinite(pz) & (pz > 0) & (pz < 10)
        pz_v = pz[valid]

        if len(pz_v) > 0:
            self.ax1.hist(pz_v, bins=50, color='darkorange',
                          edgecolor='white', alpha=0.8)
            self.ax1.set_xlabel("Photo-z")
            self.ax1.set_ylabel("Count")
            med = np.median(pz_v)
            self.ax1.set_title(f"Photo-z Distribution (N={len(pz_v)}, "
                                f"med={med:.2f})")
            sz = self._sel_val("PHOTO_Z")
            if sz is not None:
                self.ax1.axvline(sz, color='cyan', lw=2,
                                  label=f"Sel: z={sz:.3f}")
                self.ax1.legend(fontsize=8)

        # z vs MAG
        mag = self._col_values("MAG_AUTO")
        if mag is not None:
            both = valid & np.isfinite(mag) & (mag > 0) & (mag < 90)
            if np.sum(both) > 0:
                self.ax2.scatter(pz[both], mag[both], s=5, alpha=0.4,
                                  c='darkorange')
                self.ax2.set_xlabel("Photo-z")
                self.ax2.set_ylabel("MAG_AUTO")
                self.ax2.set_title("Magnitude vs Redshift")
                self.ax2.invert_yaxis()
                self._highlight_selected(
                    self.ax2, "PHOTO_Z", "MAG_AUTO")
        else:
            pz_err = self._col_values("PHOTO_Z_ERR")
            if pz_err is not None and len(pz_v) > 0:
                self.ax2.scatter(pz_v, pz_err[valid], s=5, alpha=0.4,
                                 c='darkorange')
                self.ax2.set_xlabel("Photo-z")
                self.ax2.set_ylabel("Photo-z Error")
                self.ax2.set_title("Redshift Uncertainty")

    def _plot_sed(self):
        lm = self._col_values("LOG_MASS")
        if lm is None:
            self.ax1.text(0.5, 0.5,
                          "Run SED Fitting first\n(Analysis > SED Fitting)",
                          transform=self.ax1.transAxes, ha='center',
                          fontsize=12, color='gray')
            return

        la = self._col_values("LOG_AGE")
        valid_m = np.isfinite(lm) & (lm > 0)

        # Top: SED for selected source, or mass histogram
        if self.selected_idx >= 0 and self._has_sed_data():
            self._plot_source_sed(self.selected_idx)
        elif np.sum(valid_m) > 0:
            self.ax1.hist(lm[valid_m], bins=40, color='mediumpurple',
                          edgecolor='white', alpha=0.8)
            self.ax1.set_xlabel("log(M*/M$_\\odot$)")
            self.ax1.set_ylabel("Count")
            self.ax1.set_title("Stellar Mass Distribution")
            sv = self._sel_val("LOG_MASS")
            if sv is not None and sv > 0:
                self.ax1.axvline(sv, color='red', lw=2, ls='--')

        # Bottom: SFR-Mass diagram (star-forming main sequence)
        sfr = self._col_values("SFR")
        if sfr is not None and la is not None:
            valid_s = np.isfinite(sfr) & valid_m
            if np.sum(valid_s) > 0:
                av = self._col_values("AV")
                if av is not None:
                    c = np.clip(np.where(np.isfinite(av), av, 0), 0, 3)
                    sc = self.ax2.scatter(
                        lm[valid_s], sfr[valid_s], s=8,
                        c=c[valid_s], cmap='RdYlBu_r', alpha=0.6,
                        vmin=0, vmax=3)
                    self.fig.colorbar(sc, ax=self.ax2, label="Av",
                                      shrink=0.8, pad=0.02)
                else:
                    self.ax2.scatter(lm[valid_s], sfr[valid_s], s=8,
                                     alpha=0.5, c='mediumpurple')

                # Main sequence reference (Speagle+2014 at z~1)
                ms_x = np.linspace(8, 12, 50)
                ms_y = 0.84 * ms_x - 8.3  # log SFR
                self.ax2.plot(ms_x, ms_y, 'k--', alpha=0.3, lw=1,
                               label='MS (z~1)')
                self.ax2.legend(fontsize=7)

                self.ax2.set_xlabel("log(M*/M$_\\odot$)")
                self.ax2.set_ylabel("log(SFR)")
                self.ax2.set_title("SFR-Mass Diagram")
                self._highlight_selected(
                    self.ax2, "LOG_MASS", "SFR")
        elif la is not None:
            # Fallback: Mass-Age
            valid_a = np.isfinite(la) & (la > 0)
            both = valid_m & valid_a
            if np.sum(both) > 0:
                self.ax2.scatter(la[both], lm[both], s=8, alpha=0.5,
                                 c='mediumpurple')
                self.ax2.set_xlabel("log(Age/yr)")
                self.ax2.set_ylabel("log(M*/M$_\\odot$)")
                self.ax2.set_title("Mass-Age Diagram")
                self._highlight_selected(
                    self.ax2, "LOG_AGE", "LOG_MASS")

    def _has_sed_data(self):
        if self.selected_idx < 0:
            return False
        sv = self._sel_val("LOG_MASS")
        return sv is not None and sv > 0

    def _plot_source_sed(self, row_idx):
        """Plot SED for a single source."""
        mag_cols = [c for c in self.header
                    if c.startswith("MAG_") and "ERR" not in c
                    and c not in ("MAG_PSF", "MAG_CROWD")]
        if not mag_cols:
            mag_cols = ["MAG_AUTO"]

        wavelengths, obs_mags, obs_labels = [], [], []
        for mc in mag_cols:
            idx = self._col_idx(mc)
            if idx < 0:
                continue
            try:
                val = float(self.rows[row_idx][idx])
                if val <= 0 or val >= 90:
                    continue
            except (ValueError, IndexError):
                continue

            fname = mc.replace("MAG_AUTO_", "").replace("MAG_", "")
            if fname in FILTER_WAVELENGTHS:
                wavelengths.append(FILTER_WAVELENGTHS[fname])
                obs_mags.append(val)
                obs_labels.append(fname)

        if len(wavelengths) >= 2:
            wl = np.array(wavelengths)
            om = np.array(obs_mags)
            s = np.argsort(wl)
            wl, om = wl[s], om[s]
            labels_s = [obs_labels[i] for i in s]

            self.ax1.plot(wl / 1e4, om, 'o-', color='steelblue',
                          markersize=6, label="Observed")
            for i, lbl in enumerate(labels_s):
                self.ax1.annotate(lbl, (wl[i] / 1e4, om[i]),
                                   fontsize=6, ha='center',
                                   va='bottom', rotation=45)

            # Model SED
            try:
                lm = self._sel_val("LOG_MASS")
                la = self._sel_val("LOG_AGE")
                lz = self._sel_val("LOG_Z")
                av = self._sel_val("AV")
                pz = self._sel_val("PHOTO_Z") or 0.1
                if all(v is not None for v in [lm, la, lz, av]):
                    from sed_fit.grid.generate import simple_sed_model
                    model_wl = np.linspace(
                        wl.min() * 0.8, wl.max() * 1.2, 100)
                    model_m = simple_sed_model(
                        pz, lm, la, lz, av, 9.0, model_wl)
                    self.ax1.plot(model_wl / 1e4, model_m, '-',
                                  color='tomato', alpha=0.7,
                                  label="Model")
            except Exception:
                pass

            self.ax1.invert_yaxis()
            self.ax1.set_xlabel("Wavelength (μm)")
            self.ax1.set_ylabel("AB mag")
            num_idx = self._col_idx("NUMBER")
            num = self.rows[row_idx][num_idx] if num_idx >= 0 else "?"
            self.ax1.set_title(f"SED - Source #{num}")
            self.ax1.legend(fontsize=8)
        else:
            # Text summary
            params = {}
            for pname in ["LOG_MASS", "LOG_AGE", "LOG_Z", "AV", "SFR"]:
                v = self._sel_val(pname)
                if v is not None:
                    params[pname] = v
            if params:
                text = "SED Fit Results:\n" + "\n".join(
                    f"  {k} = {v:.2f}" for k, v in params.items())
                self.ax1.text(0.1, 0.5, text,
                              transform=self.ax1.transAxes,
                              fontsize=11, family='monospace',
                              verticalalignment='center')

    def _plot_classification(self):
        has_morph = self._has_col("MORPH_TYPE")
        has_star = self._has_col("AI_STAR")

        if not has_morph and not has_star:
            self.ax1.text(0.5, 0.5,
                          "Run Classification first\n(Analysis menu)",
                          transform=self.ax1.transAxes, ha='center',
                          fontsize=12, color='gray')
            return

        # Morphology distribution
        if has_morph:
            mt = self._col_values("MORPH_TYPE", as_float=False)
            types = [t for t in mt if t and t != "-99"]
            if types:
                from collections import Counter
                counts = Counter(types)
                labels = list(counts.keys())
                values = list(counts.values())
                short = [l[:15] for l in labels]
                colors = matplotlib.colormaps['tab10'](
                    np.linspace(0, 1, len(labels)))

                bars = self.ax1.barh(range(len(labels)), values,
                                      color=colors, edgecolor='white')
                self.ax1.set_yticks(range(len(labels)))
                self.ax1.set_yticklabels(short, fontsize=8)
                self.ax1.set_xlabel("Count")
                self.ax1.set_title(
                    f"Morphology Distribution (N={len(types)})")

                si = self._col_idx("MORPH_TYPE")
                if self.selected_idx >= 0 and si >= 0:
                    sel_type = self.rows[self.selected_idx][si]
                    if sel_type in labels:
                        bi = labels.index(sel_type)
                        bars[bi].set_edgecolor('red')
                        bars[bi].set_linewidth(3)

        # Star/Galaxy
        if has_star:
            ai_star = self._col_values("AI_STAR")
            valid = np.isfinite(ai_star)
            if np.sum(valid) > 0:
                n_star = int(np.sum(ai_star[valid] > 0.5))
                n_gal = int(np.sum(ai_star[valid] <= 0.5))
                self.ax2.pie([n_star, n_gal],
                              labels=[f"Star ({n_star})",
                                      f"Galaxy ({n_gal})"],
                              colors=['cyan', 'gold'],
                              autopct='%1.0f%%', startangle=90)
                self.ax2.set_title("Star/Galaxy Classification")
        elif has_morph:
            mc = self._col_values("MORPH_CONF")
            if mc is not None:
                valid = np.isfinite(mc)
                if np.sum(valid) > 0:
                    self.ax2.hist(mc[valid], bins=30, color='teal',
                                  edgecolor='white', alpha=0.8)
                    self.ax2.set_xlabel("Confidence")
                    self.ax2.set_ylabel("Count")
                    self.ax2.set_title("Morphology Confidence")

    def _plot_structure(self):
        has_gini = self._has_col("GINI")
        has_sersic = self._has_col("SERSIC_N")
        has_cas = self._has_col("CONC")

        if not has_gini and not has_sersic and not has_cas:
            self.ax1.text(0.5, 0.5,
                          "Run CAS/Gini/M20 or Sérsic first\n"
                          "(Analysis menu)",
                          transform=self.ax1.transAxes, ha='center',
                          fontsize=12, color='gray')
            return

        # Top: Gini-M20 or C-A diagram
        if has_gini:
            gini = self._col_values("GINI")
            m20 = self._col_values("M20")
            if gini is not None and m20 is not None:
                v = np.isfinite(gini) & np.isfinite(m20) & (gini > 0)
                if np.sum(v) > 0:
                    self.ax1.scatter(m20[v], gini[v], s=8,
                                     alpha=0.5, c='teal')
                    self.ax1.set_xlabel("M$_{20}$")
                    self.ax1.set_ylabel("Gini")
                    self.ax1.set_title("Gini-M20 Diagram")
                    # Merger line (Lotz+2008)
                    x_l = np.linspace(-3, 0, 50)
                    self.ax1.plot(x_l, -0.14 * x_l + 0.33,
                                  'r--', lw=1, alpha=0.5,
                                  label="Merger (Lotz+08)")
                    self.ax1.legend(fontsize=7)
                    self._highlight_selected(self.ax1, "M20", "GINI")
        elif has_cas:
            conc = self._col_values("CONC")
            asym = self._col_values("ASYM")
            if conc is not None and asym is not None:
                v = np.isfinite(conc) & np.isfinite(asym)
                if np.sum(v) > 0:
                    self.ax1.scatter(asym[v], conc[v], s=8,
                                     alpha=0.5, c='coral')
                    self.ax1.set_xlabel("Asymmetry")
                    self.ax1.set_ylabel("Concentration")
                    self.ax1.set_title("C-A Diagram")
                    self._highlight_selected(self.ax1, "ASYM", "CONC")

        # Bottom: Sérsic n vs Re
        if has_sersic:
            sn = self._col_values("SERSIC_N")
            sre = self._col_values("SERSIC_RE")
            if sre is not None:
                v = (np.isfinite(sn) & np.isfinite(sre)
                     & (sn > 0) & (sre > 0))
                if np.sum(v) > 0:
                    self.ax2.scatter(sre[v], sn[v], s=8, alpha=0.5,
                                     c='darkorchid')
                    self.ax2.set_xlabel("R$_e$ (pix)")
                    self.ax2.set_ylabel("Sérsic n")
                    self.ax2.set_title("Sérsic Index vs R_e")
                    self.ax2.set_xscale('log')
                    self.ax2.axhline(2.5, color='gray', ls='--',
                                      alpha=0.5, label="n=2.5")
                    self.ax2.legend(fontsize=7)
                    self._highlight_selected(
                        self.ax2, "SERSIC_RE", "SERSIC_N")
            elif sn is not None:
                v = np.isfinite(sn) & (sn > 0) & (sn < 10)
                if np.sum(v) > 0:
                    self.ax2.hist(sn[v], bins=30, color='darkorchid',
                                  edgecolor='white', alpha=0.8)
                    self.ax2.set_xlabel("Sérsic n")
                    self.ax2.set_ylabel("Count")
                    self.ax2.set_title("Sérsic Index Distribution")
        elif has_cas:
            conc = self._col_values("CONC")
            if conc is not None:
                v = np.isfinite(conc)
                if np.sum(v) > 0:
                    self.ax2.hist(conc[v], bins=30, color='coral',
                                  edgecolor='white', alpha=0.8)
                    self.ax2.set_xlabel("Concentration")
                    self.ax2.set_ylabel("Count")
                    self.ax2.set_title("Concentration Distribution")

    def _plot_custom(self):
        xcol = self.plot_xcol.get()
        ycol = self.plot_ycol.get()
        ptype = self.plot_type.get()

        if ptype == "histogram":
            vals = self._col_values(xcol)
            if vals is not None:
                v = np.isfinite(vals)
                if np.sum(v) > 0:
                    self.ax1.hist(vals[v], bins=50, color='steelblue',
                                  edgecolor='white', alpha=0.8)
                    self.ax1.set_xlabel(xcol)
                    self.ax1.set_ylabel("Count")
                    self.ax1.set_title(f"Histogram: {xcol}")
                    sv = self._sel_val(xcol)
                    if sv is not None:
                        self.ax1.axvline(sv, color='red', lw=2, ls='--')
        else:
            xvals = self._col_values(xcol)
            yvals = self._col_values(ycol)
            if xvals is not None and yvals is not None:
                v = np.isfinite(xvals) & np.isfinite(yvals)
                if np.sum(v) > 0:
                    self.ax1.scatter(xvals[v], yvals[v], s=5,
                                     alpha=0.4, c='steelblue')
                    self.ax1.set_xlabel(xcol)
                    self.ax1.set_ylabel(ycol)
                    self.ax1.set_title(f"{ycol} vs {xcol}")
                    if "MAG" in ycol.upper():
                        self.ax1.invert_yaxis()
                    self._highlight_selected(self.ax1, xcol, ycol)

            # Bottom: Y histogram
            if yvals is not None and ptype != "histogram":
                v2 = np.isfinite(yvals)
                if np.sum(v2) > 0:
                    self.ax2.hist(yvals[v2], bins=40, color='lightblue',
                                  edgecolor='white', alpha=0.8)
                    self.ax2.set_xlabel(ycol)
                    self.ax2.set_ylabel("Count")

    # ================================================================
    # File Operations & Export
    # ================================================================

    def _open_catalog(self):
        path = filedialog.askopenfilename(
            title="Open Catalog",
            filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")])
        if path:
            self.load_catalog(path)

    def _open_fits(self):
        path = filedialog.askopenfilename(
            title="Open FITS File",
            filetypes=[("FITS files", "*.fits *.fit *.fz"),
                       ("All files", "*.*")])
        if path:
            self.fits_path = path
            self.fits_data = None  # reset cache
            self.fits_var.set(os.path.basename(path))
            self.status_var.set(f"FITS: {path}")
            if self.selected_idx >= 0:
                self._update_plots()

    def _export_catalog(self):
        if not self.rows:
            messagebox.showwarning("Warning", "No data to export.")
            return
        path = filedialog.asksaveasfilename(
            title="Export Catalog",
            defaultextension=".tsv",
            filetypes=[("TSV files", "*.tsv"), ("CSV files", "*.csv")])
        if not path:
            return

        sep = ',' if path.endswith('.csv') else '\t'
        with open(path, 'w') as f:
            f.write(sep.join(self.header) + '\n')
            for row in self.rows:
                f.write(sep.join(row) + '\n')
        self.status_var.set(f"Exported {len(self.rows)} sources to {path}")

    def _export_regions(self):
        """Export catalog as DS9 region file (.reg)."""
        if not self.rows:
            messagebox.showwarning("Warning", "No data to export.")
            return

        xi = self._col_idx("X_IMAGE")
        yi = self._col_idx("Y_IMAGE")
        if xi < 0 or yi < 0:
            messagebox.showerror("Error",
                                  "X_IMAGE/Y_IMAGE columns required.")
            return

        path = filedialog.asksaveasfilename(
            title="Export Regions",
            defaultextension=".reg",
            filetypes=[("DS9 Region", "*.reg")])
        if not path:
            return

        ai = self._col_idx("A_IMAGE")
        bi = self._col_idx("B_IMAGE")
        ti = self._col_idx("THETA_IMAGE")
        ni = self._col_idx("NUMBER")

        with open(path, 'w') as f:
            f.write("# Region file format: DS9 version 4.1\n")
            f.write("global color=green dashlist=8 3 width=1 "
                    "font=\"helvetica 10 normal roman\" select=1 "
                    "highlite=1 dash=0 fixed=0 edit=1 move=1 "
                    "delete=1 include=1 source=1\n")
            f.write("image\n")

            for row in self.rows:
                try:
                    x = float(row[xi])
                    y = float(row[yi])
                except (ValueError, IndexError):
                    continue

                if ai >= 0 and bi >= 0:
                    try:
                        a = float(row[ai]) * 3
                        b = float(row[bi]) * 3
                        theta = float(row[ti]) if ti >= 0 else 0
                        num = row[ni] if ni >= 0 else ""
                        f.write(f"ellipse({x:.2f},{y:.2f},"
                                f"{a:.2f},{b:.2f},{theta:.1f})")
                        if num:
                            f.write(f" # text={{{num}}}")
                        f.write("\n")
                    except (ValueError, IndexError):
                        f.write(f"circle({x:.2f},{y:.2f},5)\n")
                else:
                    f.write(f"circle({x:.2f},{y:.2f},5)\n")

        self.status_var.set(
            f"Exported {len(self.rows)} regions to {path}")

    def _export_fits_table(self):
        """Export catalog as FITS binary table (requires astropy)."""
        if not self.rows:
            messagebox.showwarning("Warning", "No data to export.")
            return

        path = filedialog.asksaveasfilename(
            title="Export FITS Table",
            defaultextension=".fits",
            filetypes=[("FITS files", "*.fits")])
        if not path:
            return

        try:
            from astropy.table import Table
            t = Table()
            for i, col in enumerate(self.header):
                vals = [row[i] if i < len(row) else "" for row in self.rows]
                # Try numeric conversion
                try:
                    t[col] = np.array([
                        float(v) if v not in ("", "-99") else np.nan
                        for v in vals])
                except ValueError:
                    t[col] = vals
            t.write(path, format='fits', overwrite=True)
            self.status_var.set(
                f"Exported FITS table to {path}")
        except ImportError:
            messagebox.showerror("Error",
                                  "astropy required for FITS export")
        except Exception as e:
            messagebox.showerror("Error", f"Export failed: {e}")

    def _on_close(self):
        self._save_settings()
        self.root.destroy()


def main():
    parser = argparse.ArgumentParser(
        description='DS9 Analysis Viewer - Standalone GUI')
    parser.add_argument('--fits', default=None,
                        help='Path to FITS image')
    parser.add_argument('--catalog', default=None,
                        help='Path to catalog TSV')
    args = parser.parse_args()

    root = tk.Tk()
    app = AnalysisGUI(root, fits_path=args.fits,
                       catalog_path=args.catalog)
    root.mainloop()


if __name__ == '__main__':
    main()
