#!/usr/bin/env python3
"""DS9 Analysis Viewer - Standalone multi-analysis GUI.

Provides an integrated analysis environment with:
- Left panel: matplotlib plots (SED, photo-z, morphology diagrams)
- Right panel: sortable catalog table
- Analysis runners: Photo-z, SED Fitting, Morphology, Star/Galaxy,
  CAS/Gini/M20, Sérsic, Bulge+Disk, PSF Phot, Custom Plot

Usage:
    ds9_analysis_gui.py [--fits image.fits] [--catalog catalog.tsv]
"""

import sys
import os
import tempfile
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

# Morphology class names (Galaxy10 DECaLS)
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

# Default columns to display
DEFAULT_COLS = [
    "NUMBER", "X_IMAGE", "Y_IMAGE", "MAG_AUTO", "FLUX_AUTO",
    "FWHM_IMAGE", "A_IMAGE", "B_IMAGE", "THETA_IMAGE",
]


class AnalysisGUI:
    """Standalone analysis viewer with plots + catalog table."""

    def __init__(self, root, fits_path=None, catalog_path=None):
        self.root = root
        self.root.title("DS9 Analysis Viewer")
        self.root.geometry("1400x800")
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

        self.fits_path = fits_path
        self.catalog_path = catalog_path

        # Data
        self.header = []
        self.rows = []          # list of lists
        self.filtered_indices = []  # indices into self.rows
        self.selected_idx = -1  # index in self.rows
        self.sort_col = None
        self.sort_asc = True

        # State
        self.mode = "overview"
        self.n_workers = max(1, os.cpu_count() - 1)

        # SED settings
        self.sed_backend = tk.StringVar(value="auto")
        self.sed_bands = tk.StringVar(value="g,r,i,z")
        self.sed_magcols = tk.StringVar(value="")

        # Custom plot settings
        self.plot_xcol = tk.StringVar(value="MAG_AUTO")
        self.plot_ycol = tk.StringVar(value="FLUX_AUTO")
        self.plot_type = tk.StringVar(value="scatter")

        self._build_ui()

        if catalog_path:
            self.load_catalog(catalog_path)

    # ================================================================
    # UI Construction
    # ================================================================

    def _build_ui(self):
        self._build_menu()
        self._build_toolbar()

        # Main paned window
        self.paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        self.paned.pack(fill=tk.BOTH, expand=True)

        self._build_plot_area()
        self._build_table_area()

        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(self.root, textvariable=self.status_var,
                  relief=tk.SUNKEN, anchor=tk.W).pack(
                      side=tk.BOTTOM, fill=tk.X)

    def _build_menu(self):
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Open Catalog...",
                              command=self._open_catalog)
        file_menu.add_command(label="Open FITS...",
                              command=self._open_fits)
        file_menu.add_separator()
        file_menu.add_command(label="Export Catalog...",
                              command=self._export_catalog)
        file_menu.add_separator()
        file_menu.add_command(label="Quit", command=self._on_close)
        menubar.add_cascade(label="File", menu=file_menu)

        # Analysis menu
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
        analysis_menu.add_command(label="Crowded Photometry",
                                  command=self.run_crowded_phot)
        menubar.add_cascade(label="Analysis", menu=analysis_menu)

        # View menu
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

        # SED backend selector
        ttk.Label(toolbar, text="SPS Backend:").pack(side=tk.LEFT, padx=2)
        ttk.Combobox(
            toolbar, textvariable=self.sed_backend, width=12,
            values=["auto", "analytic", "fsps", "bagpipes",
                    "prospector", "cigale", "dense_basis"],
            state="readonly",
        ).pack(side=tk.LEFT, padx=2)

        ttk.Separator(toolbar, orient=tk.VERTICAL).pack(
            side=tk.LEFT, fill=tk.Y, padx=5)

        # FITS path display
        ttk.Label(toolbar, text="FITS:").pack(side=tk.LEFT, padx=2)
        self.fits_var = tk.StringVar(
            value=os.path.basename(self.fits_path) if self.fits_path else "")
        ttk.Label(toolbar, textvariable=self.fits_var, width=30,
                  relief=tk.SUNKEN).pack(side=tk.LEFT, padx=2)

        # Source count
        self.count_var = tk.StringVar(value="0 sources")
        ttk.Label(toolbar, textvariable=self.count_var).pack(
            side=tk.RIGHT, padx=10)

    def _build_plot_area(self):
        plot_frame = ttk.Frame(self.paned)
        self.paned.add(plot_frame, weight=1)

        self.fig = Figure(figsize=(6, 7), dpi=96)
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)
        self.fig.tight_layout(pad=3.0)

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

        # Custom plot controls (shown when in plot mode)
        self.plot_controls = ttk.Frame(table_frame)
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
    # Data Loading
    # ================================================================

    def load_catalog(self, path):
        """Load a TSV catalog file."""
        try:
            with open(path, 'r') as f:
                lines = f.readlines()
        except Exception as e:
            messagebox.showerror("Error", f"Cannot load catalog:\n{e}")
            return

        if not lines:
            return

        # Skip comment lines
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
            # Pad to header length
            while len(vals) < len(self.header):
                vals.append("")
            self.rows.append(vals)

        self.catalog_path = path
        self.filtered_indices = list(range(len(self.rows)))
        self.selected_idx = -1

        self._rebuild_table()
        self._update_plot_col_lists()
        self.count_var.set(f"{len(self.rows)} sources")
        self.status_var.set(f"Loaded {len(self.rows)} sources from "
                            f"{os.path.basename(path)}")
        self._update_plots()

    def _rebuild_table(self):
        """Rebuild treeview with current header/data."""
        # Clear existing
        self.tree.delete(*self.tree.get_children())

        # Configure columns
        self.tree["columns"] = self.header
        for col in self.header:
            w = 90 if col != "NUMBER" else 60
            self.tree.column(col, width=w, minwidth=50, anchor=tk.E)
            self.tree.heading(
                col, text=col,
                command=lambda c=col: self._sort_by_column(c))

        # Insert rows
        for i in self.filtered_indices:
            self.tree.insert("", tk.END, iid=str(i),
                             values=self.rows[i])

    def _apply_filter(self):
        """Filter table rows by text match."""
        text = self.filter_var.get().strip().lower()
        if not text:
            self.filtered_indices = list(range(len(self.rows)))
        else:
            self.filtered_indices = []
            for i, row in enumerate(self.rows):
                if any(text in v.lower() for v in row):
                    self.filtered_indices.append(i)

        # Rebuild
        self.tree.delete(*self.tree.get_children())
        for i in self.filtered_indices:
            self.tree.insert("", tk.END, iid=str(i),
                             values=self.rows[i])
        self.count_var.set(
            f"{len(self.filtered_indices)}/{len(self.rows)} sources")

    def _sort_by_column(self, col):
        """Sort table by clicking column header."""
        if self.sort_col == col:
            self.sort_asc = not self.sort_asc
        else:
            self.sort_col = col
            self.sort_asc = True

        col_idx = self.header.index(col)

        def sort_key(i):
            v = self.rows[i][col_idx]
            try:
                return float(v)
            except (ValueError, TypeError):
                return v

        self.filtered_indices.sort(key=sort_key, reverse=not self.sort_asc)

        self.tree.delete(*self.tree.get_children())
        for i in self.filtered_indices:
            self.tree.insert("", tk.END, iid=str(i),
                             values=self.rows[i])

    def _on_select(self, event):
        """Handle row selection in table."""
        sel = self.tree.selection()
        if sel:
            self.selected_idx = int(sel[0])
            self._update_plots()

    def _col_idx(self, name):
        """Get column index by name, -1 if not found."""
        try:
            return self.header.index(name)
        except ValueError:
            return -1

    def _col_values(self, name, as_float=True):
        """Get all values for a column as numpy array."""
        idx = self._col_idx(name)
        if idx < 0:
            return None
        vals = []
        for row in self.rows:
            try:
                v = float(row[idx]) if as_float else row[idx]
                vals.append(v)
            except (ValueError, IndexError):
                vals.append(np.nan if as_float else "")
        return np.array(vals) if as_float else vals

    def _has_col(self, name):
        return name in self.header

    def _update_plot_col_lists(self):
        """Update column lists for custom plot comboboxes."""
        self.xcol_cb['values'] = self.header
        self.ycol_cb['values'] = self.header

    # ================================================================
    # Analysis Runners
    # ================================================================

    def _check_ready(self):
        """Check catalog and FITS are loaded."""
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
        """Save current catalog to temp file."""
        tmpdir = os.path.join(os.path.expanduser("~"), ".ds9")
        os.makedirs(tmpdir, exist_ok=True)
        path = os.path.join(tmpdir, f"{suffix}_gui_catalog.tsv")
        with open(path, 'w') as f:
            f.write('\t'.join(self.header) + '\n')
            for row in self.rows:
                f.write('\t'.join(row) + '\n')
        return path

    def _find_script(self, name):
        """Find a CLI script in the same directory."""
        path = os.path.join(_script_dir, name)
        if os.path.exists(path):
            return path
        return None

    def _run_script(self, script_name, extra_args, output_columns,
                    status_msg):
        """Run an analysis script and merge results into catalog.

        Runs in a background thread to keep GUI responsive.
        """
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
                cmd.extend(["--n-workers", str(self.n_workers)])

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

                # Parse and merge
                self.root.after(0, lambda: self._merge_results(
                    stdout, output_columns, status_msg.replace("...", " done")))

            except subprocess.TimeoutExpired:
                self.root.after(0, lambda: self._show_error(
                    f"{script_name} timed out"))
            except Exception as e:
                self.root.after(0, lambda: self._show_error(str(e)))

        threading.Thread(target=run, daemon=True).start()

    def _merge_results(self, tsv_data, col_names, done_msg):
        """Merge TSV result columns into current catalog."""
        lines = tsv_data.strip().split('\n')

        # Find header (skip comment lines)
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

        # Build lookup: NUMBER -> row values
        for line in lines[header_idx + 1:]:
            vals = line.strip().split('\t')
            if len(vals) > num_col:
                res_rows[vals[num_col]] = vals

        # Get column indices in result
        col_indices = {}
        for cname in col_names:
            for i, h in enumerate(res_header):
                if h.strip() == cname:
                    col_indices[cname] = i
                    break

        if not col_indices:
            self.status_var.set("No matching columns in output")
            return

        # Add new columns to header if needed
        my_num_col = self._col_idx("NUMBER")
        for cname in col_names:
            if cname not in self.header:
                self.header.append(cname)
                for row in self.rows:
                    row.append("-99")

        # Merge data
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
        self._set_mode(self.mode)  # refresh plots
        self.status_var.set(done_msg)

    def _show_error(self, msg):
        self.status_var.set("Error")
        messagebox.showerror("Analysis Error", msg)

    # Individual analysis runners

    def run_photoz(self):
        if not self._check_fits():
            return
        self._run_script(
            "ds9_photo_z.py", [],
            ["PHOTO_Z", "PHOTO_Z_ERR", "PHOTO_Z_Q68", "PHOTO_Z_OUTLIER"],
            "Running Photo-z (AI)...")
        self._set_mode("photoz")

    def run_sed_fit(self):
        if not self._check_fits():
            return
        args = ["--backend", self.sed_backend.get()]
        bands = self.sed_bands.get().strip()
        if bands:
            args.extend(["--bands", bands])
        magcols = self.sed_magcols.get().strip()
        if magcols:
            args.extend(["--mag-columns", magcols])
        self._run_script(
            "ds9_sed_fit.py", args,
            ["LOG_MASS", "LOG_MASS_ERR", "LOG_AGE", "LOG_AGE_ERR",
             "LOG_Z", "AV", "SFR"],
            f"Running SED Fitting ({self.sed_backend.get()})...")
        self._set_mode("sed")

    def run_morphology(self):
        if not self._check_fits():
            return
        self._run_script(
            "ds9_galaxy_morph.py", [],
            ["MORPH_TYPE", "MORPH_CONF"],
            "Running Galaxy Morphology (CNN)...")
        self._set_mode("classify")

    def run_starfinder(self):
        if not self._check_fits():
            return
        self._run_script(
            "ds9_star_finder.py", [],
            ["AI_STAR", "AI_STAR_CONF"],
            "Running Star/Galaxy Classification (CNN)...")
        self._set_mode("classify")

    def run_morphometry(self):
        if not self._check_fits():
            return
        self._run_script(
            "ds9_morphometry.py", [],
            ["CONC", "ASYM", "GINI", "M20", "R_PETRO"],
            "Running CAS/Gini/M20...")
        self._set_mode("structure")

    def run_sersic(self):
        if not self._check_fits():
            return
        self._run_script(
            "ds9_sersic.py", [],
            ["SERSIC_N", "SERSIC_RE", "SERSIC_MU_E", "SERSIC_CHI2"],
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
        self._run_script(
            "ds9_psf_phot.py", [],
            ["MAG_PSF", "MAGERR_PSF", "PSF_CHI2"],
            "Running PSF Photometry...")

    def run_crowded_phot(self):
        if not self._check_fits():
            return
        self._run_script(
            "ds9_crowded_phot.py", [],
            ["MAG_CROWD", "MAGERR_CROWD"],
            "Running Crowded Field Photometry...")

    # ================================================================
    # Mode & Plot Management
    # ================================================================

    def _set_mode(self, mode):
        self.mode = mode
        self.mode_var.set(mode)

        # Show/hide custom plot controls
        if mode == "plot":
            self.plot_controls.pack(side=tk.TOP, fill=tk.X, padx=2, pady=2)
        else:
            self.plot_controls.pack_forget()

        self._update_plots()

    def _update_plots(self):
        self.ax1.clear()
        self.ax2.clear()

        if not self.rows:
            self.ax1.text(0.5, 0.5, "Load a catalog to begin",
                          transform=self.ax1.transAxes, ha='center',
                          fontsize=12, color='gray')
            self.canvas.draw()
            return

        dispatch = {
            "overview":  self._plot_overview,
            "photoz":    self._plot_photoz,
            "sed":       self._plot_sed,
            "classify":  self._plot_classification,
            "structure": self._plot_structure,
            "plot":      self._plot_custom,
        }
        func = dispatch.get(self.mode, self._plot_overview)
        try:
            func()
        except Exception as e:
            self.ax1.text(0.5, 0.5, f"Plot error: {e}",
                          transform=self.ax1.transAxes, ha='center',
                          fontsize=10, color='red', wrap=True)

        self.fig.tight_layout(pad=2.5)
        self.canvas.draw()

    # ================================================================
    # Plot Functions
    # ================================================================

    def _plot_overview(self):
        """Overview: magnitude histogram + source positions."""
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
                if self.selected_idx >= 0:
                    idx = self._col_idx("MAG_AUTO")
                    try:
                        sv = float(self.rows[self.selected_idx][idx])
                        self.ax1.axvline(sv, color='red', lw=2, ls='--',
                                          label=f"Selected: {sv:.1f}")
                        self.ax1.legend(fontsize=8)
                    except (ValueError, IndexError):
                        pass
        else:
            self.ax1.text(0.5, 0.5, "No MAG_AUTO column",
                          transform=self.ax1.transAxes, ha='center')

        # Source positions
        x = self._col_values("X_IMAGE")
        y = self._col_values("Y_IMAGE")
        if x is not None and y is not None:
            valid = np.isfinite(x) & np.isfinite(y)
            self.ax2.scatter(x[valid], y[valid], s=2, alpha=0.5,
                             c='steelblue')
            self.ax2.set_xlabel("X_IMAGE")
            self.ax2.set_ylabel("Y_IMAGE")
            self.ax2.set_title("Source Positions")
            self.ax2.set_aspect('equal')
            if self.selected_idx >= 0:
                xi = self._col_idx("X_IMAGE")
                yi = self._col_idx("Y_IMAGE")
                try:
                    sx = float(self.rows[self.selected_idx][xi])
                    sy = float(self.rows[self.selected_idx][yi])
                    self.ax2.scatter([sx], [sy], s=80, c='red',
                                     marker='+', linewidths=2,
                                     zorder=10)
                except (ValueError, IndexError):
                    pass

    def _plot_photoz(self):
        """Photo-z: z distribution + color-z diagram."""
        pz = self._col_values("PHOTO_Z")
        if pz is None:
            self.ax1.text(0.5, 0.5,
                          "Run Photo-z first\n(Analysis > Photo-z)",
                          transform=self.ax1.transAxes, ha='center',
                          fontsize=12, color='gray')
            self.ax2.set_visible(False)
            return

        self.ax2.set_visible(True)
        valid = np.isfinite(pz) & (pz > 0) & (pz < 10)
        pz_valid = pz[valid]

        # Histogram
        if len(pz_valid) > 0:
            self.ax1.hist(pz_valid, bins=50, color='darkorange',
                          edgecolor='white', alpha=0.8)
            self.ax1.set_xlabel("Photo-z")
            self.ax1.set_ylabel("Count")
            self.ax1.set_title(f"Photometric Redshift Distribution "
                                f"(N={len(pz_valid)})")
            med = np.median(pz_valid)
            self.ax1.axvline(med, color='red', ls='--',
                              label=f"median = {med:.2f}")
            self.ax1.legend(fontsize=8)

            if self.selected_idx >= 0:
                idx = self._col_idx("PHOTO_Z")
                try:
                    sz = float(self.rows[self.selected_idx][idx])
                    self.ax1.axvline(sz, color='cyan', lw=2,
                                      label=f"Selected: z={sz:.3f}")
                    self.ax1.legend(fontsize=8)
                except (ValueError, IndexError):
                    pass

        # Photo-z vs MAG_AUTO
        mag = self._col_values("MAG_AUTO")
        if mag is not None:
            both = valid & np.isfinite(mag) & (mag > 0) & (mag < 90)
            if np.sum(both) > 0:
                sc = self.ax2.scatter(pz[both], mag[both], s=5, alpha=0.4,
                                      c='darkorange')
                self.ax2.set_xlabel("Photo-z")
                self.ax2.set_ylabel("MAG_AUTO")
                self.ax2.set_title("Magnitude vs Redshift")
                self.ax2.invert_yaxis()

                if self.selected_idx >= 0:
                    try:
                        sz = float(self.rows[self.selected_idx][
                            self._col_idx("PHOTO_Z")])
                        sm = float(self.rows[self.selected_idx][
                            self._col_idx("MAG_AUTO")])
                        self.ax2.scatter([sz], [sm], s=80, c='red',
                                          marker='*', zorder=10)
                    except (ValueError, IndexError):
                        pass
        else:
            pz_err = self._col_values("PHOTO_Z_ERR")
            if pz_err is not None and len(pz_valid) > 0:
                pz_err_v = pz_err[valid]
                self.ax2.scatter(pz_valid, pz_err_v, s=5, alpha=0.4,
                                 c='darkorange')
                self.ax2.set_xlabel("Photo-z")
                self.ax2.set_ylabel("Photo-z Error")
                self.ax2.set_title("Redshift Uncertainty")

    def _plot_sed(self):
        """SED Fitting: best-fit SED + mass-age diagram."""
        lm = self._col_values("LOG_MASS")
        if lm is None:
            self.ax1.text(0.5, 0.5,
                          "Run SED Fitting first\n(Analysis > SED Fitting)",
                          transform=self.ax1.transAxes, ha='center',
                          fontsize=12, color='gray')
            self.ax2.set_visible(False)
            return

        self.ax2.set_visible(True)
        la = self._col_values("LOG_AGE")
        valid_m = np.isfinite(lm) & (lm > 0)

        # Top: SED for selected source (or mass histogram)
        if self.selected_idx >= 0 and self._has_sed_data():
            self._plot_source_sed(self.selected_idx)
        else:
            # Mass histogram
            if np.sum(valid_m) > 0:
                self.ax1.hist(lm[valid_m], bins=40, color='mediumpurple',
                              edgecolor='white', alpha=0.8)
                self.ax1.set_xlabel("log(M*/M$_\\odot$)")
                self.ax1.set_ylabel("Count")
                self.ax1.set_title("Stellar Mass Distribution")
                if self.selected_idx >= 0:
                    try:
                        sv = float(self.rows[self.selected_idx][
                            self._col_idx("LOG_MASS")])
                        if np.isfinite(sv) and sv > 0:
                            self.ax1.axvline(sv, color='red', lw=2, ls='--')
                    except (ValueError, IndexError):
                        pass

        # Bottom: Mass-Age diagram
        if la is not None:
            valid_a = np.isfinite(la) & (la > 0)
            both = valid_m & valid_a
            if np.sum(both) > 0:
                av = self._col_values("AV")
                if av is not None:
                    c = np.where(np.isfinite(av), av, 0)
                    sc = self.ax2.scatter(la[both], lm[both], s=8,
                                          c=c[both], cmap='RdYlBu_r',
                                          alpha=0.6, vmin=0, vmax=3)
                    self.fig.colorbar(sc, ax=self.ax2, label="Av",
                                      shrink=0.8, pad=0.02)
                else:
                    self.ax2.scatter(la[both], lm[both], s=8,
                                     alpha=0.5, c='mediumpurple')

                self.ax2.set_xlabel("log(Age/yr)")
                self.ax2.set_ylabel("log(M*/M$_\\odot$)")
                self.ax2.set_title("Mass-Age Diagram")

                if self.selected_idx >= 0:
                    try:
                        sa = float(self.rows[self.selected_idx][
                            self._col_idx("LOG_AGE")])
                        sm = float(self.rows[self.selected_idx][
                            self._col_idx("LOG_MASS")])
                        self.ax2.scatter([sa], [sm], s=100, c='red',
                                          marker='*', zorder=10,
                                          edgecolors='black')
                    except (ValueError, IndexError):
                        pass

    def _has_sed_data(self):
        """Check if selected source has SED fit results."""
        if self.selected_idx < 0:
            return False
        idx = self._col_idx("LOG_MASS")
        if idx < 0:
            return False
        try:
            v = float(self.rows[self.selected_idx][idx])
            return np.isfinite(v) and v > 0
        except (ValueError, IndexError):
            return False

    def _plot_source_sed(self, row_idx):
        """Plot SED for a single source using best-fit parameters."""
        # Get observed magnitudes
        mag_cols = [c for c in self.header
                    if c.startswith("MAG_") and c != "MAG_PSF"
                    and c != "MAG_CROWD"]
        if not mag_cols:
            mag_cols = ["MAG_AUTO"]

        wavelengths = []
        obs_mags = []
        obs_labels = []

        for mc in mag_cols:
            idx = self._col_idx(mc)
            if idx < 0:
                continue
            try:
                mag_val = float(self.rows[row_idx][idx])
                if mag_val <= 0 or mag_val >= 90:
                    continue
            except (ValueError, IndexError):
                continue

            # Try to get wavelength from filter name
            fname = mc.replace("MAG_AUTO_", "").replace("MAG_", "")
            if fname in FILTER_WAVELENGTHS:
                wavelengths.append(FILTER_WAVELENGTHS[fname])
                obs_mags.append(mag_val)
                obs_labels.append(fname)

        if len(wavelengths) >= 2:
            wl = np.array(wavelengths)
            om = np.array(obs_mags)
            sort = np.argsort(wl)
            wl, om = wl[sort], om[sort]
            labels_sorted = [obs_labels[i] for i in sort]

            self.ax1.plot(wl / 1e4, om, 'o-', color='steelblue',
                          markersize=6, label="Observed")
            for i, lbl in enumerate(labels_sorted):
                self.ax1.annotate(lbl, (wl[i] / 1e4, om[i]),
                                   fontsize=6, ha='center',
                                   va='bottom', rotation=45)

            # Generate model SED if parameters available
            try:
                lm = float(self.rows[row_idx][self._col_idx("LOG_MASS")])
                la = float(self.rows[row_idx][self._col_idx("LOG_AGE")])
                lz = float(self.rows[row_idx][self._col_idx("LOG_Z")])
                av = float(self.rows[row_idx][self._col_idx("AV")])
                pz_idx = self._col_idx("PHOTO_Z")
                z = float(self.rows[row_idx][pz_idx]) if pz_idx >= 0 else 0.1

                from sed_fit.grid.generate import simple_sed_model
                model_wl = np.linspace(wl.min() * 0.8, wl.max() * 1.2, 100)
                model_mags = simple_sed_model(
                    z, lm, la, lz, av, 9.0, model_wl)
                self.ax1.plot(model_wl / 1e4, model_mags, '-',
                              color='tomato', alpha=0.7, label="Model")
            except Exception:
                pass

            self.ax1.invert_yaxis()
            self.ax1.set_xlabel("Wavelength (μm)")
            self.ax1.set_ylabel("AB mag")
            num = self.rows[row_idx][self._col_idx("NUMBER")] \
                if self._col_idx("NUMBER") >= 0 else "?"
            self.ax1.set_title(f"SED - Source #{num}")
            self.ax1.legend(fontsize=8)
        else:
            # Fall back to parameter summary text
            params = {}
            for pname in ["LOG_MASS", "LOG_AGE", "LOG_Z", "AV", "SFR"]:
                idx = self._col_idx(pname)
                if idx >= 0:
                    try:
                        params[pname] = float(self.rows[row_idx][idx])
                    except (ValueError, IndexError):
                        pass

            if params:
                text = "SED Fit Results:\n"
                for k, v in params.items():
                    text += f"  {k} = {v:.2f}\n"
                self.ax1.text(0.1, 0.5, text,
                              transform=self.ax1.transAxes,
                              fontsize=11, family='monospace',
                              verticalalignment='center')
                self.ax1.set_title(f"Source #{row_idx}")
            else:
                self.ax1.text(0.5, 0.5, "Select a source with SED data",
                              transform=self.ax1.transAxes, ha='center')

    def _plot_classification(self):
        """Classification: morphology type + star/galaxy."""
        has_morph = self._has_col("MORPH_TYPE")
        has_star = self._has_col("AI_STAR")

        if not has_morph and not has_star:
            self.ax1.text(0.5, 0.5,
                          "Run Classification first\n"
                          "(Analysis > Galaxy Morphology or Star/Galaxy)",
                          transform=self.ax1.transAxes, ha='center',
                          fontsize=12, color='gray')
            self.ax2.set_visible(False)
            return

        self.ax2.set_visible(True)

        # Top: Morphology distribution
        if has_morph:
            mt = self._col_values("MORPH_TYPE", as_float=False)
            types = [t for t in mt if t and t != "-99" and t != ""]
            if types:
                from collections import Counter
                counts = Counter(types)
                labels = list(counts.keys())
                values = list(counts.values())

                # Abbreviate long names
                short = [l[:15] for l in labels]
                colors = matplotlib.colormaps['tab10'](
                    np.linspace(0, 1, len(labels)))

                bars = self.ax1.barh(range(len(labels)), values,
                                      color=colors, edgecolor='white')
                self.ax1.set_yticks(range(len(labels)))
                self.ax1.set_yticklabels(short, fontsize=8)
                self.ax1.set_xlabel("Count")
                self.ax1.set_title(f"Morphology Distribution (N={len(types)})")

                if self.selected_idx >= 0:
                    idx = self._col_idx("MORPH_TYPE")
                    sel_type = self.rows[self.selected_idx][idx]
                    if sel_type in labels:
                        bar_idx = labels.index(sel_type)
                        bars[bar_idx].set_edgecolor('red')
                        bars[bar_idx].set_linewidth(3)
            else:
                self.ax1.text(0.5, 0.5, "No morphology results",
                              transform=self.ax1.transAxes, ha='center')
        else:
            self.ax1.text(0.5, 0.5, "No morphology data",
                          transform=self.ax1.transAxes, ha='center')

        # Bottom: Star/Galaxy classification
        if has_star:
            ai_star = self._col_values("AI_STAR")
            ai_conf = self._col_values("AI_STAR_CONF")
            valid = np.isfinite(ai_star)

            if np.sum(valid) > 0:
                n_star = np.sum(ai_star[valid] > 0.5)
                n_gal = np.sum(ai_star[valid] <= 0.5)
                self.ax2.pie([n_star, n_gal],
                              labels=[f"Star ({n_star})",
                                      f"Galaxy ({n_gal})"],
                              colors=['cyan', 'gold'],
                              autopct='%1.0f%%', startangle=90)
                self.ax2.set_title("Star/Galaxy Classification")
        elif has_morph:
            # Show confidence distribution instead
            mc = self._col_values("MORPH_CONF")
            if mc is not None:
                valid = np.isfinite(mc)
                if np.sum(valid) > 0:
                    self.ax2.hist(mc[valid], bins=30, color='teal',
                                  edgecolor='white', alpha=0.8)
                    self.ax2.set_xlabel("Classification Confidence")
                    self.ax2.set_ylabel("Count")
                    self.ax2.set_title("Morphology Confidence")

    def _plot_structure(self):
        """Structure: Gini-M20 + Sérsic n distribution."""
        has_gini = self._has_col("GINI")
        has_sersic = self._has_col("SERSIC_N")

        if not has_gini and not has_sersic:
            self.ax1.text(0.5, 0.5,
                          "Run CAS/Gini/M20 or Sérsic first\n"
                          "(Analysis menu)",
                          transform=self.ax1.transAxes, ha='center',
                          fontsize=12, color='gray')
            self.ax2.set_visible(False)
            return

        self.ax2.set_visible(True)

        # Top: Gini-M20 diagram
        if has_gini:
            gini = self._col_values("GINI")
            m20 = self._col_values("M20")
            if gini is not None and m20 is not None:
                valid = (np.isfinite(gini) & np.isfinite(m20)
                         & (gini > 0))
                if np.sum(valid) > 0:
                    self.ax1.scatter(m20[valid], gini[valid], s=8,
                                     alpha=0.5, c='teal')
                    self.ax1.set_xlabel("M$_{20}$")
                    self.ax1.set_ylabel("Gini")
                    self.ax1.set_title("Gini-M20 Diagram")

                    # Merger boundary (Lotz+2008)
                    x_line = np.linspace(-3, 0, 50)
                    y_line = -0.14 * x_line + 0.33
                    self.ax1.plot(x_line, y_line, 'r--', lw=1,
                                  alpha=0.5, label="Merger (Lotz+08)")
                    self.ax1.legend(fontsize=7)

                    if self.selected_idx >= 0:
                        try:
                            sg = float(self.rows[self.selected_idx][
                                self._col_idx("GINI")])
                            sm = float(self.rows[self.selected_idx][
                                self._col_idx("M20")])
                            self.ax1.scatter([sm], [sg], s=100, c='red',
                                              marker='*', zorder=10)
                        except (ValueError, IndexError):
                            pass
        elif has_sersic:
            # CAS diagram if no Gini
            conc = self._col_values("CONC")
            asym = self._col_values("ASYM")
            if conc is not None and asym is not None:
                valid = np.isfinite(conc) & np.isfinite(asym)
                if np.sum(valid) > 0:
                    self.ax1.scatter(asym[valid], conc[valid], s=8,
                                     alpha=0.5, c='coral')
                    self.ax1.set_xlabel("Asymmetry")
                    self.ax1.set_ylabel("Concentration")
                    self.ax1.set_title("C-A Diagram")

        # Bottom: Sérsic n distribution or n vs Re
        if has_sersic:
            sn = self._col_values("SERSIC_N")
            sre = self._col_values("SERSIC_RE")

            if sre is not None:
                valid = (np.isfinite(sn) & np.isfinite(sre)
                         & (sn > 0) & (sre > 0))
                if np.sum(valid) > 0:
                    self.ax2.scatter(sre[valid], sn[valid], s=8,
                                     alpha=0.5, c='darkorchid')
                    self.ax2.set_xlabel("R$_e$ (pix)")
                    self.ax2.set_ylabel("Sérsic n")
                    self.ax2.set_title("Sérsic Index vs Effective Radius")
                    self.ax2.set_xscale('log')
                    self.ax2.axhline(2.5, color='gray', ls='--',
                                      alpha=0.5, label="n=2.5")
                    self.ax2.legend(fontsize=7)

                    if self.selected_idx >= 0:
                        try:
                            ssn = float(self.rows[self.selected_idx][
                                self._col_idx("SERSIC_N")])
                            ssr = float(self.rows[self.selected_idx][
                                self._col_idx("SERSIC_RE")])
                            self.ax2.scatter([ssr], [ssn], s=100,
                                              c='red', marker='*',
                                              zorder=10)
                        except (ValueError, IndexError):
                            pass
            else:
                valid = np.isfinite(sn) & (sn > 0) & (sn < 10)
                if np.sum(valid) > 0:
                    self.ax2.hist(sn[valid], bins=30,
                                  color='darkorchid', edgecolor='white',
                                  alpha=0.8)
                    self.ax2.set_xlabel("Sérsic n")
                    self.ax2.set_ylabel("Count")
                    self.ax2.set_title("Sérsic Index Distribution")
        else:
            # Show concentration histogram if no Sérsic
            conc = self._col_values("CONC")
            if conc is not None:
                valid = np.isfinite(conc)
                if np.sum(valid) > 0:
                    self.ax2.hist(conc[valid], bins=30, color='coral',
                                  edgecolor='white', alpha=0.8)
                    self.ax2.set_xlabel("Concentration")
                    self.ax2.set_ylabel("Count")
                    self.ax2.set_title("Concentration Distribution")

    def _plot_custom(self):
        """Custom user-selected scatter or histogram plot."""
        xcol = self.plot_xcol.get()
        ycol = self.plot_ycol.get()
        ptype = self.plot_type.get()

        if ptype == "histogram":
            vals = self._col_values(xcol)
            if vals is not None:
                valid = np.isfinite(vals)
                if np.sum(valid) > 0:
                    self.ax1.hist(vals[valid], bins=50,
                                  color='steelblue', edgecolor='white',
                                  alpha=0.8)
                    self.ax1.set_xlabel(xcol)
                    self.ax1.set_ylabel("Count")
                    self.ax1.set_title(f"Histogram: {xcol}")

                    if self.selected_idx >= 0:
                        idx = self._col_idx(xcol)
                        try:
                            sv = float(self.rows[self.selected_idx][idx])
                            if np.isfinite(sv):
                                self.ax1.axvline(sv, color='red',
                                                  lw=2, ls='--')
                        except (ValueError, IndexError):
                            pass
            self.ax2.set_visible(False)
        else:
            # Scatter plot
            xvals = self._col_values(xcol)
            yvals = self._col_values(ycol)
            if xvals is not None and yvals is not None:
                valid = np.isfinite(xvals) & np.isfinite(yvals)
                if np.sum(valid) > 0:
                    self.ax1.scatter(xvals[valid], yvals[valid], s=5,
                                     alpha=0.4, c='steelblue')
                    self.ax1.set_xlabel(xcol)
                    self.ax1.set_ylabel(ycol)
                    self.ax1.set_title(f"{ycol} vs {xcol}")

                    # Invert y-axis for magnitudes
                    if "MAG" in ycol.upper():
                        self.ax1.invert_yaxis()

                    if self.selected_idx >= 0:
                        xi = self._col_idx(xcol)
                        yi = self._col_idx(ycol)
                        try:
                            sx = float(self.rows[self.selected_idx][xi])
                            sy = float(self.rows[self.selected_idx][yi])
                            self.ax1.scatter([sx], [sy], s=80, c='red',
                                              marker='*', zorder=10)
                        except (ValueError, IndexError):
                            pass

            # Bottom: X histogram
            if xvals is not None:
                valid = np.isfinite(xvals)
                if np.sum(valid) > 0:
                    self.ax2.hist(xvals[valid], bins=40,
                                  color='lightblue', edgecolor='white',
                                  alpha=0.8)
                    self.ax2.set_xlabel(xcol)
                    self.ax2.set_ylabel("Count")

    # ================================================================
    # File Operations
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
            self.fits_var.set(os.path.basename(path))
            self.status_var.set(f"FITS: {path}")

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

    def _on_close(self):
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
