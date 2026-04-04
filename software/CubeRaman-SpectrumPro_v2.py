import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from scipy.signal import savgol_filter, find_peaks
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.optimize import curve_fit
import os

# ── Palette ────────────────────────────────────────────────────────────────────
BG      = "#1e1e2e"
SURFACE = "#313244"
OVERLAY = "#45475a"
TEXT    = "#cdd6f4"
BLUE    = "#89b4fa"
RED     = "#f38ba8"
GREEN   = "#a6e3a1"
YELLOW  = "#f9e2af"
MAUVE   = "#cba6f7"
PEACH   = "#fab387"

# ═══════════════════════════════════════════════════════════════════════════════
# DSP Algorithms
# ═══════════════════════════════════════════════════════════════════════════════

def remove_spikes(y, threshold=5.0):
    """
    Modified Z-score spike removal on the first derivative.
    Returns (cleaned_y, n_spikes_removed).
    """
    dy   = np.diff(y.astype(float))
    med  = np.median(dy)
    mad  = np.median(np.abs(dy - med))
    zscore = np.abs(dy - med) / (1.4826 * mad + 1e-12)

    spike_starts = np.where(zscore > threshold)[0]
    y_clean = y.copy().astype(float)

    removed = 0
    i = 0
    while i < len(spike_starts):
        # group consecutive spike indices into one event
        s = spike_starts[i]
        e = s + 1
        while i + 1 < len(spike_starts) and spike_starts[i+1] == spike_starts[i] + 1:
            i += 1
            e = spike_starts[i] + 1
        # interpolate across the spike window [s, e]
        left  = max(0, s - 1)
        right = min(len(y) - 1, e + 1)
        if right > left:
            for j in range(s, e + 1):
                if j < len(y):
                    t = (j - left) / max(right - left, 1)
                    y_clean[j] = y_clean[left] * (1 - t) + y_clean[right] * t
        removed += (e - s + 1)
        i += 1

    return y_clean, removed


def baseline_als(y, lam=1e5, p=0.01, niter=10):
    """Asymmetric Least Squares (Eilers & Boelens 2005)."""
    L = len(y)
    D = diags([1, -2, 1], [0, 1, 2], shape=(L - 2, L))
    D = lam * D.T @ D
    w = np.ones(L)
    for _ in range(niter):
        W = diags(w, 0)
        z = spsolve(W + D, w * y)
        w = p * (y > z) + (1 - p) * (y <= z)
    return z


def baseline_arpls(y, lam=1e4, ratio=0.05, niter=50):
    """
    Asymmetrically Reweighted Penalised Least Squares (Baek et al. 2015).
    More robust than ALS — weights update automatically via logistic function
    so the p parameter is less critical.
    """
    N = len(y)
    D = diags([1, -2, 1], [0, 1, 2], shape=(N - 2, N))
    H = lam * D.T @ D
    w = np.ones(N)
    for _ in range(niter):
        W = diags(w, 0)
        z = spsolve(W + H, w * y)
        d = y - z
        d_neg = d[d < 0]
        m = d_neg.mean() if len(d_neg) > 0 else 0.0
        s = d_neg.std()  if len(d_neg) > 1 else 1.0
        w_new = 1.0 / (1.0 + np.exp(2.0 * (d - (2.0 * s - m)) / (s + 1e-12)))
        if np.linalg.norm(w_new - w) / (np.linalg.norm(w) + 1e-12) < ratio:
            break
        w = w_new
    return z


def baseline_snip(y, max_iter=40):
    """
    Statistics-sensitive Non-linear Iterative Peak-clipping (SNIP).
    Excellent for broad fluorescence backgrounds with complex curvature.
    Operates in log-log-sqrt space to linearise the background.
    """
    n = len(y)
    # Transform to linearise background
    v = np.log(np.log(np.sqrt(np.maximum(y, 0) + 1) + 1) + 1)
    for i in range(1, max_iter + 1):
        v_new = v.copy()
        lo = np.arange(i, n - i)
        v_new[lo] = np.minimum(v[lo], (v[lo - i] + v[lo + i]) / 2.0)
        v = v_new
    # Inverse transform
    return (np.exp(np.exp(v) - 1) - 1) ** 2 - 1


def pseudo_voigt(x, centre, amplitude, fwhm, eta):
    """
    Pseudo-Voigt lineshape: linear mix of Gaussian and Lorentzian.
    eta=0 → pure Gaussian, eta=1 → pure Lorentzian.
    Raman peaks are typically Lorentzian-dominated (eta 0.5–1.0).
    """
    eta  = np.clip(eta, 0.0, 1.0)
    fwhm = np.abs(fwhm) + 1e-12
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    G = np.exp(-((x - centre) ** 2) / (2.0 * sigma ** 2))
    L = 1.0 / (1.0 + ((x - centre) / (fwhm / 2.0)) ** 2)
    return amplitude * (eta * L + (1.0 - eta) * G)


def fit_peak(x, y, centre_guess, window_cm=80.0):
    """
    Fit a pseudo-Voigt to y around centre_guess ± window_cm.
    Returns dict with fit parameters and their uncertainties.
    """
    mask  = np.abs(x - centre_guess) <= window_cm
    x_fit = x[mask]
    y_fit = y[mask]
    if len(x_fit) < 6:
        raise ValueError("Too few points in fit window.")

    amp0  = y_fit.max()
    fwhm0 = window_cm * 0.3
    p0    = [centre_guess, amp0, fwhm0, 0.5]
    lo    = [x_fit[0],  0,       1.0,  0.0]
    hi    = [x_fit[-1], amp0*5,  window_cm*2, 1.0]

    popt, pcov = curve_fit(pseudo_voigt, x_fit, y_fit, p0=p0,
                           bounds=(lo, hi), maxfev=5000)
    perr = np.sqrt(np.diag(pcov))
    centre, amplitude, fwhm, eta = popt
    x_dense = np.linspace(x_fit[0], x_fit[-1], 500)
    y_fit_line = pseudo_voigt(x_dense, *popt)

    # Area under the fit curve (trapezoidal)
    area = np.trapz(y_fit_line, x_dense)

    return {
        "centre":    centre,    "centre_err":    perr[0],
        "amplitude": amplitude, "amplitude_err": perr[1],
        "fwhm":      fwhm,      "fwhm_err":      perr[2],
        "eta":       eta,       "eta_err":       perr[3],
        "area":      area,
        "x_fit":     x_dense,   "y_fit":         y_fit_line,
        "x_data":    x_fit,     "y_data":        y_fit,
    }



# ═══════════════════════════════════════════════════════════════════════════════
# Library matching
# ═══════════════════════════════════════════════════════════════════════════════

def match_peaks_to_library(query_peaks, library, tolerance_cm=30.0,
                            require_fraction=0.4, top_n=10):
    """
    Match a set of detected peak positions against a spectral library.

    Strategy:
      For each library entry that has spectral data:
        1. Find its peaks using find_peaks on the normalised spectrum.
        2. For each query peak, count how many library peaks fall within
           ±tolerance_cm  (Hungarian-style greedy nearest-neighbour).
        3. Score = matched_peaks / max(n_query, n_lib_peaks)  → Jaccard-like.
        4. Penalise entries where fewer than require_fraction of query peaks
           are matched.
      Returns top_n entries sorted by score descending.

    Parameters
    ----------
    query_peaks : array-like of float
        Detected peak positions in the query spectrum (cm⁻¹).
    library : list of dict
        Each entry must have keys: name, substance_name, x_axis, spectral_data.
    tolerance_cm : float
        Max allowed shift between query and reference peak (cm⁻¹).
    require_fraction : float
        Minimum fraction of query peaks that must be matched for an entry to
        appear in results (0 = no filter, 1 = all query peaks matched).
    top_n : int
        Number of top candidates to return.

    Returns
    -------
    list of dict with keys:
        name, substance_name, score, matched_pairs, n_query, n_lib,
        mean_error, max_error, spectral_data, x_axis
    """
    if len(query_peaks) == 0:
        return []

    q = np.asarray(query_peaks, dtype=float)
    results = []

    for entry in library:
        y_raw = entry.get("spectral_data")
        x_raw = entry.get("x_axis") or entry.get("x")
        if not y_raw or len(y_raw) < 10:
            continue

        y = np.asarray(y_raw, dtype=float)
        # Handle NaN/inf
        if not np.all(np.isfinite(y)):
            y = np.where(np.isfinite(y), y, 0.0)

        # Build x axis
        if x_raw and len(x_raw) == len(y):
            x = np.asarray(x_raw, dtype=float)
        else:
            x = np.linspace(200, 3500, len(y))

        # Normalise library spectrum 0–1
        rng = y.max() - y.min()
        if rng < 1e-12:
            continue
        yn = (y - y.min()) / rng

        # Find peaks in library spectrum (adaptive threshold)
        try:
            dx    = float(np.mean(np.diff(x))) if len(x) > 1 else 1.0
            dist  = max(1, int(15.0 / dx))   # 15 cm⁻¹ min distance
            pidx, _ = find_peaks(yn, height=0.05, distance=dist, prominence=0.03)
        except Exception:
            continue

        if len(pidx) == 0:
            continue
        lib_peaks = x[pidx]

        # Greedy nearest-neighbour matching (each lib peak used at most once)
        matched_q   = []   # query peak positions that were matched
        matched_l   = []   # corresponding library peak positions
        lib_used    = set()

        for qp in q:
            diffs = np.abs(lib_peaks - qp)
            best  = int(np.argmin(diffs))
            if diffs[best] <= tolerance_cm and best not in lib_used:
                matched_q.append(qp)
                matched_l.append(lib_peaks[best])
                lib_used.add(best)

        n_matched = len(matched_q)
        n_query   = len(q)
        n_lib     = len(lib_peaks)

        # Require minimum coverage of query peaks
        if n_query > 0 and n_matched / n_query < require_fraction:
            continue

        # Jaccard-like score: matched / union
        score = n_matched / max(n_query, n_lib)

        # Bonus: how well-centred are the matches?
        errors = [abs(a - b) for a, b in zip(matched_q, matched_l)]
        mean_err = float(np.mean(errors)) if errors else tolerance_cm
        max_err  = float(np.max(errors))  if errors else tolerance_cm

        # Slight score boost for tight matches
        tightness_bonus = max(0.0, (tolerance_cm - mean_err) / tolerance_cm) * 0.15
        score += tightness_bonus

        results.append({
            "name":          entry.get("name", ""),
            "substance_name":entry.get("substance_name", entry.get("chemicals", "")),
            "id":            entry.get("id", ""),
            "score":         round(score, 4),
            "matched_pairs": list(zip(matched_q, matched_l)),
            "n_matched":     n_matched,
            "n_query":       n_query,
            "n_lib":         n_lib,
            "mean_error":    round(mean_err, 1),
            "max_error":     round(max_err, 1),
            "spectral_data": y_raw,
            "x_axis":        list(x),
            "lib_peaks":     list(lib_peaks),
        })

    results.sort(key=lambda r: r["score"], reverse=True)
    return results[:top_n]


def cosine_similarity(query_x, query_y, lib_x, lib_y, wn_min=300, wn_max=3400, n_points=1000):
    """
    Spectral cosine similarity after interpolating both spectra onto a common
    wavenumber grid.  Returns value in [0, 1].
    """
    grid = np.linspace(wn_min, wn_max, n_points)
    try:
        q_interp = np.interp(grid, query_x, query_y, left=0, right=0)
        l_interp = np.interp(grid, lib_x,   lib_y,   left=0, right=0)
    except Exception:
        return 0.0
    # Normalise
    qn = q_interp / (np.linalg.norm(q_interp) + 1e-12)
    ln = l_interp / (np.linalg.norm(l_interp) + 1e-12)
    return float(np.clip(np.dot(qn, ln), 0.0, 1.0))


# ═══════════════════════════════════════════════════════════════════════════════
# Helper utilities
# ═══════════════════════════════════════════════════════════════════════════════

def nm_to_shift(wl, laser=532.0):
    return (1.0 / laser - 1.0 / wl) * 1e7


def norm01(y):
    mn, mx = np.nanmin(y), np.nanmax(y)
    return (y - mn) / (mx - mn + 1e-12)


def load_csv(path):
    df = pd.read_csv(path, skiprows=4)
    wl  = next((c for c in df.columns if "Wavelength" in c), None)
    avg = next((c for c in df.columns if "Averaged"   in c), None)
    if wl is None or avg is None:
        df  = pd.read_csv(path)
        wl  = next((c for c in df.columns if "Wavelength" in c or "nm" in c.lower()), None)
        avg = next((c for c in df.columns if "Average"    in c), None)
    if wl is None or avg is None:
        raise ValueError("Cannot find Wavelength / Averaged columns.")
    return df[wl].values.astype(float), df[avg].values.astype(float)


# ═══════════════════════════════════════════════════════════════════════════════
# Library match results dialog
# ═══════════════════════════════════════════════════════════════════════════════

class MatchResultsDialog(tk.Toplevel):
    """Scrollable table of top library matches with spectral overlay preview."""

    COLS = [
        ("Rank",       40),
        ("Score",      60),
        ("Cosine",     60),
        ("Name",       260),
        ("Substance",  140),
        ("Peaks",      60),
        ("Δ̄ (cm⁻¹)",  75),
        ("Max Δ",      65),
    ]

    def __init__(self, parent, results, x_q, y_q):
        super().__init__(parent)
        self.title("Library Match Results")
        self.configure(bg=BG)
        self.geometry("1000x640")
        self.resizable(True, True)
        self.results  = results
        self.x_q      = x_q
        self.y_q      = y_q
        self.parent   = parent
        self._sel_idx = 0
        self._build()
        self._populate()
        self._select_row(0)

    def _build(self):
        ttk.Label(self, text=f"Library Match Results  —  {len(self.results)} candidate(s)",
                  style="Header.TLabel").pack(padx=14, pady=(12, 4))
        ttk.Separator(self, orient="horizontal").pack(fill="x", padx=10, pady=2)

        # Paned: table left, preview right
        pw = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        pw.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

        # ── Left: results table ──
        left = ttk.Frame(pw)
        pw.add(left, weight=1)

        col_ids = [c[0] for c in self.COLS]
        self.tree = ttk.Treeview(left, columns=col_ids, show="headings",
                                  selectmode="browse", height=20)
        for cname, cw in self.COLS:
            self.tree.heading(cname, text=cname)
            self.tree.column(cname, width=cw, stretch=(cname == "Name"))
        self.tree.tag_configure("top",    background="#2a2a3e", foreground=GREEN)
        self.tree.tag_configure("good",   background=BG,        foreground=TEXT)
        self.tree.tag_configure("weak",   background=BG,        foreground=OVERLAY)

        vsb = ttk.Scrollbar(left, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        vsb.pack(side=tk.LEFT, fill=tk.Y)
        self.tree.bind("<<TreeviewSelect>>", self._on_select)

        # ── Right: spectrum overlay preview ──
        right = ttk.Frame(pw)
        pw.add(right, weight=1)

        fig = Figure(figsize=(4.5, 4.5), facecolor=BG, tight_layout=True)
        self.prev_ax = fig.add_subplot(111)
        self.prev_ax.set_facecolor("#181825")
        for sp in self.prev_ax.spines.values(): sp.set_color(OVERLAY)
        self.prev_ax.tick_params(colors=TEXT, labelsize=8)
        self.prev_ax.set_xlabel("Raman shift (cm⁻¹)", color=TEXT, fontsize=8)
        self.prev_ax.set_ylabel("Intensity (norm.)",  color=TEXT, fontsize=8)
        self.prev_canvas = FigureCanvasTkAgg(fig, master=right)
        self.prev_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Metadata text
        self.meta_var = tk.StringVar(value="")
        ttk.Label(right, textvariable=self.meta_var, background=BG, foreground=TEXT,
                  font=("Courier New", 8), justify=tk.LEFT, wraplength=420
                  ).pack(anchor=tk.W, padx=8, pady=4)

        # Buttons
        bf = ttk.Frame(self)
        bf.pack(fill=tk.X, padx=10, pady=(4, 10))
        ttk.Button(bf, text="Overlay top 3 on main plot",
                   style="Accent.TButton",
                   command=self._overlay_main).pack(side=tk.LEFT, padx=(0, 8))
        ttk.Button(bf, text="Clear overlays",
                   command=self._clear_overlays).pack(side=tk.LEFT, padx=(0, 8))
        ttk.Button(bf, text="Export results CSV",
                   command=self._export_csv).pack(side=tk.LEFT, padx=(0, 8))
        ttk.Button(bf, text="Close", command=self.destroy).pack(side=tk.RIGHT)

    def _populate(self):
        self.tree.delete(*self.tree.get_children())
        for i, r in enumerate(self.results):
            score    = f"{r['score']:.3f}"
            cos_str  = f"{r['cosine']:.3f}" if r.get("cosine") is not None else "—"
            n_str    = f"{r['n_matched']}/{r['n_query']}"
            tag      = "top" if i == 0 else ("good" if r["score"] > 0.3 else "weak")
            self.tree.insert("", tk.END, iid=str(i),
                values=(i+1, score, cos_str,
                        r["name"][:50],
                        r["substance_name"][:25],
                        n_str,
                        f"{r['mean_error']:.1f}",
                        f"{r['max_error']:.1f}"),
                tags=(tag,))
        if self.tree.get_children():
            self.tree.selection_set(self.tree.get_children()[0])

    def _on_select(self, event):
        sel = self.tree.selection()
        if sel:
            self._select_row(int(sel[0]))

    def _select_row(self, idx):
        if idx >= len(self.results):
            return
        self._sel_idx = idx
        r = self.results[idx]

        # Update preview plot
        self.prev_ax.cla()
        self.prev_ax.set_facecolor("#181825")
        self.prev_ax.set_xlabel("Raman shift (cm⁻¹)", color=TEXT, fontsize=8)
        self.prev_ax.set_ylabel("Intensity (norm.)",  color=TEXT, fontsize=8)
        self.prev_ax.tick_params(colors=TEXT, labelsize=8)

        # Query spectrum
        self.prev_ax.plot(self.x_q, self.y_q, color=BLUE,
                          linewidth=1.4, label="Query", alpha=0.9)

        # Library spectrum
        lx = r.get("x_axis")
        ly = r.get("spectral_data")
        if lx and ly:
            la = np.asarray(ly, float)
            rng = la.max() - la.min()
            if rng > 1e-12:
                la = (la - la.min()) / rng
            self.prev_ax.plot(np.asarray(lx), la, color=MAUVE,
                              linewidth=1.2, linestyle="--",
                              label=r["name"][:35], alpha=0.85)

        # Mark matched peak pairs
        for qp, lp in r.get("matched_pairs", []):
            self.prev_ax.axvline(qp, color=YELLOW, linewidth=0.7, alpha=0.5)
            self.prev_ax.axvline(lp, color=MAUVE,  linewidth=0.7, alpha=0.4, linestyle=":")
            mid = (qp + lp) / 2
            self.prev_ax.annotate("",
                xy=(lp, 0.92), xytext=(qp, 0.92),
                xycoords=("data", "axes fraction"),
                textcoords=("data", "axes fraction"),
                arrowprops=dict(arrowstyle="<->", color=PEACH, lw=0.8))

        self.prev_ax.legend(facecolor=SURFACE, edgecolor=OVERLAY,
                             labelcolor=TEXT, fontsize=7, loc="upper right")
        self.prev_canvas.draw()

        # Update metadata
        cos_txt = f"{r['cosine']:.4f}" if r.get("cosine") is not None else "n/a"
        self.meta_var.set(
            f"name:      {r['name']}\n"
            f"substance: {r['substance_name']}\n"
            f"ID:        {r.get('id','')}\n"
            f"score:     {r['score']:.4f}  cosine: {cos_txt}\n"
            f"peaks:     {r['n_matched']}/{r['n_query']} matched  "
            f"(lib has {r['n_lib']})\n"
            f"mean Δ:    {r['mean_error']:.1f} cm⁻¹  "
            f"max Δ: {r['max_error']:.1f} cm⁻¹"
        )

    def _overlay_main(self):
        """Send top 3 overlays to the main plot window."""
        colours = [MAUVE, PEACH, RED]
        self.parent._match_overlays = []
        for i, r in enumerate(self.results[:3]):
            lx = r.get("x_axis")
            ly = r.get("spectral_data")
            if lx and ly:
                la = np.asarray(ly, float)
                rng = la.max() - la.min()
                if rng > 1e-12:
                    la = (la - la.min()) / rng
                self.parent._match_overlays.append((
                    np.asarray(lx), la,
                    f"#{i+1} {r['name'][:28]}",
                    colours[i % len(colours)]))
        # Store matched pairs of top hit for span highlighting
        self.parent._last_match_pairs = self.results[0].get("matched_pairs", [])             if self.results else []
        self.parent._update_plot()

    def _clear_overlays(self):
        self.parent._match_overlays = []
        self.parent._last_match_pairs = []
        self.parent._update_plot()

    def _export_csv(self):
        import csv as _csv
        path = filedialog.asksaveasfilename(
            parent=self, defaultextension=".csv",
            filetypes=[("CSV", "*.csv")],
            initialfile="raman_matches.csv")
        if not path:
            return
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = _csv.writer(f)
            w.writerow(["Rank","Score","Cosine","Combined","Name","Substance","ID",
                        "Peaks_matched","Peaks_query","Peaks_lib",
                        "Mean_error_cm","Max_error_cm","Matched_pairs"])
            for i, r in enumerate(self.results):
                pairs_str = "; ".join(f"{a:.1f}↔{b:.1f}" for a, b in r.get("matched_pairs",[]))
                w.writerow([
                    i+1, r["score"],
                    r.get("cosine",""), r.get("combined",""),
                    r["name"], r["substance_name"], r.get("id",""),
                    r["n_matched"], r["n_query"], r["n_lib"],
                    r["mean_error"], r["max_error"], pairs_str])
        messagebox.showinfo("Exported", f"Saved to:\n{path}", parent=self)


# ═══════════════════════════════════════════════════════════════════════════════
# Peak-fit results dialog
# ═══════════════════════════════════════════════════════════════════════════════

class FitResultDialog(tk.Toplevel):
    def __init__(self, parent, result):
        super().__init__(parent)
        self.title("Peak Fit Results")
        self.configure(bg=BG)
        self.resizable(False, False)
        self.grab_set()

        c, ce   = result["centre"],    result["centre_err"]
        fw, fwe = result["fwhm"],      result["fwhm_err"]
        et, ete = result["eta"],       result["eta_err"]
        am, ame = result["amplitude"], result["amplitude_err"]
        area    = result["area"]
        char    = "Lorentzian" if et > 0.66 else ("Gaussian" if et < 0.33 else "Mixed")

        rows = [
            ("Centre",    f"{c:.2f} ± {ce:.2f} cm⁻¹"),
            ("FWHM",      f"{fw:.2f} ± {fwe:.2f} cm⁻¹"),
            ("Amplitude", f"{am:.4f} ± {ame:.4f}"),
            ("η (eta)",   f"{et:.3f} ± {ete:.3f}  →  {char}"),
            ("Area",      f"{area:.4f}  (a.u. · cm⁻¹)"),
        ]

        ttk.Label(self, text="Pseudo-Voigt Fit", style="Header.TLabel").pack(padx=16, pady=(14, 6))
        ttk.Separator(self, orient="horizontal").pack(fill="x", padx=12, pady=2)

        for lbl, val in rows:
            f = ttk.Frame(self); f.pack(fill="x", padx=16, pady=3)
            ttk.Label(f, text=lbl + ":", width=12, anchor="w").pack(side="left")
            ttk.Label(f, text=val, foreground=BLUE).pack(side="left")

        ttk.Separator(self, orient="horizontal").pack(fill="x", padx=12, pady=(6, 2))
        ttk.Label(self, text="η < 0.33 → Gaussian  |  η > 0.66 → Lorentzian",
                  foreground=OVERLAY, font=("Segoe UI", 8)).pack(padx=16, pady=(2, 4))
        ttk.Button(self, text="Close", command=self.destroy).pack(padx=16, pady=(4, 14))


# ═══════════════════════════════════════════════════════════════════════════════
# Export dialog
# ═══════════════════════════════════════════════════════════════════════════════

class ExportDialog(tk.Toplevel):
    def __init__(self, parent, default_title=""):
        super().__init__(parent)
        self.title("Export Spectrum")
        self.configure(bg=BG)
        self.resizable(False, False)
        self.grab_set()
        self.result = None

        ttk.Label(self, text="Title shown on exported image:", style="TLabel").pack(padx=16, pady=(14, 2), anchor="w")
        self.title_var = tk.StringVar(value=default_title)
        e = ttk.Entry(self, textvariable=self.title_var, width=44)
        e.pack(padx=16, pady=(0, 10), fill="x")
        e.focus_set()

        ttk.Label(self, text="Format:", style="TLabel").pack(padx=16, anchor="w")
        self.fmt_var = tk.StringVar(value="PNG")
        ff = ttk.Frame(self); ff.pack(padx=16, pady=(0, 10), anchor="w")
        ttk.Radiobutton(ff, text="PNG (700×450 px)", variable=self.fmt_var, value="PNG").pack(side="left", padx=(0, 16))
        ttk.Radiobutton(ff, text="SVG (vector)",     variable=self.fmt_var, value="SVG").pack(side="left")

        bf = ttk.Frame(self); bf.pack(padx=16, pady=(4, 14))
        ttk.Button(bf, text="Choose file & Save", style="Accent.TButton", command=self._save).pack(side="left", padx=(0, 8))
        ttk.Button(bf, text="Cancel", command=self.destroy).pack(side="left")

    def _save(self):
        fmt = self.fmt_var.get().lower()
        path = filedialog.asksaveasfilename(
            parent=self, defaultextension=f".{fmt}",
            filetypes=[(fmt.upper(), f"*.{fmt}")],
            initialfile=f"raman_spectrum.{fmt}")
        if path:
            self.result = (self.title_var.get(), fmt, path)
            self.destroy()


# ═══════════════════════════════════════════════════════════════════════════════
# Shared plot draw (canvas + export)
# ═══════════════════════════════════════════════════════════════════════════════

def _draw_spectrum(ax, x, raw, bl, proc, peaks_x, peaks_y,
                   show_raw=True, show_bl=True, small=False,
                   fit_result=None, fit_results=None):
    lw_bg   = 0.7 if small else 0.9
    lw_main = 1.4 if small else 1.6
    fs_ann  = 6.0 if small else 7.5
    ms      = 28  if small else 42

    if show_raw:
        raw_n = norm01(raw)
        ax.fill_between(x, raw_n, alpha=0.07, color=OVERLAY, zorder=1)
        ax.plot(x, raw_n, color=OVERLAY, linewidth=lw_bg, alpha=0.40,
                label="Raw (norm.)", zorder=1)

    if show_bl:
        mn, mx = np.nanmin(raw), np.nanmax(raw)
        bl_n = (bl - mn) / (mx - mn + 1e-12)
        ax.plot(x, bl_n, color=RED, linewidth=lw_bg, linestyle="--", alpha=0.50,
                label="Baseline (norm.)", zorder=2)

    ax.plot(x, proc, color=BLUE, linewidth=lw_main, label="Processed", zorder=3)

    # All-peaks fit overlays (dim green, no fill)
    if fit_results:
        for i, r in enumerate(fit_results):
            lbl = "Peak fits" if i == 0 else "_nolegend_"
            ax.plot(r["x_fit"], r["y_fit"], color=GREEN,
                    linewidth=0.9 if small else 1.1,
                    linestyle="-", alpha=0.65, label=lbl, zorder=5)
            ax.fill_between(r["x_fit"], r["y_fit"],
                            alpha=0.06, color=GREEN, zorder=4)

    # Single highlighted fit overlay (bright green + fill)
    if fit_result is not None:
        ax.plot(fit_result["x_fit"], fit_result["y_fit"],
                color=GREEN, linewidth=1.2 if small else 1.8,
                linestyle="-", alpha=0.95, label="Selected fit", zorder=7)
        ax.fill_between(fit_result["x_fit"], fit_result["y_fit"],
                        alpha=0.15, color=GREEN, zorder=6)

    if len(peaks_x) > 0:
        ax.scatter(peaks_x, peaks_y, color=YELLOW, s=ms, zorder=8, linewidths=0)
        for px in peaks_x:
            ax.axvline(px, color=YELLOW, linewidth=0.5, alpha=0.20, zorder=4)
        y_range = np.nanmax(proc) - np.nanmin(proc)
        offset  = y_range * 0.04
        for px, py in zip(peaks_x, peaks_y):
            ax.annotate(f"{px:.0f}",
                xy=(px, py), xytext=(px, py + offset),
                ha="center", va="bottom", fontsize=fs_ann, color=YELLOW,
                arrowprops=dict(arrowstyle="-", color=YELLOW, lw=0.6, alpha=0.6))


# ═══════════════════════════════════════════════════════════════════════════════
# All-peaks fit results dialog
# ═══════════════════════════════════════════════════════════════════════════════

class AllFitsDialog(tk.Toplevel):
    """Scrollable table of pseudo-Voigt results for all fitted peaks."""
    def __init__(self, parent, results):
        super().__init__(parent)
        self.title("All Peak Fit Results")
        self.configure(bg=BG)
        self.resizable(True, True)
        self.grab_set()

        ttk.Label(self, text=f"Pseudo-Voigt Fits  —  {len(results)} peaks",
                  style="Header.TLabel").pack(padx=16, pady=(14, 6))
        ttk.Separator(self, orient="horizontal").pack(fill="x", padx=12, pady=2)

        # Column headers
        cols = ["Centre (cm⁻¹)", "FWHM (cm⁻¹)", "Amplitude", "η (eta)", "Shape", "Area"]
        widths = [120, 110, 100, 80, 90, 100]

        hdr_frame = ttk.Frame(self); hdr_frame.pack(fill="x", padx=14, pady=(4, 0))
        for col, w in zip(cols, widths):
            ttk.Label(hdr_frame, text=col, foreground=BLUE, background=BG,
                      font=("Segoe UI", 8, "bold"), width=w//7,
                      anchor="w").pack(side="left", padx=2)

        ttk.Separator(self, orient="horizontal").pack(fill="x", padx=12, pady=2)

        # Scrollable rows
        scroll_frame_outer = ttk.Frame(self)
        scroll_frame_outer.pack(fill="both", expand=True, padx=14, pady=2)
        canvas = tk.Canvas(scroll_frame_outer, bg=BG, highlightthickness=0,
                           height=min(len(results) * 26 + 10, 400))
        vsb = ttk.Scrollbar(scroll_frame_outer, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        inner = ttk.Frame(canvas)
        canvas.create_window((0, 0), window=inner, anchor="nw")
        inner.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        for i, r in enumerate(results):
            eta    = r["eta"]
            shape  = "Lorentzian" if eta > 0.66 else ("Gaussian" if eta < 0.33 else "Mixed")
            colour = TEXT if i % 2 == 0 else "#b4befe"
            vals   = [
                f"{r['centre']:.2f} ± {r['centre_err']:.2f}",
                f"{r['fwhm']:.2f} ± {r['fwhm_err']:.2f}",
                f"{r['amplitude']:.4f}",
                f"{eta:.3f} ± {r['eta_err']:.3f}",
                shape,
                f"{r['area']:.4f}",
            ]
            row_f = ttk.Frame(inner); row_f.pack(fill="x", pady=1)
            for val, w in zip(vals, widths):
                ttk.Label(row_f, text=val, foreground=colour, background=BG,
                          font=("Segoe UI", 8), width=w//7,
                          anchor="w").pack(side="left", padx=2)

        ttk.Separator(self, orient="horizontal").pack(fill="x", padx=12, pady=(4, 2))

        btn_row = ttk.Frame(self); btn_row.pack(padx=14, pady=(4, 14))
        ttk.Button(btn_row, text="📋  Copy CSV", command=lambda: self._copy_csv(results, parent)).pack(side="left", padx=(0, 8))
        ttk.Button(btn_row, text="Close", command=self.destroy).pack(side="left")

    def _copy_csv(self, results, parent):
        lines = ["Centre_cm1,Centre_err,FWHM_cm1,FWHM_err,Amplitude,Amp_err,Eta,Eta_err,Area"]
        for r in results:
            lines.append(
                f"{r['centre']:.4f},{r['centre_err']:.4f},"
                f"{r['fwhm']:.4f},{r['fwhm_err']:.4f},"
                f"{r['amplitude']:.6f},{r['amplitude_err']:.6f},"
                f"{r['eta']:.4f},{r['eta_err']:.4f},"
                f"{r['area']:.6f}")
        parent.clipboard_clear()
        parent.clipboard_append("\n".join(lines))
        messagebox.showinfo("Copied", "Fit results copied to clipboard as CSV.", parent=self)


# ═══════════════════════════════════════════════════════════════════════════════
# Main application
# ═══════════════════════════════════════════════════════════════════════════════

class RamanApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Raman Spectra Processor  v3")
        self.geometry("1380x880")
        self.configure(bg=BG)

        # Data state
        self.raw_wl        = None
        self.raw_shift     = None
        self.raw_intensity = None
        self.clean_intensity = None   # after spike removal
        self.processed     = None
        self.bl_arr        = None
        self.filepath      = None
        self.peaks_x       = np.array([])
        self.peaks_y       = np.array([])
        self.fit_result    = None     # last pseudo-Voigt fit
        self.fit_results   = []       # all-peaks fit results list
        self._folder_path  = None
        self._folder_files = []
        self._click_cid    = None     # matplotlib event connection id
        self._library      = []       # loaded RamanBase library entries
        self._match_overlays = []     # (x, y, label, colour) drawn on plot
        self._last_match_pairs = []   # matched (query, lib) peak pairs for highlighting

        self._style()
        self._build_ui()

    # ── Style ──────────────────────────────────────────────────────────────────

    def _style(self):
        s = ttk.Style(self)
        s.theme_use("clam")
        s.configure("TFrame",         background=BG)
        s.configure("TLabel",         background=BG,      foreground=TEXT,  font=("Segoe UI", 9))
        s.configure("Small.TLabel",   background=BG,      foreground=OVERLAY, font=("Segoe UI", 8))
        s.configure("Header.TLabel",  background=BG,      foreground=BLUE,  font=("Segoe UI", 10, "bold"))
        s.configure("Result.TLabel",  background=SURFACE, foreground=GREEN, font=("Segoe UI", 9))
        s.configure("TScale",         background=BG,      troughcolor=SURFACE, sliderlength=15)
        s.configure("TButton",        background=SURFACE, foreground=TEXT,  font=("Segoe UI", 9), relief="flat", padding=5)
        s.map("TButton",              background=[("active", OVERLAY)])
        s.configure("Accent.TButton", background=BLUE,    foreground=BG,    font=("Segoe UI", 9, "bold"), padding=5)
        s.map("Accent.TButton",       background=[("active", "#74c7ec")])
        s.configure("Green.TButton",  background=GREEN,   foreground=BG,    font=("Segoe UI", 9, "bold"), padding=5)
        s.map("Green.TButton",        background=[("active", "#94d3a2")])
        s.configure("TCombobox",      fieldbackground=SURFACE, background=SURFACE,
                    foreground=TEXT, selectbackground=OVERLAY, selectforeground=TEXT,
                    insertcolor=TEXT, arrowcolor=TEXT)
        s.map("TCombobox",
              fieldbackground=[("readonly", SURFACE)],
              foreground=[("readonly", TEXT)],
              selectbackground=[("readonly", OVERLAY)],
              selectforeground=[("readonly", TEXT)])
        s.configure("TCheckbutton",   background=BG, foreground=TEXT)
        s.configure("TRadiobutton",   background=BG, foreground=TEXT)
        s.configure("TSeparator",     background=OVERLAY)
        s.configure("TSpinbox",       fieldbackground=SURFACE, background=SURFACE,
                    foreground=TEXT, insertcolor=TEXT, arrowcolor=TEXT)

        self.option_add("*TCombobox*Listbox.background",       SURFACE)
        self.option_add("*TCombobox*Listbox.foreground",       TEXT)
        self.option_add("*TCombobox*Listbox.selectBackground", BLUE)
        self.option_add("*TCombobox*Listbox.selectForeground", BG)
        self.option_add("*TCombobox*Listbox.relief",           "flat")
        self.option_add("*TCombobox*Listbox.font",             ("Segoe UI", 9))

    # ── UI scaffold ────────────────────────────────────────────────────────────

    def _build_ui(self):
        # Toolbar
        tb = ttk.Frame(self); tb.pack(fill="x", padx=10, pady=(10, 0))
        ttk.Button(tb, text="📂  Open File",    style="Accent.TButton", command=self._open_file).pack(side="left", padx=(0, 5))
        ttk.Button(tb, text="📁  Open Folder",  style="TButton",        command=self._open_folder).pack(side="left", padx=(0, 10))
        self.file_label = ttk.Label(tb, text="No file loaded"); self.file_label.pack(side="left", padx=4)
        for txt, cmd in [("↺  Reset", self._reset), ("📦  Batch Export", self._batch_export),
                         ("💾  Export CSV", self._export_csv), ("🖼  Save Plot", self._save_plot)]:
            ttk.Button(tb, text=txt, style="TButton", command=cmd).pack(side="right", padx=3)

        # Main area
        main = ttk.Frame(self); main.pack(fill="both", expand=True, padx=10, pady=8)

        # ── Scrollable left panel ──
        ctrl_outer = ttk.Frame(main, width=295); ctrl_outer.pack(side="left", fill="y", padx=(0, 8))
        ctrl_outer.pack_propagate(False)
        self._build_scrollable_panel(ctrl_outer)

        # ── Plot area ──
        pf = ttk.Frame(main); pf.pack(side="left", fill="both", expand=True)
        self._build_plot(pf)

    def _build_scrollable_panel(self, outer):
        canvas  = tk.Canvas(outer, bg=BG, highlightthickness=0, bd=0)
        vbar    = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vbar.set)
        vbar.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        self._ctrl_frame = ttk.Frame(canvas)
        win_id = canvas.create_window((0, 0), window=self._ctrl_frame, anchor="nw")

        def _on_frame_configure(e):
            canvas.configure(scrollregion=canvas.bbox("all"))
        def _on_canvas_configure(e):
            canvas.itemconfig(win_id, width=e.width)
        def _on_mousewheel(e):
            canvas.yview_scroll(int(-1 * (e.delta / 120)), "units")

        self._ctrl_frame.bind("<Configure>",  _on_frame_configure)
        canvas.bind("<Configure>",            _on_canvas_configure)
        canvas.bind("<MouseWheel>",           _on_mousewheel)
        self._ctrl_frame.bind("<MouseWheel>", _on_mousewheel)

        self._build_controls(self._ctrl_frame)

    def _build_controls(self, p):
        pad = {"padx": 8}

        def sep():
            ttk.Separator(p, orient="horizontal").pack(fill="x", pady=3, **pad)

        def hdr(t, colour=BLUE):
            ttk.Label(p, text=t, background=BG, foreground=colour,
                      font=("Segoe UI", 10, "bold")).pack(anchor="w", pady=(10, 1), **pad)

        def row():
            f = ttk.Frame(p); f.pack(fill="x", pady=1, **pad); return f

        def slider(from_, to, var, cb, fmt=None):
            ttk.Scale(p, from_=from_, to=to, variable=var,
                      orient="horizontal", command=cb).pack(fill="x", **pad)

        # ── INSTRUMENT ──────────────────────────────────────────────────────────
        hdr("INSTRUMENT"); sep()
        r = row(); ttk.Label(r, text="Laser (nm):").pack(side="left")
        self.laser_var = tk.DoubleVar(value=532.0)
        sb = ttk.Spinbox(r, from_=400, to=1064, increment=1, textvariable=self.laser_var,
                         width=7, command=self._recompute)
        sb.pack(side="right"); sb.bind("<Return>", lambda e: self._recompute())

        ttk.Label(p, text="Wavenumber Range (cm⁻¹):").pack(anchor="w", pady=(6, 0), **pad)
        r2 = row(); ttk.Label(r2, text="Min:").pack(side="left")
        self.xmin_var = tk.IntVar(value=650)
        self.xmin_lbl = ttk.Label(r2, text="650", width=5); self.xmin_lbl.pack(side="right")
        slider(0, 3400, self.xmin_var, self._xmin_cb)

        r3 = row(); ttk.Label(r3, text="Max:").pack(side="left")
        self.xmax_var = tk.IntVar(value=3300)
        self.xmax_lbl = ttk.Label(r3, text="3300", width=5); self.xmax_lbl.pack(side="right")
        slider(0, 3400, self.xmax_var, self._xmax_cb)

        # ── SPIKE REMOVAL ───────────────────────────────────────────────────────
        hdr("SPIKE / COSMIC RAY REMOVAL", PEACH); sep()
        self.spike_on = tk.BooleanVar(value=True)
        ttk.Checkbutton(p, text="Remove cosmic-ray spikes", variable=self.spike_on,
                        command=self._recompute).pack(anchor="w", **pad)

        r_sp = row(); ttk.Label(r_sp, text="Threshold (σ):").pack(side="left")
        self.spike_thr_var = tk.DoubleVar(value=8.5)
        self.spike_thr_lbl = ttk.Label(r_sp, text="8.5", width=4); self.spike_thr_lbl.pack(side="right")
        slider(2, 15, self.spike_thr_var, self._spike_cb)

        self.spike_info_lbl = ttk.Label(p, text="", foreground=PEACH,
                                        background=BG, font=("Segoe UI", 8))
        self.spike_info_lbl.pack(anchor="w", **pad)

        # ── BASELINE CORRECTION ─────────────────────────────────────────────────
        hdr("BASELINE CORRECTION"); sep()
        r_bl = row(); ttk.Label(r_bl, text="Method:").pack(side="left")
        self.bl_method_var = tk.StringVar(value="SNIP")
        bl_combo = ttk.Combobox(r_bl, textvariable=self.bl_method_var, width=10,
                                state="readonly", values=["ALS", "arPLS", "SNIP"])
        bl_combo.pack(side="right")
        bl_combo.bind("<<ComboboxSelected>>", self._bl_method_changed)

        # Fixed-position container — child frames swap in here, never re-packed at bottom
        self._bl_params_container = ttk.Frame(p)
        self._bl_params_container.pack(fill="x")

        # ALS / arPLS params (hidden initially)
        self.als_frame = ttk.Frame(self._bl_params_container)
        ttk.Label(self.als_frame, text="λ  smoothness  (10^n):", background=BG,
                  foreground=TEXT).pack(anchor="w", padx=8)
        self.lam_var = tk.DoubleVar(value=3.2)
        self.lam_lbl = ttk.Label(self.als_frame, text="10^3.2  ≈ 1.6k",
                                 background=BG, foreground=TEXT)
        self.lam_lbl.pack(anchor="e", padx=8)
        ttk.Scale(self.als_frame, from_=1, to=9, variable=self.lam_var,
                  orient="horizontal", command=self._lam_cb).pack(fill="x", padx=8)
        ttk.Label(self.als_frame, text="p  asymmetry:", background=BG,
                  foreground=TEXT).pack(anchor="w", padx=8, pady=(4, 0))
        self.p_var = tk.DoubleVar(value=0.01)
        self.p_lbl = ttk.Label(self.als_frame, text="0.010", background=BG, foreground=TEXT)
        self.p_lbl.pack(anchor="e", padx=8)
        ttk.Scale(self.als_frame, from_=0.001, to=0.1, variable=self.p_var,
                  orient="horizontal", command=self._p_cb).pack(fill="x", padx=8)
        self.arpls_note = ttk.Label(self.als_frame,
            text="  arPLS: p sets convergence ratio.\n  Lower = stricter (0.01–0.05 typical).",
            background=BG, foreground=OVERLAY, font=("Segoe UI", 8))

        # SNIP params (shown by default — SNIP is the default method)
        self.snip_frame = ttk.Frame(self._bl_params_container)
        r_sn = ttk.Frame(self.snip_frame); r_sn.pack(fill="x", padx=8)
        ttk.Label(r_sn, text="Iterations:", background=BG, foreground=TEXT).pack(side="left")
        self.snip_iter_var = tk.IntVar(value=11)
        self.snip_iter_lbl = ttk.Label(r_sn, text="11", width=4, background=BG, foreground=TEXT)
        self.snip_iter_lbl.pack(side="right")
        ttk.Scale(self.snip_frame, from_=5, to=100, variable=self.snip_iter_var,
                  orient="horizontal", command=self._snip_iter_cb).pack(fill="x", padx=8)
        ttk.Label(self.snip_frame,
            text="  SNIP: more iterations → wider background window.\n  5–20 narrow, 40–60 broad fluorescence.",
            background=BG, foreground=OVERLAY, font=("Segoe UI", 8)).pack(anchor="w", padx=8)
        self.snip_frame.pack(fill="x")   # visible by default

        # ── SMOOTHING ───────────────────────────────────────────────────────────
        hdr("SMOOTHING"); sep()
        r4 = row(); ttk.Label(r4, text="Method:").pack(side="left")
        self.smooth_var = tk.StringVar(value="Savitzky-Golay")
        cm = ttk.Combobox(r4, textvariable=self.smooth_var, width=14, state="readonly",
                          values=["None", "Savitzky-Golay", "Moving Average"])
        cm.pack(side="right"); cm.bind("<<ComboboxSelected>>", lambda e: self._recompute())

        r5 = row(); ttk.Label(r5, text="Window:").pack(side="left")
        self.win_var = tk.IntVar(value=7)
        self.win_lbl = ttk.Label(r5, text="7", width=3); self.win_lbl.pack(side="right")
        slider(3, 51, self.win_var, self._win_cb)

        r6 = row(); ttk.Label(r6, text="SG poly order:").pack(side="left")
        self.ord_var = tk.IntVar(value=3)
        self.ord_lbl = ttk.Label(r6, text="3", width=3); self.ord_lbl.pack(side="right")
        slider(1, 7, self.ord_var, self._ord_cb)

        # ── NORMALISATION ────────────────────────────────────────────────────────
        hdr("NORMALISATION"); sep()
        r7 = row(); ttk.Label(r7, text="Method:").pack(side="left")
        self.norm_var = tk.StringVar(value="Max = 1")
        nc = ttk.Combobox(r7, textvariable=self.norm_var, width=10, state="readonly",
                          values=["None", "Max = 1", "Area = 1", "SNV"])
        nc.pack(side="right")
        self.norm_var.trace_add("write", lambda *a: self._recompute())

        # ── PEAK DETECTION ───────────────────────────────────────────────────────
        hdr("PEAK DETECTION"); sep()
        self.peaks_on = tk.BooleanVar(value=True)
        ttk.Checkbutton(p, text="Detect & label peaks", variable=self.peaks_on,
                        command=self._recompute).pack(anchor="w", **pad)

        r8 = row(); ttk.Label(r8, text="Min height (%):").pack(side="left")
        self.ph_var = tk.IntVar(value=12)
        self.ph_lbl = ttk.Label(r8, text="12%", width=5); self.ph_lbl.pack(side="right")
        slider(1, 80, self.ph_var, self._ph_cb)

        r9 = row(); ttk.Label(r9, text="Min distance (cm⁻¹):").pack(side="left")
        self.pd_var = tk.IntVar(value=30)
        self.pd_lbl = ttk.Label(r9, text="30", width=5); self.pd_lbl.pack(side="right")
        slider(5, 300, self.pd_var, self._pd_cb)

        # ── PEAK FITTING ─────────────────────────────────────────────────────────
        hdr("PEAK FITTING  (pseudo-Voigt)", GREEN); sep()
        ttk.Label(p, text="Click a peak to fit it, or auto-fit all peaks at once.",
                  background=BG, foreground=OVERLAY, font=("Segoe UI", 8)).pack(anchor="w", **pad)

        r_fw = row(); ttk.Label(r_fw, text="Fit window (±cm⁻¹):").pack(side="left")
        self.fit_win_var = tk.IntVar(value=80)
        self.fit_win_lbl = ttk.Label(r_fw, text="80", width=4); self.fit_win_lbl.pack(side="right")
        slider(10, 200, self.fit_win_var, self._fit_win_cb)

        fit_btn_row = row()
        self.fit_btn = ttk.Button(fit_btn_row, text="🔬 Fit one peak",
                                  style="Green.TButton", command=self._start_fit_click)
        self.fit_btn.pack(side="left", fill="x", expand=True, padx=(0, 3))
        ttk.Button(fit_btn_row, text="✕", width=2,
                   command=self._clear_fit).pack(side="right")

        ttk.Button(p, text="📊  Fit all detected peaks", style="Green.TButton",
                   command=self._fit_all_peaks).pack(fill="x", padx=8, pady=(4, 2))

        self.fit_status_lbl = ttk.Label(p, text="No fit yet.", foreground=OVERLAY,
                                        background=BG, font=("Segoe UI", 8), wraplength=260)
        self.fit_status_lbl.pack(anchor="w", **pad)

        # ── LIBRARY MATCHING ─────────────────────────────────────────────────────
        hdr("LIBRARY MATCHING  (peak & spectral)", MAUVE); sep()
        ttk.Label(p, text="Load a RamanBase library JSON, then run matching\nagainst detected peaks.",
                  background=BG, foreground=OVERLAY, font=("Segoe UI", 8)).pack(anchor="w", **pad)

        lib_row = row()
        ttk.Button(lib_row, text="📂 Load library…",
                   command=self._load_library).pack(side="left", fill="x", expand=True, padx=(0,3))
        self.lib_info_lbl = ttk.Label(lib_row, text="No library", foreground=OVERLAY,
                                      background=BG, font=("Segoe UI", 8), width=12)
        self.lib_info_lbl.pack(side="right")

        r_tol = row(); ttk.Label(r_tol, text="Tolerance (±cm⁻¹):").pack(side="left")
        self.tol_var = tk.IntVar(value=30)
        self.tol_lbl = ttk.Label(r_tol, text="30", width=4); self.tol_lbl.pack(side="right")
        ttk.Scale(p, from_=5, to=100, variable=self.tol_var,
                  orient="horizontal",
                  command=lambda v: self.tol_lbl.config(text=str(int(float(v))))).pack(fill="x", **pad)

        r_frac = row(); ttk.Label(r_frac, text="Min peak coverage (%):").pack(side="left")
        self.frac_var = tk.IntVar(value=30)
        self.frac_lbl = ttk.Label(r_frac, text="30%", width=5); self.frac_lbl.pack(side="right")
        ttk.Scale(p, from_=0, to=100, variable=self.frac_var,
                  orient="horizontal",
                  command=lambda v: self.frac_lbl.config(text=f"{int(float(v))}%")).pack(fill="x", **pad)

        r_top = row(); ttk.Label(r_top, text="Results to show:").pack(side="left")
        self.top_n_var = tk.IntVar(value=10)
        ttk.Spinbox(r_top, from_=1, to=50, textvariable=self.top_n_var, width=5).pack(side="right")

        r_cos = row()
        self.use_cosine_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(r_cos, text="Also score by spectral cosine similarity",
                        variable=self.use_cosine_var).pack(side="left")

        ttk.Button(p, text="🔍  Match library", style="Accent.TButton",
                   command=self._run_library_match).pack(fill="x", padx=8, pady=(6, 2))

        self.match_status_lbl = ttk.Label(p, text="", foreground=OVERLAY,
                                          background=BG, font=("Segoe UI", 8), wraplength=260)
        self.match_status_lbl.pack(anchor="w", **pad)

        # ── DISPLAY ──────────────────────────────────────────────────────────────
        hdr("DISPLAY"); sep()
        self.show_raw = tk.BooleanVar(value=True)
        self.show_bl  = tk.BooleanVar(value=True)
        ttk.Checkbutton(p, text="Show raw (norm., background)",
                        variable=self.show_raw, command=self._update_plot).pack(anchor="w", **pad)
        ttk.Checkbutton(p, text="Show baseline (norm., background)",
                        variable=self.show_bl, command=self._update_plot).pack(anchor="w", **pad)

        ttk.Button(p, text="⚡  Apply & Replot", style="Accent.TButton",
                   command=self._recompute).pack(fill="x", padx=8, pady=(12, 4))

        self.stats_lbl = ttk.Label(p, text="", wraplength=265, justify="left",
                                   background=BG, foreground=TEXT)
        self.stats_lbl.pack(anchor="w", padx=8, pady=4)

        # bottom spacer
        ttk.Label(p, text="", background=BG).pack(pady=12)

    def _build_plot(self, parent):
        self.fig = Figure(figsize=(9, 6), facecolor=BG)
        self.ax  = self.fig.add_subplot(111)
        self._style_ax()
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        nf = ttk.Frame(parent); nf.pack(fill="x")
        nav = NavigationToolbar2Tk(self.canvas, nf)
        nav.config(background=BG); nav.update()

    def _style_ax(self):
        self.ax.set_facecolor("#181825")
        self.fig.patch.set_facecolor(BG)
        for sp in self.ax.spines.values(): sp.set_color(OVERLAY)
        self.ax.tick_params(colors=TEXT, labelsize=9)
        self.ax.xaxis.label.set_color(TEXT)
        self.ax.yaxis.label.set_color(TEXT)
        self.ax.set_xlabel("Raman Shift (cm⁻¹)", fontsize=10)
        self.ax.set_ylabel("Intensity (a.u.)", fontsize=10)
        self.ax.grid(True, color=SURFACE, linewidth=0.5, linestyle="--", alpha=0.6)

    # ── Baseline method toggle ─────────────────────────────────────────────────

    def _bl_method_changed(self, *_):
        m = self.bl_method_var.get()
        # Always forget both, then show the right one inside the fixed container
        self.als_frame.pack_forget()
        self.snip_frame.pack_forget()
        self.arpls_note.pack_forget()
        if m == "SNIP":
            self.snip_frame.pack(fill="x")
        else:
            self.als_frame.pack(fill="x")
            if m == "arPLS":
                self.arpls_note.pack(anchor="w", padx=8, pady=(2, 4))
        self._recompute()

    # ── File I/O ───────────────────────────────────────────────────────────────

    def _open_file(self):
        path = filedialog.askopenfilename(title="Open Raman CSV",
            filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if path: self._load_file(path)

    def _open_folder(self):
        folder = filedialog.askdirectory(title="Open folder of Raman CSVs")
        if not folder: return
        files = sorted(f for f in os.listdir(folder) if f.lower().endswith(".csv"))
        if not files:
            messagebox.showwarning("Empty", "No CSV files found."); return
        self._folder_path  = folder
        self._folder_files = files

        win = tk.Toplevel(self); win.title("Select file"); win.configure(bg=BG); win.grab_set()
        ttk.Label(win, text=f"Folder: {os.path.basename(folder)}",
                  style="Header.TLabel").pack(padx=12, pady=(10, 2), anchor="w")
        lb = tk.Listbox(win, bg=SURFACE, fg=TEXT, selectbackground=BLUE, selectforeground=BG,
                        font=("Segoe UI", 9), width=58, height=min(len(files), 18),
                        relief="flat", borderwidth=0)
        lb.pack(padx=12, pady=4)
        for f in files: lb.insert("end", f)
        lb.selection_set(0)

        def _open():
            idx = lb.curselection()
            if idx:
                self._load_file(os.path.join(folder, files[idx[0]]))
                win.destroy()
        ttk.Button(win, text="Open selected", style="Accent.TButton",
                   command=_open).pack(padx=12, pady=(4, 12))

    def _load_file(self, path):
        try:
            wl, intensity = load_csv(path)
            shift = nm_to_shift(wl, self.laser_var.get())
            self.raw_wl        = wl
            self.raw_shift     = shift
            self.raw_intensity = intensity
            self.filepath      = path
            self.fit_result    = None
            self.fit_results   = []
            self.file_label.config(text=os.path.basename(path))
            self._recompute()
        except Exception as e:
            messagebox.showerror("Load error", str(e))

    def _export_csv(self):
        if self.processed is None:
            messagebox.showwarning("No data", "Load a file first."); return
        path = filedialog.asksaveasfilename(defaultextension=".csv",
            filetypes=[("CSV", "*.csv")], initialfile="processed_raman.csv")
        if not path: return
        mask = self._mask()
        pd.DataFrame({
            "Raman_Shift_cm1":     self.raw_shift[mask],
            "Raw_Intensity":       self.raw_intensity[mask],
            "Cleaned_Intensity":   self.clean_intensity[mask] if self.clean_intensity is not None else np.nan,
            "Baseline":            self.bl_arr[mask],
            "Processed_Intensity": self.processed[mask],
        }).to_csv(path, index=False)
        messagebox.showinfo("Exported", f"Saved to {path}")

    def _save_plot(self):
        if self.raw_shift is None:
            messagebox.showwarning("No data", "Load a file first."); return
        default = os.path.splitext(os.path.basename(self.filepath))[0] if self.filepath else "Raman Spectrum"
        dlg = ExportDialog(self, default_title=default)
        self.wait_window(dlg)
        if dlg.result:
            title_text, fmt, path = dlg.result
            self._render_export(title_text, fmt, path)

    def _batch_export(self):
        if not self._folder_files:
            self._open_folder()
            if not self._folder_files: return

        win = tk.Toplevel(self); win.title("Batch Export"); win.configure(bg=BG); win.grab_set()
        ttk.Label(win, text="Batch Export", style="Header.TLabel").pack(padx=14, pady=(12, 4))
        ttk.Label(win, text=f"{len(self._folder_files)} files in queue").pack(padx=14, anchor="w")

        ttk.Label(win, text="Output folder:").pack(padx=14, anchor="w", pady=(8, 0))
        rf = ttk.Frame(win); rf.pack(padx=14, fill="x")
        out_var = tk.StringVar()
        ttk.Entry(rf, textvariable=out_var, width=38).pack(side="left", padx=(0, 4))
        ttk.Button(rf, text="Browse", command=lambda: out_var.set(
            filedialog.askdirectory(parent=win) or out_var.get())).pack(side="left")

        fmt_var = tk.StringVar(value="PNG")
        ff = ttk.Frame(win); ff.pack(padx=14, pady=(8, 4), anchor="w")
        for lbl, val in [("PNG", "PNG"), ("SVG", "SVG"), ("CSV (data)", "CSV")]:
            ttk.Radiobutton(ff, text=lbl, variable=fmt_var, value=val).pack(side="left", padx=(0, 12))

        prog_lbl = ttk.Label(win, text=""); prog_lbl.pack(padx=14)
        prog = ttk.Progressbar(win, length=360, mode="determinate"); prog.pack(padx=14, pady=(0, 4))
        run_btn = ttk.Button(win, text="▶  Run", style="Accent.TButton"); run_btn.pack(padx=14, pady=(4, 12))

        def _run():
            folder_out = out_var.get()
            if not folder_out or not os.path.isdir(folder_out):
                messagebox.showwarning("No folder", "Pick a valid output folder.", parent=win); return
            fmt    = fmt_var.get()
            files  = self._folder_files
            errors = []
            prog["maximum"] = len(files)
            for i, fname in enumerate(files):
                prog_lbl.config(text=f"Processing {i+1}/{len(files)}: {fname}"); win.update()
                try:
                    x, raw, bl, c = self._process_file(os.path.join(self._folder_path, fname))
                    stem = os.path.splitext(fname)[0]
                    if fmt == "CSV":
                        pd.DataFrame({"Raman_Shift_cm1": x, "Raw": raw,
                                      "Baseline": bl, "Processed": c}).to_csv(
                            os.path.join(folder_out, f"{stem}_processed.csv"), index=False)
                    else:
                        self._render_export(stem, fmt.lower(),
                            os.path.join(folder_out, f"{stem}.{fmt.lower()}"),
                            x=x, raw=raw, bl=bl, proc=c, silent=True)
                except Exception as e:
                    errors.append(f"{fname}: {e}")
                prog["value"] = i + 1
            msg = f"Done.  {len(files)-len(errors)}/{len(files)} processed."
            if errors: msg += "\n\nErrors:\n" + "\n".join(errors[:6])
            messagebox.showinfo("Batch complete", msg, parent=win)
            win.destroy()

        run_btn.config(command=_run)

    # ── Core processing pipeline ───────────────────────────────────────────────

    def _mask(self):
        return (self.raw_shift >= self.xmin_var.get()) & (self.raw_shift <= self.xmax_var.get())

    def _process_file(self, path):
        """Run full pipeline on a file. Returns (x, raw, bl, processed)."""
        wl, intensity = load_csv(path)
        x_all = nm_to_shift(wl, self.laser_var.get())
        mask  = (x_all >= self.xmin_var.get()) & (x_all <= self.xmax_var.get())
        x     = x_all[mask]
        raw   = intensity[mask].astype(float)

        # Spike removal
        y = raw.copy()
        if self.spike_on.get():
            y, _ = remove_spikes(y, threshold=self.spike_thr_var.get())

        # Baseline
        bl = self._compute_baseline(x, y)
        c  = y - bl

        # Smoothing
        c = self._smooth(c)

        # Normalise
        c = self._normalise(c, x)
        return x, raw, bl, c

    def _compute_baseline(self, x, y):
        m = self.bl_method_var.get()
        if m == "ALS":
            return baseline_als(y, lam=10**self.lam_var.get(), p=self.p_var.get())
        elif m == "arPLS":
            return baseline_arpls(y, lam=10**self.lam_var.get(), ratio=self.p_var.get())
        elif m == "SNIP":
            return baseline_snip(y, max_iter=self.snip_iter_var.get())
        return np.zeros_like(y)

    def _smooth(self, c):
        sm = self.smooth_var.get()
        if sm == "Savitzky-Golay":
            win = int(self.win_var.get()); win = win if win % 2 == 1 else win + 1
            win = max(win, self.ord_var.get() + 2)
            try: return savgol_filter(c, win, int(self.ord_var.get()))
            except: return c
        elif sm == "Moving Average":
            win = int(self.win_var.get())
            return np.convolve(c, np.ones(win) / win, mode="same")
        return c

    def _normalise(self, c, x):
        nm = self.norm_var.get()
        if nm == "Max = 1":
            m = np.max(np.abs(c)); return c / m if m > 0 else c
        elif nm == "Area = 1":
            a = np.trapz(np.abs(c), x); return c / a if a > 0 else c
        elif nm == "SNV":
            return (c - np.mean(c)) / (np.std(c) + 1e-12)
        return c

    def _recompute(self, *_):
        if self.raw_shift is None: return
        mask = self._mask()
        x    = self.raw_shift[mask]
        raw  = self.raw_intensity[mask].astype(float)

        # Spike removal
        n_spikes = 0
        y = raw.copy()
        if self.spike_on.get():
            y, n_spikes = remove_spikes(y, threshold=self.spike_thr_var.get())
        self.clean_intensity = np.full_like(self.raw_intensity, np.nan, dtype=float)
        self.clean_intensity[mask] = y

        if n_spikes > 0:
            self.spike_info_lbl.config(text=f"↳ {n_spikes} spike(s) removed")
        else:
            self.spike_info_lbl.config(text="↳ No spikes detected")

        # Baseline
        try:
            bl = self._compute_baseline(x, y)
        except Exception as e:
            bl = np.zeros_like(y)
            self.stats_lbl.config(text=f"Baseline error: {e}")

        c = y - bl
        self.bl_arr = np.zeros_like(self.raw_intensity, dtype=float)
        self.bl_arr[mask] = bl

        c = self._smooth(c)
        c = self._normalise(c, x)

        self.processed = np.full_like(self.raw_intensity, np.nan, dtype=float)
        self.processed[mask] = c

        # Peaks
        if self.peaks_on.get() and len(c) > 5:
            h    = (self.ph_var.get() / 100.0) * np.nanmax(c)
            dx   = np.mean(np.diff(x)) if len(x) > 1 else 1
            dist = max(1, int(self.pd_var.get() / dx))
            pidx, _ = find_peaks(c, height=h, distance=dist)
            self.peaks_x = x[pidx]; self.peaks_y = c[pidx]
        else:
            self.peaks_x = np.array([]); self.peaks_y = np.array([])

        snr = np.nanmax(c) / (np.std(c[int(0.85 * len(c)):]) + 1e-12)
        self.stats_lbl.config(
            text=f"Points in range: {mask.sum()}\n"
                 f"Peak intensity:  {np.nanmax(c):.4f}\n"
                 f"Peaks detected:  {len(self.peaks_x)}\n"
                 f"Est. SNR:        {snr:.1f}\n"
                 f"Range: {x.min():.0f}–{x.max():.0f} cm⁻¹")

        self._update_plot()

    def _update_plot(self, *_):
        if self.raw_shift is None: return
        mask = self._mask()
        x    = self.raw_shift[mask]
        raw  = self.raw_intensity[mask]
        bl   = self.bl_arr[mask]
        proc = self.processed[mask]

        self.ax.cla(); self._style_ax()
        _draw_spectrum(self.ax, x, raw, bl, proc,
                       self.peaks_x, self.peaks_y,
                       self.show_raw.get(), self.show_bl.get(),
                       fit_result=self.fit_result,
                       fit_results=self.fit_results if self.fit_results else None)

        # Draw library match overlays
        for lx, ly, label, colour in getattr(self, "_match_overlays", []):
            self.ax.plot(lx, ly, color=colour, linewidth=1.0, alpha=0.65,
                         linestyle="--", label=label, zorder=6)
            # Mark matched peak pairs with vertical lines
        # Mark matched peak pairs from top hit (dashed lines connecting peaks)
        if hasattr(self, "_last_match_pairs") and self._last_match_pairs:
            for qp, lp in self._last_match_pairs:
                self.ax.axvspan(min(qp, lp), max(qp, lp),
                                alpha=0.05, color=MAUVE, zorder=3)
                self.ax.axvline(qp, color=MAUVE, linewidth=0.6, alpha=0.4, zorder=4)

        title = os.path.basename(self.filepath) if self.filepath else "Raman Spectrum"
        self.ax.set_title(title, color=TEXT, fontsize=10, pad=8)
        self.ax.legend(facecolor=SURFACE, edgecolor=OVERLAY, labelcolor=TEXT, fontsize=8)
        self.canvas.draw()

    # ── Peak fitting ───────────────────────────────────────────────────────────

    # ── Library matching ──────────────────────────────────────────────────────

    def _load_library(self):
        """Load a RamanBase library JSON file."""
        path = filedialog.askopenfilename(
            title="Load RamanBase library",
            filetypes=[("JSON library", "*.json"), ("All files", "*.*")])
        if not path:
            return
        try:
            import json
            with open(path, encoding="utf-8") as f:
                data = json.load(f)
            # Accept both full library format and flat list
            entries = data.get("spectra", data.get("records", data))                 if isinstance(data, dict) else data
            if not isinstance(entries, list):
                raise ValueError("Not a valid library file")
            # Keep only entries that have spectral data
            self._library = [e for e in entries
                             if e.get("spectral_data") and len(e.get("spectral_data") or []) > 10]
            n = len(self._library)
            self.lib_info_lbl.config(text=f"{n:,} spectra", foreground=GREEN)
            self.match_status_lbl.config(
                text=f"Library loaded: {n:,} spectra with data.", foreground=TEXT)
        except Exception as e:
            messagebox.showerror("Library load error", str(e))

    def _run_library_match(self):
        """Run peak + optional cosine matching against the loaded library."""
        if not self._library:
            messagebox.showwarning("No library", "Load a library first.")
            return
        if self.processed is None:
            messagebox.showwarning("No spectrum", "Load a spectrum first.")
            return
        if len(self.peaks_x) == 0:
            messagebox.showwarning("No peaks", "No peaks detected. Lower the min-height threshold.")
            return

        tol      = self.tol_var.get()
        frac     = self.frac_var.get() / 100.0
        top_n    = self.top_n_var.get()
        use_cos  = self.use_cosine_var.get()

        self.match_status_lbl.config(text="Matching… (may take a few seconds)", foreground=YELLOW)
        self.update_idletasks()

        mask  = self._mask()
        x_q   = self.raw_shift[mask]
        y_q   = self.processed[mask]

        import threading
        def _do_match():
            try:
                results = match_peaks_to_library(
                    self.peaks_x, self._library,
                    tolerance_cm=tol,
                    require_fraction=frac,
                    top_n=top_n * 3)   # fetch more, cosine reranks

                # Optionally add cosine similarity score and re-rank
                if use_cos:
                    for r in results:
                        lx = r["x_axis"]
                        ly = r["spectral_data"]
                        if lx and ly:
                            ly_n = np.asarray(ly, float)
                            rng = ly_n.max() - ly_n.min()
                            if rng > 1e-12:
                                ly_n = (ly_n - ly_n.min()) / rng
                            cos = cosine_similarity(x_q, y_q, np.asarray(lx), ly_n)
                            r["cosine"] = round(cos, 4)
                            # Combined score: 60% peak match + 40% cosine
                            r["combined"] = round(0.6 * r["score"] + 0.4 * cos, 4)
                        else:
                            r["cosine"]   = 0.0
                            r["combined"] = round(0.6 * r["score"], 4)
                    results.sort(key=lambda r: r["combined"], reverse=True)
                else:
                    for r in results:
                        r["cosine"]   = None
                        r["combined"] = r["score"]

                results = results[:top_n]
                self.after(0, lambda: self._show_match_results(results, x_q, y_q, tol))
            except Exception as e:
                self.after(0, lambda: (
                    self.match_status_lbl.config(
                        text=f"Match error: {e}", foreground=RED),
                    messagebox.showerror("Match error", str(e))))

        threading.Thread(target=_do_match, daemon=True).start()

    def _show_match_results(self, results, x_q, y_q, tol):
        """Show match results in a dialog and optionally overlay top hit on plot."""
        n = len(results)
        if n == 0:
            self.match_status_lbl.config(
                text="No matches found. Try increasing tolerance or reducing min coverage.",
                foreground=RED)
            return

        top = results[0]
        score_key = "combined"
        self.match_status_lbl.config(
            text=f"Top hit: {top['name'][:40]}  score={top[score_key]:.3f}  "
                 f"({top['n_matched']}/{top['n_query']} peaks, "
                 f"Δ̄={top['mean_error']:.1f} cm⁻¹)",
            foreground=GREEN)

        # Update plot overlays with top 3 library spectra
        self._match_overlays = []
        overlay_colours = [MAUVE, PEACH, RED]
        for i, r in enumerate(results[:3]):
            lx = r.get("x_axis")
            ly = r.get("spectral_data")
            if lx and ly:
                ly_a = np.asarray(ly, float)
                rng = ly_a.max() - ly_a.min()
                if rng > 1e-12:
                    ly_a = (ly_a - ly_a.min()) / rng
                self._match_overlays.append((
                    np.asarray(lx), ly_a,
                    f"#{i+1} {r['name'][:30]}",
                    overlay_colours[i % len(overlay_colours)]
                ))
        self._update_plot()

        # Open results dialog
        MatchResultsDialog(self, results, x_q, y_q)

    def _clear_match_overlays(self):
        self._match_overlays = []
        self._update_plot()

    # ── Peak fitting ──────────────────────────────────────────────────────────

    def _start_fit_click(self):
        """Enable click-to-fit mode."""
        if self.processed is None:
            messagebox.showwarning("No data", "Load a file first."); return
        if len(self.peaks_x) == 0:
            messagebox.showwarning("No peaks", "No peaks detected. Lower the min-height threshold."); return
        self.fit_status_lbl.config(text="Click on a peak in the plot…", foreground=GREEN)
        # Disconnect previous handler if any
        if self._click_cid is not None:
            self.canvas.mpl_disconnect(self._click_cid)
        self._click_cid = self.canvas.mpl_connect("button_press_event", self._on_plot_click)

    def _on_plot_click(self, event):
        if event.inaxes != self.ax or event.xdata is None:
            return
        self.canvas.mpl_disconnect(self._click_cid)
        self._click_cid = None

        click_x = event.xdata
        # Snap to nearest detected peak
        if len(self.peaks_x) > 0:
            nearest_idx = np.argmin(np.abs(self.peaks_x - click_x))
            centre = self.peaks_x[nearest_idx]
        else:
            centre = click_x

        mask = self._mask()
        x    = self.raw_shift[mask]
        proc = self.processed[mask]

        try:
            result = fit_peak(x, proc, centre_guess=centre,
                              window_cm=self.fit_win_var.get())
            self.fit_result = result
            self.fit_status_lbl.config(
                text=f"Fit: {result['centre']:.1f} cm⁻¹  "
                     f"FWHM {result['fwhm']:.1f} cm⁻¹  "
                     f"η={result['eta']:.2f}",
                foreground=GREEN)
            self._update_plot()
            FitResultDialog(self, result)
        except Exception as e:
            self.fit_status_lbl.config(text=f"Fit failed: {e}", foreground=RED)

    def _fit_all_peaks(self):
        """Auto-fit pseudo-Voigt to every detected peak and show results table."""
        if self.processed is None:
            messagebox.showwarning("No data", "Load a file first."); return
        if len(self.peaks_x) == 0:
            messagebox.showwarning("No peaks", "No peaks detected. Lower the min-height threshold."); return

        mask = self._mask()
        x    = self.raw_shift[mask]
        proc = self.processed[mask]
        win  = self.fit_win_var.get()

        results  = []
        failures = []
        for centre in self.peaks_x:
            try:
                r = fit_peak(x, proc, centre_guess=centre, window_cm=win)
                results.append(r)
            except Exception as e:
                failures.append(f"{centre:.0f} cm⁻¹: {e}")

        self.fit_results = results
        self.fit_result  = None   # clear single-peak highlight when doing all

        n = len(results)
        fail_txt = f"  ({len(failures)} failed)" if failures else ""
        self.fit_status_lbl.config(
            text=f"{n} peak(s) fitted{fail_txt}. See table →",
            foreground=GREEN)

        self._update_plot()
        if results:
            AllFitsDialog(self, results)

    def _clear_fit(self):
        self.fit_result  = None
        self.fit_results = []
        self.fit_status_lbl.config(text="Fits cleared.", foreground=OVERLAY)
        if self._click_cid is not None:
            self.canvas.mpl_disconnect(self._click_cid)
            self._click_cid = None
        self._update_plot()

    # ── Export renderer ────────────────────────────────────────────────────────

    def _render_export(self, title_text, fmt, path,
                        x=None, raw=None, bl=None, proc=None, silent=False):
        if x is None:
            mask = self._mask()
            x    = self.raw_shift[mask]
            raw  = self.raw_intensity[mask]
            bl   = self.bl_arr[mask]
            proc = self.processed[mask]

        px, py = np.array([]), np.array([])
        if self.peaks_on.get() and len(proc) > 5:
            h    = (self.ph_var.get() / 100.0) * np.nanmax(proc)
            dx   = np.mean(np.diff(x)) if len(x) > 1 else 1
            dist = max(1, int(self.pd_var.get() / dx))
            pidx, _ = find_peaks(proc, height=h, distance=dist)
            px, py  = x[pidx], proc[pidx]

        DPI = 100
        fig = Figure(figsize=(7, 4.5), dpi=DPI, facecolor=BG)
        ax  = fig.add_subplot(111)
        ax.set_facecolor("#181825")
        for sp in ax.spines.values(): sp.set_color(OVERLAY)
        ax.tick_params(colors=TEXT, labelsize=8)
        ax.xaxis.label.set_color(TEXT); ax.yaxis.label.set_color(TEXT)
        ax.set_xlabel("Raman Shift (cm⁻¹)", fontsize=9)
        ax.set_ylabel("Intensity (a.u.)", fontsize=9)
        ax.grid(True, color=SURFACE, linewidth=0.5, linestyle="--", alpha=0.6)

        _draw_spectrum(ax, x, raw, bl, proc, px, py,
                       self.show_raw.get(), self.show_bl.get(),
                       small=True, fit_result=self.fit_result,
                       fit_results=self.fit_results if self.fit_results else None)
        ax.set_title(title_text, color=TEXT, fontsize=10, pad=8)
        ax.legend(facecolor=SURFACE, edgecolor=OVERLAY, labelcolor=TEXT, fontsize=7)
        fig.tight_layout()

        if fmt == "svg":
            from matplotlib.backends.backend_svg import FigureCanvasSVG
            FigureCanvasSVG(fig).print_figure(path, bbox_inches="tight", facecolor=BG)
        else:
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            FigureCanvasAgg(fig).print_figure(path, dpi=DPI, bbox_inches="tight", facecolor=BG)

        if not silent:
            messagebox.showinfo("Saved", f"Exported to:\n{path}")

    # ── Slider callbacks ───────────────────────────────────────────────────────

    def _xmin_cb(self, v):
        val = int(float(v)); self.xmin_lbl.config(text=str(val))
        if val >= self.xmax_var.get(): self.xmin_var.set(self.xmax_var.get() - 50)
        self._recompute()

    def _xmax_cb(self, v):
        val = int(float(v)); self.xmax_lbl.config(text=str(val))
        if val <= self.xmin_var.get(): self.xmax_var.set(self.xmin_var.get() + 50)
        self._recompute()

    def _spike_cb(self, v):
        self.spike_thr_lbl.config(text=f"{float(v):.1f}"); self._recompute()

    def _lam_cb(self, v):
        val = float(v)
        eq  = 10 ** val
        rd  = f"≈ {eq/1000:.1f}k" if eq >= 1000 else f"≈ {eq:.0f}"
        self.lam_lbl.config(text=f"10^{val:.1f}  {rd}"); self._recompute()

    def _p_cb(self, v):
        self.p_lbl.config(text=f"{float(v):.3f}"); self._recompute()

    def _snip_iter_cb(self, v):
        self.snip_iter_lbl.config(text=str(int(float(v)))); self._recompute()

    def _win_cb(self, v):
        val = int(float(v)); val = val if val % 2 == 1 else val + 1
        self.win_lbl.config(text=str(val)); self._recompute()

    def _ord_cb(self, v):
        self.ord_lbl.config(text=str(int(float(v)))); self._recompute()

    def _ph_cb(self, v):
        self.ph_lbl.config(text=f"{int(float(v))}%"); self._recompute()

    def _pd_cb(self, v):
        self.pd_lbl.config(text=str(int(float(v)))); self._recompute()

    def _fit_win_cb(self, v):
        self.fit_win_lbl.config(text=str(int(float(v))))

    def _reset(self):
        self.xmin_var.set(650);          self.xmin_lbl.config(text="650")
        self.xmax_var.set(3300);         self.xmax_lbl.config(text="3300")
        self.spike_on.set(True)
        self.spike_thr_var.set(8.5);     self.spike_thr_lbl.config(text="8.5")
        self.bl_method_var.set("SNIP");  self._bl_method_changed()
        self.lam_var.set(3.2);           self.lam_lbl.config(text="10^3.2  ≈ 1.6k")
        self.p_var.set(0.01);            self.p_lbl.config(text="0.010")
        self.snip_iter_var.set(11);      self.snip_iter_lbl.config(text="11")
        self.smooth_var.set("Savitzky-Golay")
        self.norm_var.set("Max = 1")
        self.win_var.set(7);            self.win_lbl.config(text="7")
        self.ord_var.set(3);             self.ord_lbl.config(text="3")
        self.peaks_on.set(True)
        self.ph_var.set(12);             self.ph_lbl.config(text="12%")
        self.pd_var.set(30);             self.pd_lbl.config(text="30")
        self.fit_win_var.set(80);        self.fit_win_lbl.config(text="80")
        self.show_raw.set(True)
        self.show_bl.set(True)
        self.fit_result  = None
        self.fit_results = []
        self.fit_status_lbl.config(text="No fit yet.", foreground=OVERLAY)
        self._recompute()


if __name__ == "__main__":
    app = RamanApp()
    app.mainloop()
