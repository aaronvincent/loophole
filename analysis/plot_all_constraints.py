#!/usr/bin/env python3
"""
Combined-constraint tension plots for the Loophole scan: relic + ALL the
Stage-2 constraints (Delta-Neff BBN/CMB + SN1987A + neutron-star + 5th force) on
one figure, with TRANSPARENT excluded overlays so overlapping exclusions blend.

This is the single constraint plotter: all constraints are always computed.  The
sidecar is loaded via constraint_io.load, which reads the unified
``<base>_postprocessed.hdf5`` (built by postprocess.py -- relic + every offline
likelihood in one self-contained file).  It carries

    LogLike_relic, LogLike_neff_bbn, LogLike_neff_cmb,
    LogLike_sn1987a, LogLike_neutronstar, LogLike_stellar, LogLike_fifthforce,

and the combined likelihood is

    LogLike = relic + neff_bbn + neff_cmb + sn1987a + neutronstar + stellar
              + fifthforce.

For every requested (x, y) plane it draws the 2-sigma survivor region (relic +
ALL constraints, flat blue fill = what survives everything), then overlays, for
each constraint that actually bites, the bins that are relic-allowed at 2 sigma
but pushed out of 2 sigma by THAT constraint alone:

    { bin relic-allowed }  AND NOT  { bin (relic + C) allowed }

Each constraint is a translucent fill (alpha 0.7): colour = family of probe, hatch =
which member (see CONSTRAINTS); where several overlap the alphas stack.  Bins
killed only by the SUM of constraints (no single one alone) are drawn as every
bound in an inclusion-minimal violated set there, outlined by a dotted line (see
make_plane).  The survivor fill is labelled "Omega h^2 = 0.12, unconstrained"
along its straightest edge (place_survivor_label).

Usage
-----
    python plot_all_constraints.py [BASE_OR_SIDECAR] [GROUP] [PLANES]

    BASE_OR_SIDECAR  the run's raw scan OR its sidecar; the script strips a
                     trailing _postprocessed and locates the
                     sidecar.  (default: runs/heavyScalar_3p.hdf5)
    GROUP            hdf5 group (default: "data")
    PLANES           comma-separated "<x>:<y>" planes (scan_axes.KNOWN names).
                     Default: auto.

    python plot_all_constraints.py --legend PLOTS_DIR [PLOTS_DIR ...]

    rebuild only legend.pdf in each <run>_allconstraints dir from the
    legend.json a full run leaves there (seconds; no data read).

Env knobs:  LOOP_PLOT_ALPHA (0.7 overlay opacity), LOOP_PLOT_PNG (0 --
PDF only), LOOP_PL_GRID / LOOP_PL_BWFRAC / LOOP_PL_KPTS / LOOP_PL_MINAREA
(region engine; see the config block), LOOP_SIGMA_YMIN / LOOP_SIGMA_YMAX
(log10 cm^2: widen the sigma_chiN window to at least this range; unset = the
relic-allowed percentiles only).
"""

import json
import math
import os
import re
import sys

import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib

import scan_axes
import pl_regions
import constraint_io

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import patheffects
from matplotlib.colors import to_rgba
from matplotlib.ticker import FixedLocator, FuncFormatter

import plot_style
plot_style.use_latex_fonts()
# The panels are shrunk to ~0.44\textwidth in the multi-panel LaTeX grids, so the
# default sizes are illegible there.  Push the base sizes up (env-tunable via
# LOOP_PLOT_FONTSCALE) so axis labels/ticks stay readable at grid scale.
_FS = float(os.environ.get("LOOP_PLOT_FONTSCALE", "1.0"))
matplotlib.rcParams.update({
    "font.size":       17 * _FS,
    "axes.labelsize":  24 * _FS,
    "axes.titlesize":  24 * _FS,
    "xtick.labelsize": 19 * _FS,
    "ytick.labelsize": 19 * _FS,
    "legend.fontsize": 19 * _FS,
})

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))

def vertex_label(interaction):
    """LaTeX vertex label for the in-plot title, e.g. 'heavyScalarSP' -> '$S \\otimes P$'.

    The interaction suffix orders the two factors as (phi-chi-chi, phi-SM); plain
    scalar (no suffix) is S(x)S.  Falls back to the raw string for axialVector etc.
    """
    s = interaction
    for pre in ("heavy", "light"):
        if s.startswith(pre):
            s = s[len(pre):]
            break
    if s.lower().startswith("axialvector"):
        return r"$A \otimes A$"
    if s.lower().startswith("scalar"):
        suf = s[len("scalar"):] or "SS"          # '' -> scalar(x)scalar
        if len(suf) == 2 and all(c in "SP" for c in suf):
            return rf"${suf[0]} \otimes {suf[1]}$"
    return interaction


alpha =float(os.environ.get("LOOP_PLOT_ALPHA", "0.7"))
save_png = os.environ.get("LOOP_PLOT_PNG", "0") not in ("0", "", "false", "False")
# Density level-set engine knobs (see pl_regions.point_density).  GRID is the fine
# raster resolution.  BW_FRAC is the gaussian bandwidth as a FRACTION of the axis
# span; pl_regions wants it in fine cells, so BW is derived.  Anchoring the
# bandwidth to the axis (not the raster) makes GRID a pure rendering knob: raising
# it traces a smoother boundary without shrinking the kernel, which would silently
# cut the samples-per-kernel that KPTS thresholds and pull the region inward.
# Small BW_FRAC keeps sharp features; raise KPTS to reject more fringe/outlier
# points; MIN_AREA prunes speckle islands/pinholes.
GRID = int(os.environ.get("LOOP_PL_GRID", "600"))
BW_FRAC = float(os.environ.get("LOOP_PL_BWFRAC", "0.00625"))
BW = BW_FRAC * GRID
KPTS = float(os.environ.get("LOOP_PL_KPTS", "2.0"))
MIN_AREA = float(os.environ.get("LOOP_PL_MINAREA", "0.0008"))
# LOOP_PL_XENV = x-bin width in dex: replace the density engine by the per-column
# sample envelope (pl_regions.xenvelope_region), for scans whose whole plane is
# one thin band (MFV, fixed m_phi).  Draws it at its true width; wrong for 2-D
# regions with holes.  Unset / 0 = density engine.
XENV = float(os.environ.get("LOOP_PL_XENV", "0") or 0)

CL_2SIG = -3.090

# The 2-sigma survivor region is a single flat, opaque fill (navy), outlined by a
# darker navy envelope (SURVIVOR_EDGE).  No constraint family is blue, so blue
# always means "allowed".
SURVIVOR_BLUE = "#08519c"
SURVIVOR_EDGE = "#08306b"

# Constraint registry: key -> (colour, label, hatch).  Source/datasets live in
# constraint_io now; this only carries plot styling.
#
# Colour is the FAMILY of probe; within it, one member is solid and every other
# member has its own hatch, unique across all constraints:
#   stars            yellow      SN1987A (solid), neutron star //, stellar cooling \\
#   lab              teal        meson decays (solid), 5th force --
#                                (meson decays = K -> pi chi chi, K -> pi phi and
#                                pi0/eta/eta' -> invisible, merged into one region
#                                by constraint_io.MERGED)
#   early universe   vermillion  Delta Neff CMB (solid), CMB <sigma v> .., Delta Neff BBN oo
#   DM astrophysics  wine        bullet cluster __ (horizontal dashes), INTEGRAL **
# So a solid fill is identified by its colour, a hatched one by its pattern alone.
# The wine family has no solid member: a solid bullet fill was too close to the
# starred INTEGRAL one.  pl_regions.DASHES is drawn by hand (not a matplotlib hatch).
# xx and ++ are avoided: overlapping // and \\ (or || and --) already draw them.
# Ten translucent fills cannot all be told apart by colour alone.  Four family
# hues can: drawn at LOOP_PLOT_ALPHA = 0.7 over white, next to the opaque navy
# survivor fill and the joint grey, every pair is >= 15.3 apart in OKLab x100
# for normal vision and >= 10.2 under protan/deutan simulation.  Hatches are drawn
# in a darker shade of the family colour and survive greyscale.
CONSTRAINTS = {
    "neff_bbn":    ("#D55E00", r"$\Delta N_\mathrm{eff}$ (BBN)", "oo"),
    "neff_cmb":    ("#D55E00", r"$\Delta N_\mathrm{eff}$ (CMB)", None),
    "sn1987a":     ("#F0E442", "SN1987A", None),
    "neutronstar": ("#F0E442", "Neutron star", "//"),
    "stellar":     ("#F0E442", "Stellar cooling", "\\\\"),
    "fifthforce":  ("#009988", "5th force", "--"),
    "bullet":      ("#882255", "Bullet cluster", pl_regions.DASHES),
    "meson":       ("#009988", "Meson decays", None),
    "cmb":         ("#D55E00", r"CMB ($p_\mathrm{ann}$)", ".."),
    "integral":    ("#882255", "INTEGRAL", "**"),
}


def hatch_color(color, shade=0.55):
    """Hatch lines: the fill colour darkened, so a hatch reads as the same family."""
    r, g, b, _ = to_rgba(color)
    return (r * shade, g * shade, b * shade, 1.0)
# Joint kills (relic-allowed, killed by the full likelihood, but by no single
# constraint alone) are drawn as the constraints that CONTRIBUTE there -- those
# that remove some of the cell's relic-allowed points -- in their usual styles,
# with the joint region outlined by a dotted line (JOINT_COLOR).  A plane cell is "allowed by C" if ANY of its points
# (profiled over the hidden axes) passes relic + C, so with hard cuts a joint cell
# is one where every point fails some bound but different points fail different
# ones -- typically SN1987A takes the light-mphi points and K -> pi phi the heavy.
# Only the soft terms (relic, Delta Neff) can still kill by summing penalties.
JOINT_COLOR = "#333333"
SPLIT = None   # survivor split; set below from LOOP_SURV_SPLIT (unset in --legend mode)
JOINT_LABEL = "Several bounds jointly"


def write_legend(active, include_combined, out_dir):
    """Write the shared legend (identical across all planes) to its own file.

    Entries: survives-all, one
    translucent (hatched) patch per active constraint, and the dotted joint-kill
    outline when any plane produced that category.  Rendered into a
    legend-only figure cropped to the legend's bounding box.  Patch labels are
    bare constraint names: a filled patch always means "excluded by".

    Those two inputs are all the legend depends on, so they are saved beside it
    as legend.json; ``plot_all_constraints.py --legend <dir>...`` rebuilds the
    legend from that alone, without loading the sidecar or redrawing a plane.
    """
    handles = ([Patch(fc=SURVIVOR_BLUE, label=r"Survives all ($2\sigma$)")] if SPLIT is None else
               [Patch(fc=col, label=rf"Survives all ($2\sigma$): {name}")
                for name, col in zip(SPLIT["names"], SPLIT["colors"])])
    handler_map = {}
    for n in active:
        c, lab, h = CONSTRAINTS[n]
        dashes = h == pl_regions.DASHES
        p = Patch(fc=to_rgba(c, alpha), hatch=None if dashes else h,
                  ec=hatch_color(c) if h and not dashes else "none", lw=0, label=lab)
        if dashes:
            handler_map[p] = pl_regions.HandlerDashes(hatch_color(c))
        handles.append(p)
    if include_combined:
        handles.append(Line2D([], [], color=JOINT_COLOR, lw=1.3, ls=":", label=JOINT_LABEL))
    fig = plt.figure(figsize=(11, 1.6))
    leg = fig.legend(handles=handles, loc="center", ncol=3, frameon=True,
                     framealpha=1.0, handler_map=handler_map)
    fig.canvas.draw()
    bbox = leg.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    p = os.path.join(out_dir, "legend.pdf")
    fig.savefig(p, bbox_inches=bbox.expanded(1.05, 1.15))
    if save_png:
        fig.savefig(p.replace(".pdf", ".png"), dpi=200,
                    bbox_inches=bbox.expanded(1.05, 1.15))
    plt.close(fig)
    with open(os.path.join(out_dir, "legend.json"), "w") as f:
        json.dump({"active": list(active), "joint": bool(include_combined)}, f)
    print(f"  wrote {os.path.basename(p)}{' (+ .png)' if save_png else ''}")


# Legend-only mode: relabel/restyle without touching the data.
if len(sys.argv) > 1 and sys.argv[1] == "--legend":
    for d in sys.argv[2:]:
        with open(os.path.join(d, "legend.json")) as f:
            spec = json.load(f)
        print(d)
        write_legend(spec["active"], spec["joint"], d)
    sys.exit(0)

if len(sys.argv) < 2:
    sys.exit("usage: plot_all_constraints.py SCAN.hdf5 [group] [planes]")
in_arg = sys.argv[1]
group_name = sys.argv[2] if len(sys.argv) > 2 else "data"
planes_arg = sys.argv[3] if len(sys.argv) > 3 else None

# Locate + load the constraint sidecar <base>_postprocessed.hdf5.  Accepts the
# raw scan or the sidecar path as the argument.
try:
    C = constraint_io.load(in_arg, group_name)
except FileNotFoundError as e:
    sys.exit(str(e))

run_tag = os.path.basename(C.base)
# Optional survivor split.  LOOP_SURV_SPLIT names an .npz with per-row int
# `labels` aligned to the sidecar, and `names` / `colors` per label; label 0 is
# the main survivor region.  Each component gets its OWN density region, so
# features the smoothing would otherwise merge stay apart.  LOOP_PLOT_SUFFIX
# names a separate output dir, so a split figure never replaces the standard one.
SPLIT = None
if os.environ.get("LOOP_SURV_SPLIT"):
    _z = np.load(os.environ["LOOP_SURV_SPLIT"])
    SPLIT = dict(labels=_z["labels"], names=[str(v) for v in _z["names"]],
                 colors=[str(v) for v in _z["colors"]])
out_dir = os.path.join(HERE, "plots", run_tag + "_allconstraints"
                       + os.environ.get("LOOP_PLOT_SUFFIX", ""))
title_tag = run_tag.replace("_", " ")

interaction = run_tag.split("_")[0]
overlay_title = vertex_label(interaction)
if "mphi10" in run_tag.split("_"):
    overlay_title += "\n" + r"$m_\phi = 10$ GeV"

# -------------------------------------------------------------------------
# Loaded arrays (full-length / unfiltered; constraint_io applied the valid mask
# to relic+total finiteness only).
# -------------------------------------------------------------------------
axis_vals = C.axis_raw
avail = C.avail
valid = C.valid
relic_LL = C.relic_LL
comp = C.comp
total_LL = C.total_LL
interaction = C.interaction
print(f"Reading: {C.source}  ({C.kind})")
print(f"interaction:  {interaction}")
print(f"Valid points: {valid.sum():,} of {valid.size:,}")
if SPLIT is not None and SPLIT["labels"].size != valid.size:
    sys.exit(f"LOOP_SURV_SPLIT has {SPLIT['labels'].size:,} labels for {valid.size:,} rows")

# Which constraints actually bite (any penalty on a valid point).  The two-sided
# Neff Gaussians penalize ~every point a little, so report the strongly-penalized
# fraction (< -0.5) -- the 2-sigma kill attribution below is what the figure uses.
active = C.active
print(f"Active constraints: {active or '(none)'}")
for n in active:
    print(f"  {n}: {(comp[n][valid] < -0.5).mean():.3%} of valid points "
          f"strongly penalized (<-0.5)")

# Axis machinery.
AX = {k: (avail[k][1], avail[k][2]) for k in avail}
derived = scan_axes.derived_available(avail)
for dk, (_srcs, dmode, dlabel) in derived.items():
    axis_vals[dk] = scan_axes.compute_derived(dk, axis_vals)
    AX[dk] = (dmode, dlabel)
log_axis = {k: scan_axes.transform(axis_vals[k], AX[k][0]) for k in AX}

varied = scan_axes.varied_axes({k: axis_vals[k][valid] for k in avail})
have_sigma = "sigma_chiN" in avail
if planes_arg:
    planes = planes_arg.split(",")
    for pl in planes:
        xk, yk = pl.split(":")
        assert xk in AX and yk in AX, f"plane {pl!r}: axis missing; have {list(AX)}"
else:
    planes = [f"{x}:{y}" for x, y in scan_axes.default_pairs(varied, have_sigma)]
    for dk in derived:
        if scan_axes._is_varied(axis_vals[dk][valid]):
            for c in ("gx",) + tuple(m for m in scan_axes.MASSES if m in varied):
                if c in AX:
                    planes.append(f"{c}:{dk}")

BOUNDS = {k: scan_axes.auto_bounds(log_axis[k][valid]) for k in AX}
print(f"Varied axes:  {varied}")
print(f"Planes:       {planes}")

relic_max = relic_LL[valid].max()
full_max = total_LL[valid].max()
full_c_LL = {n: relic_LL + comp[n] for n in active}             # relic + one constraint
full_c_max = {n: full_c_LL[n][valid].max() for n in active}
# Per point: which constraints push it out of 2 sigma on their own (bit i ->
# active[i]), and how many.  Used to attribute joint-only kills.
viol_bits = np.zeros(valid.shape, dtype=np.uint16)
for i, n in enumerate(active):
    viol_bits |= (full_c_LL[n] < full_c_max[n] + CL_2SIG).astype(np.uint16) << i
n_viol = np.zeros(valid.shape, dtype=np.uint8)
for i in range(len(active)):
    n_viol += ((viol_bits >> i) & 1).astype(np.uint8)
os.makedirs(out_dir, exist_ok=True)


# -------------------------------------------------------------------------
# Log-axis presentation.  Data are plotted in log10 (the region engine works in
# log space), so a "log" axis is really a LINEAR axis of log10 values.  We dress
# it as a proper log plot: major ticks 10^n, decade minor ticks, and a label
# with the log10(...) wrapper stripped (the 10^n ticks carry the scaling now).
# -------------------------------------------------------------------------
def delog_label(label):
    """'$\\log_{10}(m_\\chi/\\mathrm{GeV})$' -> '$m_\\chi/\\mathrm{GeV}$'."""
    inner = label.strip("$")
    m = re.match(r"\\log_\{10\}\s*\((.*)\)$", inner) or \
        re.match(r"\\log_\{10\}\s*(.*)$", inner)
    return f"${m.group(1)}$" if m else label


def style_log_axis(mpl_axis, lo, hi, target=7):
    """Dress a (linear-in-log10) axis as a log plot: 10^n major ticks at a sane
    decade step, with minor ticks below.  ``sigma_chiN`` can span ~50 decades, so
    a tick on every decade is illegible -- pick a nice step (1,2,5,10,...) that
    keeps ~``target`` labelled ticks.  Sub-decade (2..9) minors only when the step
    is one decade; otherwise minor ticks mark each unlabelled decade."""
    raw = (hi - lo) / target
    step = next((s for s in (1, 2, 5, 10, 20, 50, 100) if raw <= s), 100)
    start = step * math.ceil(lo / step)
    majors = list(range(start, math.floor(hi) + 1, step))
    mpl_axis.set_major_locator(FixedLocator(majors))
    mpl_axis.set_major_formatter(
        FuncFormatter(lambda v, _pos: r"$10^{%d}$" % int(round(v))))
    if step == 1:
        minors = [n + math.log10(k)
                  for n in range(math.floor(lo) - 1, math.ceil(hi) + 1)
                  for k in range(2, 10)]
    else:
        minors = list(range(math.floor(lo), math.ceil(hi) + 1))
    mpl_axis.set_minor_locator(FixedLocator([m for m in minors if lo <= m <= hi]))


# -------------------------------------------------------------------------
# In-region label for the survivor fill
# -------------------------------------------------------------------------
SURVIVOR_LABEL = r"$\Omega h^2 = 0.12$, unconstrained"


def place_survivor_label(ax, mask, xf, yf, segs, text, fontsize, avoid=(),
                         color="white"):
    """Write `text` in `mask` -- the survivor region (white text), or for the
    fallback the empty space around it (navy text) -- parallel to the straightest
    stretch of the survivor boundary.

    Works in DISPLAY coordinates, so the angle is the one the reader sees.  A
    window as long as the label (plus margin) slides along every boundary polyline;
    a line is fitted to it and the window is scored by its worst perpendicular
    deviation.  The label sits just inside the region along the fitted line's
    inward normal, and a position is kept only if the whole rotated text box lies
    in `mask` and clears the `avoid` artists.  Returns the Text, or None if nothing
    fits (e.g. a survivor region too thin for the label)."""
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    probe = ax.text(0, 0, text, fontsize=fontsize, color=color)
    bb = probe.get_window_extent(rend)
    probe.remove()
    tw, th = bb.width, bb.height
    to_disp, to_data = ax.transData.transform, ax.transData.inverted().transform
    dx, dy = xf[1] - xf[0], yf[1] - yf[0]

    def inside(pts_disp):
        d = to_data(pts_disp)
        ix = np.round((d[:, 0] - xf[0]) / dx).astype(int)
        iy = np.round((d[:, 1] - yf[0]) / dy).astype(int)
        ok = (ix >= 0) & (ix < mask.shape[1]) & (iy >= 0) & (iy < mask.shape[0])
        return ok.all() and mask[iy, ix].all()

    avoid_bbs = [a.get_window_extent(rend) for a in avoid]
    gu, gv = np.meshgrid(np.linspace(-0.5, 0.5, 25), np.linspace(-0.5, 0.5, 7))
    box = np.column_stack([gu.ravel(), gv.ravel()])        # unit text box, sampled
    L = 1.15 * tw                                           # window length along the edge
    # Pass 1 sits the label a fixed distance in from the edge.  Pass 2, only when
    # that fits nowhere, slides it across a thin band along the normal (up to the
    # band's midline and beyond) and takes the nearest offset that fits.
    for offsets in (None, np.arange(0.6, 6.0, 0.1) * th):
        best = _best_label_window(segs, to_disp, inside, avoid_bbs, box, tw, th, L, offsets)
        if best is not None:
            break
    if best is None:
        return None
    x, y = to_data(best[1][None, :])[0]
    return ax.text(x, y, text, fontsize=fontsize, color=color, rotation=best[2],
                   rotation_mode="anchor", ha="center", va="center", zorder=6)


def _best_label_window(segs, to_disp, inside, avoid_bbs, box, tw, th, L, offsets):
    """Straightest boundary window (see place_survivor_label) whose label box
    fits: (deviation, centre in display coords, angle), or None.  `offsets` None
    means the fixed pass-1 distance from the edge."""
    best = None
    for seg in segs:
        P = to_disp(seg)
        arc = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(P, axis=0).T))])
        if arc[-1] < L:
            continue
        for a0 in np.arange(0, arc[-1] - L, max(tw / 20, 2.0)):
            w = P[(arc >= a0) & (arc <= a0 + L)]
            if len(w) < 4:
                continue
            c = w.mean(axis=0)
            _, _, vt = np.linalg.svd(w - c)
            t = vt[0] if vt[0][0] >= 0 else -vt[0]            # along-edge, left to right
            nrm = np.array([-t[1], t[0]])
            dev = np.abs((w - c) @ nrm).max()
            ang = np.degrees(np.arctan2(t[1], t[0]))
            if abs(ang) > 60:
                continue
            if best is not None and dev >= best[0]:
                continue
            ds = [0.5 * th + 0.35 * th + dev] if offsets is None else offsets
            for side in (1, -1):                              # inward = the side in the mask
                for dd in ds:
                    ctr = c + side * nrm * dd
                    pts = ctr + (box * [tw * 1.04, th * 1.2]) @ np.column_stack([t, nrm]).T
                    if not inside(pts):
                        continue
                    if any(b.contains(*q) for b in avoid_bbs for q in pts):
                        continue
                    best = (dev, ctr, ang)
                    break
    return best


# -------------------------------------------------------------------------
# Per-plane profiling + plotting
# -------------------------------------------------------------------------
def _env_float(name):
    v = os.environ.get(name, "").strip()
    return float(v) if v else None


# log10(cm^2) limits the sigma_chiN window is widened to (never narrowed).
SIGMA_YMIN = _env_float("LOOP_SIGMA_YMIN")
SIGMA_YMAX = _env_float("LOOP_SIGMA_YMAX")


def make_plane(x_key, y_key):
    lx, ly = log_axis[x_key], log_axis[y_key]
    fin = valid & np.isfinite(lx) & np.isfinite(ly)

    xb, yb = BOUNDS[x_key], BOUNDS[y_key]
    for key, bb, lv in ((x_key, "x", lx), (y_key, "y", ly)):
        if key == "sigma_chiN":
            sel = lv[fin & (relic_LL > CL_2SIG)]
            if sel.size:
                lo = np.floor(np.nanpercentile(sel, 0.1)) - 1
                hi = np.ceil(np.nanpercentile(sel, 99.9)) + 1
                # Opt-in widening for overlays (e.g. direct-detection curves)
                # that sit above/below the relic-allowed band: only extends.
                if SIGMA_YMIN is not None:
                    lo = min(lo, SIGMA_YMIN)
                if SIGMA_YMAX is not None:
                    hi = max(hi, SIGMA_YMAX)
                if bb == "x":
                    xb = (lo, hi)
                else:
                    yb = (lo, hi)

    xlabel = delog_label(AX[x_key][1]) if AX[x_key][0] == "log" else AX[x_key][1]
    ylabel = delog_label(AX[y_key][1]) if AX[y_key][0] == "log" else AX[y_key][1]

    # --- Density level-set regions from the cut-passing samples ---------------
    # A location (x,y) is allowed iff SOME sample there passes the 2-sigma cut,
    # which per sample is exactly  LL >= LL_max + CL_2SIG  (the binned profile-
    # max threshold, but with no binning and no max statistic).  Each region is a
    # level set of the SMOOTHED DENSITY of its cut-passing points: weighting by
    # sample count (not a binary hull) lets lone outliers wash out below the
    # threshold without a heavy blur, so smooth edges stay smooth and sharp real
    # features (dense ridge/fingers) are kept.  Set differences/intersections of
    # the rasterised region masks give the same kill categories as before; the
    # density field + threshold are returned so the boundary LINE traces the
    # smooth field, matching the fill edge.
    relic_pt = fin & (relic_LL >= relic_max + CL_2SIG)
    full_pt = fin & (total_LL >= full_max + CL_2SIG)

    def region(pt):
        if XENV:
            return pl_regions.xenvelope_region(lx[pt], ly[pt], xb, yb, grid_res=GRID,
                                               xbin=XENV)
        return pl_regions.density_region(lx[pt], ly[pt], xb, yb, grid_res=GRID,
                                         bw=BW, kpts=KPTS, min_area_frac=MIN_AREA)

    xf, yf, relic_ok, relic_dens, relic_thr = region(relic_pt)
    _, _, full_ok, full_dens, full_thr = region(full_pt)
    comp_reg = None
    if SPLIT is not None:
        comp_reg = [region(full_pt & (SPLIT["labels"] == k)) for k in range(len(SPLIT["names"]))]
        full_ok = np.logical_or.reduce([r[2] for r in comp_reg])
        full_dens, full_thr = comp_reg[0][3], comp_reg[0][4]   # main part: edge + label
    allowed_c = {n: region(fin & (full_c_LL[n] >= full_c_max[n] + CL_2SIG))[2]
                 for n in active}

    fig, ax = plt.subplots(figsize=(7.8, 5.8))
    # Survivor region (relic + ALL constraints, allowed at 2 sigma), flat fill.
    # No interior colour gradient: the profiled Lambda = ln(L/Lmax) is
    # uninformative here (the relic likelihood is a thin Gaussian, so Lambda
    # collapses from 0 to CL_2SIG across a hairline and any ramp reads as an
    # abrupt band).  The survivor edge is the envelope drawn below.
    if comp_reg is None:
        pl_regions.fill_region(ax, full_ok, xf, yf, SURVIVOR_BLUE, zorder=1)
    else:
        for r, col in zip(comp_reg, SPLIT["colors"]):
            pl_regions.fill_region(ax, r[2], xf, yf, col, zorder=1)

    # Translucent per-constraint excluded overlays (relic-allowed but killed by C
    # alone).  Each is one colour at fixed alpha; overlapping regions stack alpha
    # and darken, so joint exclusions read directly.
    killed = relic_ok & ~full_ok                # relic-allowed but killed by the SUM
    kc = {n: ~allowed_c[n] for n in active}
    attributed = np.zeros_like(killed)
    for n in active:
        attributed |= killed & kc[n]
    # Joint-only kills: excluded by the COMBINED likelihood but by no single
    # constraint alone -> grey (else they'd render white and read as "not
    # relic-viable", which they are not).  Drawn first, under the colours.
    joint_only = killed & ~attributed
    # In a joint-only cell every constraint alone lets SOME point through, but
    # the sum kills them all.  Each killed point has a violated set E(p) (the
    # bounds that exclude it on their own).  Draw, in each cell, every bound that
    # belongs to an INCLUSION-MINIMAL E(p) present there -- e.g. {SN1987A} from
    # points only SN1987A excludes and {K -> pi phi, neutron star} from points that
    # pass SN1987A.  Since every single bound lets some point through, a joint cell
    # always gets at least two bounds, and a bound in a superset (a point that
    # SN1987A AND something else exclude) adds nothing new.  (Any bound removing
    # some point is too weak a test: projected over the hidden axes, nearly every
    # bound removes something in every cell.)  Each violated set's cells are a
    # density region like everything else, taken in order of set size.
    joint_bits = np.zeros(joint_only.shape, dtype=np.uint16)
    if joint_only.any():
        killed_pts = relic_pt & (n_viol > 0)
        codes = np.unique(viol_bits[killed_pts])
        present = {}                                    # violated set -> its cells
        for b in sorted(codes.tolist(), key=lambda c: bin(c).count("1")):
            cells = joint_only & region(relic_pt & (viol_bits == b))[2]
            if not cells.any():
                continue
            blocked = np.zeros_like(cells)
            for b2, cells2 in present.items():
                if (b2 & b) == b2:                      # strict subset already present
                    blocked |= cells2
            joint_bits[cells & ~blocked] |= b
            present[b] = cells
    joint_c = {n: joint_only & (((joint_bits >> i) & 1) == 1) for i, n in enumerate(active)}
    n_joint = np.zeros(joint_only.shape, dtype=np.uint8)
    for n in active:
        n_joint += joint_c[n]
    tallies = []
    covered = np.zeros_like(joint_only)
    for n in active:
        killed_c = killed & kc[n]
        joint_c_n = joint_c[n]
        covered |= joint_c_n
        tallies.append(f"{CONSTRAINTS[n][1]} {int(killed_c.sum())}+{int(joint_c_n.sum())}j")
        pl_regions.fill_region(ax, killed_c | joint_c_n, xf, yf, CONSTRAINTS[n][0],
                               alpha=alpha, zorder=2.5, hatch=CONSTRAINTS[n][2],
                               hatch_color=hatch_color(CONSTRAINTS[n][0]))
    # Dotted outline of the joint-only region (lightly smoothed so the raster
    # mask's staircase does not show).
    if joint_only.any():
        cs = ax.contour(xf, yf, gaussian_filter(joint_only.astype(float), 1.5), levels=[0.5],
                        colors=[JOINT_COLOR], linestyles=[":"], linewidths=[1.4], zorder=4.5)
        cs.set_path_effects([patheffects.withStroke(linewidth=3.2, foreground="white")])

    # Survivor envelope (solid), the smooth threshold contour of its density field.
    edges = ([(full_dens, full_thr, SURVIVOR_EDGE)] if comp_reg is None else
             [(r[3], r[4], hatch_color(col)) for r, col in zip(comp_reg, SPLIT["colors"])])
    for dens_k, thr_k, edge_k in edges:
        for seg in pl_regions.field_boundaries(dens_k, xf, yf, thr_k):
            ax.plot(seg[:, 0], seg[:, 1], color=edge_k, lw=1.4, zorder=4)

    ax.set(xlabel=xlabel, ylabel=ylabel, xlim=xb, ylim=yb)
    if AX[x_key][0] == "log":
        style_log_axis(ax.xaxis, *xb)
    if AX[y_key][0] == "log":
        style_log_axis(ax.yaxis, *yb)
    title = ax.text(0.03, 0.97, overlay_title, transform=ax.transAxes,
                    va="top", ha="left", fontsize=26 * _FS, zorder=6,
                    bbox=dict(boxstyle="round", facecolor="white",
                              alpha=0.7, edgecolor="none"))
    # Legend lives in a single standalone file (write_legend), not on each plot.
    fig.tight_layout()
    if full_ok.any():
        segs = pl_regions.field_boundaries(full_dens, xf, yf, full_thr)
        # Inside the blue in white; if the region is too thin or ragged for that,
        # just outside it in navy, over empty space only (never on a fill).
        if (place_survivor_label(ax, full_ok, xf, yf, segs, SURVIVOR_LABEL,
                                 16 * _FS, avoid=(title,)) is None
                and place_survivor_label(ax, ~relic_ok, xf, yf, segs, SURVIVOR_LABEL,
                                         16 * _FS, avoid=(title,),
                                         color=SURVIVOR_EDGE) is None):
            print(f"  [{x_key}:{y_key}] survivor label: no straight edge with room -- omitted")
    p = os.path.join(out_dir, f"all_overlay__{x_key}__{y_key}.pdf")
    fig.savefig(p, bbox_inches="tight")
    if save_png:
        fig.savefig(p.replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [{x_key}:{y_key}] wrote {os.path.basename(p)}"
          f"{' (+ .png)' if save_png else ''};  killed: {', '.join(tallies)}"
          f", joint-only {int(joint_only.sum())} (cells by #bounds drawn {np.bincount(n_joint[joint_only], minlength=4).tolist()})"
          f" / {int(relic_ok.sum())} relic-allowed"
          f"  (fine-grid cells)")
    return bool(joint_only.any())


any_joint = False
for pl in planes:
    any_joint |= make_plane(*pl.split(":"))
write_legend(active, any_joint, out_dir)
