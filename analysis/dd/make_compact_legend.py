#!/usr/bin/env python3
"""
Compact legend for ONE (m_chi, sigma_chiN) panel, exactly one panel wide.

Entries are read off the cached figure itself, not off legend.json: that file
lists every constraint active anywhere in the run, whereas this legend carries
only what this panel draws (e.g. INTEGRAL is active in the MFV S(x)S run but
kills nothing on this plane).  Constraint styling -- colours, hatches, alpha,
the darkened hatch shade -- is parsed from plot_all_constraints.py with ``ast``
so it cannot drift from the figures; the overlay entries come from the same
decorate() that draws them.

Width: the page is made exactly as wide as the panel PDF (measured, in pt), and
the font is the panel's own legend size.  Placed at the same \\includegraphics
width as the panel, the text therefore comes out the same size as the panel's
tick labels.  Height is minimised by picking the most columns that fit, tight
handle/column spacing, no frame, and short labels (``--long`` keeps the
production labels).

``--ddproj`` builds it for a plot_sigma_mx_ddproj.py panel instead (existing SI
limit + projections, or the SD set for A(x)A), from that driver's cache, and
writes ``legend_compact.pdf`` beside it; ``--nolegend-panel`` also re-renders
that panel without its inline legend (``..._nolegend.pdf``), the pair to use
when the legend is placed separately.

    ./make_compact_legend.py [RUN_TAG] [--long] [--ddproj [--nolegend-panel]]
"""

import ast
import os
import pickle
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba, to_hex
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import Bbox

HERE = os.path.dirname(os.path.abspath(__file__))
PLOTTING = os.path.abspath(os.path.join(HERE, ".."))
sys.path[:0] = [HERE, PLOTTING]
import plot_style
import plot_sigma_mx_dd as drv

PAC = os.path.join(PLOTTING, "plot_all_constraints.py")
TAG_DEFAULT = "heavyScalar_relic_mfv_merged"


def pac_constants():
    """Module-level literals of plot_all_constraints.py, without running it."""
    want = {"CONSTRAINTS", "SURVIVOR_BLUE", "SURVIVOR_EDGE", "JOINT_COLOR"}
    out = {}
    # Bullet's style is pl_regions.DASHES (hand-drawn dashes), not a literal.
    import pl_regions
    src = open(PAC).read().replace("pl_regions.DASHES", repr(pl_regions.DASHES))
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 \
                and getattr(node.targets[0], "id", None) in want:
            out[node.targets[0].id] = ast.literal_eval(node.value)
    return out


def hatch_color(color, shade=0.55):          # as plot_all_constraints.hatch_color
    r, g, b, _ = to_rgba(color)
    return (r * shade, g * shade, b * shade, 1.0)


# Short labels.  Only ever lossless on a panel: "(CMB)" is dropped from Delta
# N_eff only when the BBN term is not also drawn (see drawn_entries).
SHORT = {
    "survivor":    "Survives all",
    "neff_cmb":    r"$\Delta N_\mathrm{eff}$ (CMB)",
    "neff_bbn":    r"$\Delta N_\mathrm{eff}$ (BBN)",
    "sn1987a":     "SN1987A",
    "neutronstar": "Neutron star",
    "stellar":     "Stellar",
    "fifthforce":  "5th force",
    "bullet":      "Bullet cluster",
    "meson":       "Meson decays",
    "cmb":         r"CMB ($p_\mathrm{ann}$)",
    "integral":    "INTEGRAL",
    "joint":       "Joint",
    "dd":          "Direct detection",
    "XQC":         "XQC",
    "MWSatellites": "MW satellites",
}


def drawn_entries(fig, K):
    """Constraint keys + survivor/joint flags present on the cached panel.

    Each constraint has a unique (colour, hatch) pair -- that is how the figure
    tells them apart -- so the filled contour sets map back to keys exactly.
    """
    import pl_regions
    # A DASHES fill is drawn as a plain (unhatched) contourf plus clipped lines.
    by_style = {(to_hex(c), None if h == pl_regions.DASHES else (h or None)): k
                for k, (c, _l, h) in K["CONSTRAINTS"].items()}
    keys, survivor, joint = [], False, False
    for coll in fig.axes[0].collections:
        fc, ec = coll.get_facecolor(), coll.get_edgecolor()
        if len(fc):
            hatch = (getattr(coll, "hatches", None) or [None])[0]
            hx = to_hex(fc[0])
            if hx == to_hex(K["SURVIVOR_BLUE"]) and hatch is None:
                survivor = True
            elif (hx, hatch) in by_style:
                keys.append(by_style[(hx, hatch)])
        elif len(ec) and to_hex(ec[0]) == to_hex(K["JOINT_COLOR"]):
            joint = True
    return [k for k in K["CONSTRAINTS"] if k in keys], survivor, joint


def overlay_entries(tag):
    """(name, handle) for each overlay curve decorate() actually draws."""
    fig = pickle.load(open(drv.cache_path_for(tag), "rb"))
    old = drv.SHOW_LEGEND
    drv.SHOW_LEGEND = True
    try:
        drv.decorate(fig.axes[0])
    finally:
        drv.SHOW_LEGEND = old
    leg = fig.axes[0].get_legend()
    labels = [t.get_text() for t in leg.get_texts()] if leg else []
    plt.close(fig)
    import si_limits as S
    out = []
    for lab in labels:
        if lab.startswith("SI direct detection"):
            out.append(("dd", (Patch(fc="none", ec=S.DD_COLOR, hatch="xx", lw=0),
                               Line2D([], [], color=S.DD_COLOR, lw=2.2))))
        elif lab == S.CATALOG["XQC"][2]:
            out.append(("XQC", (Patch(fc="none", ec=S.ASTRO_COLOR, hatch="..", lw=0),
                                Line2D([], [], color=S.ASTRO_COLOR, lw=1.6, ls="--"))))
        elif lab == S.CATALOG["MWSatellites"][2]:
            out.append(("MWSatellites", (Patch(fc="none", ec=S.ASTRO_COLOR, hatch="\\\\", lw=0),
                                         Line2D([], [], color=S.ASTRO_COLOR, lw=1.8, ls="-."))))
    return out


def overlay_entries_ddproj(tag):
    """(label, handle) for each projection-overlay item drawn on the panel."""
    import plot_sigma_mx_ddproj as dp
    fig = pickle.load(open(os.path.join(dp.CACHE_DIR, f"{tag}__mx__sigma_chiN.fig.pkl"), "rb"))
    dp.set_mode(tag)
    old = dp.SHOW_LEGEND
    dp.SHOW_LEGEND = True
    try:
        dp.decorate(fig.axes[0])
    finally:
        dp.SHOW_LEGEND = old
    leg = fig.axes[0].get_legend()
    out = []
    if leg is not None:
        for h, t in zip(leg.legend_handles, leg.get_texts()):
            out.append((t.get_text(), h))
    plt.close(fig)
    return out


def nolegend_panel(tag):
    """Re-render the ddproj panel (and zoom) from cache with no inline legend."""
    import plot_sigma_mx_ddproj as dp
    fig = pickle.load(open(os.path.join(dp.CACHE_DIR, f"{tag}__mx__sigma_chiN.fig.pkl"), "rb"))
    dp.set_mode(tag)
    old = dp.SHOW_LEGEND
    dp.SHOW_LEGEND = False
    try:
        dp.decorate(fig.axes[0])
    finally:
        dp.SHOW_LEGEND = old
    out = os.path.join(dp.OUT_ROOT, f"{tag}_ddproj", f"{drv.TARGET}_nolegend.pdf")
    fig.savefig(out, bbox_inches="tight")
    print(f"  wrote {os.path.relpath(out, HERE)}")
    plt.close(fig)


def out_root(sub):
    """The ddproj driver writes beside the production plots, the prototype here."""
    if sub == "_ddproj":
        import plot_sigma_mx_ddproj as dp
        return dp.OUT_ROOT
    return drv.OUT_ROOT


def panel_width_in(tag, sub="_dd"):
    pdf = os.path.join(out_root(sub), f"{tag}{sub}", f"{drv.TARGET}.pdf")
    info = subprocess.run(["pdfinfo", pdf], capture_output=True, text=True, check=True).stdout
    w_pt = float(next(l for l in info.splitlines() if l.startswith("Page size")).split()[2])
    return w_pt / 72.0


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    long_labels = "--long" in sys.argv[1:]
    ddproj = "--ddproj" in sys.argv[1:]
    tag = args[0] if args else TAG_DEFAULT
    sub = "_ddproj" if ddproj else "_dd"
    if ddproj and "--nolegend-panel" in sys.argv[1:]:
        nolegend_panel(tag)

    # Same text setup as the panels: CM fonts and plot_all_constraints' sizes.
    plot_style.use_latex_fonts()
    fs = 19 * float(os.environ.get("LOOP_PLOT_FONTSCALE", "1.0"))
    alpha = float(os.environ.get("LOOP_PLOT_ALPHA", "0.7"))

    K = pac_constants()
    if ddproj:
        import plot_sigma_mx_ddproj as dp
        fig0 = pickle.load(open(os.path.join(dp.CACHE_DIR, f"{tag}__mx__sigma_chiN.fig.pkl"), "rb"))
    else:
        fig0 = pickle.load(open(drv.cache_path_for(tag), "rb"))
    keys, survivor, joint = drawn_entries(fig0, K)
    plt.close(fig0)

    def label(k, prod):
        if long_labels:
            return prod
        if k == "neff_cmb" and "neff_bbn" not in keys:
            return r"$\Delta N_\mathrm{eff}$"
        return SHORT.get(k, prod)

    handles, labels = [], []
    if survivor:
        handles.append(Patch(fc=K["SURVIVOR_BLUE"], lw=0))
        labels.append(label("survivor", r"Survives all ($2\sigma$)"))
    import pl_regions
    dash_map = {}
    for k in keys:
        c, prod, h = K["CONSTRAINTS"][k]
        if h == pl_regions.DASHES:
            p = Patch(fc=to_rgba(c, alpha), lw=0)
            dash_map[p] = pl_regions.HandlerDashes(hatch_color(c))
        else:
            p = Patch(fc=to_rgba(c, alpha), hatch=h, ec=hatch_color(c) if h else "none", lw=0)
        handles.append(p)
        labels.append(label(k, prod))
    if joint:
        handles.append(Line2D([], [], color=K["JOINT_COLOR"], lw=1.4, ls=":"))
        labels.append(label("joint", "Several bounds jointly"))
    if ddproj:
        for lab, h in overlay_entries_ddproj(tag):
            handles.append(h)
            labels.append(lab)
    else:
        for k, h in overlay_entries(tag):
            handles.append(h)
            labels.append(label(k, {"dd": "SI direct detection (best limit)",
                                    "XQC": "XQC (overburden)",
                                    "MWSatellites": "MW satellites"}[k]))

    W = panel_width_in(tag, sub)
    kw = dict(frameon=False, fontsize=fs, handlelength=1.3, handleheight=0.75,
              handletextpad=0.4, columnspacing=0.9, labelspacing=0.18,
              borderpad=0.0, borderaxespad=0.0,
              handler_map={tuple: HandlerTuple(ndivide=1, pad=0), **dash_map})  # overlay, as drawn

    def build(ncol, expand):
        fig = plt.figure(figsize=(W, 4))
        leg = fig.legend(handles, labels, ncol=ncol, loc="upper left",
                         bbox_to_anchor=(0, 1, 1, 0) if expand else (0, 1),
                         mode="expand" if expand else None, **kw)
        fig.canvas.draw()
        bb = leg.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        return fig, bb

    # Most columns whose natural width fits the panel -> fewest rows.
    ncol = 1
    for n in range(len(handles), 0, -1):
        f, bb = build(n, expand=False)
        plt.close(f)
        if bb.width <= W:
            ncol = n
            break
    fig, bb = build(ncol, expand=True)
    pad = 0.02                                            # in: keeps descenders
    box = Bbox([[0, bb.y0 - pad], [W, bb.y1 + pad]])
    out = os.path.join(out_root(sub), f"{tag}{sub}",
                       "legend_compact_long.pdf" if long_labels else "legend_compact.pdf")
    fig.savefig(out, bbox_inches=box)
    plt.close(fig)
    rows = -(-len(handles) // ncol)
    print(f"  wrote {os.path.relpath(out, HERE)}: {len(handles)} entries, "
          f"{ncol} cols x {rows} rows, {W * 72:.1f} x {box.height * 72:.1f} pt "
          f"(panel is {W * 72:.1f} pt wide)")


if __name__ == "__main__":
    main()
