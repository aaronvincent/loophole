#!/usr/bin/env python3
"""
The combined-constraint (m_chi, sigma_chiN) panel for the HEAVY-MEDIATOR S(x)S
scans (m_phi = 10 GeV), with

  * the existing spin-independent direct-detection limit, drawn as the closed
    region it is (floor = best limit, ceiling = overburden), hatched, and
  * the projected sensitivities of HeRALD, QUEST-DMC and SuperCDMS SNOLAB,
    drawn as outlines only (a projection excludes nothing yet).

All four come from ``data/projections/`` (read-only; see projections.py).

An axialVector run (A(x)A) gets the SPIN-DEPENDENT set instead: the existing
SD-proton limit, hatched above its curve, and the QUEST-DMC neutron-only
projection.  The panel's sigma_chiN is the chi-PROTON cross section, so the
neutron curve is not comparable to it point by point: sd_reach() computes
sigma_chin = sigma_chiN (a_n/a_p)^2 for every survivor and the panel outlines the
survivors whose sigma_chin lies above the QUEST-DMC curve.  They
assume a contact interaction, which is why only the m_phi = 10 GeV runs are
drawn: at these m_chi the momentum transfer is <= 0.1 GeV << m_phi.

Machinery is shared with plot_sigma_mx_dd.py: plot_all_constraints
runs under ``runpy`` behind a wrapped ``Figure.savefig``, which draws this
overlay onto the finished axes and redirects every output into
``analysis/plots/<tag>_ddproj/``.  The main figures are not touched.  The sigma window
is widened with ``LOOP_SIGMA_YMAX`` so the region above the relic band -- where
the existing limit actually removes samples -- is inside the frame.

Usage
-----
    python analysis/dd/plot_sigma_mx_ddproj.py RUN_OR_SIDECAR ...    # scan file(s)
    python analysis/dd/plot_sigma_mx_ddproj.py --from-cache [TAG]    # restyle only (~2 s)

Env knobs
---------
    DDP_SIGMA_YMAX   top of the sigma window, log10 cm^2 (default -26)
    DDP_LEGEND       1 (default) inline legend on the full panel (not the zoom);
                     DDP_LEGEND_LOC / DDP_LEGEND_ANCHOR (upper left, "0.005,0.79")
    DDP_ZOOM         1 (default) also write *_zoom.pdf; DDP_ZOOM_BOX
                     "mlo,mhi,slo,shi" in log10 (default "-4,-0.6,-47,-30")
    DDP_PROJ         comma list of projections to draw (default all three)
    plus every LOOP_* knob plot_all_constraints understands.
"""

import os
import pickle
import runpy
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plot_sigma_mx_dd as base          # noqa: E402  (also sets sys.path for plotting/)
import projections as P                  # noqa: E402

REPO = base.REPO
PLOTTING = base.PLOTTING
DEFAULT_RUNS = [
    os.path.join(REPO, "runs", t, "samples", f"{t}_postprocessed.hdf5")
    for t in ("heavyScalar_relic_mphi10", "heavyScalar_relic_mphi10_mfv")
]

SIGMA_YMAX = os.environ.get("DDP_SIGMA_YMAX", "-26")
SHOW_LEGEND = os.environ.get("DDP_LEGEND", "1") not in ("0", "", "false", "False")
# Upper-left, under the S(x)S / m_phi title: the one corner the relic band, the
# existing region and the projection floors all leave empty.
LEGEND_LOC = os.environ.get("DDP_LEGEND_LOC", "upper left")
LEGEND_ANCHOR = tuple(float(v) for v in
                      os.environ.get("DDP_LEGEND_ANCHOR", "0.005,0.79").split(","))
SHOW_ZOOM = os.environ.get("DDP_ZOOM", "1") not in ("0", "", "false", "False")
ZOOM_BOX = tuple(float(v) for v in
                 os.environ.get("DDP_ZOOM_BOX", "-4,-0.6,-47,-30").split(","))
PROJ = [p for p in os.environ.get("DDP_PROJ", ",".join(P.PROJECTIONS)).split(",") if p]

OUT_ROOT = os.path.join(PLOTTING, "plots")      # beside <tag>_allconstraints/
# DDP_OUT_SUFFIX: write to <tag>_ddproj<suffix>/ instead (variants beside the main output).
OUT_SUFFIX = os.environ.get("DDP_OUT_SUFFIX", "")
# Both HeRALD curves drawn: say which is which.
if "HeRALD_c4" in PROJ:
    P.CATALOG["HeRALD"] = P.CATALOG["HeRALD"][:2] + (r"HeRALD 1 kg$\cdot$day, 40 eV (proj.)",)

# Spin-dependent mode (axialVector runs), set per run by set_mode().
SD = False
REACH = None                 # (log10 mx, log10 sigma_chiN) of survivors QUEST-DMC SD-n reaches
DQ_P = np.array([0.847, -0.407, -0.035])   # Delta_q^p (u, d, s), FLAG 2024, as observables.jl
DQ_N = DQ_P[[1, 0, 2]]                     # isospin: Delta_u^n = Delta_d^p, ...
REACH_COLOR = "#000000"
SD_STYLE = {"QUESTDMC_SDn": dict(color="#555555", ls=(0, (1.5, 2)), lw=2.2)}
CACHE_DIR = os.path.join(HERE, "cache_proj")

# Existing limit: black, hatched (it excludes).  Projections: outlines only, one
# dash pattern each, achromatic as in plot_sigma_mx_dd.py -- the
# coloured fills in these panels are terms in the likelihood, and none of these is.
EXIST_COLOR = "#000000"
PROJ_STYLE = {
    "HeRALD":    dict(color="#000000", ls=(0, (6, 3)),          lw=1.8),
    "HeRALD_c4": dict(color="#000000", ls=(0, (3, 1.5, 1, 1.5)), lw=1.4),
    "QUESTDMC":  dict(color="#555555", ls=(0, (1.5, 2)),        lw=2.2),
    "SuperCDMS": dict(color="#000000", ls=(0, (7, 2, 1.5, 2)), lw=1.8),
}


def _visible(x, y, xlim, ylim):
    return not (np.nanmin(x) > xlim[1] or np.nanmax(x) < xlim[0]
                or np.nanmin(y) > ylim[1] or np.nanmax(y) < ylim[0])


def decorate(ax):
    """Draw the existing-limit region and the projection outlines."""
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    base.PANEL["xlim"], base.PANEL["ylim"] = xlim, ylim
    if SD:
        return decorate_sd(ax, xlim, ylim)
    handles = []

    px, py = P.polygon_log("existing")
    if _visible(px, py, xlim, ylim):
        ax.fill(px, py, facecolor="none", edgecolor=EXIST_COLOR, hatch="xx",
                linewidth=0.0, alpha=0.6, zorder=5.0)
        lm, lo, _ = P.region_log("existing")
        # Solid line on the floor (the bound); the ceiling and sides are the
        # overburden / mass reach, drawn thin.
        ax.plot(px, py, color=EXIST_COLOR, lw=0.8, zorder=5.1)
        ax.plot(lm, lo, color=EXIST_COLOR, lw=2.2, zorder=5.2)
        handles.append(Patch(facecolor="none", edgecolor=EXIST_COLOR, hatch="xx",
                             lw=1.5, label=P.CATALOG["existing"][2]))

    for n in PROJ:
        st = PROJ_STYLE[n]
        if P.CATALOG[n][1] == "region":
            x, y = P.polygon_log(n)
        else:
            x, y, _ = P.region_log(n)
        if not _visible(x, y, xlim, ylim):
            continue
        ax.plot(x, y, zorder=5.3, solid_capstyle="round", **st)
        handles.append(Line2D([], [], label=P.CATALOG[n][2], **st))

    if SHOW_LEGEND and handles:
        leg = ax.legend(handles=handles, loc=LEGEND_LOC, bbox_to_anchor=LEGEND_ANCHOR,
                        fontsize=11, framealpha=0.9, facecolor="white",
                        edgecolor="0.6", handlelength=2.6)
        leg.set_zorder(8)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def sd_handles():
    return ([Patch(facecolor="none", edgecolor=EXIST_COLOR, hatch="xx", lw=1.5,
                   label=P.CATALOG["existing_SD"][2])]
            + [Line2D([], [], label=P.CATALOG[n][2], **SD_STYLE[n]) for n in P.SD_PROJECTIONS]
            + [Line2D([], [], color=REACH_COLOR, lw=1.6, ls=(0, (4, 1.5, 1, 1.5)),
                      label="QUEST-DMC SD-$n$ reach:\n" r"allowed points with $\sigma_{\chi n}$ above it")])


def decorate_sd(ax, xlim, ylim):
    """A(x)A: existing SD-p limit (hatched above the curve), the QUEST-DMC SD-n
    projection, and the outline of the survivors that projection reaches."""
    lm, lo, _ = P.region_log("existing_SD")
    if _visible(lm, lo, xlim, ylim):
        top = np.full(lm.size, ylim[1])
        ax.fill_between(lm, lo, top, facecolor="none", edgecolor=EXIST_COLOR, hatch="xx",
                        linewidth=0.0, alpha=0.6, zorder=5.0)
        ax.plot(lm, lo, color=EXIST_COLOR, lw=2.2, zorder=5.2)
    for n in P.SD_PROJECTIONS:
        x, y, _ = P.region_log(n)
        if _visible(x, y, xlim, ylim):
            ax.plot(x, y, zorder=5.3, solid_capstyle="round", **SD_STYLE[n])
    if REACH is not None and len(REACH[0]) >= 5:
        import pl_regions
        xf, yf, ok, dens, thr = pl_regions.density_region(
            REACH[0], REACH[1], xlim, ylim, grid_res=600, bw=0.00625 * 600, kpts=2.0)
        for seg in pl_regions.field_boundaries(dens, xf, yf, thr):
            ax.plot(seg[:, 0], seg[:, 1], color=REACH_COLOR, lw=1.6,
                    ls=(0, (4, 1.5, 1, 1.5)), zorder=5.4)
    if SHOW_LEGEND:
        leg = ax.legend(handles=sd_handles(), loc=LEGEND_LOC, bbox_to_anchor=LEGEND_ANCHOR,
                        fontsize=10, framealpha=0.9, facecolor="white",
                        edgecolor="0.6", handlelength=2.6)
        leg.set_zorder(8)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def sd_reach(in_arg, group="data"):
    """Survivors of an A(x)A run whose chi-NEUTRON cross section lies above the
    QUEST-DMC SD-n curve, as (log10 mx, log10 sigma_chiN, log10 sigma_chin).
    sigma_chin = sigma_chiN * (a_n / a_p)^2: the kinematic factor is common to
    both nucleons, only the spin charge differs.  Survivors use the plotter's
    own cut (LogLike >= max - 3.090)."""
    import h5py
    with h5py.File(in_arg, "r") as h:
        g = h[group]
        v = g["LogLike_isvalid"][:].astype(bool)
        LL = g["LogLike"][:]
        mx, sp = g["input::mx"][:], g["output::sigma_chiN"][:]
        q = np.column_stack([g["input::gu"][:], g["input::gd"][:], g["input::gs"][:]])
    fin = v & np.isfinite(LL) & np.isfinite(sp) & (sp > 0)
    surv = fin & (LL >= LL[fin].max() - 3.090)
    ap, an = q[surv] @ DQ_P, q[surv] @ DQ_N
    with np.errstate(divide="ignore", invalid="ignore"):
        sn = sp[surv] * (an / ap) ** 2
    lx, lyp, lyn = np.log10(mx[surv]), np.log10(sp[surv]), np.log10(sn)
    f = P.floor_at("QUESTDMC_SDn", lx)
    k = np.isfinite(f) & np.isfinite(lyn) & (lyn >= f)
    return lx[k], lyp[k], lyn[k], dict(n_surv=int(surv.sum()), ratio=np.log10((an / ap) ** 2))


def set_mode(tag, in_arg=None):
    """SD overlay for axialVector runs; the reach is cached beside the figure."""
    global SD, REACH
    SD = tag.startswith("axialVector")
    REACH = None
    if not SD:
        return
    cp = os.path.join(CACHE_DIR, f"{tag}__sd_reach.npz")
    if in_arg is not None:
        lx, ly, lyn, info = sd_reach(in_arg)
        os.makedirs(CACHE_DIR, exist_ok=True)
        np.savez(cp, lx=lx, ly=ly, lyn=lyn, n_surv=info["n_surv"])
        r = info["ratio"][np.isfinite(info["ratio"])]
        print(f"  SD: log10(sigma_n/sigma_p) over survivors: median {np.median(r):.2f}, "
              f"5-95% [{np.percentile(r, 5):.2f}, {np.percentile(r, 95):.2f}]")
    elif not os.path.exists(cp):
        raise SystemExit(f"no SD reach cache at {cp} -- run once without --from-cache")
    z = np.load(cp)
    REACH = (z["lx"], z["ly"])


def write_legend(out_dir):
    import matplotlib.pyplot as plt
    if SD:
        handles = sd_handles()
    else:
        handles = [Patch(facecolor="none", edgecolor=EXIST_COLOR, hatch="xx",
                         label=P.CATALOG["existing"][2])]
        handles += [Line2D([], [], label=P.CATALOG[n][2], **PROJ_STYLE[n]) for n in PROJ]
    fig = plt.figure(figsize=(10, 1.2))
    leg = fig.legend(handles=handles, loc="center", ncol=len(handles), frameon=True,
                     framealpha=1.0, handlelength=2.6)
    fig.canvas.draw()
    bbox = leg.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    p = os.path.join(out_dir, "legend_ddproj.pdf")
    fig.savefig(p, bbox_inches=bbox.expanded(1.05, 1.15))
    plt.close(fig)
    print(f"  wrote {os.path.basename(p)}")


def report(g, out_dir):
    """How many samples each region contains -- counted on SAMPLES, not on the
    drawn (smoothed) regions.  Writes ``ddproj_overlap.txt``."""
    if SD:
        return report_sd(g, out_dir)
    lx, ly = g["log_axis"]["mx"], g["log_axis"]["sigma_chiN"]
    fin = g["valid"] & np.isfinite(lx) & np.isfinite(ly)
    surv = fin & (g["total_LL"] >= g["full_max"] + g["CL_2SIG"])
    relic = fin & (g["relic_LL"] >= g["relic_max"] + g["CL_2SIG"])
    ex = P.inside("existing", lx, ly)
    lines = [f"run: {g['run_tag']}",
             f"m_chi range of survivors: {10**lx[surv].min():.3g} -- {10**lx[surv].max():.3g} GeV",
             f"relic-allowed (2 sigma): {int(relic.sum()):,}   survives all scan constraints: {int(surv.sum()):,}",
             "",
             f"inside the EXISTING SI region: relic-allowed {int((relic & ex).sum()):,}, "
             f"survivors {int((surv & ex).sum()):,}",
             "",
             "projected reach (survivors not already excluded by the existing limit):"]
    base_s = surv & ~ex
    for n in PROJ:
        inn = P.inside(n, lx, ly)
        k = base_s & inn
        lines.append(f"  {P.CATALOG[n][2]:28s} {int(k.sum()):>10,}  "
                     f"({100 * k.sum() / max(base_s.sum(), 1):.1f}%)"
                     + (f"   m_chi {10**lx[k].min():.3g} -- {10**lx[k].max():.3g} GeV"
                        if k.any() else ""))
    anyp = np.zeros_like(base_s)
    for n in PROJ:
        anyp |= P.inside(n, lx, ly)
    lines.append(f"  {'any projection':28s} {int((base_s & anyp).sum()):>10,}  "
                 f"({100 * (base_s & anyp).sum() / max(base_s.sum(), 1):.1f}%)")
    ylo, yhi = base.PANEL.get("ylim", (np.nan, np.nan))
    off = int((surv & (ly > yhi)).sum())
    if off:
        lines += ["", f"NOTE: {off:,} survivors sit above the panel top (10^{yhi:.0f})."]
    p = os.path.join(out_dir, "ddproj_overlap.txt")
    with open(p, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join("  " + l for l in lines if l))


def report_sd(g, out_dir):
    """A(x)A counts: survivors inside the existing SD-p limit, and those whose
    sigma_chin the QUEST-DMC SD-n projection reaches.  Writes ddproj_overlap.txt."""
    lx, ly = g["log_axis"]["mx"], g["log_axis"]["sigma_chiN"]
    fin = g["valid"] & np.isfinite(lx) & np.isfinite(ly)
    surv = fin & (g["total_LL"] >= g["full_max"] + g["CL_2SIG"])
    f = P.floor_at("existing_SD", lx)
    with np.errstate(invalid="ignore"):
        ex = np.isfinite(f) & (ly >= f)
    lines = [f"run: {g['run_tag']}",
             f"m_chi range of survivors: {10**lx[surv].min():.3g} -- {10**lx[surv].max():.3g} GeV",
             f"survives all scan constraints: {int(surv.sum()):,}",
             f"existing SD-p limit covers m_chi >= {10**P.region_log('existing_SD')[0].min():.3g} GeV;"
             f" survivors above it: {int((surv & ex).sum()):,}",
             f"QUEST-DMC SD-n (proj.) reaches {len(REACH[0]):,} survivors (sigma_chin above the curve)"
             + (f", m_chi {10**REACH[0].min():.3g} -- {10**REACH[0].max():.3g} GeV" if len(REACH[0]) else "")]
    p = os.path.join(out_dir, "ddproj_overlap.txt")
    with open(p, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join("  " + l for l in lines))


_base_save_zoom = base.save_zoom


def save_zoom(fig, path, saver=None):
    """The base zoom, with the inline legend hidden: at zoom scale it
    lands on the region the zoom is for (legend_ddproj.pdf carries it)."""
    leg = fig.axes[0].get_legend()
    if leg is not None:
        leg.set_visible(False)
    try:
        _base_save_zoom(fig, path, saver=saver)
    finally:
        if leg is not None:
            leg.set_visible(True)


def configure_base():
    """Point the plot_sigma_mx_dd wrapper at this overlay and these output dirs."""
    base.decorate = decorate
    base.SHOW_ZOOM = SHOW_ZOOM
    base.ZOOM_BOX = ZOOM_BOX
    base.LEGEND_LOC = LEGEND_LOC
    base.save_zoom = save_zoom


def run_one(in_arg, group="data"):
    tag = base.run_tag_of(in_arg)
    out_dir = os.path.join(OUT_ROOT, f"{tag}_ddproj{OUT_SUFFIX}")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(CACHE_DIR, exist_ok=True)
    cache_path = os.path.join(CACHE_DIR, f"{tag}__mx__sigma_chiN.fig.pkl")
    print(f"run:    {in_arg}\noutput: {out_dir}")
    set_mode(tag, in_arg)
    os.environ["LOOP_SIGMA_YMAX"] = SIGMA_YMAX
    import matplotlib.figure
    orig = matplotlib.figure.Figure.savefig
    import builtins
    orig_open = builtins.open
    try:
        base.install_wrapper(out_dir, cache_path,
                             os.path.join(PLOTTING, "plots", f"{tag}_allconstraints"))
        sys.argv = [os.path.join(PLOTTING, "plot_all_constraints.py"), in_arg, group,
                    base.PLANE]
        g = runpy.run_path(sys.argv[0], run_name="__main__")
    finally:
        matplotlib.figure.Figure.savefig = orig      # un-wrap before the next run
        builtins.open = orig_open
    g["run_tag"] = tag
    write_legend(out_dir)
    report(g, out_dir)


def main():
    configure_base()
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if "--from-cache" in sys.argv[1:]:
        tags = args or [base.run_tag_of(r) for r in DEFAULT_RUNS]
        for tag in tags:
            cp = os.path.join(CACHE_DIR, f"{tag}__mx__sigma_chiN.fig.pkl")
            if not os.path.exists(cp):
                raise SystemExit(f"no cache at {cp} -- run once without --from-cache")
            with open(cp, "rb") as f:
                fig = pickle.load(f)
            set_mode(tag)
            out_dir = os.path.join(OUT_ROOT, f"{tag}_ddproj{OUT_SUFFIX}")
            os.makedirs(out_dir, exist_ok=True)
            decorate(fig.axes[0])
            p = os.path.join(out_dir, f"{base.TARGET}.pdf")
            fig.savefig(p, bbox_inches="tight")
            print(f"  wrote {os.path.relpath(p, HERE)} (from cache)")
            if SHOW_ZOOM:
                save_zoom(fig, p[:-4] + "_zoom.pdf")
            write_legend(out_dir)
        return
    for r in (args or DEFAULT_RUNS):
        if not os.path.exists(r):
            raise SystemExit(f"{r} not found")
        run_one(r)


if __name__ == "__main__":
    main()
