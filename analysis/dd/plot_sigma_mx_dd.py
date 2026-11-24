#!/usr/bin/env python3
"""
The combined-constraint (m_chi, sigma_chiN) panel with the published
spin-independent direct-detection limits laid over it.

Default run is the MFV S(x)S merged scan
(``merged/heavyScalar_relic_mfv_merged``).

How it reuses the main plotter
------------------------------------
It does NOT copy ``plot_all_constraints.py``.  It executes that script with
``runpy`` (PLANES restricted to ``mx:sigma_chiN``) behind two wrappers:

  * ``Figure.savefig`` is wrapped -- the direct-detection overlay is drawn onto
    the finished axes just before the file is written, and the output path is
    redirected into ``analysis/dd/plots/<tag>_dd/`` so the main
    figures in ``analysis/plots/<tag>_allconstraints/`` are never touched;
  * the undecorated figure is pickled to ``cache/`` on the way past, so
    ``--from-cache`` can restyle the overlay in a second without re-reading the
    6 GB sidecar.

So every region, colour and label in the panel underneath is exactly what the
main figure shows; the only addition is the overlay.

Usage
-----
    ./plot_sigma_mx_dd.py [RUN_OR_SIDECAR] [GROUP]
    ./plot_sigma_mx_dd.py --from-cache [RUN_TAG]      # restyle only, no data

Env knobs
---------
    DD_MODE          envelope (default) | individual | off
    DD_FILL          band (default) | full | none -- how much of the excluded
                     side of each curve is hatched (see si_limits.draw_upper)
    DD_CEILING       log10 sigma the SI band is hatched up to (default -30, an
                     approximate overburden ceiling); "none" = thin strip
    DD_ASTRO         1 (default) draw the MW-satellite bound, 0 to hide
    DD_BANDS         1 (default) draw the XQC overburden contour, 0 to hide
    DD_LEGEND        1 (default) inline legend for the overlay only
    DD_ZOOM          1 (default) also write a *_zoom.pdf cropped to the corner
                     where the SI bound and the survivor region actually meet
    DD_ZOOM_BOX      "mlo,mhi,slo,shi" in log10 (default "-1.3,0.3,-45,-30")
    DD_LEGEND_LOC    matplotlib legend loc (default "lower right")
    DD_COLOR / DD_ASTRO_COLOR   overlay colours (see si_limits)
    plus every LOOP_* knob plot_all_constraints understands.
"""

import os
import pickle
import runpy
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib.figure import Figure
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
PLOTTING = os.path.abspath(os.path.join(HERE, ".."))
REPO = os.path.abspath(os.path.join(PLOTTING, ".."))
sys.path.insert(0, PLOTTING)   # constraint_io / scan_axes / pl_regions / plot_style
sys.path.insert(0, HERE)

import si_limits

PLANE = "mx:sigma_chiN"
TARGET = "all_overlay__mx__sigma_chiN"
DEFAULT_RUN = os.path.join(
    REPO, "merged", "heavyScalar_relic_mfv_merged", "samples",
    "heavyScalar_relic_mfv_merged_postprocessed.hdf5")

MODE = os.environ.get("DD_MODE", "envelope")
FILL = os.environ.get("DD_FILL", "band")          # band | full | none
# Approximate overburden ceiling for the SI band, log10(cm^2); "none" = strip.
_ceil = os.environ.get("DD_CEILING", "-30")
CEILING = None if _ceil.lower() in ("none", "") else float(_ceil)
SHOW_ASTRO = os.environ.get("DD_ASTRO", "1") not in ("0", "", "false", "False")
SHOW_BANDS = os.environ.get("DD_BANDS", "1") not in ("0", "", "false", "False")
SHOW_LEGEND = os.environ.get("DD_LEGEND", "1") not in ("0", "", "false", "False")
LEGEND_LOC = os.environ.get("DD_LEGEND_LOC", "lower right")
SHOW_ZOOM = os.environ.get("DD_ZOOM", "1") not in ("0", "", "false", "False")
ZOOM_BOX = tuple(float(v) for v in
                 os.environ.get("DD_ZOOM_BOX", "-1.3,0.3,-45,-30").split(","))

OUT_ROOT = os.path.join(HERE, "plots")
CACHE_DIR = os.path.join(HERE, "cache")


# ---------------------------------------------------------------------------
# Caching the undecorated figure
# ---------------------------------------------------------------------------
def make_picklable(fig):
    """Replace unpicklable tick formatters so ``pickle.dump(fig)`` succeeds."""
    from matplotlib.ticker import FuncFormatter
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            fmt = axis.get_major_formatter()
            if isinstance(fmt, FuncFormatter):
                axis.set_major_formatter(
                    FuncFormatter(si_limits.PowerOfTenFormatter()))


# ---------------------------------------------------------------------------
# The overlay
# ---------------------------------------------------------------------------
# What the finished panel actually shows.  make_plane recomputes the sigma
# window locally (percentiles of the relic-allowed samples) rather than using
# the scan_axes BOUNDS entry, so the only reliable source is the axes itself.
PANEL = {}


def decorate(ax):
    """Draw the direct-detection overlay on a finished (m_chi, sigma) panel.

    The panel's axes are log10 values on linear axes (plot_all_constraints
    dresses them with 10^n ticks), which is exactly what si_limits returns, so
    no transform is needed.  The x/y limits are left alone: the overlay must not
    be allowed to rescale the physics panel.
    """
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    PANEL["xlim"], PANEL["ylim"] = xlim, ylim
    handles = []

    if MODE == "envelope":
        lm, ls = si_limits.envelope_log(mlim=xlim)
        handles += si_limits.draw_upper(
            ax, lm, ls, ylim, label="SI direct detection (best limit)", fill=FILL,
            ceiling=CEILING)
    elif MODE == "individual":
        # turbo restricted to its yellow-to-red tail.  Its blue end would put
        # exclusion curves in the colour this figure reserves for "allowed",
        # and its green middle disappears into the teal constraint fills.
        cmap = matplotlib.colormaps["turbo"]
        names = [n for n in si_limits.DEFAULT_ENVELOPE if n in si_limits.available()]
        for i, n in enumerate(names):
            cm, cs = si_limits.curve_log(n)
            keep = (cm >= xlim[0]) & (cm <= xlim[1])
            if keep.sum() < 2:
                continue
            handles += si_limits.draw_upper(
                ax, cm[keep], cs[keep], ylim,
                color=cmap(0.55 + 0.45 * i / max(len(names) - 1, 1)),
                label=si_limits.CATALOG[n][2], fill="none", lw=1.4)
    elif MODE != "off":
        raise SystemExit(f"DD_MODE={MODE!r}: expected envelope | individual | off")

    if SHOW_BANDS:
        for n in si_limits.BANDS:
            if n in si_limits.available():
                handles += si_limits.draw_band(
                    ax, n, color=si_limits.ASTRO_COLOR,
                    label=si_limits.CATALOG[n][2])

    if SHOW_ASTRO:
        for n in si_limits.COSMO_ASTRO:
            if n not in si_limits.available():
                continue
            cm, cs = si_limits.curve_log(n)
            keep = (cm >= xlim[0] - 1) & (cm <= xlim[1] + 1)
            handles += si_limits.draw_upper(
                ax, cm[keep], cs[keep], ylim, color=si_limits.ASTRO_COLOR,
                label=si_limits.CATALOG[n][2], hatch="\\\\", ls_style="-.", lw=1.8,
                fill=FILL)

    if SHOW_LEGEND and handles:
        # The main panels carry no legend (it is a separate legend.pdf).
        # Here the new curves are named in place; the
        # shared legend.pdf is still written beside it, unchanged.
        # DD_MODE=individual is a diagnostic ("which experiment sets the bound
        # here?") and can produce a dozen entries -- shrink and column them so
        # the box does not swallow the panel.
        many = len(handles) > 6
        leg = ax.legend(handles=handles, loc=LEGEND_LOC,
                        fontsize=8 if many else 11, ncol=2 if many else 1,
                        framealpha=0.9, facecolor="white", edgecolor="0.6")
        leg.set_zorder(8)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def save_zoom(fig, path, saver=None):
    """Re-save the decorated panel cropped to ZOOM_BOX.

    The whole panel is ~60 decades tall, so the only place the SI bound can meet
    the survivor region -- the m_chi ~ 0.1 GeV corner -- is a few pixels there.
    Everything is already drawn in data coordinates, so this is just new limits.
    """
    ax = fig.axes[0]
    keep = (ax.get_xlim(), ax.get_ylim())
    leg = ax.get_legend()
    saver = saver or (lambda f, pth, **kw: f.savefig(pth, **kw))
    # The survivor label sits in DATA coordinates out at m_chi ~ 1e-5 and is not
    # clipped to the axes, so bbox_inches="tight" would pad the page out to
    # reach it.  Hide any data-space text that falls outside the zoom box.
    hidden = []
    for t in ax.texts:
        if t.get_transform() is not ax.transData:
            continue
        tx, ty = t.get_position()
        if not (ZOOM_BOX[0] <= tx <= ZOOM_BOX[1] and ZOOM_BOX[2] <= ty <= ZOOM_BOX[3]):
            hidden.append(t)
    try:
        ax.set_xlim(ZOOM_BOX[0], ZOOM_BOX[1])
        ax.set_ylim(ZOOM_BOX[2], ZOOM_BOX[3])
        for t in hidden:
            t.set_visible(False)
        if leg is not None:
            leg.set_loc("upper right")
        saver(fig, path, bbox_inches="tight")
        print(f"  wrote {os.path.basename(path)}")
    finally:
        ax.set_xlim(*keep[0])
        ax.set_ylim(*keep[1])
        for t in hidden:
            t.set_visible(True)
        if leg is not None:
            leg.set_loc(LEGEND_LOC)


def write_dd_legend(out_dir):
    """Standalone legend for the overlay entries, to sit next to legend.pdf."""
    import matplotlib.pyplot as plt
    handles = []
    if MODE == "envelope":
        handles.append(Patch(facecolor="none", edgecolor=si_limits.DD_COLOR,
                             hatch="xx", label="SI direct detection (best limit)"))
    if SHOW_BANDS and "XQC" in si_limits.available():
        handles.append(Patch(facecolor="none", edgecolor=si_limits.ASTRO_COLOR,
                             hatch="..", label=si_limits.CATALOG["XQC"][2]))
    if SHOW_ASTRO and "MWSatellites" in si_limits.available():
        handles.append(Patch(facecolor="none", edgecolor=si_limits.ASTRO_COLOR,
                             hatch="\\\\", label=si_limits.CATALOG["MWSatellites"][2]))
    if not handles:
        return
    fig = plt.figure(figsize=(9, 1.2))
    leg = fig.legend(handles=handles, loc="center", ncol=3, frameon=True, framealpha=1.0)
    fig.canvas.draw()
    bbox = leg.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    p = os.path.join(out_dir, "legend_dd.pdf")
    fig.savefig(p, bbox_inches=bbox.expanded(1.05, 1.15))
    plt.close(fig)
    print(f"  wrote {os.path.basename(p)}")


# ---------------------------------------------------------------------------
# Does the overlay actually bite?
# ---------------------------------------------------------------------------
def report_overlap(g, out_dir):
    """Quantify the overlay against the samples, not against the drawn region.

    ``runpy`` hands back plot_all_constraints' module globals, so the same
    arrays the panel was built from are still in hand.  A drawn region is a
    smoothed density level set (see plotting/README.md -- it is deliberately
    fatter than the samples); asking the SAMPLES whether any of them sit above
    the SI envelope is the statement that can be quoted.

    Writes ``dd_overlap.txt`` beside the figure and returns the summary lines.
    """
    lx, ly = g["log_axis"]["mx"], g["log_axis"]["sigma_chiN"]
    fin = g["valid"] & np.isfinite(lx) & np.isfinite(ly)
    surv = fin & (g["total_LL"] >= g["full_max"] + g["CL_2SIG"])
    relic = fin & (g["relic_LL"] >= g["relic_max"] + g["CL_2SIG"])

    lm, lsig = si_limits.envelope_log(n=8000)
    good = np.isfinite(lsig)
    mlo, mhi = lm[good].min(), lm[good].max()
    env = np.interp(lx, lm[good], lsig[good], left=np.nan, right=np.nan)
    covered = np.isfinite(env)

    lines = [
        f"run:                  {g['run_tag']}",
        f"SI envelope members:  {', '.join(si_limits.DEFAULT_ENVELOPE)}",
        f"envelope mass range:  {10**mlo:.4g} -- {10**mhi:.4g} GeV",
        "",
    ]
    for name, sel in (("relic-allowed (2 sigma)", relic),
                      ("survives everything (2 sigma)", surv)):
        n_tot = int(sel.sum())
        n_cov = int((sel & covered).sum())
        n_above = int((sel & covered & (ly > env)).sum())
        lines.append(f"{name}: {n_tot:,} samples")
        lines.append(f"    within the envelope's mass range: {n_cov:,}"
                     f"  (m_chi max = {10**lx[sel].max():.4g} GeV)")
        lines.append(f"    ALSO above the SI envelope (excluded by DD): {n_above:,}"
                     f"  ({(n_above / n_cov * 100) if n_cov else 0:.2f}% of covered)")
        if n_cov and n_above < n_cov:
            gap = (env - ly)[sel & covered & (ly <= env)]
            lines.append(f"    closest surviving approach to the bound: "
                         f"{np.nanmin(gap):.2f} dex below it")
        lines.append("")

    # The panel picks its sigma window from the 0.1/99.9 percentiles of the
    # relic-allowed samples (plot_all_constraints.make_plane), which for a thin
    # freeze-in locus can sit BELOW the DD bound -- the excluded samples are
    # then off the top of the frame and the region engine never rasterised them,
    # so the figure cannot show the exclusion even though it is real.
    ylo, yhi = PANEL.get("ylim", g["BOUNDS"]["sigma_chiN"])
    off = int((surv & covered & (ly > env) & (ly > yhi)).sum())
    n_above = int((surv & covered & (ly > env)).sum())
    if off:
        lines.append(
            f"WARNING: {off} of the {n_above} DD-excluded surviving samples sit "
            f"ABOVE the panel's sigma window (top = 10^{yhi:.0f} cm^2), so the figure "
            f"cannot show them.  The window comes from the relic-allowed "
            f"percentiles in plot_all_constraints.make_plane; widen it "
            f"(LOOP_SIGMA_YMAX) to include the DD curve.")
        lines.append("")

    p = os.path.join(out_dir, "dd_overlap.txt")
    with open(p, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join("  " + l for l in lines if l))
    print(f"  wrote {os.path.basename(p)}")
    return lines


# ---------------------------------------------------------------------------
# Running plot_all_constraints.py under the savefig wrapper
# ---------------------------------------------------------------------------
def cache_path_for(tag):
    return os.path.join(CACHE_DIR, f"{tag}__{PLANE.replace(':', '__')}.fig.pkl")


def run_tag_of(path):
    base = os.path.basename(path)
    for suf in ("_postprocessed.hdf5", ".hdf5"):
        if base.endswith(suf):
            return base[: -len(suf)]
    return base


def install_wrapper(out_dir, cache_path, prod_dir):
    """Wrap Figure.savefig: decorate the target plane, redirect every path.

    Also wraps builtins.open for writes into ``prod_dir``: write_legend saves
    legend.json with a plain open(), which savefig never sees, so without this
    this script would overwrite the main run's legend.json.
    """
    import builtins
    orig = Figure.savefig
    state = {"decorated": False}
    orig_open = builtins.open
    prod_dir = os.path.realpath(prod_dir)

    def guarded_open(file, mode="r", *a, **kw):
        if isinstance(file, (str, os.PathLike)) and any(c in mode for c in "wax+") \
                and os.path.dirname(os.path.realpath(file)) == prod_dir:
            file = os.path.join(out_dir, os.path.basename(file))
        return orig_open(file, mode, *a, **kw)

    builtins.open = guarded_open

    def savefig(self, fname, *a, **kw):
        name = os.path.basename(str(fname))
        dest = os.path.join(out_dir, name)
        first = TARGET in name and not state["decorated"]
        if first:
            state["decorated"] = True
            if cache_path:
                try:
                    make_picklable(self)
                    with open(cache_path, "wb") as f:
                        pickle.dump(self, f)
                    print(f"  cached undecorated figure -> {os.path.relpath(cache_path, HERE)}")
                except Exception as e:                       # pickling a figure is best-effort
                    print(f"  (figure not cacheable: {e})")
            decorate(self.axes[0])
        out = orig(self, dest, *a, **kw)
        if first and SHOW_ZOOM and dest.endswith(".pdf"):
            # orig, not self.savefig: the wrapper would redirect it a second time.
            save_zoom(self, dest[:-4] + "_zoom.pdf",
                      saver=lambda f, pth, **kw2: orig(f, pth, **kw2))
        return out

    Figure.savefig = savefig
    return state


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    from_cache = "--from-cache" in sys.argv[1:]

    if from_cache:
        tag = args[0] if args else run_tag_of(DEFAULT_RUN)
        cache_path = cache_path_for(tag)
        if not os.path.exists(cache_path):
            raise SystemExit(f"no cache at {cache_path} -- run once without --from-cache")
        with open(cache_path, "rb") as f:
            fig = pickle.load(f)
        out_dir = os.path.join(OUT_ROOT, f"{tag}_dd")
        os.makedirs(out_dir, exist_ok=True)
        decorate(fig.axes[0])
        p = os.path.join(out_dir, f"{TARGET}.pdf")
        fig.savefig(p, bbox_inches="tight")
        print(f"  wrote {os.path.relpath(p, HERE)} (from cache)")
        if SHOW_ZOOM:
            save_zoom(fig, p[:-4] + "_zoom.pdf")
        write_dd_legend(out_dir)
        return

    in_arg = args[0] if args else DEFAULT_RUN
    group = args[1] if len(args) > 1 else "data"
    if not os.path.exists(in_arg):
        raise SystemExit(f"{in_arg} not found")
    tag = run_tag_of(in_arg)
    out_dir = os.path.join(OUT_ROOT, f"{tag}_dd")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(CACHE_DIR, exist_ok=True)
    cache_path = cache_path_for(tag)

    print(f"run:       {in_arg}")
    print(f"plane:     {PLANE}")
    print(f"overlay:   mode={MODE} astro={SHOW_ASTRO} bands={SHOW_BANDS}")
    print(f"limits:    {len(si_limits.available())} curves in {si_limits.LIMITS_DIR}")
    print(f"output:    {out_dir}")

    install_wrapper(out_dir, cache_path,
                    os.path.join(PLOTTING, "plots", f"{tag}_allconstraints"))
    sys.argv = [os.path.join(PLOTTING, "plot_all_constraints.py"), in_arg, group, PLANE]
    g = runpy.run_path(sys.argv[0], run_name="__main__")
    write_dd_legend(out_dir)
    g["run_tag"] = tag
    try:
        report_overlap(g, out_dir)
    except Exception as e:
        print(f"  (overlap report skipped: {type(e).__name__}: {e})")


if __name__ == "__main__":
    main()
