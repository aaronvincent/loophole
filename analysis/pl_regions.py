#!/usr/bin/env python3
"""
Shared profile-likelihood region/boundary helpers for the Loophole
constraint plot (plot_all_constraints.py).

Why this module exists
----------------------
The constraint plot reduces the Diver scan to a 2D allowed/excluded region.  Two
earlier engines each failed in a characteristic way on a sparse, irregular
optimizer sample:

  * binned ``max`` profile + threshold -- "choppy in smooth places" (the per-bin
    max is a high-variance statistic) and "missing real features" (a coarse grid
    can't resolve a thin ridge finger);
  * binary alpha-shape of the cut-passing points -- because the region indicator
    is BINARY, every lone outlier sample gets full weight and the boundary chases
    it (wobble); suppressing that needs a heavy blur, which then rounds genuine
    features (over-smoothing).  Both complaints at once, unfixable by one blur.

The ENGINE here is a smoothed point-DENSITY level set, which breaks that
wobble-vs-feature tension by weighting the boundary by how many samples are at a
location, not merely whether any is.  Key fact: the profiled threshold
``profile_max(LL) - LL_max >= CL`` is exactly the PER-SAMPLE test
``LL[i] >= LL_max + CL`` -- a location is allowed iff SOME sample there passes.
So the region is the support of the cut-passing point cloud; we estimate its edge
robustly:

  1. ``point_density``  -- histogram the cut-passing points, gaussian-smooth the
                           COUNTS (small bandwidth), and define the level-set
                           threshold as the smoothed density of ``kpts`` points
                           per kernel.  Lone outliers (count 1) wash out below
                           threshold; dense real features (the ridge, thin
                           fingers) stay above it.  The true edge is a steep
                           density cliff -> the contour sits there SMOOTHLY (no
                           outlier chasing) yet keeps SHARP features (small
                           bandwidth, because the threshold -- not a big blur --
                           is what rejects speckle).
  2. ``density_region`` -- threshold the density -> boolean mask (for fills and
                           the set-difference kill categories), speckle pruned.
  3. ``field_boundaries``- trace the density field at the threshold -> smooth
                           boundary polylines for ax.plot (a level set of a smooth
                           field, so already smooth; no spline needed).
  4. ``fill_region``    -- contourf a category mask (smooth filled region).

``profile_grid`` / ``region_mask`` are binned-field helpers, not used by the
plotters.

All functions are pure (no global state); plotting helpers take an Axes.
Dependencies: numpy, scipy, contourpy (vendored by matplotlib) -- NOT skimage.
"""

import matplotlib
import matplotlib.legend_handler
import numpy as np
from scipy.stats import binned_statistic_2d
from scipy.ndimage import gaussian_filter, binary_fill_holes, label
from contourpy import contour_generator, LineType


# -------------------------------------------------------------------------
# 1. Per-bin max profile (extracted verbatim from the three scripts)
# -------------------------------------------------------------------------
def profile_grid(lx, ly, weight, fin, xedges, yedges):
    """Coarse (ny, nx) max-profiled grid; un-sampled bins are -inf.

    Identical reduction to the old inline ``profile_max``: within each (x,y)
    bin take the maximum of ``weight`` over all (profiled-out) points that fall
    in it.  ``fin`` is the validity-&-finite mask shared across all weights.
    """
    m = fin & np.isfinite(weight)
    stat, _, _, _ = binned_statistic_2d(
        lx[m], ly[m], weight[m], statistic="max", bins=[xedges, yedges])
    return np.where(np.isfinite(stat), stat, -np.inf).T   # -> (y, x)


# -------------------------------------------------------------------------
# 2. Clean boolean mask (speckle pruning)
# -------------------------------------------------------------------------
def prune_mask(m, min_area_frac=0.002, fill_holes=True):
    """Drop sub-min-area connected components and fill sub-min-area pinholes.

    Connected components smaller than ``min_area_frac`` of the grid are removed
    (the diamond islands), and interior pinholes smaller than the same area are
    filled.  Genuine large regions/voids are preserved, so the true outer
    envelope and any real disconnected island survive.
    """
    m = np.asarray(m, bool)
    if not m.any():
        return m
    min_px = max(1, int(round(min_area_frac * m.size)))

    lab, n = label(m)
    if n:
        sizes = np.bincount(lab.ravel())
        sizes[0] = 0                                   # background
        m = (sizes >= min_px)[lab]

    if fill_holes and m.any():
        holes = binary_fill_holes(m) & ~m
        labh, nh = label(holes)
        if nh:
            hsz = np.bincount(labh.ravel())
            hsz[0] = 0
            small = (hsz < min_px) & (hsz > 0)
            m = m | small[labh]
    return m


def region_mask(Zf, thr, min_area_frac=0.002, fill_holes=True):
    """Boolean allowed-mask {Zf >= thr} with speckle removed (see prune_mask)."""
    return prune_mask(np.isfinite(Zf) & (Zf >= thr), min_area_frac, fill_holes)


# -------------------------------------------------------------------------
# 2b. Point-density level-set engine (binning-free regions, robust edges)
# -------------------------------------------------------------------------
def point_density(px, py, xb, yb, grid_res=400, bw=2.5, kpts=2.0):
    """Smoothed cut-passing-sample density on a fine grid + the level threshold.

    Histograms the points onto a ``grid_res``^2 grid and gaussian-smooths the
    COUNTS with bandwidth ``bw`` (in fine cells).  The level-set threshold is the
    smoothed density of ``kpts`` points concentrated in one kernel,
    ``thr = kpts / (2*pi*bw^2)`` -- so the boundary marks "local passing density
    >= ~kpts samples".  Returns ``(xf, yf, dens, thr)``: 1D fine axes, the
    density field of shape (len(yf), len(xf)), and the scalar threshold.

    Weighting by COUNT (not a binary occupancy) is what removes the
    wobble-vs-feature trade-off: a lone outlier carries weight 1 and falls below
    ``thr`` without needing a big blur, so ``bw`` can stay small and keep sharp
    real features, while the populated bulk forms a steep density cliff whose
    threshold crossing sits at the true edge and barely moves laterally (smooth).
    """
    xf = np.linspace(xb[0], xb[1], grid_res)
    yf = np.linspace(yb[0], yb[1], grid_res)
    thr = kpts / (2.0 * np.pi * max(bw, 1e-6) ** 2)
    dens = np.zeros((grid_res, grid_res))

    xr = (xb[1] - xb[0]) or 1.0
    yr = (yb[1] - yb[0]) or 1.0
    u = (np.asarray(px) - xb[0]) / xr
    v = (np.asarray(py) - yb[0]) / yr
    ok = (np.isfinite(u) & np.isfinite(v)
          & (u >= 0) & (u <= 1) & (v >= 0) & (v <= 1))
    if ok.sum() < 1:
        return xf, yf, dens, thr

    iu = np.clip((u[ok] * (grid_res - 1)).astype(int), 0, grid_res - 1)
    iv = np.clip((v[ok] * (grid_res - 1)).astype(int), 0, grid_res - 1)
    cnt = np.zeros((grid_res, grid_res))
    np.add.at(cnt, (iv, iu), 1.0)                      # (y, x)
    dens = gaussian_filter(cnt, bw)
    return xf, yf, dens, thr


def xenvelope_region(px, py, xb, yb, grid_res=400, xbin=0.02, xgap=0.15):
    """Region for a THIN 1-D band: in every fine x column, everything between the
    lowest and highest sample within +-xbin/2 (x units) of it.

    The density engine cannot draw such a band at its true width: its kernel
    fattens the band by several bandwidths, and shrinking the kernel drops the
    sparse stretches of the band below threshold.  This traces the samples
    themselves.  Only valid where the region is one band that is single-valued
    in x at each x (an S-shaped curve is fine; a hole or two branches in one
    column is filled in).  Empty stretches up to ``xgap`` (x units) are bridged.
    Same return as ``density_region``: ``(xf, yf, mask,
    field, thr)``, the field being the mask lightly smoothed with thr = 0.5, so
    boundary tracing works unchanged.
    """
    xf = np.linspace(xb[0], xb[1], grid_res)
    yf = np.linspace(yb[0], yb[1], grid_res)
    mask = np.zeros((grid_res, grid_res), dtype=bool)
    px, py = np.asarray(px, float), np.asarray(py, float)
    ok = np.isfinite(px) & np.isfinite(py)
    if ok.any():
        dx = xf[1] - xf[0]
        col = np.clip(np.round((px[ok] - xb[0]) / dx).astype(int), 0, grid_res - 1)
        lo = np.full(grid_res, np.inf)
        hi = np.full(grid_res, -np.inf)
        np.minimum.at(lo, col, py[ok])
        np.maximum.at(hi, col, py[ok])
        from scipy.ndimage import minimum_filter1d, maximum_filter1d
        w = max(1, int(round(xbin / dx)) | 1)            # odd window, in columns
        lo = minimum_filter1d(lo, w, mode="constant", cval=np.inf)
        hi = maximum_filter1d(hi, w, mode="constant", cval=-np.inf)
        # Where the band is steep its samples are sparse in x: interpolate across
        # empty stretches up to `xgap` wide, then make neighbouring columns
        # overlap in y, so a near-vertical stretch is drawn continuous rather
        # than as a staircase of separate dashes.
        filled = np.flatnonzero(np.isfinite(lo))
        for a, b in zip(filled[:-1], filled[1:]):
            if 1 < b - a <= max(1, int(round(xgap / dx))):
                t = (np.arange(a + 1, b) - a) / (b - a)
                lo[a + 1:b] = lo[a] + t * (lo[b] - lo[a])
                hi[a + 1:b] = hi[a] + t * (hi[b] - hi[a])
        fin = np.isfinite(lo)
        both = fin[:-1] & fin[1:]
        lo2, hi2 = lo.copy(), hi.copy()
        j = np.flatnonzero(both)
        hi2[j] = np.maximum(hi2[j], lo[j + 1])          # rising: reach the next column
        lo2[j + 1] = np.minimum(lo2[j + 1], hi[j])
        hi2[j + 1] = np.maximum(hi2[j + 1], lo[j])      # falling
        lo2[j] = np.minimum(lo2[j], hi[j + 1])
        lo, hi = lo2, hi2
        dy = yf[1] - yf[0]
        for j in np.flatnonzero(np.isfinite(lo)):
            mask[(yf >= lo[j] - dy) & (yf <= hi[j] + dy), j] = True
    field = gaussian_filter(mask.astype(float), 0.7)
    return xf, yf, mask, field, 0.5


def density_region(px, py, xb, yb, grid_res=400, bw=2.5, kpts=2.0,
                   min_area_frac=0.0008):
    """Density level-set region: ``(xf, yf, mask, dens, thr)``.

    ``mask = dens >= thr`` (speckle-pruned) is the boolean region used for fills
    and the set-difference kill categories; ``dens, thr`` are returned so the
    boundary LINE can be traced from the smooth field (see field_boundaries),
    keeping line and fill edges consistent.
    """
    xf, yf, dens, thr = point_density(px, py, xb, yb, grid_res, bw, kpts)
    mask = dens >= thr
    if min_area_frac > 0:
        mask = prune_mask(mask, min_area_frac)
    return xf, yf, mask, dens, thr


# -------------------------------------------------------------------------
# 4. Smooth boundary polylines (level set of the density field)
# -------------------------------------------------------------------------
def _contour_segments(field, xf, yf, level):
    """Marching-squares polylines of {field >= level} (list of (N,2) arrays)."""
    cg = contour_generator(xf, yf, np.asarray(field, float),
                           line_type=LineType.Separate)
    return cg.lines(level)


def field_boundaries(field, xf, yf, level, min_len=8):
    """Smooth boundary polylines: the ``level`` contour of a smooth scalar field.

    The density field is already gaussian-smooth, so its threshold crossing is a
    smooth curve -- marching squares returns it directly with no spline needed.
    Loops shorter than ``min_len`` vertices (residual specks) are dropped.
    Returns a list of (M, 2) arrays in data coordinates, ready for ``ax.plot``.
    """
    out = []
    for seg in _contour_segments(field, xf, yf, level):
        seg = np.asarray(seg)
        if len(seg) < min_len:
            continue
        out.append(seg)
    return out


# -------------------------------------------------------------------------
# 5. Smooth filled category region
# -------------------------------------------------------------------------
# Horizontal dashes, staggered row to row.  matplotlib has no such hatch (and
# custom hatch characters stop working in 3.11), so fill_region draws them as
# dashed lines clipped to the fill, and HandlerDashes draws the legend swatch.
DASHES = "__"
_DASH_PITCH = 6.0             # row spacing [pt], as for the built-in '--'
_DASH_ON_OFF = (21.6, 14.4)   # dash / gap [pt]; odd rows shifted by half a period


def _dash_rows(ax, clip, color, zorder):
    """Staggered horizontal dashes across `ax`, clipped to the filled contour `clip`."""
    from matplotlib.collections import LineCollection
    from matplotlib.path import Path
    paths = [p for p in clip.get_paths() if len(p.vertices)]
    if not paths:
        return
    h_pt = ax.get_position().height * ax.figure.get_figheight() * 72.0
    n = max(1, int(round(h_pt / _DASH_PITCH)))
    ys = (np.arange(n) + 0.5) / n
    period = sum(_DASH_ON_OFF)
    lw = matplotlib.rcParams["hatch.linewidth"]
    lc = LineCollection([[(0.0, y), (1.0, y)] for y in ys], colors=[color],
                        linewidths=lw, transform=ax.transAxes, zorder=zorder,
                        linestyles=[((i % 2) * period / 2 / lw,
                                     tuple(d / lw for d in _DASH_ON_OFF))
                                    for i in range(n)])
    lc.set_clip_path(Path.make_compound_path(*paths), ax.transData)
    ax.add_collection(lc, autolim=False)


class HandlerDashes(matplotlib.legend_handler.HandlerPatch):
    """Legend swatch for DASHES: the patch plus two staggered dash rows."""

    def __init__(self, color, **kw):
        super().__init__(**kw)
        self.dash_color = color

    def create_artists(self, legend, orig, xd, yd, w, h, fontsize, trans):
        from matplotlib.lines import Line2D
        arts = super().create_artists(legend, orig, xd, yd, w, h, fontsize, trans)
        lw = matplotlib.rcParams["hatch.linewidth"]
        on, off = w * 0.3, w * 0.2
        for k, fy in enumerate((0.3, 0.75)):
            y = -yd + fy * h
            arts.append(Line2D([-xd, -xd + w], [y, y], color=self.dash_color, lw=lw,
                               ls=((on + off) / 2 * k / lw, (on / lw, off / lw)),
                               dash_capstyle="butt", transform=trans))
        return arts


def fill_region(ax, mask, xf, yf, color, alpha=1.0, zorder=2, hatch=None,
                hatch_color=None):
    """Draw a boolean category mask as a smooth filled contour (not imshow).
    `hatch` (a matplotlib pattern such as '//', or DASHES) overlays lines in
    `hatch_color`; matplotlib reads the hatch colour from rcParams when the
    artist is created."""
    if not np.asarray(mask).any():
        return None
    if not hatch or hatch == DASHES:
        cs = ax.contourf(xf, yf, mask.astype(float), levels=[0.5, 1.5],
                         colors=[color], alpha=alpha, zorder=zorder)
        if hatch == DASHES:
            _dash_rows(ax, cs, hatch_color or color, zorder)
        return cs
    with matplotlib.rc_context({"hatch.color": hatch_color or color}):
        return ax.contourf(xf, yf, mask.astype(float), levels=[0.5, 1.5],
                           colors=[color], alpha=alpha, zorder=zorder,
                           hatches=[hatch])
