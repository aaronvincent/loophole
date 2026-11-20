#!/usr/bin/env python3
"""
Spin-independent direct-detection PROJECTIONS (and the combined existing limit)
from ``data/projections/``, as closed regions in (log10 m_chi, log10 sigma).

The files are read-only reference data (see data/projections/readme.txt):

    column 1: log10(m / MeV)        column 2: log10(sigma / cm^2)

and apply to the heavy-mediator S(x)S case (contact limit).

Most of them are DIGITISED HATCHED REGIONS, not curves: the points trace the
region's outline and then zig-zag between its floor and ceiling (the hatch
strokes), in whatever order the digitiser visited them.  They cannot be drawn
as a polyline.  :func:`region_log` recovers the region as a floor and a ceiling
on a common mass grid by binning the points in log mass and taking the min / max
in each bin -- every hatch stroke puts one point on each edge, and the outline
dots fill the gaps where there are no strokes.

SuperCDMS is an ordinary curve (excluded above), kind ``"curve"``.
"""

import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
PROJ_DIR = os.environ.get("DD_PROJ_DIR", os.path.join(REPO, "data", "projections"))

# name -> (filename, kind, label)
#   kind "region": digitised hatched region, floor + ceiling
#   kind "curve" : excluded above the curve (no ceiling given)
CATALOG = {
    "existing": ("Combined Existing Spin Independent Limits.csv", "region",
                 "Existing SI limits"),
    # HeRALD Fig. 9 curves: 1 = 1 kg-day, 40 eV threshold; 2 = 1 kg-yr, 10 eV;
    # 3 = 10 kg-yr, 0.1 eV; 4 = 100 kg-yr, 1 meV.  Upper edge = 100 m overburden.
    # Curve 1 is HeRALD's recommendation (2026-10-26: closest to the TESSERACT
    # requirement).  Floor digitised from the figure; ceiling taken from the
    # earlier curve-4 file ("1 MeV 100 m.csv" -- the "MeV" there is the 1 meV
    # threshold).
    "HeRALD":   ("Herald Spin Independent curve 1 100 m.csv", "region",
                 "HeRALD (proj.)"),
    # The earlier HeRALD file: Fig. 9 curve 4 (100 kg-yr, 1 meV; reaches m ~ 1 MeV).
    # Not drawn by default; DDP_PROJ=HeRALD,HeRALD_c4,... shows both.
    "HeRALD_c4": ("Herald Spin Independent 1 MeV 100 m.csv", "region",
                  r"HeRALD 100 kg$\cdot$yr, 1 meV (proj.)"),
    "QUESTDMC": ("QUESTDMC Spin Independent Projection.csv", "region",
                 "QUEST-DMC (proj.)"),
    "SuperCDMS": ("SuperCDMS Spin Independent Projection.csv", "curve",
                  "SuperCDMS SNOLAB (proj.)"),
    # Spin-dependent, for the axial-vector (A(x)A) scans.  Both are plain curves
    # (excluded above).  The existing limit is read as SD-PROTON, which is what
    # the scans' sigma_chiN is (a_N built from the proton charges Delta_q^p);
    # QUEST-DMC is NEUTRON-only and bounds sigma_chin instead (see
    # plot_sigma_mx_ddproj.sd_reach).
    "existing_SD": ("Spin Dependent Limits.csv", "curve", "Existing SD-$p$ limits"),
    "QUESTDMC_SDn": ("QUEST DMC Neutron Only Spin Dependent.csv", "curve",
                     r"QUEST-DMC SD-$n$ (proj.; bounds $\sigma_{\chi n}$)"),
}
PROJECTIONS = ["HeRALD", "QUESTDMC", "SuperCDMS"]
SD_PROJECTIONS = ["QUESTDMC_SDn"]


def raw_log(name):
    """``(log10 m [GeV], log10 sigma [cm^2])`` in file order."""
    fn = CATALOG[name][0]
    d = np.loadtxt(os.path.join(PROJ_DIR, fn), delimiter=",")
    return d[:, 0] - 3.0, d[:, 1]            # MeV -> GeV


# The digitised "Combined Existing" floor sits ~0.1 dex above the published
# XENON1T-Migdal curve at m_chi = 0.1--0.13 GeV (where the survivors of the
# m_phi = 10 GeV scans sit within 0.02 dex of the bound).  So by default the
# existing floor is the pointwise best of the digitised file and the published
# curves the sibling prototype uses (si_limits.envelope_log); the file alone
# supplies the reach below 0.104 GeV and the ceiling.  DDP_EXIST_ENVELOPE=0: file only.
EXIST_WITH_ENVELOPE = os.environ.get("DDP_EXIST_ENVELOPE", "1") not in ("0", "", "false", "False")
# The existing region is extended to log10(m / GeV) = DDP_EXIST_EXTEND_TO with a
# flat ceiling (see region_log); "none" keeps the digitised mass reach only.
_ext = os.environ.get("DDP_EXIST_EXTEND_TO", "3").strip().lower()
EXIST_EXTEND_TO = -np.inf if _ext in ("none", "") else float(_ext)


def _edge_points(x, y, dx):
    """Per-bin (x_at_min, min) and (x_at_max, max), for bins spanning both edges."""
    edges = np.arange(x.min(), x.max() + dx, dx)
    idx = np.clip(np.digitize(x, edges) - 1, 0, edges.size - 2)
    fl, ce = [], []
    for b in range(edges.size - 1):
        sel = np.flatnonzero(idx == b)
        if sel.size < 2:
            continue
        i_lo, i_hi = sel[np.argmin(y[sel])], sel[np.argmax(y[sel])]
        if y[i_hi] - y[i_lo] < 0.3:
            continue
        fl.append((x[i_lo], y[i_lo]))
        ce.append((x[i_hi], y[i_hi]))
    return np.array(fl), np.array(ce)


def region_log(name, dx=0.03):
    """Recover ``(lm, floor, ceil)`` for one entry.

    ``region``: points are binned in log mass with width ``dx`` (dex); a bin
    contributes only if it holds points on BOTH edges (spread > 0.3 dex), so a
    bin that caught just the tip of the outline does not pinch the region.  Each
    bin's lowest point is a floor vertex and its highest a ceiling vertex, AT
    THEIR OWN MASSES (not the bin centre), so steps in the floor stay sharp;
    both are then interpolated onto one grid.
    ``curve``: sorted in mass, ``ceil`` is NaN.
    """
    x, y = raw_log(name)
    if CATALOG[name][1] == "curve":
        o = np.argsort(x)
        return x[o], y[o], np.full(x.size, np.nan)
    fl, ce = _edge_points(x, y, dx)
    lm = np.unique(np.r_[fl[:, 0], ce[:, 0]])
    lm = np.unique(np.r_[x.min(), lm, x.max()])
    lo = np.interp(lm, fl[:, 0], fl[:, 1])
    hi = np.interp(lm, ce[:, 0], ce[:, 1])
    # A bin that caught a vertical stretch of the outline but no hatch top reads
    # as a notch in the ceiling; a running max over +-dx removes it.  The floor
    # -- the sensitivity itself -- is not smoothed.
    hi = np.array([hi[np.abs(lm - m) <= dx].max() for m in lm])
    if name == "existing" and EXIST_WITH_ENVELOPE:
        import si_limits
        grid = np.linspace(lm.min(), lm.max(), 1500)
        lo_g, hi_g = np.interp(grid, lm, lo), np.interp(grid, lm, hi)
        em, es = si_limits.envelope_log(n=4000, mlim=(lm.min(), lm.max()))
        env = np.interp(grid, em, es, left=np.nan, right=np.nan)
        env = np.where(np.isfinite(env), env, np.inf)
        lm, lo, hi = grid, np.minimum(lo_g, env), hi_g
    if name == "existing" and EXIST_EXTEND_TO > lm.max():
        # The digitised region stops at ~1.9 GeV.  NAIVE extension to the right:
        # the ceiling (overburden) held flat at its last value; the floor never
        # weaker than its last digitised value, and the published envelope where
        # that is stronger.
        ext = np.linspace(lm.max(), EXIST_EXTEND_TO, 400)[1:]
        lo_ext = np.full(ext.size, lo[-1])
        if EXIST_WITH_ENVELOPE:
            import si_limits
            em, es = si_limits.envelope_log(n=2000, mlim=(ext.min(), ext.max()))
            if em.size:
                env = np.interp(ext, em, es, left=np.nan, right=np.nan)
                lo_ext = np.fmin(lo_ext, env)
        lm, lo, hi = np.r_[lm, ext], np.r_[lo, lo_ext], np.r_[hi, np.full(ext.size, hi[-1])]
    return lm, lo, hi


def polygon_log(name, **kw):
    """Closed outline ``(x, y)`` of a ``region`` entry (floor out, ceiling back)."""
    lm, lo, hi = region_log(name, **kw)
    return np.r_[lm, lm[::-1], lm[:1]], np.r_[lo, hi[::-1], lo[:1]]


def floor_at(name, lm_query):
    """The region's floor (the sensitivity) interpolated at ``lm_query``;
    NaN outside its mass support."""
    lm, lo, _ = region_log(name)
    return np.interp(lm_query, lm, lo, left=np.nan, right=np.nan)


def ceil_at(name, lm_query):
    """The region's ceiling at ``lm_query`` (+inf for a plain curve)."""
    lm, _, hi = region_log(name)
    if not np.isfinite(hi).any():
        out = np.where((lm_query >= lm.min()) & (lm_query <= lm.max()), np.inf, np.nan)
        return out
    return np.interp(lm_query, lm, hi, left=np.nan, right=np.nan)


def inside(name, lx, ly):
    """Boolean: sample points (log10 m, log10 sigma) inside the region."""
    f, c = floor_at(name, lx), ceil_at(name, lx)
    with np.errstate(invalid="ignore"):
        return np.isfinite(f) & (ly >= f) & (ly <= c)
